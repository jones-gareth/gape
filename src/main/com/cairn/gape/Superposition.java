package com.cairn.gape;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Logger;

import com.cairn.common.utils.CommonUtils;
import com.cairn.common.utils.Timer;
import com.cairn.gape.chromosome.BinaryAndIntegerChromosome;
import com.cairn.gape.chromosome.SuperpositionChromosome;
import com.cairn.gape.chromosome.SuperpositionCrossover;
import com.cairn.gape.chromosome.SuperpositionMutation;
import com.cairn.gape.feature.AcceptorAtomFeature;
import com.cairn.gape.feature.AromaticRingFeature;
import com.cairn.gape.feature.DonorHydrogenFeature;
import com.cairn.gape.feature.Feature;
import com.cairn.gape.feature.Feature.FeatureType;
import com.cairn.gape.feature.FeatureMapping;
import com.cairn.gape.feature.FeatureOverlay;
import com.cairn.gape.feature.HydrogenBondingType;
import com.cairn.gape.feature.UserFeatureSet;
import com.cairn.gape.ga.BaseSupervisor;
import com.cairn.gape.ga.IslandModel;
import com.cairn.gape.ga.Population;
import com.cairn.gape.ga.ProgressReport;
import com.cairn.gape.molecule.GaMolecule;
import com.cairn.gape.molecule.TorsionalDistributions;
import com.cairn.gape.utils.GaussianList;
import com.cairn.molecule.IncrementalSuperpositionComparison;
import com.cairn.molecule.Molecule;

/**
 * Superposition is the supervisor class for GAPE ({@link BaseSupervisor}). Runs
 * the GAPE (Genetic Algorithm for Pharmacophore Elucidation) Algorithm.
 * 
 * @author Gareth Jones
 * @version 2.0.0
 */
public class Superposition extends BaseSupervisor {

	private static Logger logger = Logger.getLogger(Superposition.class);
	static {
		// logger.setLevel(Level.DEBUG);
	}
	protected volatile GaMolecule baseMolecule, fittingMolecule;
	protected volatile List<GaMolecule> molecules;

	protected volatile Population pop;

	protected volatile SuperpositionChromosome best;
	protected final List<SuperpositionChromosome> solutions = new ArrayList<>();

	protected volatile TorsionalDistributions tordist;

	protected volatile int baseMoleculeNo, fittingMoleculeNo, integerEntryPoints[],
			binaryEntryPoints[], integerStringLength, integerStringRanges[],
			binaryStringLength, nRuns, nMapped;

	protected volatile boolean[] allowNulls;

	private volatile int nRebuilds, nIterations, popSize, nicheSize, crossoverWt,
			mutateWt, migrateWt, nIslands, reportInterval,
			nFeatures = Feature.N_BUILTIN_FEATURES, nUserFeatures = 0;

	protected volatile double acceptorAtomWt, aromaticRingWt, donorHydrogenWt, volumeWt,
			geometricWeightCriteria = .5, finishFittingRadius;

	private volatile double pharmFactor = 1.0, nicheDistance, confWt, startFittingRadius,
			selectPressure, constraintWt = 100, nichingOff;

	protected volatile boolean scalePharmacophore = false, incLonePairs,
			useFeatureClustering;

	protected volatile FeatureOverlay featureOverlay;

	private volatile boolean compareAll = false, singleMatchOnly = false,
			useNiches = false, nichesOn = false, scaleFittingOn = false,
			useTordistFile = true, solvate = true, addHydrogens = true,
			groupResults = false, scaleFitting, dumpCurrentBest, useActivities,
			doRigidFit = true, guessGAParameters = false;

	Molecule.FileType fileFormat = Molecule.FileType.MOL2;

	public static final int MAX_MOLS = 100;

	protected final Map<FeatureType, UserFeatureSet> userFeatures = new LinkedHashMap<>();

	static final boolean TEST_REFERENCE = false;
	private boolean canCompareToXtal, foundXtal;

	private static final String VERSION = "2.0.0";

	protected Timer totalTimer, runTimer;

	/**
	 * Constructor
	 */
	public Superposition() {
		super();
	}

	/**
	 * Initializes from the the default settings
	 */
	public void initFromResource() {
		initFromResource("/com/cairn/gape/utils/superposition.conf");
	}

	/**
	 * Runs GAPE. Requires a configuration file as an argument. Additional
	 * optional files contain input molecular structures. Just calls
	 * {@link #fitMolecules}.
	 * 
	 * @param args
	 *            filename arguments.
	 */
	public static void main(String args[]) {
		Superposition sp = null;
		try {
			if (args.length < 1) {
				System.out.println("Usage: " + Superposition.class.getName()
						+ " -version | <conf file> [molecule files..]");
				System.exit(0);
			}
			if (args.length == 1 && args[0].equals("-version")) {
				System.out.println("version = " + VERSION);
				System.exit(1);
			}

			sp = new Superposition();
			System.out.print("\n\nGAPE (c)\n");
			System.out.print("Gareth Jones\n\n");

			String molFiles[] = null;
			if (args.length > 1) {
				molFiles = new String[args.length - 1];
				for (int i = 1; i < args.length; i++)
					molFiles[i - 1] = args[i];
			}

			// add progress listener
			// sp.addProgressListener(new TestProgressListener());
			sp.fitMolecules(args[0], molFiles);
			sp.infoMessageLogger.finish();

		} catch (Exception ex) {
			String stackTrace = CommonUtils.getStackTrace(ex);
			try {
				sp.infoMessageln("Exception:");
				sp.infoMessageln(stackTrace);
			} catch (Exception e) {
				;
			}
			System.err.println(Superposition.class.getName() + ": Exception: "
					+ stackTrace);
			logger.error("Fatal exception: ", ex);
		}
	}

	/**
	 * Loads Torsional Distribution and hydrogen bonding definition files.
	 */
	protected void commonMolSetup() {
		if (useTordistFile) {
			tordist = new TorsionalDistributions(infoMessageLogger);
			if (hasKey("tordist_file")) {
				String tordistFile = fileName(getStringValue("tordist_file"));
				tordist.readTordistFile(tordistFile);
			} else
				tordist.readTordistFile();
		}

		if (hasKey("hydrogen_bonding_file")) {
			String hydrogenBondingFile = fileName(getStringValue("hydrogen_bonding_file"));
			HydrogenBondingType.loadParameters(hydrogenBondingFile);
		} else
			HydrogenBondingType.loadParameters();

		// TODO the code for user feature clustering was here- it has to be in
		// setup molecules as all features cannot be found until the molecules
		// are setup- in the old code this was effectively turned off- need to
		// check implications of this.

	}

	/**
	 * Initializes molecules for superposition. Removes lone pairs. Checks to
	 * see if any molecules are rigid or fixed. Assigns any activities. Matches
	 * torsional distributions. Calls {@link #findBaseMolecule()} and
	 * {@link GaMolecule#setup()}. Finds binary string length and entry points.
	 */
	public void setupMolecules(List<GaMolecule> m) {
		molecules = m;

		commonMolSetup();

		if (molecules.size() > MAX_MOLS) {
			infoMessageln("Maximum number of molecules that can be fitted is " + MAX_MOLS);
			infoMessageln("Fitting the first " + MAX_MOLS + " molecules");
			molecules = new ArrayList<>(molecules.subList(0, MAX_MOLS));
		}
		if (molecules.size() < 2) {
			throw new RuntimeException("GAPE needs at least 2 molecules: only "
					+ molecules.size() + " present");
		}

		String rigidMol = null;
		if (hasKey("rigid_mol"))
			rigidMol = getStringValue("rigid_mol");
		boolean allMoleculesRigid = false;
		if (hasKey("all_molecules_rigid") && getBooleanValue("all_molecules_rigid")) {
			allMoleculesRigid = true;
		}

		for (GaMolecule molecule : molecules) {
			// Remove pharmacophore information
			infoMessageln("\nProcessing " + molecule.getName());
			molecule.removeLonePairs();
			if (rigidMol != null && molecule.getName().equals(rigidMol)) {
				molecule.setRigid(true);
			}
			if (allMoleculesRigid) {
				molecule.setRigid(allMoleculesRigid);
			}

			if (molecule.isRigid())
				infoMessageln(molecule.getName() + " is rigid");
			if (molecule.isFixed())
				infoMessageln(molecule.getName() + " is fixed");

			molecule.setProblem(this);
			molecule.setup();

			if (useActivities) {
				molecule.setActivity(getActivity(molecule.getName()));
				infoMessageln("pKa activity is " + molecule.getActivity());
			}

			if (tordist != null)
				tordist.moleculeMatchTordist(molecule);
		}

		findBaseMolecule();

		int nBits = 0;
		binaryEntryPoints = new int[molecules.size()];
		for (int i = 0; i < molecules.size(); i++) {
			binaryEntryPoints[i] = nBits;
			GaMolecule molecule = molecules.get(i);
			nBits += molecule.getConformationalBitLength();
			if (molecule.isRelaxMolecule() && molecule != fittingMolecule)
				nBits += 6 * 8;
		}
		binaryStringLength = nBits;

		setupFeatureClustering();
	}

	/**
	 * Initializes the data structures required for 3D feature clustering
	 */
	protected void setupFeatureClustering() {
		if (useFeatureClustering) {
			// use feature clustering
			if (molecules.size() > 2) {
				featureOverlay = new FeatureOverlay(this);
				infoMessageln("Using feature clustering");
			} else {
				useFeatureClustering = false;
				infoMessageln("Only two compounds - turning off feature clustering");
			}
		}
	}

	/**
	 * Determine the length of the integer string portion.
	 */
	protected void determineIntegerStringLength() {
		integerStringLength = fittingMolecule.getNMappingFeatures() * (nMapped);
	}

	/**
	 * Finds the base molecule and fitting molecule for the overlay. Strategies
	 * for finding the base molecule include maximum features, most activity or
	 * minimum rotatable bonds. Any rigid or fixed molecule will automatically
	 * be chosen as the base molecule.
	 */
	public void findBaseMolecule() {
		baseMolecule = molecules.get(0);
		fittingMolecule = molecules.get(0);

		boolean rigidBase = false, multipleRigid = false;
		for (GaMolecule molecule : molecules) {
			if (molecule.isRigid()) {
				if (rigidBase) {
					multipleRigid = true;
					baseMolecule = molecules.get(0);
					fittingMolecule = molecules.get(0);
					rigidBase = false;
					break;
				}
				infoMessageln("Base molecule is rigid");
				baseMolecule = fittingMolecule = molecule;
				rigidBase = true;
			}
		}

		// I thought for feature clustering it would be best to use the least
		// flexible molecule as the fitting molecule, but this doesn't seem to
		// work in practice.
		boolean featureClusteringUseLeastFlexible = false;
		String baseMoleculeSelection = getStringValue("base_molecule_selection")
				.toLowerCase();

		// Base molecule selection
		if (rigidBase) {
			;
		}

		else if (multipleRigid) {
			// in the case of multiple rigids use max features whn min rotatable
			// bonds is selected
			if (baseMoleculeSelection.equals("min_rotatable_bonds")
					|| baseMoleculeSelection.equals("max_features")) {
				infoMessageln("Base molecule selected on " + "maximum features");
				baseMolecule = mostFeaturedRigidMolecule();
			} else if (baseMoleculeSelection.equals("max_activity")) {
				infoMessageln("Base molecule selected on " + "maximum activity");
				if (!useActivities)
					throw new RuntimeException("Base molecule selection method is "
							+ "maxmimum activity, but use_activities is false");
				baseMolecule = mostActiveRigidMolecule();
			} else
				throw new RuntimeException("Unknown selection strategy "
						+ baseMoleculeSelection);
		}

		else if (featureClusteringUseLeastFlexible && useFeatureClustering) {
			infoMessageln("Using feature clustering - base molecule is least flexible");
			if (MultiMolSuperposition.class.isAssignableFrom(this.getClass()))
				throw new RuntimeException(
						"Can't use least flexible base molecule criteria for multi conformer superposition");
			baseMolecule = leastFlexibleMolecule();
		}

		else if (hasKey("base_molecule")) {
			String name = getStringValue("base_molecule");
			infoMessageln("Base molecule is specified in configuration file");
			Optional<? extends GaMolecule> result = molecules.stream()
					.filter(m -> m.getName().equals(name)).findFirst();
			if (result.isPresent()) {
				baseMolecule = result.get();
			} else {
				throw new RuntimeException("Can't find base molecule " + name);
			}
		}

		else {
			if (baseMoleculeSelection.equals("min_rotatable_bonds")) {
				infoMessageln("Base molecule selected " + "on minimum rotatable bonds");
				if (MultiMolSuperposition.class.isAssignableFrom(this.getClass()))
					throw new RuntimeException(
							"Can't use least flexible base molecule criteria for multi conformer superposition");
				baseMolecule = leastFlexibleMolecule();
			} else if (baseMoleculeSelection.equals("max_features")) {
				infoMessageln("Base molecule selected on " + "maximum features");
				baseMolecule = mostFeaturedMolecule();
			} else if (baseMoleculeSelection.equals("max_activity")) {
				infoMessageln("Base molecule selected on " + "maximum activity");
				if (!useActivities)
					throw new RuntimeException("Base molecule selection method is "
							+ "maxmimum activity, but use_activities is false");
				baseMolecule = mostActiveMolecule();
			} else
				throw new RuntimeException("Unknown selection strategy "
						+ baseMoleculeSelection);
		}

		String fittingMoleculeSelection = null;
		if (hasKey("fitting_molecule_selection")) {
			fittingMoleculeSelection = getStringValue("fitting_molecule_selection")
					.toLowerCase();
		}

		// fitting molecule selection
		if (rigidBase) {
			;
		}

		else if (multipleRigid) {
			if (fittingMoleculeSelection != null) {
				if (fittingMoleculeSelection.equals("min_rotatable_bonds")
						|| fittingMoleculeSelection.equals("max_features")) {
					infoMessageln("Fitting molecule selected on " + "maximum features");
					fittingMolecule = mostFeaturedRigidMolecule();
				} else if (fittingMoleculeSelection.equals("base_molecule")) {
					infoMessageln("Fitting molecule is base molecule");
					fittingMolecule = baseMolecule;
				} else
					throw new RuntimeException("Unknown selection strategy "
							+ fittingMoleculeSelection);
			} else
				fittingMolecule = baseMolecule;
		}

		else if (featureClusteringUseLeastFlexible && useFeatureClustering) {
			infoMessageln("Using feature clustering - fitting molecule is least flexible");
			fittingMolecule = leastFlexibleMolecule();
		}

		else {
			if (fittingMoleculeSelection != null) {
				if (fittingMoleculeSelection.equals("min_rotatable_bonds")) {
					infoMessageln("Fitting molecule selected "
							+ "on minimum rotatable bonds");
					if (MultiMolSuperposition.class.isAssignableFrom(this.getClass()))
						throw new RuntimeException(
								"Can't use least flexible fitting molecule criteria for multi conformer superposition");
					fittingMolecule = leastFlexibleMolecule();
				} else if (fittingMoleculeSelection.equals("max_features")) {
					infoMessageln("Fitting molecule selected on " + "maximum features");
					fittingMolecule = mostFeaturedMolecule();
				} else if (fittingMoleculeSelection.equals("base_molecule")) {
					infoMessageln("Fitting molecule is base molecule");
					fittingMolecule = baseMolecule;
				} else
					throw new RuntimeException("Unknown selection strategy "
							+ fittingMoleculeSelection);
			} else
				fittingMolecule = baseMolecule;
		}

		infoMessageln("Base molecule is " + baseMolecule.getName());
		infoMessageln("Fitting molecule is " + fittingMolecule.getName() + "\n");

		nMapped = (int) molecules.stream()
				.filter(m -> !(m.isFixed() || m == fittingMolecule)).count();

		integerEntryPoints = new int[nMapped];
		determineIntegerStringLength();
		int no = 0;
		for (int i = 0; i < molecules.size(); i++) {
			GaMolecule molecule = molecules.get(i);
			if (molecule == baseMolecule) {
				baseMoleculeNo = i;
			}
			if (molecule == fittingMolecule) {
				fittingMoleculeNo = i;
				continue;
			}
			if (molecule.isFixed())
				continue;
			integerEntryPoints[no] = no * fittingMolecule.getNMappingFeatures();
			no++;
		}
		integerStringRanges = new int[integerStringLength];
		int pos = 0;
		Map<FeatureType, FeatureMapping> fittingFeatures = fittingMolecule
				.getFeatureMappings();
		for (GaMolecule molecule : molecules) {
			if (molecule == fittingMolecule)
				continue;
			if (molecule.isFixed())
				continue;
			for (Entry<FeatureType, FeatureMapping> entry : fittingFeatures.entrySet()) {
				if (!entry.getValue().isMapping())
					continue;
				FeatureMapping otherFeatures = molecule
						.getFeatureMappings(entry.getKey());
				for (int k = 0; k < entry.getValue().getNFeatures(); k++) {
					integerStringRanges[pos] = otherFeatures.getNFeatures();
					pos++;
				}
			}
		}

		setWeights();
	}

	/**
	 * Finds and returns the least flexible molecule.
	 */
	public GaMolecule leastFlexibleMolecule() {
		return molecules
				.stream()
				.min(Comparator.comparing(m -> m.getnRotatableBonds()
						+ m.getNFreeCorners())).get();
	}

	/**
	 * Finds and returns the most featured molecule.
	 */
	public GaMolecule mostFeaturedMolecule() {
		return molecules.stream().max(Comparator.comparing(GaMolecule::getNFeatures))
				.get();
	}

	/**
	 * Finds and returns the most active molecule
	 */
	public GaMolecule mostActiveMolecule() {
		return molecules.stream().max(Comparator.comparing(GaMolecule::getActivity))
				.get();
	}

	/**
	 * Finds and returns the most featured rigid molecule.
	 */
	public GaMolecule mostFeaturedRigidMolecule() {
		return molecules.stream().filter(GaMolecule::isRigid)
				.max(Comparator.comparing(GaMolecule::getNFeatures)).get();
	}

	/**
	 * Finds and returns the most active rigid molecule
	 */
	public GaMolecule mostActiveRigidMolecule() {
		return molecules.stream().filter(GaMolecule::isRigid)
				.max(Comparator.comparing(GaMolecule::getActivity)).get();
	}

	/**
	 * Reads in run-time setting from the configuration file. Prints out
	 * compiled settings and runtime settings.
	 * 
	 * @param configFile
	 *            Configuration <name.
	 */
	@Override
	public void init(String configFile) {
		super.init(configFile);

		if (infoMessageLogger.getLogLevel() > 1) {
			infoMessageln("\nCompiled Settings:");
			infoMessageln("Gaussian operations use exp lookup \t"
					+ GaussianList.USE_EXP_LOOKUP);
		}

		compareAll = getBooleanValue("compare_all");
		singleMatchOnly = getBooleanValue("single_match_only");
		// if (hasKey("use_feature_clustering"))
		useFeatureClustering = getBooleanValue("use_feature_clustering");

		useNiches = getBooleanValue("use_niches");
		nicheDistance = getDoubleValue("niche_distance");
		nichingOff = getDoubleValue("niching_off");

		volumeWt = getDoubleValue("volume_wt");
		confWt = getDoubleValue("conf_wt");
		donorHydrogenWt = getDoubleValue("donor_hydrogen_wt");
		acceptorAtomWt = getDoubleValue("acceptor_atom_wt");
		aromaticRingWt = getDoubleValue("aromatic_ring_wt");
		if (hasKey("constraint_wt"))
			constraintWt = getDoubleValue("constraint_wt");
		scalePharmacophore = getBooleanValue("scale_pharmacophore");
		pharmFactor = getDoubleValue("pharm_factor");
		geometricWeightCriteria = getDoubleValue("geometric_weight_criteria");
		useTordistFile = getBooleanValue("use_tordist_file");

		if (hasKey("add_hydrogens"))
			addHydrogens = getBooleanValue("add_hydrogens");
		if (hasKey("solvate"))
			solvate = getBooleanValue("solvate");

		Molecule.setAddHydrogensFlag(addHydrogens);
		Molecule.setSolvateFlag(solvate);
		Molecule.setChargeFlag(true);
		Molecule.setAssignTypesFlag(true);

		nIterations = getIntValue("n_iterations");
		popSize = getIntValue("popsize");
		selectPressure = getDoubleValue("select_pressure");
		nicheSize = getIntValue("niche_size");
		crossoverWt = getIntValue("crossover_wt");
		mutateWt = getIntValue("mutate_wt");
		migrateWt = getIntValue("migrate_wt");
		nIslands = getIntValue("n_islands");
		guessGAParameters = getBooleanValue("guess_ga_parameters");
		nRuns = getIntValue("n_runs");
		reportInterval = getIntValue("report_interval");
		dumpCurrentBest = getBooleanValue("dump_current_best");
		groupResults = getBooleanValue("group_results");
		incLonePairs = getBooleanValue("inc_lone_pairs");
		if (hasKey("do_rigid_fit"))
			doRigidFit = getBooleanValue("do_rigid_fit");

		if (hasKey("use_activities"))
			useActivities = getBooleanValue("use_activities");

		for (int i = 1; i <= 10; i++) {
			String key = "user_features_file_" + i;
			if (hasKey(key)) {
				UserFeatureSet set = new UserFeatureSet(this);
				String parameterFile = fileName(getStringValue(key));
				set.loadParameterFile(parameterFile);
			}
		}

		double hBondLen = getDoubleValue("h_bond_len");
		double solventVolOk = getDoubleValue("solvent_vol_ok");
		double atomSolventRadius = getDoubleValue("atom_solvent_radius");
		double chargeFactor = getDoubleValue("charge_factor");
		double matchFactor = getDoubleValue("match_factor");
		double normalLen = getDoubleValue("normal_len");
		double maximumGaussianScore = getDoubleValue("maximum_gaussian_score");
		boolean scoreDonorAtoms = getBooleanValue("score_donor_atoms");
		double maxDonorDonorAngle = getDoubleValue("max_donor_donor_angle");
		double minDonorDonorAngle = getDoubleValue("min_donor_donor_angle");
		boolean scaleLonePairs = getBooleanValue("scale_lone_pairs");

		double maxForwardAcceptorAngle = getDoubleValue("max_forward_acceptor_angle");
		double minForwardAcceptorAngle = getDoubleValue("min_forward_acceptor_angle");
		double maxLonePairLonePairAngle = getDoubleValue("max_lone_pair_lone_pair_angle");
		double minLonePairLonePairAngle = getDoubleValue("min_lone_pair_lone_pair_angle");
		double maxPlaneLonePairAngle = getDoubleValue("max_plane_lone_pair_angle");
		double minPlaneLonePairAngle = getDoubleValue("min_plane_lone_pair_angle");
		double maxPlanePlaneAngle = getDoubleValue("max_plane_plane_angle");
		double minPlanePlaneAngle = getDoubleValue("min_plane_plane_angle");

		Feature.setMaximumGaussianScore(maximumGaussianScore);
		Feature.setSolvationAlpha(atomSolventRadius);
		DonorHydrogenFeature.setHBondLen(hBondLen);
		Feature.setSolventVolOk(solventVolOk);
		DonorHydrogenFeature.setChargeFactor(chargeFactor);
		DonorHydrogenFeature.setMatchFactor(matchFactor);
		DonorHydrogenFeature.setScoreDonorAtoms(scoreDonorAtoms);
		DonorHydrogenFeature.setMaxDonorDonorAngle(maxDonorDonorAngle * Math.PI / 180);
		DonorHydrogenFeature.setMinDonorDonorAngle(minDonorDonorAngle * Math.PI / 180);

		AcceptorAtomFeature.setHBondLen(hBondLen);
		AcceptorAtomFeature.setChargeFactor(chargeFactor);
		AcceptorAtomFeature.setScaleLonePairs(scaleLonePairs);

		AcceptorAtomFeature.setMaxForwardAcceptorAngle(maxForwardAcceptorAngle * Math.PI
				/ 180);
		AcceptorAtomFeature.setMinForwardAcceptorAngle(minForwardAcceptorAngle * Math.PI
				/ 180);
		AcceptorAtomFeature.setMaxLonePairLonePairAngle(maxLonePairLonePairAngle
				* Math.PI / 180);
		AcceptorAtomFeature.setMinLonePairLonePairAngle(minLonePairLonePairAngle
				* Math.PI / 180);
		AcceptorAtomFeature.setMaxPlaneLonePairAngle(maxPlaneLonePairAngle * Math.PI
				/ 180);
		AcceptorAtomFeature.setMinPlaneLonePairAngle(minPlaneLonePairAngle * Math.PI
				/ 180);
		AcceptorAtomFeature.setMaxPlanePlaneAngle(maxPlanePlaneAngle * Math.PI / 180);
		AcceptorAtomFeature.setMinPlanePlaneAngle(minPlanePlaneAngle * Math.PI / 180);

		AcceptorAtomFeature.setMatchFactor(matchFactor);
		AromaticRingFeature.setNormalLen(normalLen);

		startFittingRadius = getDoubleValue("start_fitting_radius");
		finishFittingRadius = getDoubleValue("finish_fitting_radius");
		nRebuilds = getIntValue("n_rebuilds");
		scaleFitting = getBooleanValue("scale_fitting");
		useNiches = getBooleanValue("use_niches");

	}

	/**
	 * @return the name of the chromosome class used by the superposition
	 *         algorithm
	 */
	protected String chromosomeClassName() {
		return "com.cairn.gape.chromosome.SuperpositionChromosome";
	}

	/**
	 * Performs a single GA superposition run. See {@link Population} and
	 * {@link IslandModel}.
	 */
	public void run(int runNo) {

		if (scaleFitting)
			Feature.setRadius(startFittingRadius);
		else
			Feature.setRadius(finishFittingRadius);
		infoMessageln("Setting fitting radius to " + Feature.getRadius());
		nichesOn = useNiches;
		scaleFittingOn = scaleFitting;

		canCompareToXtal = TEST_REFERENCE && (new File("xtal.sdf")).exists()
				&& groupResults;
		foundXtal = false;

		// pop = new LinkedPopLinearSel(this);
		pop = new IslandModel(this, nIslands, migrateWt);
		pop.setSelectPressure(selectPressure);
		pop.setSize(popSize);
		pop.setNicheSize(nicheSize);
		pop.create(chromosomeClassName(), this);

		if (logger.isDebugEnabled())
			logger.debug(pop.popInfo());
		best = (SuperpositionChromosome) pop.getBest();

		pop.addOperator(new SuperpositionCrossover(this), crossoverWt);
		pop.addOperator(new SuperpositionMutation(this), mutateWt);

		double d = nIterations * nichingOff;
		int oppOff = (int) d;
		infoMessageln("Turning off niching and fitting radius scaling at " + oppOff);
		int rebuildInterval = oppOff / nRebuilds;

		for (int i = 0; i < nIterations; i++) {

			if (cancel) {
				throw new RuntimeException("Program terminated at user request");
			}

			// Turn off niching after 80% to allow convergence
			if (nichesOn && i >= oppOff) {
				infoMessageln("OP " + i + ": turning off niching");
				nichesOn = false;
			}

			// anneal fitting point radius
			if (scaleFittingOn && i > 0 && i % rebuildInterval == 0) {
				if (i >= oppOff) {
					Feature.setRadius(finishFittingRadius);
					scaleFittingOn = false;
				} else {
					double frac = ((double) (oppOff - i)) / ((double) oppOff);
					double radius = finishFittingRadius
							+ (startFittingRadius - finishFittingRadius) * frac;
					Feature.setRadius(radius);
				}
				infoMessageln("Setting fitting radius to "
						+ nf.format(Feature.getRadius()));
				pop.rebuild();
			}

			// report progress to any listeners
			if (CollectionUtils.isNotEmpty(progressListeners) && i % 5000 == 0) {
				double percent = ((double) (i + nIterations * runNo))
						/ ((double) (nIterations * nRuns));
				String message = "Run " + String.valueOf(runNo + 1) + ", Operation "
						+ String.valueOf(i);
				percent *= 100;
				tellListeners(new ProgressReport(percent, message));
			}

			pop.iterate();
			if (reportInterval > 0 && i > 1 && (i + 1) % reportInterval == 0)
				infoMessage(pop.info());

			SuperpositionChromosome testBest = (SuperpositionChromosome) pop.getBest();
			if (testBest != best) {
				best = testBest;
				if (dumpCurrentBest)
					outputSolution(best, "Current Best", false);
				if (canCompareToXtal) {
					boolean test = compareToXtal();
					if (test) {
						foundXtal = true;
						infoMessageln("Matches xtal structure operation " + i);
					}
				}
			}
		}

		if (logger.isDebugEnabled()) {
			logger.debug(pop.popInfo());
		}

		if (dumpCurrentBest) {
			for (GaMolecule molecule : molecules) {
				String n = "Current Best" + molecule.getName() + ".mol2";
				n = fileName(n.replace(' ', '_'));
				File f = new File(n);
				if (f.exists())
					f.delete();
			}
		}

		if (canCompareToXtal) {
			boolean test = compareToXtal();
			if (test) {
				foundXtal = true;
				infoMessageln("Solution matches xtal structure");
			} else if (foundXtal)
				infoMessageln("WARNING: found xtal in run but solution does not match xtal");
		}

		File f = new File(fileName("Current_Best_Pharm.mol2"));
		if (f.exists())
			f.delete();
		f = new File(fileName("Current_Best.mol2"));
		if (f.exists())
			f.delete();
	}

	/**
	 * Compare the current overlay to a crystal structure template or reference
	 * file
	 * 
	 * @return
	 * 
	 */
	public boolean compareToXtal() {
		String file = outputSolution(best, "Current Best", false);
		IncrementalSuperpositionComparison superpositonComparison = new IncrementalSuperpositionComparison(
				true, false, false);

		superpositonComparison.loadMolecules(file, "xtal.sdf");
		boolean test = superpositonComparison.buildOverlays();
		return test;
	}

	/**
	 * cleans the current directory of GAPE output files. Removes the most
	 * common output files, but depending on options may not remove all files.
	 * 
	 * Doesn't remove seed file as that is created before this routine is called
	 */
	public void cleanSuperposition() {
		File files[] = (new File(StringUtils.isNotEmpty(directory) ? directory : "."))
				.listFiles();
		for (File file : files) {
			String name = file.getName();
			if ((name.startsWith("GA_") && name.endsWith(".mol2"))
					|| (name.startsWith("GA_") && name.endsWith(".sdf"))
					|| (name.startsWith("Randomized_") && name.endsWith(".mol2"))
					|| (name.startsWith("Current_") && name.endsWith(".mol2"))) {
				if (file.exists())
					file.delete();
			}
		}
	}

	/**
	 * Copies the current best chromosome, rebuilds then returns it.
	 */
	public SuperpositionChromosome saveBest() {
		SuperpositionChromosome b = (SuperpositionChromosome) BinaryAndIntegerChromosome
				.getFreedChromosome();
		pop.getBest().copyGene(b);
		b.rebuild();
		return b;
	}

	/**
	 * Loads in the molecules from a file;
	 * 
	 * @param molFiles
	 * @return
	 * 
	 */
	protected List<GaMolecule> loadMolecules(String[] molFiles) {
		List<GaMolecule> m = GaMolecule.loadFiles(molFiles, Molecule.FileType.UNKNOWN,
				Molecule.Source.FILE, infoMessageLogger, true);
		fileFormat = Molecule.getType(molFiles[0]);
		return m;
	}

	/**
	 * Runs GAPE. Loads up the files. Repeatedly runs the algorithm saving the
	 * best solutions. Sorts solutions. Output to ga.log and structure files
	 * 
	 * @param configFile
	 *            configuration file
	 * @param molFiles
	 *            array of structure files
	 */
	public void fitMolecules(String configFile, String molFiles[]) {

		final File statusFile = new File("gape.complete");
		if (statusFile.exists()) {
			statusFile.delete();
		}

		totalTimer = new Timer();
		runTimer = new Timer();
		Timer cpuTimer = new Timer();

		init(configFile);
		cleanSuperposition();

		infoMessage("\n\nGAPE (c) \n");
		infoMessage("Gareth Jones gjones@cairn.com\n\n");

		if (molFiles == null) {
			if (!hasKey("molecules"))
				throw new RuntimeException("No structure files [supply molecule files "
						+ "in argument list or configuration file]");
			molFiles = new String[] { getValue("molecules") };
		}

		for (int i = 0; i < molFiles.length; i++)
			molFiles[i] = fileName(molFiles[i]);
		List<GaMolecule> m = loadMolecules(molFiles);

		setupMolecules(m);

		// now we have the number of molecules guess ga parameters if
		// required.
		if (guessGAParameters) {
			infoMessageln("Guessing GA Parameters:");
			nIslands = molecules.size() + 1;
			nIterations = molecules.size() * 15000;
			popSize = 100;
			infoMessageln("No Islands " + nIslands + " Popsize " + popSize
					+ " no iterations " + nIterations);
		}

		solutions.clear();

		for (int j = 0; j < nRuns; j++) {
			runTimer.reset();
			cpuTimer.reset();

			run(j);

			best = saveBest();
			int no = j + 1;
			infoMessage("\n\nSolution " + no + " " + best.fitnessInfo());
			outputSolution(best, "GA " + StringUtils.leftPad(String.valueOf(no), 3, '0'),
					true);
			if (!useFeatureClustering)
				infoMessage(best.mappingInfo());

			solutions.add(best);
			infoMessageln();
			infoMessageln();
			pop.freePop();

			runTimer.interval();
			infoMessageln("Run " + no + " Thread Execution Time " + runTimer.info()
					+ "\n\n");
			runTimer.reset();

			cpuTimer.interval();
			infoMessageln("Run " + no + " CPU Execution Time " + cpuTimer.info() + "\n\n");
			cpuTimer.reset();
		}

		if (nRuns > 1) {
			Map<SuperpositionChromosome, Integer> initialPositions = new HashMap<>();
			for (int i = 0; i < solutions.size(); i++) {
				initialPositions.put(solutions.get(i), i);
			}
			Collections.sort(solutions, chromosomeComparator());
			infoMessageln();
			RigidFit rigidFit = null;
			if (doRigidFit)
				infoMessageln("Aligning solutions to best with rigid body fitting");
			infoMessageln("Ranked order of solutions:");
			int rank = 0;
			for (SuperpositionChromosome solution : solutions) {
				solution.rebuild();
				if (doRigidFit) {
					if (rank == 0)
						rigidFit = new RigidFit(solution);
					else
						rigidFit.fit(solution);
				}
				rank++;
				int no = initialPositions.get(solution) + 1;
				infoMessage("Rank " + rank + " solution " + no + " "
						+ solution.fitnessInfo());
				outputSolution(solution,
						"GA ranked " + StringUtils.leftPad(String.valueOf(rank), 3, '0'),
						false);
			}
		}

		totalTimer.interval();
		infoMessageln("Total Thread Execution Time " + totalTimer.info());

		try (Writer statusWriter = new FileWriter(statusFile)) {
			statusWriter.write("Difgape search completed "
					+ SimpleDateFormat.getInstance().format(new Date()) + "\n");
		} catch (IOException e) {
			logger.error("IOException saving status file", e);
			throw new RuntimeException(e);
		}

	}

	/**
	 * Writes the overlay encoded in a chromosome out to a molecule file or
	 * files (depends on the class variable groupResults). Returns the fileName
	 * if group results is set.
	 * 
	 * @param prefix
	 *            used to generate filnames and molecule names within files
	 * @param log
	 *            set to log file output
	 */
	public String outputSolution(SuperpositionChromosome c, String prefix, boolean log) {
		if (SuperpositionChromosome.getCurrentChromosome() != c)
			c.rebuild();
		try {
			String fileName = null;
			// addPharmInfo should be included in c.rebuild
			if (scalePharmacophore)
				c.addPharmInfo();
			FileWriter outFile = null;
			if (groupResults) {
				String file = prefix
						+ (fileFormat == Molecule.FileType.SDF ? ".sdf" : ".mol2");
				file = file.replace(' ', '_');
				file = fileName(file);
				outFile = new FileWriter(file);
				if (log)
					infoMessageln("Writing solution to " + file);
				fileName = file;
			}
			for (GaMolecule molecule : molecules) {
				// for (int i = 0; i < nMolecules; i++) {
				if (!groupResults) {
					String file = prefix + " " + molecule.getName()
							+ (fileFormat == Molecule.FileType.SDF ? ".sdf" : ".mol2");
					file = file.replace(' ', '_');
					file = fileName(file);
					outFile = new FileWriter(file);
					if (log)
						infoMessageln("Writing solution to " + file);
				}
				String name = molecule.getName();
				if (molecule == baseMolecule)
					molecule.setName(prefix + " " + name + " Base Molecule");
				else if (molecule == fittingMolecule)
					molecule.setName(prefix + " " + name + " Fitting Molecule");
				else
					molecule.setName(prefix + " " + name);

				// add conformational energy to sd fields
				if (confWt > .0) {
					molecule.removeSdfField("CONFORMATIONAL_ENERGY");
					molecule.addSdfField("CONFORMATIONAL_ENERGY",
							String.valueOf(molecule.getConformationalEnergy()));
				}
				// need a write routine based on file format
				if (fileFormat == Molecule.FileType.SDF)
					molecule.writeSdfMol(outFile, "GA Superposition Molecule");
				else
					molecule.writeSybylMol2(outFile, "GA Superposition Molecule",
							incLonePairs);

				molecule.setName(name);
				if (!groupResults)
					outFile.close();
			}
			if (scalePharmacophore) {
				if (!groupResults) {
					String file = prefix + " Pharm"
							+ (fileFormat == Molecule.FileType.SDF ? ".sdf" : ".mol2");
					file = file.replace(' ', '_');
					file = fileName(file);
					outFile = new FileWriter(file);
				}
				c.outputPharmMol(outFile, prefix + " Pharmacophore");
				if (!groupResults)
					outFile.close();
			}
			if (groupResults)
				outFile.close();
			return fileName;
		}

		catch (IOException ex) {
			throw new RuntimeException("writeMol2File IO exception: " + ex);
		}
	}

	/**
	 * Creates a comparator for comparing superpositon chromosomes
	 * 
	 * @param o1
	 *            {@link SuperpositionChromosome} class instance
	 * @param o2
	 *            {@link SuperpositionChromosome} class instance
	 */
	public Comparator<SuperpositionChromosome> chromosomeComparator() {
		return (c1, c2) -> Double.compare(c1.getFitness(), c2.getFitness());
	}

	/**
	 * Sets activity weights
	 */
	public void setWeights() {
		if (!useActivities) {
			for (GaMolecule molecule : molecules) {
				molecule.setWeight(1);
			}
			return;
		}

		if (!compareAll) {
			double totalAct = molecules
					.stream()
					.filter(m -> m != baseMolecule)
					.mapToDouble(
							m -> (m.getActivity() + baseMolecule.getActivity()) / 2.0)
					.sum();
			double avgAct = totalAct / (((double) molecules.size()) - 1);
			for (GaMolecule molecule : molecules) {
				molecule.setWeight(molecule.getActivity() / avgAct);
			}
		}

		else {
			double totalAct = 0;
			int nMolecules = molecules.size();
			for (int i = 0; i < nMolecules; i++) {
				for (int j = i + 1; j < nMolecules; j++) {
					totalAct += (molecules.get(i).getActivity() + molecules.get(j)
							.getActivity()) / 2.0;
				}
			}
			int cnt = (nMolecules * (nMolecules - 1)) / 2;
			double avgAct = totalAct / cnt;
			for (GaMolecule molecule : molecules) {
				molecule.setWeight(molecule.getActivity() / avgAct);
			}
		}

		for (GaMolecule molecule : molecules) {
			infoMessageln("Molecule " + molecule.getName() + " weight "
					+ nf.format(molecule.getWeight()));
		}
		infoMessageln();
	}

	/**
	 * Registers a new user feature set and gives it a feature set id.
	 * 
	 * @params userFeatureSet UserFeatureSet to be registered
	 */
	@Override
	public FeatureType registerNextFeatureSetNo(UserFeatureSet userFeatureSet) {
		FeatureType featureType = Feature.userFeatureType(nUserFeatures);
		nUserFeatures++;
		userFeatures.put(featureType, userFeatureSet);
		assert (userFeatures.size() == nUserFeatures);
		nFeatures++;
		return featureType;
	}

	/**
	 * Returns total number of feature types defined
	 */
	@Override
	public int getnFeatureTypes() {
		return nFeatures;
	}

	/**
	 * Returns total number of user feature sets defined
	 */
	@Override
	public int getnUserFeatureTypes() {
		return nUserFeatures;
	}

	/**
	 * Returns a user feature set
	 * 
	 * @param no
	 *            User feature set number
	 */
	@Override
	public UserFeatureSet getUserFeatureSet(FeatureType featureType) {
		return userFeatures.get(featureType);
	}

	public List<? extends GaMolecule> getMolecules() {
		return molecules;
	}

	public GaMolecule getBaseMolecule() {
		return baseMolecule;
	}

	public GaMolecule getFittingMolecule() {
		return fittingMolecule;
	}

	public int getBaseMoleculeNo() {
		return baseMoleculeNo;
	}

	public int getFittingMoleculeNo() {
		return fittingMoleculeNo;
	}

	public int getBinaryStringLength() {
		return binaryStringLength;
	}

	public int getIntegerStringLength() {
		return integerStringLength;
	}

	public int[] getIntegerStringRanges() {
		return integerStringRanges;
	}

	public int[] getBinaryEntryPoints() {
		return binaryEntryPoints;
	}

	public int getNMolecules() {
		return molecules.size();
	}

	public int getNChromosomes() {
		return popSize * nIslands;
	}

	public int getNRuns() {
		return nRuns;
	}

	public int[] getIntegerEntryPoints() {
		return integerEntryPoints;
	}

	public double getPharmFactor() {
		return pharmFactor;
	}

	public boolean isCompareAll() {
		return compareAll;
	}

	public boolean isScaleFitting() {
		return scaleFitting;
	}

	public boolean isScaleFittingOn() {
		return scaleFittingOn;
	}

	public boolean isScalePharmacophore() {
		return scalePharmacophore;
	}

	public Molecule.FileType getFileFormat() {
		return fileFormat;
	}

	public double getNicheDistance() {
		return nicheDistance;
	}

	public double getGeometricWeightCriteria() {
		return geometricWeightCriteria;
	}

	public double getAcceptorAtomWt() {
		return acceptorAtomWt;
	}

	public double getConstraintWt() {
		return constraintWt;
	}

	public double getDonorHydrogenWt() {
		return donorHydrogenWt;
	}

	public int getNicheSize() {
		return nicheSize;
	}

	public boolean isNichesOn() {
		return nichesOn;
	}

	public double getAromaticRingWt() {
		return aromaticRingWt;
	}

	public boolean isSingleMatchOnly() {
		return singleMatchOnly;
	}

	public double getConfWt() {
		return confWt;
	}

	/**
	 * @param confWt
	 *            the confWt to set
	 */
	public void setConfWt(double confWt) {
		this.confWt = confWt;
	}

	public double getVolumeWt() {
		return volumeWt;
	}

	public double getFinishFittingRadius() {
		return finishFittingRadius;
	}

	public double getStartFittingRadius() {
		return startFittingRadius;
	}

	public FeatureOverlay getFeatureOverlay() {
		return featureOverlay;
	}

	public boolean isUseFeatureClustering() {
		return useFeatureClustering;
	}

	/**
	 * @return the useActivities
	 */
	public boolean isUseActivities() {
		return useActivities;
	}

	/**
	 * @return the allowNulls
	 */
	public boolean[] getAllowNulls() {
		return allowNulls;
	}

}
