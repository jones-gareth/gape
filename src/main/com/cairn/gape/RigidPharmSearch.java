package com.cairn.gape;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

import org.apache.log4j.Logger;

import com.cairn.common.utils.CommonUtils;
import com.cairn.common.utils.Coord;
import com.cairn.common.utils.Timer;
import com.cairn.gape.feature.Feature;
import com.cairn.gape.feature.Feature.FeatureType;
import com.cairn.gape.feature.HydrogenBondingType;
import com.cairn.gape.feature.LonePairAddition;
import com.cairn.gape.feature.Pharmacophore;
import com.cairn.gape.feature.UserFeature;
import com.cairn.gape.ga.BaseSupervisor;
import com.cairn.gape.molecule.GaMolecule;
import com.cairn.gape.utils.Clique;
import com.cairn.molecule.Atom;
import com.cairn.molecule.AtomType;
import com.cairn.molecule.Molecule;
import com.cairn.molecule.Molecule.MolReadException;

/**
 * 
 * RigidPharmSearch is a supervisor class ({@link BaseSupervisor}) that uses
 * clique detection to match rigid conformations against a GAPE pharmacophore.
 * The conformations can then be scored against the query molecule that contains
 * the pharmacophore using the GAPE scoring function.
 * 
 * @author Gareth Jones
 * @version 0.1
 */
public class RigidPharmSearch extends Superposition {

	private static Logger logger;
	static {
		logger = Logger.getLogger(RigidPharmSearch.class);
		// logger.setLevel(Level.DEBUG);
	}
	private volatile List<GaMolecule> queries;

	private volatile Pharmacophore queryPharmacophore;

	private volatile int correspondenceMatrixSize, targetFileNo, moleculeNo,
			checkPoint = -1, savePoint = -1, skipMolecules = -1, nSavedHits = 0,
			nCurrentHits = 0, reportMols = 1000, maxSavedHits = 200,
			nConformerMatches = 1;

	private volatile BufferedReader targetReader;

	private volatile List<Feature> queryFeatures, targetFeatures;

	private volatile GaMolecule target, previousTarget;

	// tolerance is a ratio of difference to distance
	private volatile double tolerance = 0.2, maxScore = -Double.MAX_VALUE,
			minScore = Double.MAX_VALUE, proximityDistance = 5.0, targetScore,
			featureScore, volumeIntegral, queryDistanceMatrix[][],
			targetDistanceMatrix[][];

	private volatile Molecule.FileType targetFileType = Molecule.FileType.MOL2;

	private volatile boolean correspondenceMatrix[][], openEyeConformers, skipping,
			noPharmacophore, useProximity = true, textSummary;

	private final boolean baseMoleculeOnly = false;

	private volatile boolean conformersPrepared = true;

	private final Clique clique = new Clique();

	private volatile String queryFile;

	private final String targetFiles[];

	private volatile Writer molErrorWriter, textSummaryWriter;

	private final ArrayList<ConformerSummarySolution> conformerSummarySolutions = new ArrayList<ConformerSummarySolution>();

	// The summary file can contain all hits or just those in the matches
	// structure file
	private static final boolean SUMMARIZE_ALL_HITS = false;

	/**
	 * Constructor.
	 * 
	 * @param queryFile
	 *            Structure file containing query compounds. Must contain GAPE
	 *            pharmacophore definition.
	 * @param targetFiles
	 *            Structure files of targets
	 */
	public RigidPharmSearch(String configFile, String queryFile, String targetFiles[]) {
		super();

		cleanRigidPharmSearch();

		System.out.print("\n\nGAPE(GRIPS)\n");
		System.out.print("Gareth Jones\n\n");

		setupInfoOutput();

		this.queryFile = queryFile;
		this.targetFiles = targetFiles;

		init(configFile);

		long mem = Runtime.getRuntime().maxMemory();
		if (mem == Long.MAX_VALUE)
			infoMessageln(0, "Maximum memory Available");
		else {
			double meg = mem / (1024 * 1024);
			infoMessageln(0, meg + " MB memory avaliable");
		}

		// load in queries
		if (queryFile.toLowerCase().endsWith(".pharm")) {

			// Pharmacophore definition file
			queryPharmacophore = new Pharmacophore(queryFile);

		} else {

			// Structural query/queries

			queries = GaMolecule.loadFiles(new String[] { queryFile },
					Molecule.FileType.UNKNOWN, Molecule.Source.FILE, infoMessageLogger,
					true);
			queries = queries.stream()
					.filter(m -> !m.getName().contains("Pharmacophore"))
					.collect(Collectors.toList());
		}

		// Open mol error file
		try {
			molErrorWriter = new FileWriter("failed.mol2");
		} catch (IOException ex) {
			infoMessageln("failed to open molecule error file");
			System.exit(0);
		}

	}

	/**
	 * cleans the current directory of GRIPS output files. Removes the most
	 * comon output files, but depending on options may not remove all files.
	 */
	private void cleanRigidPharmSearch() {
		String files[] = (new File(".")).list();
		for (int i = 0; i < files.length; i++) {
			String name = files[i];
			if ((name.startsWith("check_") && name.endsWith(".mol2"))
					|| (name.startsWith("check_") && name.endsWith(".sdf"))
					|| (name.startsWith("Matches") && name.endsWith(".mol2"))
					|| (name.startsWith("Matches") && name.endsWith(".sdf"))
					|| (name.startsWith("GRIPS") && name.endsWith(".mol2"))
					|| (name.startsWith("GRIPS") && name.endsWith(".sdf"))
					|| (name.startsWith("failed") && name.endsWith(".sdf"))
					|| name.equals(INFO_FILE) || name.equals(SEED_FILE)
					|| name.equals("grips_summary.txt")) {
				File f = new File(name);
				if (f.exists())
					f.delete();
			}
		}
	}

	/**
	 * Runs the program
	 * 
	 * @param args
	 *            . Program arguments, at least three: configuration file, query
	 *            structure file and target structure file(s).
	 */
	public static void main(String args[]) {

		if (args.length < 3) {
			System.err
					.println("Usage RigidPharmSearch <configFile> <query> <target_files..>");
			System.err.println();
			System.err.println("      configFile:   parameter file");
			System.err.println("      query:        structure file containing "
					+ "GAPE overlay or single molecule");
			System.err.println("                    or pharmacophore (.pharm) file ");
			System.err
					.println("      target_files: one or more files of target conformations");
			System.exit(0);
		}

		String configFile = args[0];
		String queryFile = args[1];
		String targetFiles[] = new String[args.length - 2];
		for (int i = 0; i < args.length - 2; i++)
			targetFiles[i] = args[i + 2];

		try {
			RigidPharmSearch search = new RigidPharmSearch(configFile, queryFile,
					targetFiles);
			search.searchMolecules();
		} catch (Exception ex) {
			ex.printStackTrace(System.err);
			System.err.println("RigidPharmSearch: Exception: " + ex);
		}
	}

	/**
	 * Reads in run-time setting from the configuration file. Prints out
	 * compiled settings and runtime settings.
	 * 
	 * @param configFile
	 *            Configuration filename.
	 */
	@Override
	public void init(String configFile) {
		super.init(configFile);

		if (hasKey("rigid_pharm_tolerance"))
			tolerance = getDoubleValue("rigid_pharm_tolerance");
		if (hasKey("max_saved_hits"))
			maxSavedHits = getIntValue("max_saved_hits");
		if (hasKey("report_mols"))
			reportMols = getIntValue("report_mols");
		if (hasKey("proximity_distance"))
			proximityDistance = getDoubleValue("proximity_distance");
		if (hasKey("use_proximity")) {
			useProximity = getBooleanValue("use_proximity");
		}
		if (hasKey("openeye_conformers")) {
			openEyeConformers = getBooleanValue("openeye_conformers");
		}
		if (hasKey("conformers_prepared")) {
			conformersPrepared = getBooleanValue("conformers_prepared");
		}
		if (hasKey("n_conformer_matches")) {
			nConformerMatches = getIntValue("n_conformer_matches");
			if (!openEyeConformers)
				nConformerMatches = -1;
		}
		if (hasKey("no_pharmacophore"))
			noPharmacophore = getBooleanValue("no_pharmacophore");
		if (hasKey("check_point"))
			checkPoint = getIntValue("check_point");
		if (hasKey("save_point"))
			savePoint = getIntValue("save_point");
		if (hasKey("skip_molecules"))
			skipMolecules = getIntValue("skip_molecules");
		if (hasKey("text_summary"))
			textSummary = getBooleanValue("text_summary");
		// base molecule only does not make sense for grips 2.0
		// if (hasKey("base_molecule_only"))
		// baseMoleculeOnly = getBooleanValue("base_molecule_only");

		infoMessageln();
		infoMessageln("Saving up to " + maxSavedHits + " hits");
		infoMessageln("Clique detection tolerance is " + tolerance);
		if (useProximity)
			infoMessageln("Using proximity distance of " + proximityDistance);
		if (openEyeConformers)
			infoMessageln("Using openeye conformer names");
		if (conformersPrepared)
			infoMessageln("Conformers already prepared");
		if (nConformerMatches > 0)
			infoMessageln("Maximum number of matches per structure is "
					+ nConformerMatches);
		if (noPharmacophore)
			infoMessageln("Ignoring any pharmacophore in input structures");
		if (baseMoleculeOnly)
			infoMessageln("Matching only on base molecule features");

		if (textSummary) {
			try {
				infoMessageln("Writing summary to grips_summary.txt");
				textSummaryWriter = new BufferedWriter(
						new FileWriter("grips_summary.txt"));
				textSummaryWriter.write("STRUCTURE\tRIGID_PHARMACOPHORE_SCORE\t"
						+ "FEATURE_SCORE\tVOLUME_SCORE\n");
				textSummaryWriter.flush();
			} catch (IOException ex) {
				throw new RuntimeException("can't open grips_summary.txt");
			}
		}

		queryFile = fileName(queryFile);
		for (int i = 0; i < targetFiles.length; i++)
			targetFiles[i] = fileName(targetFiles[i]);
	}

	/**
	 * Searches a set of rigid compounds against a GAPE solution or base
	 * molecule.
	 * 
	 */
	private void searchMolecules() {

		totalTimer = new Timer();
		runTimer = new Timer();

		Feature.setRadius(finishFittingRadius);
		setupQuery();

		if (conformersPrepared) {
			Molecule.setAddHydrogensFlag(false);
			Molecule.setSolvateFlag(false);
		}

		while (true) {

			if (!getNextTarget()) {
				break;
			}

			try {
				matchMolecule();
				if (!skipping)
					previousTarget = target;
			} catch (Exception ex) {
				String err = "RigidPharmSearch: failed to process " + target.getName();
				System.err.println(err);
				infoMessageln(0, err);
				infoMessageln(0, "Molecule written to failed.mol2");
				String stack = CommonUtils.getStackTrace(ex);
				infoMessageln(0, stack);
				writeFailedTarget();
			}

			// if (moleculeNo == 1000)
			// break;
		}

		infoMessageln("Searched " + moleculeNo + " conformations.");

		if (best != null) {
			infoMessageln(0, "Found " + nCurrentHits + " hits");
			infoMessageln(0, "Best score " + nf.format(maxScore));
			infoMessageln(0, "Worst score " + nf.format(minScore));
			infoMessageln(0, "Writing " + nSavedHits + " hits");

			writeSolutions("Matches");
		} else
			infoMessageln(0, "No Matches found!");

		try {
			molErrorWriter.close();
		} catch (IOException ex) {
			System.err.println("Error closing molecule error file");
		}

		if (SUMMARIZE_ALL_HITS && textSummary) {
			try {
				// write out any summary data
				if (nConformerMatches > 0)
					writeSummarySolutions();
				textSummaryWriter.close();
			} catch (IOException ex) {
				logger.error("Error closing text_summary.txt");
			}
		}

		totalTimer.interval();
		infoMessageln(0, "Run Time " + totalTimer.info());
	}

	/**
	 * Writes out the current set of matches
	 * 
	 */
	private void writeSolutions(String prefix) {
		String suffix = targetFileType == Molecule.FileType.MOL2 ? ".mol2" : ".sdf";
		String file = prefix + suffix;

		if (best == null) {
			infoMessageln(0, "No solutions to write to " + file);
			return;
		} else {
			infoMessageln(0, "Writing current solutions to " + file);
		}

		try {
			FileWriter out = new FileWriter(file);

			// Replace by non-recursive function
			// best.writeSolutions(out);
			for (Solution curr = best; curr != null; curr = curr.next) {
				curr.writeSolution(out);
				if (!SUMMARIZE_ALL_HITS && textSummary) {
					String summary = curr.solutionMolecule.getName() + "\t"
							+ nf.format(curr.targetScore) + "\t"
							+ nf.format(curr.featureScore) + "\t"
							+ nf.format(curr.volumeIntegral) + "\n";
					logger.debug(summary);
					textSummaryWriter.write(summary);

				}
			}
			out.close();
			if (!SUMMARIZE_ALL_HITS && textSummary)
				textSummaryWriter.close();
		} catch (IOException ex) {
			System.err
					.println("PharmSearch: IOException writing to " + file + " : " + ex);
		}

	}

	/**
	 * Writes out a target to the file of failed molecules
	 * 
	 */
	private void writeFailedTarget() {
		try {
			target.writeSybylMol2(molErrorWriter, null);
			molErrorWriter.flush();
		} catch (IOException x) {
			System.err.println("Error writing to molecule error file");
		}
	}

	/**
	 * Reads in the next target from the target files. Returns false when there
	 * are no more targets to be read.
	 */
	private boolean getNextTarget() {
		target = new GaMolecule();
		target.setInfoMessageLogger(infoMessageLogger);

		if (targetReader == null) {

			// No more files
			if (targetFileNo == targetFiles.length)
				return false;

			String targetFile = targetFiles[targetFileNo++];
			targetFileType = Molecule.getType(targetFile);
			infoMessageln(0, "Reading from " + targetFile);

			// open reader
			try {
				if (targetFile.toUpperCase().endsWith(".GZ")) {
					GZIPInputStream zipStream = new GZIPInputStream(new FileInputStream(
							targetFile));
					targetReader = new BufferedReader(new InputStreamReader(zipStream));
				} else {
					targetReader = new BufferedReader(new FileReader(targetFile));
				}
			} catch (FileNotFoundException ex) {
				throw new RuntimeException("File " + targetFile + " not found");
			} catch (IOException ex) {
				throw new RuntimeException("IO Exception " + ex);
			}
		}

		try {

			// skip molecules
			skipping = false;
			if (moleculeNo < skipMolecules) {
				skipping = true;
				if (targetFileType == Molecule.FileType.MOL2) {
					int molCnt = 0;
					while (true) {
						targetReader.mark(4096);
						String val = targetReader.readLine();
						if (val.startsWith("@<TRIPOS>MOLECULE")) {
							molCnt++;
							if (molCnt == 2)
								targetReader.reset();
							break;
						}
					}
				} else {
					while (true) {
						String val = targetReader.readLine();
						if (val.startsWith("$$$$"))
							break;
					}
				}
			}

			// actually read in a molecule!!
			else if (targetFileType == Molecule.FileType.MOL2)
				target.loadSybylMol2(targetReader);
			else {
				target.loadSdfMol(targetReader);
				// check target contains hydrogens
				if (!Molecule.getAddHydrogensFlag() && target.countHydrogens() == 0) {
					String err = "SD Molecule " + target.getName()
							+ " has no hydrogens\n"
							+ "GRIPS requires that SD files contain all hydrogens";
					infoMessageln(err);
					System.err.println(err);
					// Try to read next file
					return getNextTarget();

				}
			}

		} catch (IOException | MolReadException ex) {
			try {
				targetReader.close();
			} catch (IOException e) {
				throw new RuntimeException("can't close reader");
			}
			targetReader = null;
			// Try to read next file
			return getNextTarget();
		}

		return true;
	}

	/**
	 * Sees if the current target can be matched against the pharmacophore.
	 * Returns false if there are no more molecules to be read.
	 */
	private void matchMolecule() {
		int logLevel = infoMessageLogger.getLogLevel();

		moleculeNo++;
		if (skipping) {
			if (moleculeNo % reportMols == 0)
				infoMessageln(0, "Skipped " + moleculeNo + " confs");
			return;
		}

		if (target.getName().endsWith("Pharmacophore"))
			return;

		infoMessageln(2, "Searching Molecule " + target.getName());

		if (moleculeNo % reportMols == 0) {
			runTimer.interval();
			totalTimer.interval();
			String info = "Searched " + moleculeNo + " confs " + nCurrentHits
					+ " hits [interval " + runTimer.total() + " total "
					+ totalTimer.total() + "]";
			infoMessageln(1, info);
			runTimer.reset();
			if (textSummary) {
				try {
					textSummaryWriter.flush();
				} catch (IOException e) {
					throw new RuntimeException("can't flush grips_summary.txt");
				}
			}
		}

		if (checkPoint > 0 && moleculeNo % checkPoint == 0) {
			String file = "check_point_" + moleculeNo;
			writeSolutions(file);
		}

		if (savePoint > 0 && moleculeNo % savePoint == 0) {
			String file = "Matches_" + moleculeNo;
			writeSolutions(file);
			infoMessageln(0, "Deleting saved hits");
			best = worst = null;
		}

		// try to reuse structure if openEyeConformers is set
		if (openEyeConformers && sameAsPrevious() && !Molecule.getAddHydrogensFlag()
				&& !Molecule.getSolvateFlag()) {
			copyPrevious();
		}

		else {

			logger.debug("Initializing Target");
			target.init();

			logger.debug("Setting Up Target");
			setupMolecule(target);
			targetFeatures = target.getAllMappingFeatures();
		}

		createTargetDistanceMatrix();
		createCorrespondenceMatrix();

		clique.setEdges(correspondenceMatrix, correspondenceMatrixSize);
		clique.setSetlim(3);
		clique.setVerbose(false);
		if (clique.solve() < 3) {
			infoMessageln(2, "No Match");
			return;
		}

		if (logLevel > 2)
			infoMessageln(2, clique.bestInfo());
		findMatch();

		if (!fitTarget()) {
			infoMessageln(2, "Can't fit target");
			return;
		}

		if (useProximity && !checkProximity())
			return;

		targetScore = scoreFeatureMatch() + scoreVolumeMatch();
		if (Double.isNaN(targetScore))
			throw new RuntimeException("target score is NaN");
		if (logLevel > 1)
			infoMessageln("Molecule Score " + nf.format(targetScore));
		if (logger.isDebugEnabled()) {
			logger.debug(target.getName() + " Score " + nf.format(targetScore));
			logger.debug(target.getName() + " Feature Score " + nf.format(featureScore));
			logger.debug(target.getName() + " Volume Score " + nf.format(volumeIntegral));
		}

		if (targetScore > maxScore)
			maxScore = targetScore;
		if (targetScore < minScore)
			minScore = targetScore;
		nCurrentHits++;

		addSolution();

		if (SUMMARIZE_ALL_HITS && textSummary) {
			// write out summary data- note that we restrict the number of
			// conformers in the summary.
			ConformerSummarySolution conformerSummarySolution = new ConformerSummarySolution(
					target.getName(), targetScore, featureScore, volumeIntegral);

			if (nConformerMatches <= 0)
				conformerSummarySolution.print();

			else {
				// if we've moved onto a new structure output the best n
				// conformers
				if (conformerSummarySolutions.size() > 0
						&& !sameStructure(target.getName(),
								conformerSummarySolutions.get(0).name))
					writeSummarySolutions();
				conformerSummarySolutions.add(conformerSummarySolution);
			}
		}

		// target.writeSybylMol2File
		// ("match_"+moleculeNo+".mol2", null);
	}

	/**
	 * Writes summary solutions. We need to restrict the number of solutions for
	 * a structure to the number of conformer matches.
	 * 
	 */
	private void writeSummarySolutions() {
		if (conformerSummarySolutions.size() == 0)
			return;
		if (conformerSummarySolutions.size() == 1) {
			conformerSummarySolutions.get(0).print();
			conformerSummarySolutions.clear();
			return;
		}
		int no = conformerSummarySolutions.size() > nConformerMatches ? nConformerMatches
				: conformerSummarySolutions.size();
		List<ConformerSummarySolution> sortedSolutions = new ArrayList<>(
				conformerSummarySolutions.subList(0, no));
		Collections.sort(sortedSolutions);
		for (ConformerSummarySolution solution : sortedSolutions) {
			solution.print();
		}
		conformerSummarySolutions.clear();
	}

	/**
	 * Sets up the query. Calls {@link setupMolecule} for all queries. Finds the
	 * pharmacophore definition.
	 * 
	 */
	private void setupQuery() {

		commonMolSetup();

		if (queryPharmacophore != null) {

			// Get query features from pharmacophore
			queryFeatures = queryPharmacophore.getAllFeatures();

		} else {

			// Get query from structure files.
			// molecules = new GaMolecule[queries.length + 1];
			// nMolecules = queries.length + 1;
			// binaryEntryPoints = new int[nMolecules];
			//
			GaMolecule baseMolecule = null;

			for (GaMolecule query : queries) {
				setupMolecule(query);
				query.getAtomicGaussians();

				if (query.getName().endsWith("Base Molecule")) {
					baseMolecule = query;
				}
			}

			if (noPharmacophore && queries.size() > 1)
				throw new RuntimeException(
						"No pharmacophore set and multiple query structures\n");

			// if noPharmacophore is set we match against the first query
			// structure and use all features in it (also do this if we only
			// have one structure and no pharmacophore).

			if (noPharmacophore
					|| (queries.size() == 1 && queries.get(0).getNPharmDescriptions() == 0)) {
				infoMessageln(0, "No pharmacophore set: using all " + "features from "
						+ queries.get(0).getName());
				defaultPharmacophore();

			} else {

				// Otherwise use pharmacophore features defined in the input
				// structures

				// Pharmacophore pharmacophore = new Pharmacophore();
				// pharmacophore.createFromMoleculePharmDescriptions(queries);

				if (baseMoleculeOnly) {
					// just use base molecule features
					if (baseMolecule == null)
						throw new RuntimeException("Unable to find base molecule");
					baseMolecule.findPharmFeatures();
					queryFeatures = baseMolecule.getPharmFeatures();
					infoMessageln(0, "Using base molecule pharmacophore in "
							+ baseMolecule.getName());
				}

				else {
					// use all features

					queryFeatures = queries.stream().flatMap(q -> {
						q.findPharmFeatures();
						return q.getPharmFeatures().stream();
					}).collect(Collectors.toList());

					infoMessageln(0,
							"Using pharmacophore features found in all molecules");
				}

			}
		}

		if (queryFeatures.size() < 3)
			throw new RuntimeException("Pharmacophore should have at least three points!");

		if (infoMessageLogger.getLogLevel() > 0) {
			infoMessageln("Pharmacophore");
			for (Feature queryFeature : queryFeatures) {
				infoMessageln(queryFeature.info());
			}
		}

		createQueryDistanceMatrix();

		// Plenty for space as multiple target features may map to a
		// single pharmacophore point
		maxMatches = queryFeatures.size() * 5;
		matches = new Pair[maxMatches];
		for (int i = 0; i < maxMatches; i++)
			matches[i] = new Pair();

		xPoints = new double[maxMatches][];
		yPoints = new double[maxMatches][];
	}

	/**
	 * Default Pharmacophore for a query. Uses all features except hydrophobic
	 * atoms. Valid only for a single molecule.
	 * 
	 */
	private void defaultPharmacophore() {
		GaMolecule query = queries.get(0);
		if (queries.size() > 1)
			throw new RuntimeException(
					"Default pharmacophore valid for single query only!");

		queryFeatures = new ArrayList<Feature>();
		for (Feature f : query.getAllFeatures()) {
			if (f.getFeatureType() != FeatureType.HYDROPHOBIC_ATOM)
				queryFeatures.add(f);
		}

		// add pharmacophore descriptions to query
		for (Feature queryFeature : queryFeatures)
			query.addPharmDescription(queryFeature.featureLabel());

		// save query with pharmacophore info
		Molecule.FileType queryFileType = Molecule.getType(queryFile);
		String outFile = "GRIPS_query.sdf";
		if (queryFileType == Molecule.FileType.MOL2)
			outFile = "GRIPS_query.mol2";
		query.write(outFile,
				"Grips query with default pharmacophore for a single molecule");

	}

	/**
	 * Sets up a molecule. Stripped down version of {@link GaMolecule#setup}.
	 * 
	 * @param molecule
	 *            Structure to initialize.
	 */
	private void setupMolecule(GaMolecule molecule) {

		molecule.setRigid(true);
		molecule.setFixed(true);
		molecule.setProblem(this);

		infoMessageln(3, "\nFinding Dean and Mills Donors and Acceptors");
		HydrogenBondingType.searchMolecule(molecule, infoMessageLogger);

		LonePairAddition.addLonePairs(molecule);
		infoMessageln(2, "Finding Features");
		molecule.findFeatures();

		setupFeatureClustering();

		for (double[] coord : molecule.getCoords()) {
			double[] ref = new double[4];
			Coord.copy(coord, ref);
			molecule.getReference().add(ref);
		}

		for (int i = 0; i < molecule.getNFeatures(); i++)
			molecule.getAllFeature(i).calculateCoordinate();

	}

	/**
	 * Creates an inter-feature distance matrix for the query.
	 */
	private void createQueryDistanceMatrix() {
		int nQueryFeatures = queryFeatures.size();
		queryDistanceMatrix = new double[nQueryFeatures][nQueryFeatures];
		for (int i = 0; i < nQueryFeatures; i++) {
			Feature f1 = queryFeatures.get(i);
			for (int j = i + 1; j < nQueryFeatures; j++) {
				Feature f2 = queryFeatures.get(j);
				double distance = Coord.distance(f1.getSavedCoordinate(),
						f2.getSavedCoordinate());
				queryDistanceMatrix[i][j] = queryDistanceMatrix[j][i] = distance;
			}
		}
	}

	private int _targetDistanceMatrixSize;

	/**
	 * Creates an inter-feature distance maxtix for the target.
	 */
	private void createTargetDistanceMatrix() {
		int nTargetFeatures = targetFeatures.size();
		if (nTargetFeatures > _targetDistanceMatrixSize) {
			targetDistanceMatrix = new double[nTargetFeatures][nTargetFeatures];
			_targetDistanceMatrixSize = nTargetFeatures;
		}
		for (int i = 0; i < nTargetFeatures; i++) {
			Feature f1 = targetFeatures.get(i);
			for (int j = i + 1; j < nTargetFeatures; j++) {
				Feature f2 = targetFeatures.get(j);
				double distance = Coord.distance(f1.getSavedCoordinate(),
						f2.getSavedCoordinate());
				targetDistanceMatrix[i][j] = targetDistanceMatrix[j][i] = distance;
			}
		}
	}

	/**
	 * Class for representing a correspondance node between query and target.
	 */
	private class Pair {
		int queryNo;

		int targetNo;

		Feature queryFeature;

		Feature targetFeature;

		/**
		 * Given an index into the correspondence matrix find the target and
		 * query feature numbers.
		 * 
		 * @param index
		 *            Row or column index into correspondence matrix.
		 * @param pair
		 *            Class for holding feature numbers.
		 */
		void getCorrespondanceIndex(int index) {
			int nTargetFeatures = targetFeatures.size();
			queryNo = index / nTargetFeatures;
			targetNo = index % nTargetFeatures;
		}

		/**
		 * Checks that the corresponding features are of the same type.
		 */
		boolean sameType() {
			if (queryFeatures.get(queryNo).getFeatureType() == targetFeatures.get(
					targetNo).getFeatureType())
				return true;
			return false;
		}
	}

	private final Pair column = new Pair(), row = new Pair();

	private int _allocatedCorrepondenceMatrixSize;

	/**
	 * Forms the correspondence matrix
	 */
	private void createCorrespondenceMatrix() {
		correspondenceMatrixSize = queryFeatures.size() * targetFeatures.size();
		if (correspondenceMatrixSize > _allocatedCorrepondenceMatrixSize) {
			correspondenceMatrix = new boolean[correspondenceMatrixSize][correspondenceMatrixSize];
			_allocatedCorrepondenceMatrixSize = correspondenceMatrixSize;
		}

		for (int i = 0; i < correspondenceMatrixSize; i++) {
			row.getCorrespondanceIndex(i);
			for (int j = i + 1; j < correspondenceMatrixSize; j++) {
				column.getCorrespondanceIndex(j);

				if (!row.sameType() || !column.sameType()) {
					correspondenceMatrix[i][j] = correspondenceMatrix[j][i] = false;
					continue;
				}

				double queryDistance = queryDistanceMatrix[row.queryNo][column.queryNo];
				double targetDistance = targetDistanceMatrix[row.targetNo][column.targetNo];
				if (Math.abs(queryDistance - targetDistance) / queryDistance < tolerance)
					correspondenceMatrix[i][j] = correspondenceMatrix[j][i] = true;
				else
					correspondenceMatrix[i][j] = correspondenceMatrix[j][i] = false;
			}
		}
	}

	private Pair matches[];

	private int nMatches, maxMatches;

	private double xPoints[][], yPoints[][];

	/**
	 * Find the feature match associated with the best clique.
	 */
	private void findMatch() {
		nMatches = clique.getBestsize();
		int logLevel = infoMessageLogger.getLogLevel();
		for (int i = 0; i < nMatches; i++) {
			matches[i].getCorrespondanceIndex(clique.getBestset(i));
			matches[i].queryFeature = queryFeatures.get(matches[i].queryNo);
			matches[i].targetFeature = targetFeatures.get(matches[i].targetNo);

			if (logLevel > 2)
				infoMessageln(matches[i].queryFeature.info() + " --> "
						+ matches[i].targetFeature.info());
		}
	}

	private final double trans[][] = new double[4][4];

	/**
	 * Fits the target onto the matched clique
	 */
	private boolean fitTarget() {
		for (int i = 0; i < nMatches; i++) {
			xPoints[i] = matches[i].targetFeature.getSavedCoordinate();
			yPoints[i] = matches[i].queryFeature.getSavedCoordinate();
		}
		for (int i = nMatches; i < maxMatches; i++)
			xPoints[i] = yPoints[i] = null;

		Coord.leastSquaresFit(xPoints, yPoints, trans);

		int nAtoms = target.getnAtoms();
		for (int i = 0; i < nAtoms; i++) {
			if (target.getCoord(i)[3] > 1.0001 || target.getCoord(i)[3] < 0.99999)
				throw new RuntimeException("target coordinate error");
			Coord.transPointInPlace(trans, target.getCoord(i));
		}
		return true;
	}

	/**
	 * Returns the feature score from the clique match.
	 */
	private double scoreFeatureMatch() {
		featureScore = .0;

		int logLevel = infoMessageLogger.getLogLevel();
		boolean debugLog = logger.isDebugEnabled();
		List<String> pharmDescriptions = new ArrayList<String>();

		for (int i = 0; i < nMatches; i++) {
			Feature queryFeature = matches[i].queryFeature;
			Feature targetFeature = matches[i].targetFeature;
			targetFeature.calculateCoordinate();
			double score = queryFeature.score(targetFeature);
			if (Double.isNaN(score)) {
				String queryLabel = queryFeature.atomLabel();
				String targetLabel = targetFeature.atomLabel();
				throw new RuntimeException("scoreFeatureMatch: NAN error: query: "
						+ queryLabel + " target: " + targetLabel);
			}
			if (logLevel > 2)
				infoMessageln(queryFeature.info() + " -> " + targetFeature.info()
						+ " score " + nf.format(score));
			if (debugLog)
				logger.debug(queryFeature.info() + " -> " + targetFeature.info()
						+ " score " + nf.format(score));
			double weight = .0;
			FeatureType featureType = queryFeature.getFeatureType();
			switch (featureType) {
			case HYDROPHOBIC_ATOM:
				break;
			case ACCEPTOR_ATOM:
				weight = acceptorAtomWt;
				break;
			case DONOR_INTERACTION_POINT:
				weight = donorHydrogenWt;
				break;
			case AROMATIC_RING:
				weight = aromaticRingWt;
				break;
			default:
				UserFeature userFeature = (UserFeature) queryFeature;
				weight = userFeature.getUserFeatureSet().getFeatureWeight();
				break;
			}

			featureScore += weight * score;

			// add pharmacophore info
			if (queryFeature.getGeometricScore() > geometricWeightCriteria)
				pharmDescriptions.add(targetFeature.featureLabel());
		}

		target.setPharmDescriptions(pharmDescriptions);

		return featureScore;
	}

	/**
	 * Returns the volume score for overlay with query molecule(s)
	 */
	private double scoreVolumeMatch() {
		int nQueries = queries.size();
		if (nQueries == 0 || volumeWt == 0)
			return 0;
		// Only get the atomic gaussians if we need them
		target.getAtomicGaussians();
		volumeIntegral = .0;
		for (int i = 0; i < nQueries; i++)
			volumeIntegral += target.gaussianIntegral(queries.get(i));
		volumeIntegral = volumeIntegral / nQueries;
		if (infoMessageLogger.getLogLevel() > 2)
			infoMessageln("volume integral " + nf.format(volumeIntegral));
		volumeIntegral = volumeIntegral * volumeWt;
		return volumeIntegral;
	}

	private Solution best, worst;

	/**
	 * Stores the current target into a linked list of solutions. The linked
	 * list retains only the top hits.
	 */
	private void addSolution() {
		logger.debug("called addSolution score " + targetScore);

		if (nSavedHits >= maxSavedHits && targetScore <= worst.targetScore)
			return;
		Solution s = null;

		if (openEyeConformers) {
			// Here we need to copy the target as we may reuse the
			// structure.
			GaMolecule sol = new GaMolecule();
			sol.setName(target.getName());
			sol.setPharmDescriptions(target.getPharmDescriptions());
			target.copyMinimum(sol);
			s = new Solution(sol, targetScore, featureScore, volumeIntegral);
		} else {
			s = new Solution(target, targetScore, featureScore, volumeIntegral);
		}

		// only store nConformerMatches conformer hits of the same molecule
		if (nConformerMatches >= 1) {
			int nMatches = 0;
			Solution other = null;
			for (Solution curr = best; curr != null; curr = curr.next) {
				// look for solutions with the same conformer
				if (sameStructure(curr.solutionMolecule.getName(),
						s.solutionMolecule.getName())) {
					other = curr;
					nMatches++;
					if (nMatches == nConformerMatches)
						break;
				}
			}

			if (nMatches == nConformerMatches) {
				// If the number of matching conformers is the same as
				// nConformerMatches and this match is better than the current
				// worst, remove the current worst.
				if (s.targetScore > other.targetScore) {
					if (best == other)
						best = other.next;
					if (worst == other)
						worst = other.prev;
					other.remove();
					nSavedHits--;
				} else {
					// this solution is worse than the one(s) already present so
					// don't add it
					return;
				}
			}
		}

		if (best == null && worst == null) {
			best = worst = s;
		} else {
			// Replace by non recursive function
			// best.addSolution(s);
			boolean added = false;
			for (Solution curr = best; curr != null; curr = curr.next) {
				if (s.targetScore > curr.targetScore) {
					curr.addBefore(s);
					added = true;
					break;
				}
			}
			if (!added)
				worst.addAfter(s);
		}
		if (targetScore > best.targetScore)
			best = s;
		if (targetScore <= worst.targetScore)
			worst = s;
		if (nSavedHits >= maxSavedHits) {
			if (worst.prev == null)
				System.err.println("Linked list error");
			worst = worst.prev;
			worst.next = null;
		} else
			nSavedHits++;
		if (logger.isDebugEnabled())
			logger.debug("added solution " + target.getName() + " score " + targetScore
					+ " list length " + solutionListLength());
	}

	/**
	 * Checks that all atoms in that target are within proximityDistance of a
	 * query atom.
	 */
	private boolean checkProximity() {
		double tol = proximityDistance * proximityDistance;

		// Skip check if there is no query structure
		if (queries.size() == 0 && queryPharmacophore != null)
			return true;

		for (Atom targetAtom : target.getAtoms()) {
			if (targetAtom.getAtomType() == AtomType.Type.LP)
				continue;
			boolean ok = false;

			search: for (GaMolecule query : queries) {

				for (Atom queryAtom : query.getAtoms()) {
					if (queryAtom.getAtomType() == AtomType.Type.LP)
						continue;
					if (Coord.sqrDistance(query.getCoord(queryAtom.getNo()),
							target.getCoord(targetAtom.getNo())) < tol) {
						ok = true;
						break search;
					}

				}
			}

			if (!ok)
				return false;
		}
		return true;
	}

	/**
	 * Uses the molecule name to see if the current conformer is the same as the
	 * previous one. Conformer names are like conformer_1, conformer_2 or
	 * conformer_1_rconf, conformer_2_rconf etc
	 */
	private boolean sameAsPrevious() {
		if (previousTarget == null)
			return false;

		String prev = previousTarget.getName();
		String name = target.getName();

		return sameStructure(name, prev);
	}

	// These strings are patterns for conformer matching
	// This patern considers "chiral" enumerations to be the same structure-
	// took it out
	// "_\\d+_\\d+\\z"
	static private String openEyeStrings[] = new String[] { "_\\d+_rconf_\\d+\\z",
			"_\\d+_rconf\\z", "_\\d+\\z" };

	static private Pattern openEyePatterns[];
	static {
		openEyePatterns = new Pattern[openEyeStrings.length];
		for (int i = 0; i < openEyePatterns.length; i++) {
			openEyePatterns[i] = Pattern.compile(openEyeStrings[i]);
		}
	}

	/**
	 * @param name
	 * @return The openye baser name for the conformer.
	 * 
	 */
	static private String findBaseName(String name) {
		for (int i = 0; i < openEyePatterns.length; i++) {
			Matcher m = openEyePatterns[i].matcher(name);
			if (m.find())
				return m.replaceAll("");
		}
		return name;
	}

	/**
	 * Compares molecule names to determine if they are two conformers of the
	 * same compound. Uses Openeye conformer naming system (conformer_1,
	 * conformer_2 or conformer_1_rconf, conformer_2_rconf, conformer_1_1,
	 * conformer_1_2 etc)
	 * 
	 * @param name
	 * @param otherName
	 * @return
	 */
	private boolean sameStructure(String name, String otherName) {

		name = findBaseName(name);
		otherName = findBaseName(otherName);
		if (logger.isDebugEnabled())
			logger.debug("Previous " + otherName + " target " + name);

		if (otherName.equals(name)) {
			if (logger.isDebugEnabled())
				logger.debug(name + " is the same molecule as " + otherName);
			return true;
		}
		return false;
	}

	/**
	 * Copies coordinates from the current molecule to the previous one and sets
	 * the previous target as the current target. This won't work if we do
	 * things that add or remove atoms such as adding hydrogens or solvate
	 * atoms.
	 */
	private void copyPrevious() {
		int nAtoms = target.getnAtoms();
		for (int i = 0; i < nAtoms; i++)
			previousTarget.getCoords().set(i, target.getCoord(i));
		previousTarget.setName(target.getName());

		target = previousTarget;

		// LonePairAddition.updateLonePairs(target);
		target.removeLonePairs();
		LonePairAddition.addLonePairs(target);

		for (int i = 0; i < nAtoms; i++) {
			Coord.copy(target.getCoord(i), target.getReference().get(i));
		}
		for (int i = 0; i < target.getNFeatures(); i++) {
			target.getAllFeature(i).calculateCoordinate();
		}
	}

	/**
	 * Return (non-recursively) the length of the solution list.
	 * 
	 */
	private int solutionListLength() {
		int len = 0;
		for (Solution curr = best; curr != null; curr = curr.next)
			len++;
		return len;
	}

	/**
	 * Linked list class to store solutions. Ordered by GAPE score. Had to
	 * rewrite use non-recursive functions to handle large hitlists.
	 */
	private class Solution {
		private final GaMolecule solutionMolecule;

		private final double featureScore, targetScore, volumeIntegral;

		private Solution next, prev;

		/**
		 * Constructor
		 * 
		 * @param solutionMolecule
		 *            Hit molecule.
		 * @param featureScore
		 *            Associated Gape score.
		 */
		private Solution(GaMolecule solutionMolecule, double targetScore,
				double featureScore, double volumeIntegral) {
			this.solutionMolecule = solutionMolecule;
			this.featureScore = featureScore;
			this.targetScore = targetScore;
			this.volumeIntegral = volumeIntegral;
		}

		/**
		 * Adds a new node before this node
		 * 
		 * @param s
		 *            The new node.
		 */
		private void addBefore(Solution s) {
			s.prev = prev;
			s.next = this;
			if (prev != null)
				prev.next = s;
			prev = s;
		}

		/**
		 * Adds a new node after this node
		 * 
		 * @param s
		 *            The new node.
		 */
		private void addAfter(Solution s) {
			s.prev = this;
			s.next = next;
			if (next != null)
				next.prev = s;
			next = s;
		}

		/**
		 * removes this node from the linked list.
		 */
		private void remove() {
			if (prev != null && next != null) {
				prev.next = next;
				next.prev = prev;
			} else if (prev != null)
				prev.next = null;
			else if (next != null)
				next.prev = null;

		}

		/**
		 * Writes a solution out. The structure format is the same as the target
		 * file format. The score is in the MOL2 comment block or as an SD
		 * field.
		 * 
		 * @param out
		 *            Writer to the structure file
		 */
		private void writeSolution(Writer out) throws IOException {
			logger.debug("Writing solution for " + solutionMolecule.getName());
			if (targetFileType == Molecule.FileType.MOL2) {
				String comments = "RIGID_PHARMACOPHORE_SCORE " + nf.format(targetScore)
						+ " FEATURE_SCORE " + nf.format(featureScore) + " VOLUME_SCORE "
						+ nf.format(volumeIntegral);
				solutionMolecule.writeSybylMol2(out, comments);
			} else {
				solutionMolecule.addSdfField("RIGID_PHARMACOPHORE_SCORE",
						nf.format(targetScore));
				solutionMolecule.addSdfField("FEATURE_SCORE", nf.format(featureScore));
				solutionMolecule.addSdfField("VOLUME_SCORE", nf.format(volumeIntegral));
				solutionMolecule.writeSdfMol(out, null);
			}
		}
	}

	/**
	 * This class is used to store summary data for conformations. If we have a
	 * structure we only want to output the best n conformations matches to the
	 * summary writer. We store summary solutions in this class until we've
	 * processed the structure.
	 * 
	 */
	private class ConformerSummarySolution implements
			Comparable<ConformerSummarySolution> {
		private final String name;

		private final double targetScore, featureScore, volumeIntegral;

		private ConformerSummarySolution(String _name, double _targetScore,
				double _featureScore, double _volumeIntegral) {
			name = _name;
			targetScore = _targetScore;
			featureScore = _featureScore;
			volumeIntegral = _volumeIntegral;
		}

		/*
		 * Used to sort the solutions by store so we only write out the
		 * best.(non-Javadoc)
		 * 
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@Override
		public int compareTo(ConformerSummarySolution other) {
			if (targetScore > other.targetScore)
				return -1;
			else if (targetScore < other.targetScore)
				return 1;
			else
				return 0;
		}

		private void print() {

			String summary = name + "\t" + nf.format(targetScore) + "\t"
					+ nf.format(featureScore) + "\t" + nf.format(volumeIntegral) + "\n";
			logger.debug(summary);
			try {
				textSummaryWriter.write(summary);
			} catch (IOException ex) {
				throw new RuntimeException("Can't write to grips_summary.txt");
			}

		}

	}

}
