package com.cairn.gape.feature;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;
import com.cairn.gape.molecule.GaMolecule;
import com.cairn.gape.utils.InfoMessageLogger;
import com.cairn.molecule.Atom;
import com.cairn.molecule.Molecule;

/**
 * Class to represent Donor Hydrogens
 * 
 * @author Gareth Jones
 * 
 */
public class DonorHydrogenFeature extends Feature {

	private static Logger logger;
	private static boolean logDebug;

	static {
		logger = Logger.getLogger(DonorHydrogenFeature.class);
		// logger.setLevel(Level.DEBUG);
		logDebug = logger.isDebugEnabled();
	}
	private volatile Atom donor;

	private volatile double donorCoord[];

	private volatile HydrogenBondingType hydrogenBondingType;

	/**
	 * Length of a hydrogen bond. May be changed.
	 */
	private static final ThreadLocal<Double> hBondLen = new ThreadLocal<Double>() {
		@Override
		protected Double initialValue() {
			return 2.9;
		}
	};

	/**
	 * Pairwise feature score is scaled by this if both features have negative
	 * charge. May be changed.
	 */
	private static final ThreadLocal<Double> chargeFactor = new ThreadLocal<Double>() {
		@Override
		protected Double initialValue() {
			return 2.0;
		}
	};

	/**
	 * Pairwise feature score is scaled by this if both features have the same
	 * type. May be changed.
	 */
	private static final ThreadLocal<Double> matchFactor = new ThreadLocal<Double>() {
		@Override
		protected Double initialValue() {
			return 1.0;
		}
	};

	private static final ThreadLocal<Double> maxDonorDonorAngle = new ThreadLocal<Double>() {
		@Override
		protected Double initialValue() {
			return 120.0 * Math.PI / 180;
		}
	};
	private static final ThreadLocal<Double> minDonorDonorAngle = new ThreadLocal<Double>() {
		@Override
		protected Double initialValue() {
			return 90.0 * Math.PI / 180;
		}
	};

	private static final ThreadLocal<Boolean> scoreDonorAtoms = new ThreadLocal<Boolean>() {
		@Override
		protected Boolean initialValue() {
			return true;
		}
	};

	private boolean charged;

	/**
	 * Constructor
	 */
	private DonorHydrogenFeature() {
		featureType = FeatureType.DONOR_INTERACTION_POINT;
		featureSetName = "DONOR_HYDROGEN";
		atomFeature = false;
	}

	/**
	 * Creates a feature from the donor hydrogen atom
	 * 
	 * @param m
	 * @param a
	 * @param no
	 */
	public DonorHydrogenFeature(GaMolecule m, Atom a, int no) {
		this();
		molecule = m;
		atom = a;
		donor = atom.getNotDummyNeighbours().get(0);
		hydrogenBondingType = donor.getDonorType();
		if (donor.getPartialCharge() > 0.24)
			charged = true;
		// else if (donor.getChargeFromBondOrder() > 0.25)
		// charged = true;
		featureSetNo = no;
	}

	/**
	 * Creates a donor hydrogen feature from a string pharmacophore definition.
	 * 
	 * @param featureName
	 * @param options
	 * @param FeatureDefinion
	 * @param geometry
	 * @param _moleculeNo
	 */
	public DonorHydrogenFeature(String featureName, HashMap<String, String> options,
			String FeatureDefinion, PharmFeatureGeometry geometry, int _moleculeNo) {
		this();

		hydrogenBondingType = HydrogenBondingType.getTypeFromName(featureName);

		if (options.containsKey("charged") && options.get("charged").equals("yes"))
			charged = true;
		virtual = true;
		moleculeNo = _moleculeNo;

		pharmFeatureGeometry = geometry;
		VectorPharmFeatureGeometry g = (VectorPharmFeatureGeometry) geometry;
		donorCoord = g.points[0];
		coordinate = g.points[1];

	}

	/**
	 * Searches a molecules and finds all donor hydrogens.
	 * 
	 * @param m
	 * @return list of donor hydrogen features.
	 */
	public static List<Feature> findFeatures(GaMolecule m) {
		InfoMessageLogger infoMessageLogger = m.getInfoMessageLogger();

		ArrayList<Feature> v = new ArrayList<Feature>();
		infoMessageLogger.infoMessage(2, "Looking for Donor Hydrogens: ");
		int no = 0;
		for (Atom atom : m.getAtoms()) {
			if (atom.isDonorHydrogen()) {
				DonorHydrogenFeature feature = new DonorHydrogenFeature(m, atom, no);
				v.add(feature);

				if (infoMessageLogger.getLogLevel() > 2) {
					infoMessageLogger.infoMessage(String.valueOf(atom.getNo() + 1));
					if (feature.charged)
						infoMessageLogger.infoMessage("[charged]");
					infoMessageLogger.infoMessage(" ");
				}

				no++;
			}
		}

		infoMessageLogger.infoMessageln(2, "");
		return v;
	}

	static private ThreadLocal<double[]> bondVec = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[4];
		}
	};

	/*
	 * Determines the fitting point for the donor hydrogen. (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Feature#getCoordinate()
	 */
	@Override
	public double[] calculateCoordinate() {
		if (virtual)
			return coordinate;

		donorCoord = molecule.getCoord(donor.getNo());
		double hydrogenCoord[] = molecule.getCoord(atom.getNo());
		Coord.subtract(hydrogenCoord, donorCoord, bondVec.get());
		Coord.setLength(bondVec.get(), hBondLen.get());
		Coord.add(donorCoord, bondVec.get(), coordinate);

		return coordinate;
	}

	/*
	 * Descriptive string for the feature. (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Feature#info()
	 */
	@Override
	public String info() {
		if (virtual)
			return "Donor H [virtual]";
		int no = atom.getNo() + 1;
		return "Donor H [" + no + "]";
	}

	/*
	 * Label to be used in descriptions and structure file comments
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Feature#pharmLabel()
	 */
	@Override
	public String pharmLabel() {
		String name = hydrogenBondingType.getName().toUpperCase().replace(' ', '_');
		String options = "";
		if (charged)
			options += "charged = yes";
		if (!options.equals(""))
			options = " [" + options + "]";
		String label = virtual ? "VIRTUAL" : "ATOM_" + String.valueOf(atom.getNo() + 1);
		return featureSetName + " " + name + " " + label + options;
	}

	static private final ThreadLocal<double[]> midPoint = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[4];
		}
	};
	static private final ThreadLocal<double[]> vec1 = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[4];
		}
	};
	static private final ThreadLocal<double[]> vec2 = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[4];
		}
	};

	/*
	 * Returns similarity between two donor hydrogens. (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Feature#score(com.cairn.gape.Feature)
	 */
	@Override
	public double score(Feature f) {
		DonorHydrogenFeature other = (DonorHydrogenFeature) f;

		// fitting point Gaussian
		double vol = Feature.score(Coord.sqrDistance(coordinate, other.coordinate));
		Coord.midPoint(coordinate, other.coordinate, midPoint.get());

		double correct = .0;
		if (!virtual) {
			// Corrections to make sure fitting point is solvent accessible
			double molVol = Feature.solvationPenalty(midPoint.get(), molecule, atom);
			double otherMolVol = Feature.solvationPenalty(midPoint.get(), other.molecule,
					other.atom);

			if (molVol < .0)
				molVol = .0;
			if (otherMolVol < .0)
				otherMolVol = .0;
			correct = molVol > otherMolVol ? molVol : otherMolVol;

			if (logDebug) {
				logger.debug(info() + " " + other.info() + " fitting vol " + vol
						+ " molA " + molVol + " molB " + otherMolVol);
			}
		}
		double fitScore = vol - correct;
		if (fitScore < 0)
			fitScore = 0;

		geometricScore = fitScore;

		double score;
		if (scoreDonorAtoms.get()) {
			// now check that donor atoms are reasonably positioned.
			if (fitScore == 0)
				return .0;
			Coord.subtract(donorCoord, midPoint.get(), vec1.get());
			Coord.subtract(other.donorCoord, midPoint.get(), vec2.get());
			double angle = Coord.angle(vec1.get(), vec2.get());
			if (logDebug)
				logger.debug("Angle " + angle + " max " + maxDonorDonorAngle + " min "
						+ minDonorDonorAngle);
			double donorScore = .0;
			if (angle > maxDonorDonorAngle.get())
				return .0;
			else if (angle < minDonorDonorAngle.get())
				donorScore = 1.0;
			else {
				donorScore = 1 - (angle - minDonorDonorAngle.get())
						/ (maxDonorDonorAngle.get() - minDonorDonorAngle.get());
			}

			if (logDebug)
				logger.debug(" donor score " + donorScore);
			score = fitScore * donorScore;
		} else {
			score = fitScore;
			if (fitScore > maximumGaussianScore.get())
				score = maximumGaussianScore.get();
		}

		if (logDebug)
			logger.debug(" geometric score " + geometricScore);

		if (score < .0)
			return .0;

		// Mills and Dean probabilities for types
		double prob = hydrogenBondingType.getProbability()
				+ other.hydrogenBondingType.getProbability();
		score *= prob;

		// increase for matching types
		if (hydrogenBondingType.getId() == other.hydrogenBondingType.getId())
			score *= matchFactor.get();

		// increase if both donors are charged
		if (charged && other.charged)
			score *= chargeFactor.get();

		if (logDebug)
			logger.debug(" type score " + score);
		return score;
	}

	/**
	 * Writes out donor hydrogen fitting points to a mol2 file. Used for
	 * debugging.
	 * 
	 * @param mol
	 * @param file
	 */
	static private void donorFittingPoints2Mol2(GaMolecule mol, String file) {
		List<Feature> features = mol.getFeatureMappings(
				FeatureType.DONOR_INTERACTION_POINT).getFeatures();

		List<Atom> atoms = new ArrayList<>();
		List<double[]> coords = new ArrayList<>();
		int no = 0;
		for (Feature feature : features) {
			atoms.add(new Atom(no, String.valueOf(no), "Du"));
			coords.add(feature.calculateCoordinate());
			no++;
			atoms.add(new Atom(no, String.valueOf(no), "Du"));
			coords.add(((DonorHydrogenFeature) feature).donorCoord);
			no++;
		}
		Molecule vpMol = new Molecule("Donor Hydrogen Fitting Points", atoms, coords);
		logger.info(file);
		vpMol.writeSybylMol2File(file, "Donor Interaction Points");
	}

	/**
	 * Main routine for testing and debugging
	 * 
	 * @param args
	 */
	public static void main(String args[]) {
		GaMolecule molA = new GaMolecule();
		GaMolecule molB = new GaMolecule();

		HydrogenBondingType.loadParameters();
		molA.loadFile(args[0]);
		HydrogenBondingType.searchMolecule(molA);
		molA.findFeatures();
		molA.getAtomicGaussians();
		donorFittingPoints2Mol2(molA, "donorVPsA.mol2");
		List<Feature> features = molA.getFeatureMappings(
				FeatureType.DONOR_INTERACTION_POINT).getFeatures();
		features.stream().forEach(f -> f.calculateCoordinate());

		molB.loadFile(args[1]);
		HydrogenBondingType.searchMolecule(molB);
		molB.findFeatures();
		molB.getAtomicGaussians();
		donorFittingPoints2Mol2(molA, "donorVPsB.mol2");
		features = molB.getFeatureMappings(FeatureType.DONOR_INTERACTION_POINT)
				.getFeatures();
		features.stream().forEach(f -> f.calculateCoordinate());

		double score = molA.donorHydrogenScore(molB);
		logger.warn("score " + score);

	}

	/*
	 * Returns the appropriate pharmacophore feature geometry definition for
	 * this donor hydrogen.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Feature#getPharmFeatureGeometry()
	 */
	@Override
	public PharmFeatureGeometry getPharmFeatureGeometry() {
		if (virtual)
			return pharmFeatureGeometry;

		double point1[] = molecule.getCoord(donor.getNo());
		double point2[] = coordinate;
		pharmFeatureGeometry = new VectorPharmFeatureGeometry(point1, point2);

		return pharmFeatureGeometry;
	}

	public static double getChargeFactor() {
		return chargeFactor.get();
	}

	public static void setChargeFactor(double chargeFactor) {
		DonorHydrogenFeature.chargeFactor.set(chargeFactor);
	}

	public static double getHBondLen() {
		return hBondLen.get();
	}

	public static void setHBondLen(double bondLen) {
		hBondLen.set(bondLen);
	}

	public static double getMatchFactor() {
		return matchFactor.get();
	}

	public static void setMatchFactor(double matchFactor) {
		DonorHydrogenFeature.matchFactor.set(matchFactor);
	}

	public static double getMaxDonorDonorAngle() {
		return maxDonorDonorAngle.get();
	}

	public static void setMaxDonorDonorAngle(double maxDonorDonorAngle) {
		DonorHydrogenFeature.maxDonorDonorAngle.set(maxDonorDonorAngle);
	}

	public static double getMinDonorDonorAngle() {
		return minDonorDonorAngle.get();
	}

	public static void setMinDonorDonorAngle(double minDonorDonorAngle) {
		DonorHydrogenFeature.minDonorDonorAngle.set(minDonorDonorAngle);
	}

	public static boolean isScoreDonorAtoms() {
		return scoreDonorAtoms.get();
	}

	public static void setScoreDonorAtoms(boolean scoreDonorAtoms) {
		DonorHydrogenFeature.scoreDonorAtoms.set(scoreDonorAtoms);
	}

}
