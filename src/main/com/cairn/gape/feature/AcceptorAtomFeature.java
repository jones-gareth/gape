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
 * Class for representing an atom-centered hydrogen bond acceptor.
 * 
 * @author Gareth Jones
 * 
 */
public class AcceptorAtomFeature extends Feature {

	private static Logger logger = Logger.getLogger(AcceptorAtomFeature.class);
	private static boolean logDebug;
	static {
		// logger.setLevel(Level.DEBUG);
		logDebug = logger.isDebugEnabled();
	}

	/**
	 * We can use a scoring function based on lone-pair positions or a new more
	 * complex function that takes account of differing acceptor geometries, The
	 * default is to use the second. In the acceptor geometry model acceptors
	 * with three sp3 lone pairs are represented by one lone pair in the forward
	 * direction.
	 */
	public static final boolean USE_ACCEPTOR_GEOMETRY = true;

	/**
	 * New improved Mills and Dean classification
	 */
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

	private final double lonePairCoords[][] = new double[4][];

	private int nLonePairs;

	// Coordinate for storing a point that should be solvent accessible
	private double solvationPoint[] = new double[4];

	// Normal coordinate- for planar acceptors
	private final double normal[] = new double[4];

	/**
	 * For scoring functions based on lone pair overlap we can scale totals
	 * based on the number of lone-pairs
	 */
	private static final ThreadLocal<Boolean> scaleLonePairs = new ThreadLocal<Boolean>() {
		@Override
		protected Boolean initialValue() {
			return true;
		}
	};

	// The default is to use Mills and Dean types in preference to GASP/GOLD
	// I've now stopped using the GASP/GOLD types completely

	// Acceptor geometry e.g. LP, PLANE
	private HydrogenBondingType.AcceptorGeometry geometry;

	/**
	 * Set for a charged acceptor
	 */
	private boolean charged;

	private AcceptorAtomFeature() {
		featureType = FeatureType.ACCEPTOR_ATOM;
		featureSetName = "ACCEPTOR_ATOM";
		atomFeature = true;
	}

	/**
	 * Constructor. Created a Feature from an acceptor atom.
	 * 
	 * @param m
	 * @param a
	 * @param no
	 */
	public AcceptorAtomFeature(GaMolecule m, Atom a, int no) {
		this();
		molecule = m;
		atom = a;
		hydrogenBondingType = atom.getAcceptorType();
		geometry = hydrogenBondingType.getGeometry();

		if (atom.getPartialCharge() < -0.24)
			charged = true;
		// else if (atom.getChargeFromBondOrder() < -0.24)
		// charged = true;
		featureSetNo = no;

		nLonePairs = atom.getnLonePairs();
		logger.debug("N Lone Pairs " + nLonePairs);
	}

	/**
	 * Used for constructing a virtual feature from a string description.
	 * 
	 * @param featureName
	 * @param options
	 * @param FeatureDefinion
	 * @param geometry
	 * @param _moleculeNo
	 */
	public AcceptorAtomFeature(String featureName, HashMap<String, String> options,
			String FeatureDefinion, PharmFeatureGeometry _pharmFeatureGeometry,
			int _moleculeNo) {
		this();

		hydrogenBondingType = HydrogenBondingType.getTypeFromName(featureName);
		geometry = hydrogenBondingType.getGeometry();

		virtual = true;
		if (options.containsKey("charged") && options.get("charged").equals("yes"))
			charged = true;
		virtual = true;
		moleculeNo = _moleculeNo;

		pharmFeatureGeometry = _pharmFeatureGeometry;

		if (geometry == HydrogenBondingType.AcceptorGeometry.DIR) {
			if (pharmFeatureGeometry instanceof VectorPharmFeatureGeometry) {
				VectorPharmFeatureGeometry g = (VectorPharmFeatureGeometry) pharmFeatureGeometry;
				coordinate = g.points[0];
				lonePairCoords[0] = g.points[1];
				nLonePairs = 1;
			} else if (pharmFeatureGeometry instanceof MultiVectorPharmFeatureGeometry) {
				MultiVectorPharmFeatureGeometry g = (MultiVectorPharmFeatureGeometry) pharmFeatureGeometry;
				coordinate = g.getCenter();
				nLonePairs = g.getNVectors();
				for (int i = 0; i < nLonePairs; i++)
					lonePairCoords[i] = g.getEnd(i);
			}
		}

		else if (geometry == HydrogenBondingType.AcceptorGeometry.PLANE) {
			ArcPharmFeatureGeometry g = (ArcPharmFeatureGeometry) pharmFeatureGeometry;
			coordinate = g.points[0];
			lonePairCoords[0] = g.points[1];
			lonePairCoords[1] = g.points[2];
			nLonePairs = 2;

			if (USE_ACCEPTOR_GEOMETRY) {
				if (geometry == HydrogenBondingType.AcceptorGeometry.PLANE)
					getPlaneNormal();
			}
		}

		else if (geometry == HydrogenBondingType.AcceptorGeometry.CONE) {
			ConePharmFeatureGeometry g = (ConePharmFeatureGeometry) pharmFeatureGeometry;
			coordinate = g.points[0];
			double[] mid = new double[4], diff = new double[4];
			mid[3] = 1.0;
			diff[3] = 1.0;
			Coord.midPoint(g.points[1], g.points[2], mid);
			Coord.subtract(mid, coordinate, diff);
			double d = Coord.distance(coordinate, g.points[1]);
			Coord.setLength(diff, d);
			Coord.add(coordinate, diff, mid);
			lonePairCoords[0] = mid;
			nLonePairs = 1;
		}

		else if (geometry == HydrogenBondingType.AcceptorGeometry.AG_NONE) {
			SpherePharmFeatureGeometry g = (SpherePharmFeatureGeometry) pharmFeatureGeometry;
			coordinate = g.point;
			nLonePairs = 0;
		}

		else
			throw new RuntimeException("Unable to create acceptor from definition "
					+ FeatureDefinion);

	}

	/**
	 * Finds all acceptors in a molecule and returns them in an ArrayList.
	 * 
	 * @param m
	 * @return
	 */
	public static ArrayList<Feature> findFeatures(GaMolecule m) {
		InfoMessageLogger infoMessageLogger = m.getInfoMessageLogger();
		ArrayList<Feature> v = new ArrayList<Feature>();

		infoMessageLogger.infoMessage(2, "Looking for Acceptors: ");
		int no = 0;
		for (Atom atom : m.getAtoms()) {
			// if ((!useMillsAndDeanTypes && atom.donorAcceptorType != null &&
			// atom.donorAcceptorType.acceptor)
			// || (useMillsAndDeanTypes && atom.acceptorType != null)) {
			if (atom.getAcceptorType() != null) {
				AcceptorAtomFeature feature = new AcceptorAtomFeature(m, atom, no);
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Feature#<()
	 * 
	 * Gets the feature coordinate, which is the acceptor atom.
	 */
	@Override
	public double[] calculateCoordinate() {
		coordinate = super.calculateCoordinate();

		for (int i = 0; i < nLonePairs; i++)
			lonePairCoords[i] = molecule.getCoord(atom.getLonePairs().get(i).getNo());

		if (USE_ACCEPTOR_GEOMETRY) {
			getSolvationPoint();
			if (geometry == HydrogenBondingType.AcceptorGeometry.PLANE)
				getPlaneNormal();
		}

		return coordinate;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Feature#info()
	 */
	@Override
	public String info() {
		if (virtual)
			return "Acceptor [virtual]";
		int no = atom.getNo() + 1;
		return "Acceptor [" + no + "]";
	}

	/*
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
		return featureSetName + " " + name + " " + label + " " + options;
	}

	static private final ThreadLocal<double[]> midPoint = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[4];
		}
	};

	/*
	 * Scores the overlap between two acceptors. This function contains old code
	 * that uses overlap lone pair positions in addition to acceptor atom
	 * overlap.
	 * 
	 * It's currently defined to call the new preferred
	 * acceptorGeometryScoreFeature.
	 * 
	 * @see com.cairn.gape.Feature#score(com.cairn.gape.Feature)
	 * 
	 * @see #acceptorGeometryScore
	 */
	@Override
	public double score(Feature f) {

		if (USE_ACCEPTOR_GEOMETRY)
			return acceptorGeometryScore(f);

		AcceptorAtomFeature other = (AcceptorAtomFeature) f;

		geometricScore = .0;

		double accScore = Feature.score(Coord.sqrDistance(coordinate, other.coordinate));
		if (accScore > maximumGaussianScore.get())
			accScore = maximumGaussianScore.get();

		double lpScore = .0;

		for (int i = 0; i < nLonePairs; i++) {
			for (int j = 0; j < other.nLonePairs; j++) {

				double vol = Feature.score(Coord.sqrDistance(lonePairCoords[i],
						other.lonePairCoords[j]));
				Coord.midPoint(lonePairCoords[i], other.lonePairCoords[j], midPoint.get());

				// Corrections to make sure fitting point is solvent accessible
				double molVol = virtual ? 0 : Feature.solvationPenalty(midPoint.get(),
						molecule, atom);
				double otherMolVol = virtual ? 0 : Feature.solvationPenalty(
						midPoint.get(), other.molecule, other.atom);

				if (molVol < .0)
					molVol = .0;
				if (otherMolVol < .0)
					otherMolVol = .0;
				double correct = molVol > otherMolVol ? molVol : otherMolVol;

				double score = vol - correct;
				if (score < .0)
					score = .0;
				if (logDebug)
					logger.debug("Lp " + atom.getLonePairs().get(i).info() + " Other Lp "
							+ other.atom.getLonePairs().get(j).info() + " score " + score
							+ " vol " + vol + " molVol " + molVol + " otherMolVol "
							+ otherMolVol);

				if (score > maximumGaussianScore.get())
					score = maximumGaussianScore.get();

				// lp score shouldn't be greater that the acceptor score
				if (score > accScore)
					score = accScore;

				lpScore += score;
			}
		}

		if (scaleLonePairs.get())
			lpScore = lpScore / nLonePairs;

		// double score = (lpScore+accScore)/2.0;
		double score = lpScore;

		geometricScore = score;

		double prob = hydrogenBondingType.getProbability()
				+ other.hydrogenBondingType.getProbability();
		score *= prob;
		if (hydrogenBondingType.getId() == other.hydrogenBondingType.getId())
			score *= matchFactor.get();

		if (charged && other.charged)
			score *= chargeFactor.get();

		if (logDebug) {
			logger.debug(info() + " " + other.info() + " score " + score
					+ " gemoetric score " + geometricScore + " accScore " + accScore
					+ " lpScore " + lpScore);
		}

		return score;
	}

	/**
	 * Writes out the acceptor fitting points to a mol2 file. Used for
	 * debugging.
	 * 
	 * @param mol
	 * @param file
	 */
	private static void acceptorFittingPoints2Mol2(GaMolecule mol, String file) {
		List<Feature> features = mol.getFeatureMappings(FeatureType.ACCEPTOR_ATOM)
				.getFeatures();

		List<Atom> atoms = new ArrayList<>();
		List<double[]> coords = new ArrayList<>();
		for (int i = 0; i < features.size(); i++) {
			int no = features.get(i).atom.getNo() + 1;
			atoms.add(new Atom(i, String.valueOf(no), "Du"));
			coords.add(features.get(i).calculateCoordinate());
		}
		Molecule vpMol = new Molecule("Acceptor Fitting Points", atoms, coords);
		System.out.println(file);
		vpMol.writeSybylMol2File(file, "Acceptor Interaction Points");
	}

	static private final ThreadLocal<double[]> midPoint2 = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[4];
		}
	};

	static private final ThreadLocal<double[]> diff = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[4];
		}
	};

	/**
	 * Determines a point that should be solvent accessible if an acceptor with
	 * one or two lone pairs is able to accept a hydrogen bond. Note that in the
	 * acceptor geometry model an acceptor with three sp3 lone pairs has a
	 * single lone pair in the forward direction.
	 */
	private void getSolvationPoint() {
		if (nLonePairs == 1)
			solvationPoint = lonePairCoords[0];
		else if (nLonePairs == 2) {
			Coord.midPoint(lonePairCoords[0], lonePairCoords[1], midPoint2.get());
			Coord.subtract(midPoint2.get(), coordinate, diff.get());
			Coord.setLength(diff.get(), hBondLen.get());
			Coord.add(coordinate, diff.get(), solvationPoint);
		} else
			System.out.println("getSolvationPoint: molecule " + molecule.getName()
					+ " atom " + atom.info() + " has " + nLonePairs + " lone pairs");
	}

	/**
	 * Gets the normal for two lone pairs.
	 */
	private void getPlaneNormal() {
		Coord.subtract(lonePairCoords[0], coordinate, vec1.get());
		Coord.subtract(lonePairCoords[1], coordinate, vec2.get());
		Coord.vectorProduct(vec1.get(), vec2.get(), normal);
	}

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

	/**
	 * Determines the angle between vectors (p1-p2) and (p3-p4).
	 * 
	 * @param p1
	 * @param p2
	 * @param p3
	 * @param p4
	 * @return
	 */
	private static double angle(double p1[], double p2[], double p3[], double p4[]) {
		Coord.subtract(p1, p2, vec1.get());
		Coord.subtract(p3, p4, vec2.get());
		return Coord.angle(vec1.get(), vec2.get());
	}

	/**
	 * Scoring function for an angle constraint. Returns 1 id angle < min, 0 if
	 * angle > max and linear interpolates otherwise.
	 * 
	 * @param angle
	 * @param max
	 * @param min
	 * @return
	 */
	private static double scoreAngle(double angle, double max, double min) {
		if (angle > max)
			return .0;
		else if (angle < min)
			return 1.0;
		else {
			return 1 - (angle - min) / (max - min);
		}
	}

	/**
	 * Main scoring function.
	 * 
	 * The acceptor score includes the following:
	 * 
	 * Gaussian overlay of acceptor centers Solvent accessibility of acceptors
	 * Acceptor geometry compatibility Check both acceptors pointing in the same
	 * direction Mills and Dean acceptor probabilities Scaling for charge
	 * Scaling for matching types
	 * 
	 * @param f
	 * @return
	 * 
	 * @see #lonePairLonePairScore(AcceptorAtomFeature)
	 * @see #planeLonePairScore(AcceptorAtomFeature)
	 * @see #planePlaneScore(AcceptorAtomFeature)
	 */
	private double acceptorGeometryScore(Feature f) {

		AcceptorAtomFeature other = (AcceptorAtomFeature) f;
		geometricScore = .0;

		// Score for acceptor overlay
		double accScore = Feature.score(Coord.sqrDistance(coordinate, other.coordinate));
		if (accScore > maximumGaussianScore.get())
			accScore = maximumGaussianScore.get();

		// Corrections to make sure fitting point is solvent accessible
		double correct = .0;
		if (!virtual) {
			Coord.midPoint(solvationPoint, other.solvationPoint, midPoint.get());
			double molVol = Feature.solvationPenalty(midPoint.get(), molecule, atom);
			double otherMolVol = Feature.solvationPenalty(midPoint.get(), other.molecule,
					other.atom);

			correct = molVol > otherMolVol ? molVol : otherMolVol;
			if (correct < .0)
				correct = .0;
			accScore -= correct;
			if (accScore <= 0)
				return 0;
		}

		// Both acceptors pointing in the same direction?
		double forwardScore = forwardScore(other);
		if (forwardScore <= .0)
			return .0;

		// Compatable geometries?
		double geometryScore = .0;
		HydrogenBondingType.AcceptorGeometry otherGeometry = other.geometry;

		if (geometry == HydrogenBondingType.AcceptorGeometry.AG_NONE
				|| otherGeometry == HydrogenBondingType.AcceptorGeometry.AG_NONE
				|| geometry == HydrogenBondingType.AcceptorGeometry.CONE
				|| otherGeometry == HydrogenBondingType.AcceptorGeometry.CONE)
			// For acceptors that accept in a cone or have no directionality
			// assume that the general forward constraint is sufficient.
			geometryScore = 1.0;

		else if (geometry == HydrogenBondingType.AcceptorGeometry.DIR
				&& otherGeometry == HydrogenBondingType.AcceptorGeometry.DIR)
			geometryScore = lonePairLonePairScore(other);

		else if (geometry == HydrogenBondingType.AcceptorGeometry.PLANE
				&& otherGeometry == HydrogenBondingType.AcceptorGeometry.PLANE)
			geometryScore = planePlaneScore(other);

		else if (geometry == HydrogenBondingType.AcceptorGeometry.PLANE
				&& otherGeometry == HydrogenBondingType.AcceptorGeometry.DIR)
			geometryScore = planeLonePairScore(other);

		else if (geometry == HydrogenBondingType.AcceptorGeometry.DIR
				&& otherGeometry == HydrogenBondingType.AcceptorGeometry.PLANE)
			geometryScore = other.planeLonePairScore(this);

		assert !Double.isNaN(geometryScore) : "geometry score is NaN";

		double score = accScore * forwardScore * geometryScore;
		geometricScore = score;

		// Acceptor strengths/probabilities

		// Type matching scale up
		double prob = hydrogenBondingType.getProbability()
				+ other.hydrogenBondingType.getProbability();
		score *= prob;
		if (hydrogenBondingType.getId() == other.hydrogenBondingType.getId())
			score *= matchFactor.get();

		// charge scale up
		if (charged && other.charged)
			score *= chargeFactor.get();

		if (logger.isDebugEnabled()) {
			logger.debug(info() + " " + other.info() + " score " + score
					+ " geometric score " + geometricScore + " geometry score "
					+ geometryScore + " forward score " + forwardScore
					+ " solvation correction " + correct + " accScore " + accScore);
		}

		return score;
	}

	/**
	 * Test routine for checking scoring function
	 * 
	 * @param args
	 */
	public static void main(String args[]) {
		GaMolecule molA = new GaMolecule();
		GaMolecule molB = new GaMolecule();

		HydrogenBondingType.loadParameters();
		molA.loadFile(args[0]);
		// DonorAcceptorType.findMoleculeDonorAcceptors(molA);
		HydrogenBondingType.searchMolecule(molA);
		LonePairAddition.addLonePairs(molA, 2.9);
		molA.findFeatures();
		molA.getAtomicGaussians();
		acceptorFittingPoints2Mol2(molA, "acceptorsA.mol2");
		List<Feature> features = molA.getFeatureMappings(FeatureType.ACCEPTOR_ATOM)
				.getFeatures();
		features.stream().forEach(f -> f.calculateCoordinate());

		molB.loadFile(args[1]);
		// DonorAcceptorType.findMoleculeDonorAcceptors(molB);
		HydrogenBondingType.searchMolecule(molB);
		LonePairAddition.addLonePairs(molB);
		molB.findFeatures();
		molB.getAtomicGaussians();
		acceptorFittingPoints2Mol2(molB, "acceptorsB.mol2");
		features = molB.getFeatureMappings(FeatureType.ACCEPTOR_ATOM).getFeatures();
		features.stream().forEach(f -> f.calculateCoordinate());

		double score = molA.acceptorAtomScore(molB);
		System.out.println("score " + score);

	}

	/**
	 * Simple label for the acceptor
	 * 
	 * @return
	 */
	public String label() {
		int no = featureSetNo + 1;
		return "ACCEPTOR_" + no;
	}

	// Check both acceptors are facing forward
	private static final ThreadLocal<Double> maxForwardAcceptorAngle = new ThreadLocal<Double>() {
		@Override
		protected Double initialValue() {
			return 90.0 * Math.PI / 180;
		}
	};
	private static final ThreadLocal<Double> minForwardAcceptorAngle = new ThreadLocal<Double>() {
		@Override
		protected Double initialValue() {
			return 60.0 * Math.PI / 180;
		}
	};

	/**
	 * Using the solvation point and angle constraints check that both acceptors
	 * are able to accept from roughly the same direction.
	 * 
	 * @param other
	 * @return
	 */
	private double forwardScore(AcceptorAtomFeature other) {
		double angle = angle(solvationPoint, coordinate, other.solvationPoint,
				other.coordinate);
		return scoreAngle(angle, maxForwardAcceptorAngle.get(),
				minForwardAcceptorAngle.get());
	}

	// Lone pairs -- Lone Pairs directionality score

	private static final ThreadLocal<Double> maxLonePairLonePairAngle = new ThreadLocal<Double>() {
		@Override
		protected Double initialValue() {
			return 60.0 * Math.PI / 180;
		}
	};
	private static final ThreadLocal<Double> minLonePairLonePairAngle = new ThreadLocal<Double>() {
		@Override
		protected Double initialValue() {
			return 30.0 * Math.PI / 180;
		}
	};

	/**
	 * For two acceptors that each accept along a single lone pair, use angle
	 * constraints to check that they have compatible orientation.
	 * 
	 * @param other
	 * @return
	 */
	private double lonePairLonePairScore(AcceptorAtomFeature other) {
		double score = 0;
		for (int i = 0; i < nLonePairs; i++) {
			double minAngle = Math.PI;
			for (int j = 0; j < other.nLonePairs; j++) {
				double angle = angle(lonePairCoords[i], coordinate,
						other.lonePairCoords[j], other.coordinate);
				if (minAngle > angle)
					minAngle = angle;
			}
			if (logger.isDebugEnabled())
				logger.debug("LP LP angle " + String.valueOf(minAngle * 180.0 / Math.PI));
			score += scoreAngle(minAngle, maxLonePairLonePairAngle.get(),
					minLonePairLonePairAngle.get());
		}

		if (scaleLonePairs.get())
			score = score / nLonePairs;
		return score;
	}

	// Planar -- Planar acceptor score
	private static final ThreadLocal<Double> maxPlanePlaneAngle = new ThreadLocal<Double>() {
		@Override
		protected Double initialValue() {
			return 60.0 * Math.PI / 180;
		}
	};
	private static final ThreadLocal<Double> minPlanePlaneAngle = new ThreadLocal<Double>() {
		@Override
		protected Double initialValue() {
			return 30.0 * Math.PI / 180;
		}
	};

	/**
	 * For two acceptors that each accept along a plane of two lone pairs, use
	 * angle constraints to check that they have compatible orientation.
	 * 
	 * @param other
	 * @return
	 */
	private double planePlaneScore(AcceptorAtomFeature other) {

		double angle = Coord.angle(normal, other.normal);
		if (logDebug)
			logger.debug("Plane Plane angle " + String.valueOf(angle * 180.0 / Math.PI));
		if (angle > Math.PI / 2.0)
			angle = Math.PI - angle;
		if (logDebug)
			logger.debug("Plane Plane angle " + String.valueOf(angle * 180.0 / Math.PI));
		return scoreAngle(angle, maxPlanePlaneAngle.get(), minPlanePlaneAngle.get());
	}

	// Plane -- Lone pair acceptor score. Acceptor other has the LonePairs.
	private static final ThreadLocal<Double> maxPlaneLonePairAngle = new ThreadLocal<Double>() {
		@Override
		protected Double initialValue() {
			return 60.0 * Math.PI / 180;
		}
	};
	private static final ThreadLocal<Double> minPlaneLonePairAngle = new ThreadLocal<Double>() {
		@Override
		protected Double initialValue() {
			return 30.0 * Math.PI / 180;
		}
	};

	static private final ThreadLocal<double[]> vec3 = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[4];
		}
	};

	/**
	 * For an acceptors that accepts along a single lone pair and another that
	 * accepts in the plane of two lone pairs, use angle constraints to check
	 * that they have compatible orientation.
	 * 
	 * @param other
	 * @return
	 */
	private double planeLonePairScore(AcceptorAtomFeature other) {
		double score = 0;
		for (int i = 0; i < other.nLonePairs; i++) {
			Coord.subtract(other.lonePairCoords[i], other.coordinate, vec3.get());
			double angle = Coord.angle(normal, vec3.get());
			angle = angle - Math.PI / 2.0;
			if (angle < 0)
				angle = -angle;
			if (logger.isDebugEnabled())
				logger.debug("Plane LP angle " + String.valueOf(angle * 180.0 / Math.PI));
			score += scoreAngle(angle, maxPlaneLonePairAngle.get(),
					minPlaneLonePairAngle.get());
		}
		if (scaleLonePairs.get())
			score = score / other.nLonePairs;
		return score;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Feature#getPharmFeatureGeometry()
	 * 
	 * Return the appropriate pharmacophore feature definition for this
	 * acceptor.
	 */
	@Override
	public PharmFeatureGeometry getPharmFeatureGeometry() {
		// Assume that we're using ACCEPTOR_GEOMETRY
		if (!USE_ACCEPTOR_GEOMETRY)
			return null;

		double[] point1 = null, point2 = null;
		switch (geometry) {

		case DIR:
			point1 = coordinate;
			if (nLonePairs == 1) {
				point2 = lonePairCoords[0];
				return new VectorPharmFeatureGeometry(point1, point2);
			} else {
				return new MultiVectorPharmFeatureGeometry(nLonePairs, point1,
						lonePairCoords);
			}

		case PLANE:
			point1 = coordinate;
			point2 = lonePairCoords[0];
			double point3[] = lonePairCoords[1];
			return new ArcPharmFeatureGeometry(point1, point2, point3);

		case CONE:
			// Cone acceptor is defined by a single lone pair , but for the cone
			// feature we want two opposite points on the top of the cone. We
			// set the cone angle to be that for lone pairs in an sp3 setting
			// and hi-jack the LonePairAddition code to generate two points.
			point1 = coordinate;
			double forward[] = lonePairCoords[0];
			double diff[] = new double[4];
			double backward[] = new double[4];
			point2 = new double[4];
			point3 = new double[4];
			diff[3] = backward[3] = point2[3] = point3[3] = 1.0;
			Coord.subtract(point1, forward, diff);
			Coord.add(point1, diff, backward);

			// in sp3 system angle between lp is 109.47
			double angle = 109.47 * Math.PI / 180.0;
			// This is the size of the cone edge- so that cone length is the
			// same
			// as the lone pair distance.
			double dist = Coord.mag(diff);
			LonePairAddition.addTwoPairsRandomlyToTrigonal(point1, backward, point2,
					point3, angle / 2, dist);
			return new ConePharmFeatureGeometry(point1, point2, point3);

		case AG_NONE:
			point1 = coordinate;
			return new SpherePharmFeatureGeometry(point1, hBondLen.get());

		}
		return null;
	}

	public static double getChargeFactor() {
		return chargeFactor.get();
	}

	public static void setChargeFactor(double chargeFactor) {
		AcceptorAtomFeature.chargeFactor.set(chargeFactor);
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
		AcceptorAtomFeature.matchFactor.set(matchFactor);
	}

	public static double getMaxForwardAcceptorAngle() {
		return maxForwardAcceptorAngle.get();
	}

	public static void setMaxForwardAcceptorAngle(double maxForwardAcceptorAngle) {
		AcceptorAtomFeature.maxForwardAcceptorAngle.set(maxForwardAcceptorAngle);
	}

	public static double getMaxLonePairLonePairAngle() {
		return maxLonePairLonePairAngle.get();
	}

	public static void setMaxLonePairLonePairAngle(double maxLonePairLonePairAngle) {
		AcceptorAtomFeature.maxLonePairLonePairAngle.set(maxLonePairLonePairAngle);
	}

	public static double getMaxPlaneLonePairAngle() {
		return maxPlaneLonePairAngle.get();
	}

	public static void setMaxPlaneLonePairAngle(double maxPlaneLonePairAngle) {
		AcceptorAtomFeature.maxPlaneLonePairAngle.set(maxPlaneLonePairAngle);
	}

	public static double getMaxPlanePlaneAngle() {
		return maxPlanePlaneAngle.get();
	}

	public static void setMaxPlanePlaneAngle(double maxPlanePlaneAngle) {
		AcceptorAtomFeature.maxPlanePlaneAngle.set(maxPlanePlaneAngle);
	}

	public static double getMinForwardAcceptorAngle() {
		return minForwardAcceptorAngle.get();
	}

	public static void setMinForwardAcceptorAngle(double minForwardAcceptorAngle) {
		AcceptorAtomFeature.minForwardAcceptorAngle.set(minForwardAcceptorAngle);
	}

	public static double getMinLonePairLonePairAngle() {
		return minLonePairLonePairAngle.get();
	}

	public static void setMinLonePairLonePairAngle(double minLonePairLonePairAngle) {
		AcceptorAtomFeature.minLonePairLonePairAngle.set(minLonePairLonePairAngle);
	}

	public static double getMinPlaneLonePairAngle() {
		return minPlaneLonePairAngle.get();
	}

	public static void setMinPlaneLonePairAngle(double minPlaneLonePairAngle) {
		AcceptorAtomFeature.minPlaneLonePairAngle.set(minPlaneLonePairAngle);
	}

	public static double getMinPlanePlaneAngle() {
		return minPlanePlaneAngle.get();
	}

	public static void setMinPlanePlaneAngle(double minPlanePlaneAngle) {
		AcceptorAtomFeature.minPlanePlaneAngle.set(minPlanePlaneAngle);
	}

	public static boolean isScaleLonePairs() {
		return scaleLonePairs.get();
	}

	public static void setScaleLonePairs(boolean scaleLonePairs) {
		AcceptorAtomFeature.scaleLonePairs.set(scaleLonePairs);
	}

}
