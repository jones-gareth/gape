package com.cairn.gape.feature;

import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.util.FastMath;
import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;
import com.cairn.gape.molecule.GaMolecule;
import com.cairn.gape.utils.GaussianList;
import com.cairn.molecule.Atom;
import com.cairn.molecule.Molecule;

/**
 * This class is used to generally describe pharmacophore features. It's
 * sub-classed for actual features.
 * 
 * Features are represented by Gaussians and the feature overlay scores
 * determined by Gaussian integrals. If the feature contains geometric
 * constraints such as angles these are generally represented by bounded ranges
 * with linear drop-offs.
 * 
 * @author Gareth Jones
 * 
 * @see com.cairn.gape.feature.HydrophobicFeature
 * @see com.cairn.gape.feature.AcceptorAtomFeature
 * @see com.cairn.gape.feature.HydrophobicFeature
 * @see com.cairn.gape.feature.UserFeature
 */
public abstract class Feature {
	private static Logger logger;
	static {
		logger = Logger.getLogger(Feature.class);
		// logger.setLevel(Level.DEBUG);
	}

	// The ordering here is important as we also want to use these
	// as indices in the molecule featureMapping array
	public enum FeatureType {
		HYDROPHOBIC_ATOM, DONOR_INTERACTION_POINT, ACCEPTOR_ATOM, AROMATIC_RING,

		USER_FEATURES1, USER_FEATURES2, USER_FEATURES3, USER_FEATURES4, USER_FEATURES5,

		USER_FEATURES6, USER_FEATURES7, USER_FEATURES8, USER_FEATURES9, USER_FEATURES10
	};

	/**
	 * This field corrsponds to one of the final int above - for example
	 * HYDROPHOBIC_ATOM
	 */
	protected volatile FeatureType featureType;

	/**
	 * This is the feature number- e.g. HYDROPHOBIC_ATOM number 2
	 */
	protected volatile int featureSetNo;

	/**
	 * A descriptive name for the feature
	 */
	protected volatile String featureSetName;

	/**
	 * Set this if the feature is atom-centered, e.g. an acceptor
	 */
	protected boolean atomFeature;

	/**
	 * This is set if created the feature from a string- i.e it's a
	 * pharmacophore feature without any associated molecule.
	 */
	protected volatile boolean virtual;

	/**
	 * Atom if feature is atom centered
	 */
	protected volatile Atom atom;

	/**
	 * Molecule containing this feature
	 */
	protected volatile GaMolecule molecule;

	/**
	 * For virtual features we can store the molecule itself, but we can see if
	 * pharmacophore features come from the same molecule.
	 */
	protected volatile int moleculeNo;

	/**
	 * Set this if we've stored a mapping
	 */
	protected volatile boolean mapped;

	/**
	 * If this feature (in the base molecule) is mapped to a feature in another
	 * molecule, we can store it here.
	 */
	protected volatile Feature otherFeature;

	// sqr distance to other feature
	protected volatile double sqrDist;

	// These variables are used to track the best mappings so far- for
	// resolving one-to-one mappings.
	protected volatile double bestScore, bestGeometricScore;

	/**
	 * If this feature (in another molecule) is mappped to a feature in the base
	 * molecule, we can store the base molecule feature here.
	 */
	protected volatile Feature baseFeature;

	/**
	 * A feature normally has an point which can be considered it's interaction
	 * point. E.g. for acceptors the atom, for donors the virtual acceptor on
	 * the donor-hydrogen line and for rings the ring center. We can use this
	 * point in chromosome encoding for mapping features and place Gaussians on
	 * this point for scoring feature overlap.
	 */
	protected volatile double coordinate[] = new double[4];

	/**
	 * Set this true if we consider this feature to be part of a pharmacophore.
	 */
	protected volatile boolean pharmPoint;

	/**
	 * number of molecules significantly contribution to this feature if it is a
	 * pharmacophore point.
	 * 
	 */
	protected volatile int nMatched;

	/**
	 * Set this true if we can use this feature in chromosome encoding.
	 */
	private volatile boolean mappingFeature = true;

	static protected ThreadLocal<Double> maximumGaussianScore = new ThreadLocal<Double>() {
		@Override
		protected Double initialValue() {
			return 2.0;
		}
	};
	static protected ThreadLocal<Double> solventVolOk = new ThreadLocal<Double>() {
		@Override
		protected Double initialValue() {
			return .0;
		}
	};

	/**
	 * The score function must set this variable
	 */
	protected volatile double geometricScore;

	public static final int N_BUILTIN_FEATURES = 4, MAX_FEATURES = 14;

	/**
	 * The current number of user-defined feature sets.
	 */
	// public int nUserFeatureSets = 0;

	/**
	 * Gaussian parameters for the core built in features.
	 */
	private static final ThreadLocal<Double> alpha = new ThreadLocal<Double>(),
			gaussianN = new ThreadLocal<Double>(), radius = new ThreadLocal<Double>();
	// fitting point Gaussian: radius of 2.0
	static {
		setRadius(2.0);
	}

	/**
	 * Geometry of feature
	 */
	protected volatile PharmFeatureGeometry pharmFeatureGeometry;

	protected Feature() {
	}

	/**
	 * Each subclass must implement the score method which returns the
	 * similarity between two features.
	 * 
	 * @param f
	 * @return
	 */
	abstract public double score(Feature f);

	/**
	 * Returns a label for the feature to use in mol2 or sdf pharmacophore
	 * files.
	 * 
	 * @return
	 */
	public String atomLabel() {
		int no = featureSetNo + 1;
		return featureSetName + "_" + no;
	}

	/**
	 * Generates a label to use to describe the pharmacophore and feature
	 * geometry in the descriptions that appear in structure files (either in
	 * sd-field or a comment section in mol2 files). Similar to pharmLabel, but
	 * includes feature geometry definition
	 * 
	 * @return
	 * 
	 * @see #pharmLabel()
	 */
	public String featureLabel() {
		String pharmLabel = pharmLabel();
		PharmFeatureGeometry featureGeometry = getPharmFeatureGeometry();
		String featureLabel = pharmLabel + " " + featureGeometry.summary();
		if (isPharmPoint())
			featureLabel = featureLabel + " N_MATCHED=" + getnMatched();
		return featureLabel;
	}

	/**
	 * Generates a label to use to describe the pharmacophore in the
	 * descriptions that appear in structure files (either in sd-field or a
	 * comment section in mol2 files).
	 * 
	 * @return
	 */
	abstract public String pharmLabel();

	/**
	 * Returns an information/description string for this pharmacophore.
	 * 
	 * @return
	 */
	abstract public String info();

	/**
	 * Returns the geometry associated with this feature
	 * 
	 * @return
	 */
	abstract public PharmFeatureGeometry getPharmFeatureGeometry();

	/**
	 * Initializes a feature by building it from the string generated by
	 * featureLabel. There should be a sub-class string constructor which uses
	 * this method and the constructor should set virtual to true.
	 * 
	 * @see #featureLabel()
	 */
	// abstract void loadFromString(String s) ;
	/**
	 * Print the information string.
	 * 
	 * @see #info
	 */
	@Deprecated
	public void printInfo() {
		System.out.println(info());
	}

	/**
	 * Set the maximum Gaussian overlay score that you can have for acceptors.
	 * Gaussian overlay varies between 0 and 1. It the Gaussian overlay score is
	 * greater than s, it will be set to s.
	 * 
	 * @param s
	 */
	public static void setMaximumGaussianScore(double s) {
		maximumGaussianScore.set(s);
	}

	/**
	 * Sets the gaussian radius to r. Adjusts the other gausian parameters so
	 * that the overlap is normalized. This is used only by the built-in
	 * features and not the user-defined features, which have their own defined
	 * radius.
	 * 
	 * @param r
	 */
	public static void setRadius(double r) {
		double alpha = getAlpha(r);
		Feature.alpha.set(alpha);
		gaussianN.set(FastMath.pow((2.0 * alpha) / Math.PI, 0.75));
		radius.set(r);
		logger.debug("Setting alpha " + alpha + " gaussianN " + gaussianN.get()
				+ " radius " + r);
	}

	/**
	 * Determines the integral between two gaussians. The square distance
	 * between gaussian centers is given. The gaussian parameter for alpha is
	 * taken for the class variable, so this method only works for the class
	 * features.
	 * 
	 * @param sqrDistance
	 * @return
	 */
	public static double score(double sqrDistance) {
		return score(sqrDistance, alpha.get());
	}

	/**
	 * Determines the integral between two gaussians. The square distance
	 * between gaussian centers is given as is the alpha parameter for the
	 * gaussians (both gaussians have the same alpha parameter).
	 * 
	 * @param sqrDistance
	 * @param a
	 * @return
	 */
	public static double score(double sqrDistance, double a) {
		double val = -0.5 * a * sqrDistance;
		if (GaussianList.USE_EXP_LOOKUP)
			return GaussianList.lookupExp(val);
		else
			return FastMath.exp(val);
	}

	/**
	 * Determine a gaussian alpha parameter from a "radius". It's more intuitive
	 * to specify a radius rather that alpha.
	 * 
	 * @param r
	 * @return
	 */
	public static double getAlpha(double r) {
		return (-2.0 * Math.log(0.5)) / (r * r);
	}

	/**
	 * The alpha used in a gaussian to determine if a point is solvent
	 * accessible.
	 */
	protected static double solvationAlpha = getAlpha(1.0);

	/**
	 * Set the alpha used in a gaussian to determine if a point is solvent
	 * accessible.
	 * 
	 * @param r
	 */
	public static void setSolvationAlpha(double r) {
		solvationAlpha = getAlpha(r);
	}

	/**
	 * Returns a penalty if point is not solvent accessible. Used to determine
	 * if acceptor lone pairs or donor hydrogen fitting points are available. A
	 * gaussian of alpha solvationAlpha is placed at point and we check for
	 * overlap with any atoms in mol. Overlap with a specific atom is ignored
	 * (this atom being the acceptor or donor hydrogen).
	 * 
	 * @param point
	 * @param mol
	 * @param atom
	 * @return
	 */
	static public double solvationPenalty(double point[], Molecule mol, Atom atom) {
		double penalty = .0;

		for (Atom atom2 : mol.getAtoms()) {
			if (!atom2.isNotDummy())
				continue;
			if (atom2 == atom)
				continue;
			double sqrDist = Coord.sqrDistance(point, mol.getCoord(atom2.getNo()));
			double val = -0.5 * solvationAlpha * sqrDist;
			if (GaussianList.USE_EXP_LOOKUP)
				penalty += GaussianList.lookupExp(val);
			else
				penalty += FastMath.exp(val);
		}
		return penalty - solventVolOk.get();
	}

	/**
	 * Returns the atom-centered coordinate. As a side effect copies the
	 * coordinate from the molecule to this class- if we do this other methods
	 * are free to manipulate the molecule coordinates without affecting the
	 * feature coordinate. For features without atom-centered coordinates
	 * (rings, donor hydrogens and some user-features) you need to override this
	 * method.
	 * 
	 * @return
	 */
	public double[] calculateCoordinate() {
		// if (atomFeature)
		Coord.copy(molecule.getCoord(atom.getNo()), coordinate);
		// feature coordinates not implemented
		// else
		// Coord.copy(molecule.featureCoords[coordinateNo], coordinate);
		return coordinate;
	}

	/**
	 * As getCoordinate, but returns the saved class copy of the coordinate.
	 * 
	 * @return
	 */
	public double[] getSavedCoordinate() {
		return coordinate;
	}

	/**
	 * Returns the square distance between this feature and the other feature
	 * mapped to it.
	 */
	public double calculateSqrDist() {
		sqrDist = Coord.sqrDistance(coordinate, otherFeature.coordinate);
		return sqrDist;
	}

	/**
	 * Prints out information about any mapping of this feature to another
	 * feature
	 */
	@Deprecated
	public void printMapping() {
		printInfo();
		if (mapped) {
			System.out.print(" --> ");
			otherFeature.printInfo();
			System.out.println(" sqrDist " + sqrDist);
		} else {
			System.out.println(":Not Mapped");
		}
	}

	public String mappingInfo() {
		StringBuilder sb = new StringBuilder();
		sb.append(info());
		if (mapped) {
			sb.append(" --> ").append(otherFeature.info()).append(" sqrDist ")
					.append(sqrDist);
		} else {
			sb.append(":Not Mapped");
		}
		return sb.toString();

	}

	private static Pattern optionMatch = Pattern.compile("\\[.*\\]"),
			geometryMatch = Pattern.compile("(\\{.*\\}) N_MATCHED=\\d+");

	/**
	 * Recreates a feature from a pharmacophore description string.
	 * 
	 * @param str
	 * @param moleculeNo
	 * @return
	 * 
	 * @see PharmFeatureGeometry#parseString(String)
	 * @see AcceptorAtomFeature#AcceptorAtomFeature(String, HashMap, String,
	 *      PharmFeatureGeometry, int)
	 * @see DonorHydrogenFeature#DonorHydrogenFeature(String, HashMap, String,
	 *      PharmFeatureGeometry, int)
	 * @see AromaticRingFeature#AromaticRingFeature(HashMap, String,
	 *      PharmFeatureGeometry, int)
	 */
	public static Feature featureFromPharmDescription(String str, int moleculeNo) {
		String orig = str;
		HashMap<String, String> options = new HashMap<String, String>();
		String geometryStr, featureSetName, featureName = "", featureDefinition = "";

		logger.debug("Processing " + str);
		// Find geometry string
		Matcher m = geometryMatch.matcher(str);
		if (m.find()) {
			geometryStr = m.group(1);
			str = m.replaceFirst("");
			str = str.trim();
		} else {
			throw new RuntimeException("can't find geometry string in " + str);
		}

		// Find and remove options string if present
		m = optionMatch.matcher(str);
		if (m.find()) {
			String optStr = m.group(0);
			optStr = optStr.substring(1, optStr.length() - 1);
			String opts[] = optStr.split(",");
			for (String s : opts) {
				String vals[] = s.split(" = ");
				String k = vals[0], v = vals[1];
				options.put(k, v);
				logger.debug("adding option " + k + " = " + v);
			}
			str = m.replaceFirst("");
		}

		String vals[] = str.split(" ");
		featureSetName = vals[0];
		if (vals.length == 3) {
			featureName = vals[1];
			featureDefinition = vals[2];
		} else {
			featureDefinition = vals[1];
		}

		PharmFeatureGeometry geometry = PharmFeatureGeometry.parseString(geometryStr);

		logger.debug("Feature: Set Name " + featureSetName + " Name " + featureName
				+ " Definition " + featureDefinition + " Geometry " + geometryStr);

		if (featureSetName.equals("DONOR_HYDROGEN")) {
			return new DonorHydrogenFeature(featureName, options, featureDefinition,
					geometry, moleculeNo);
		} else if (featureSetName.equals("ACCEPTOR_ATOM")) {
			return new AcceptorAtomFeature(featureName, options, featureDefinition,
					geometry, moleculeNo);

		} else if (featureSetName.equals("AROMATIC_RING")) {
			return new AromaticRingFeature(options, featureDefinition, geometry,
					moleculeNo);
		} else {
			throw new RuntimeException("Can't load feature: " + orig);
		}

	}

	/**
	 * Returns the user feature type for the feature number. The feature number
	 * is indexed by 0 (so 1 will return USER_FEATURES2).
	 * 
	 * @param featureNo
	 * @return
	 */
	public static FeatureType userFeatureType(int featureNo) {
		return FeatureType.valueOf("USER_FEATURES" + String.valueOf(featureNo + 1));
	}

	/**
	 * Returns the index of a user feature- where USER_FEATURES1 has an index of
	 * 0.
	 * 
	 * @param featureType
	 * @return
	 */
	public static int userFeatureNo(FeatureType featureType) {
		return featureType.ordinal() - FeatureType.USER_FEATURES1.ordinal();
	}

	public boolean isMappingFeature() {
		return mappingFeature;
	}

	public static double getSolventVolOk() {
		return solventVolOk.get();
	}

	public static void setSolventVolOk(double solventVolOk) {
		Feature.solventVolOk.set(solventVolOk);
	}

	public boolean isMapped() {
		return mapped;
	}

	public void setMapped(boolean mapped) {
		this.mapped = mapped;
	}

	public Feature getOtherFeature() {
		return otherFeature;
	}

	public void setOtherFeature(Feature otherFeature) {
		this.otherFeature = otherFeature;
	}

	public boolean isPharmPoint() {
		return pharmPoint;
	}

	public void setPharmPoint(boolean pharmPoint) {
		this.pharmPoint = pharmPoint;
	}

	public void setMappingFeature(boolean mapping) {
		this.mappingFeature = mapping;
	}

	public void setPharmFeatureGeometry(PharmFeatureGeometry pharmFeatureGeometry) {
		this.pharmFeatureGeometry = pharmFeatureGeometry;
	}

	public FeatureType getFeatureType() {
		return featureType;
	}

	public Feature getBaseFeature() {
		return baseFeature;
	}

	public void setBaseFeature(Feature baseFeature) {
		this.baseFeature = baseFeature;
	}

	public double getBestGeometricScore() {
		return bestGeometricScore;
	}

	public void setBestGeometricScore(double bestGeometricScore) {
		this.bestGeometricScore = bestGeometricScore;
	}

	public double getBestScore() {
		return bestScore;
	}

	public void setBestScore(double bestScore) {
		this.bestScore = bestScore;
	}

	public double getGeometricScore() {
		return geometricScore;
	}

	public void setGeometricScore(double geometricScore) {
		this.geometricScore = geometricScore;
	}

	public GaMolecule getMolecule() {
		return molecule;
	}

	public double getSqrDist() {
		return sqrDist;
	}

	public void setSqrDist(double sqrDist) {
		this.sqrDist = sqrDist;
	}

	public static double getRadius() {
		return radius.get();
	}

	/**
	 * @return the nMatched
	 */
	public int getnMatched() {
		return nMatched;
	}

	/**
	 * @param nMatched
	 *            the nMatched to set
	 */
	public void setnMatched(int nMatched) {
		this.nMatched = nMatched;
	}

	/**
	 * @return the atom
	 */
	public Atom getAtom() {
		return atom;
	}

}
