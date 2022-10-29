package com.cairn.gape.feature;

import com.cairn.common.utils.Coord;
import com.cairn.gape.molecule.GaMolecule;
import com.cairn.molecule.Atom;

/**
 * Class for handling user-defined features
 */
public class UserFeature extends Feature {
	private volatile Atom[] featureAtoms;
	private volatile UserFeatureDefinition userFeatureType;
	private volatile UserFeatureSet userFeatureSet;

	static final boolean DEBUG = false;

	/**
	 * Constructor. The feature is an atom feature if the user feature is atom
	 * centered.
	 *
	 * @param m
	 *            Molecule Containing feature
	 * @param userFeatureType
	 *            Defining class for feature
	 * @param featureAtoms
	 *            All the atoms in the feature
	 */
	public UserFeature(GaMolecule m, UserFeatureDefinition userFeatureType,
			Atom featureAtoms[]) {
		this.userFeatureType = userFeatureType;
		this.featureAtoms = featureAtoms;
		atom = featureAtoms[0];
		userFeatureSet = userFeatureType.getUserFeatureSet();
		featureType = userFeatureSet.getFeatureType();
		featureSetName = userFeatureSet.getFeatureSetName().toUpperCase()
				.replace(' ', '_');
		molecule = m;
		if (userFeatureSet.isAtomCentered()) {
			atomFeature = true;
			atom = featureAtoms[0];
		}
		if (userFeatureSet.getFeatureWeight() <= 0)
			setMappingFeature(false);
	}

	/**
	 * The feature coordinate is first atom in the feature for atom centered
	 * features and the centroid of all feature atoms otherwise.
	 */
	@Override
	public double[] calculateCoordinate() {
		if (atomFeature)
			coordinate = super.calculateCoordinate();
		else {
			double x = 0, y = 0, z = 0;
			for (int i = 0; i < featureAtoms.length; i++) {
				double[] point = molecule.getCoord(featureAtoms[i].getNo());
				x += point[0];
				y += point[1];
				z += point[2];
			}
			coordinate[0] = x / featureAtoms.length;
			coordinate[1] = y / featureAtoms.length;
			coordinate[2] = z / featureAtoms.length;
			coordinate[3] = 1;
		}

		return coordinate;
	}

	/**
	 * Returns an informative string
	 */
	@Override
	public String info() {
		String rtn = userFeatureType.getName() + " [";
		for (int i = 0; i < featureAtoms.length; i++) {
			rtn += String.valueOf(featureAtoms[i].getNo() + 1);
			if (i != featureAtoms.length - 1)
				rtn += ",";
		}
		rtn += "]";
		return rtn;
	}

	/**
	 * Returns a string for identifying the pharmacophore feature
	 */
	@Override
	public String pharmLabel() {
		String label = featureSetName + " "
				+ userFeatureType.getName().toUpperCase().replace(' ', '_');
		if (userFeatureSet.isAtomCentered()) {
			label += " ATOM_" + String.valueOf(atom.getNo() + 1);
		} else {
			String alist = "";
			for (int i = 0; i < featureAtoms.length; i++) {
				alist += String.valueOf(featureAtoms[i].getNo() + 1);
				if (i < featureAtoms.length - 1)
					alist += "_";
			}
			label += " " + " ATOMS_" + alist;
		}
		return label;
	}

	/**
	 * Score this feature's overlap with another feature (in another molecule).
	 * Calculates the features gaussian overlap and scaled with sum of the
	 * feature weights.
	 *
	 * @param f
	 *            The other feature.
	 */
	@Override
	public double score(Feature f) {

		UserFeature other = (UserFeature) f;

		geometricScore = Feature.score(Coord.sqrDistance(coordinate, other.coordinate),
				userFeatureSet.getAlpha());
		double score = geometricScore
				* (userFeatureType.getWeight() + other.userFeatureType.getWeight());
		if (DEBUG) {
			System.out.println("Feature score " + info() + " " + other.info());
			System.out.println("Geom " + geometricScore + " score " + score + " dist "
					+ Coord.sqrDistance(coordinate, other.coordinate) + " wt "
					+ userFeatureType.getWeight() + " other wt "
					+ other.userFeatureType.getWeight());
		}

		return score;
	}

	/**
	 * main method for development and debugging
	 */
	public static void main(String args[]) {
		// GaMolecule molA = new GaMolecule();
		// GaMolecule molB = new GaMolecule();

		// try {
		// }
		// catch (GaException ex) {
		// ex.printStackTrace(System.err);
		// System.err.println("AcceptorAtomFeature: GaException: "+ex);
		// }

	}

	public String label() {
		int no = featureSetNo + 1;
		// int id = featureSet - Feature.USER_FEATURES1+1;
		return "USER_FEATURE_" + "_" + no;
	}

	@Override
	public PharmFeatureGeometry getPharmFeatureGeometry() {
		return new SpherePharmFeatureGeometry(coordinate, userFeatureSet.getRadius());
	}

	public Atom[] getFeatureAtoms() {
		return featureAtoms;
	}

	public UserFeatureSet getUserFeatureSet() {
		return userFeatureSet;
	}

	public UserFeatureDefinition getUserFeatureType() {
		return userFeatureType;
	}

}
