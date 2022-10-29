package com.cairn.gape.feature;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import com.cairn.gape.molecule.GaMolecule;
import com.cairn.molecule.Atom;
import com.cairn.molecule.MolPattern;
import com.cairn.molecule.PatternMatch;

/**
 * Class for defining user-defined feature types and finding them in input
 * structures.
 */
public class UserFeatureDefinition {
	private final UserFeatureSet userFeatureSet;

	private final String sln, name;

	private final double weight;

	private final int matchWeight;

	private final boolean atomCentered;

	private final MolPattern pattern;

	private volatile UserFeaturePatternMatch match;

	private static final boolean DEBUG = false;

	// these two may need to be thread local
	private volatile ArrayList<UserFeature> featureMatches;

	private volatile GaMolecule currentMolecule;

	/**
	 * Constructor
	 * 
	 * @param userFeatureSet
	 *            Equivalent feature set
	 * @param name
	 *            Feature type name
	 * @param weight
	 *            Feature weight
	 * @param sln
	 *            Pattern for feature
	 */
	public UserFeatureDefinition(UserFeatureSet userFeatureSet, String name,
			double weight, String sln) {

		this.userFeatureSet = userFeatureSet;
		this.name = name;
		this.weight = weight;
		this.sln = sln;
		atomCentered = userFeatureSet.isAtomCentered();

		pattern = MolPattern.generateMolPattern(sln);

		match = new UserFeaturePatternMatch(pattern);
		matchWeight = pattern.getPatternMol().getnAtoms();
	}

	/**
	 * Returns an informational string about the feature type
	 */
	public String info() {
		return "Feature " + name + " " + sln + " [wt:" + weight + "]";
	}

	/**
	 * Seaches molecule mol for features in this set. Matches for atom centered
	 * features are added to the Atom structure- matches that are more specific
	 * (have more atoms) replace less specific matches. If the feature is not
	 * atom centered, the list for matched features is kept in featureMatches
	 * and returned.
	 * 
	 * @param mol
	 *            Molecule to be searched
	 */
	public List<UserFeature> findUserFeatureTypes(GaMolecule mol) {
		if (!atomCentered)
			featureMatches = new ArrayList<UserFeature>();
		else
			featureMatches = null;
		currentMolecule = mol;
		match.match(mol);
		return featureMatches;
	}

	/**
	 * Checks feature match to see if it's already been found. This occurs if
	 * the fragment has symmetry (eg benzene can match six ways). If it's not
	 * present featureMatches is updated.
	 * 
	 * @param atomlist
	 *            Array of matching atoms.
	 */
	private void addMatch(Atom atomList[]) {

		// sort list by atom id so that the atom list is now unique.

		if (DEBUG)
			System.out.println("Found Match " + formatAtomList(atomList));
		Arrays.sort(atomList, new Comparator<Atom>() {
			@Override
			public int compare(Atom atom1, Atom atom2) {
				if (atom1.getNo() > atom2.getNo())
					return +1;
				if (atom1.getNo() < atom2.getNo())
					return -1;
				return 0;
			}

		});
		if (DEBUG)
			System.out.println("Sorted match " + formatAtomList(atomList));

		for (int i = 0; i < featureMatches.size(); i++) {
			Atom otherList[] = featureMatches.get(i).getFeatureAtoms();
			boolean match = true;
			for (int j = 0; j < atomList.length; j++) {
				if (atomList[j] != otherList[j]) {
					match = false;
					break;
				}
			}
			if (match) {
				if (DEBUG)
					System.out.println("Already Present");
				return;
			}
		}

		if (DEBUG)
			System.out.println("New Match");

		featureMatches.add(new UserFeature(currentMolecule, UserFeatureDefinition.this,
				atomList));
	}

	/**
	 * Returns a String of atom list numbers
	 */
	public String formatAtomList(Atom[] list) {
		String rtn = "[";
		for (int i = 0; i < list.length; i++) {
			int no = list[i].getNo() + 1;
			rtn += String.valueOf(no);
			if (i < list.length - 1)
				rtn += ",";
		}
		rtn += "]";
		return rtn;
	}

	/**
	 * Class for matching feature sln to structure
	 */
	private class UserFeaturePatternMatch extends PatternMatch {
		UserFeaturePatternMatch(MolPattern p) {
			super(p);
		}

		/**
		 * Found a match. The match is in queryMatches array
		 */
		@Override
		public void process() {
			int match[] = queryMatches[nMatches - 1];
			int nQueryAtoms = query.getnAtoms();
			Atom aMatch[] = new Atom[nQueryAtoms];
			for (int i = 0; i < nQueryAtoms; i++) {
				aMatch[i] = target.getAtom(match[i]);
			}
			if (atomCentered) {
				Atom featureAtom = target.getAtom(match[0]);
				if (featureAtom.userFeature == null
						|| featureAtom.userFeature.getUserFeatureType().matchWeight < matchWeight) {
					featureAtom.userFeature = new UserFeature(currentMolecule,
							UserFeatureDefinition.this, aMatch);
				}
			} else {
				addMatch(aMatch);
			}
		}
	}

	public String getName() {
		return name;
	}

	public String getSln() {
		return sln;
	}

	public UserFeatureSet getUserFeatureSet() {
		return userFeatureSet;
	}

	public double getWeight() {
		return weight;
	}

}