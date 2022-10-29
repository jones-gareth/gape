package com.cairn.molecule;

import org.apache.log4j.Logger;

import com.cairn.molecule.Atom.AtomFragmentType;

/**
 * Searches molecules for fragments that may be used in defining torsional
 * distributions.
 * 
 * @author Gareth Jones
 * 
 */
public class Fragments {
	Molecule molecule;

	// Array of common fragments- fragment ids are defined in the Atom class
	static Fragment[] fragments;

	static {
		fragments = new Fragment[] {
				new Fragment(AtomFragmentType.ATM_BENZENE, "Benzene", "C[1]:C:C:C:C:C:@1"),
				new Fragment(AtomFragmentType.ATM_RIBOSE, "Ribose", "O[1]-C-C-C-C-@1"),
				new Fragment(AtomFragmentType.ATM_ADENINE, "Adenine",
						"C[1]:N:C[3]:C:N:C:N:C:@3:NH:@1"),
				new Fragment(AtomFragmentType.ATM_URACIL, "Uracil",
						"C[1]-C-NH-C(=O)-NH-C=@1"),
				new Fragment(AtomFragmentType.ATM_CYTOSINE, "Cytosine",
						"C[1]-C=N-C(=O)-NH=@1") };
	}

	/**
	 * This class just provides static methods
	 */
	private Fragments() {
	}

	/**
	 * Searches a molecule for fragments and marks all atoms in any fragments
	 * found
	 * 
	 * @param molecule
	 */
	public static void searchMolecule(Molecule molecule) {
		for (int i = 0; i < fragments.length; i++) {
			fragments[i].match(molecule);
		}
	}

	/**
	 * @param frag
	 * @return a fragment name given a fragment id
	 */
	public static String idToFragName(AtomFragmentType frag) {
		if (frag == AtomFragmentType.ATM_NONE)
			return "none";
		for (int i = 0; i < fragments.length; i++) {
			if (frag == fragments[i].fragType)
				return fragments[i].name;
		}
		System.out.println("idToFragName: unknown ID" + frag);
		return null;
	}

	/**
	 * @param name
	 * @return a fragment id given a fragment name
	 */
	public static AtomFragmentType fragNameToId(String name) {
		for (int i = 0; i < fragments.length; i++) {
			if (fragments[i].name.equalsIgnoreCase(name))
				return fragments[i].fragType;
		}
		System.out.println("fragNameToID: unknown Fragment" + name);
		return AtomFragmentType.ATM_NONE;
	}

	/**
	 * Class to represent a single fragment- includes matching functions.
	 * 
	 * @author Gareth Jones
	 * @see com.cairn.molecule.PatternMatch
	 */
	private static class Fragment extends PatternMatch {

		// private final String sln;

		private final String name;

		private final AtomFragmentType fragType;

		private static final Logger logger = Logger.getLogger(Fragment.class);

		/**
		 * Constructor
		 * 
		 * @param _fragType
		 * @param _name
		 * @param _sln
		 */
		Fragment(AtomFragmentType _fragType, String _name, String _sln) {
			super(_sln);
			// sln = _sln;
			name = _name;
			fragType = _fragType;
		}

		/*
		 * Records a match between this fragment and a molecule. Each matching
		 * atom in the molecule has a flag set indicating that it matches this
		 * fragment. (non-Javadoc)
		 * 
		 * @see com.cairn.molecule.PatternMatch#process()
		 */
		@Override
		public void process() {
			int match[] = queryMatches[nMatches - 1];
			logger.debug("Match " + name);
			for (int i = 0; i < query.getnAtoms(); i++) {
				Atom qAtom = query.getAtom(i);
				Atom tAtom = target.getAtom(match[i]);
				tAtom.addFragment(fragType);
				logger.debug(qAtom.info() + " --> " + tAtom.info());
			}
		}
	}
}
