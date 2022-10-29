package com.cairn.molecule;

import java.util.ArrayList;
import java.util.List;

/**
 * A class that extends PatternMatch to search an SLN or Molpattern against a
 * Molecule. Generates an ArrayList of atoms, one atom for each match of the
 * pattern in the molecule. Each atom is the target atom which matches the first
 * atom of the query SLN.
 * 
 * @author Gareth Jones
 * @see com.cairn.molecule.MolPattern
 * @see com.cairn.molecule.PatternMatch
 */
public class AtomPatternMatch extends PatternMatch {
	private final ArrayList<Atom> matches;

	static final boolean DEBUG = false;

	/**
	 * Constructor, takes a MolPattern
	 * 
	 * @param p
	 */
	public AtomPatternMatch(MolPattern p) {
		super(p);
		matches = new ArrayList<Atom>();
	}

	/**
	 * Returns the list of target atom which match the first atom of the query.
	 * 
	 * @return
	 */
	public ArrayList<Atom> getMatches() {
		return matches;
	}

	/*
	 * Finds the first target atom in the match and adds it to the array of
	 * matching atoms.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.molecule.PatternMatch#process()
	 */
	@Override
	public void process() {
		int match[] = queryMatches[nMatches - 1];
		Atom tAtom = target.getAtom(match[0]);
		if (!matches.contains(tAtom))
			matches.add(tAtom);
		if (DEBUG) {
			System.out.println("Match " + nMatches);
			System.out.println("Match " + tAtom.info());
		}
	}

	/**
	 * Static method. Supply an sln and a molecule and you get an ArrayList of
	 * atoms in the molecule which match the first atom in the SLN.
	 * 
	 * @param sln
	 * @param mol
	 * @return
	 */
	public static ArrayList<Atom> matchSln(String sln, Molecule mol) {
		MolPattern p = MolPattern.generateMolPattern(sln);
		return matchPattern(p, mol);
	}

	/**
	 * Static method. Supply a Molpattern and a molecule and you get an
	 * ArrayList of atoms in the molecule which match the first atom in the
	 * pattern
	 * 
	 * @param p
	 * @param mol
	 * @return
	 */
	public synchronized static ArrayList<Atom> matchPattern(MolPattern p, Molecule mol) {
		AtomPatternMatch atomPatternMatch = new AtomPatternMatch(p);
		atomPatternMatch.match(mol);
		return atomPatternMatch.getMatches();
	}

	/*
	 * Makes sure matches list is clear before maching molecule.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.molecule.PatternMatch#match(com.cairn.molecule.Molecule
	 * )
	 */
	@Override
	public synchronized void match(Molecule t) {
		matches.clear();
		super.match(t);
	}

	/**
	 * Test method. Supply an sln and molecule file and list the matching atoms.
	 * 
	 * @param args
	 */
	public static void main(String[] args) {

		if (args.length != 2) {
			System.err.println("Usage: AtomPatternMatch <sln_query> <structure_file>");
			System.exit(0);
		}

		String sln = args[0];
		String molFile = args[1];

		MolPattern pattern = MolPattern.generateMolPattern(sln);
		List<Molecule> mols = Molecule.loadFiles(new String[] { molFile });
		for (Molecule mol : mols) {
			System.out.println("Molecule " + mol.getName());
			ArrayList<Atom> atomMatches = AtomPatternMatch.matchPattern(pattern, mol);
			for (Atom atom : atomMatches) {
				System.out.println(atom.info());
			}
		}
	}

}
