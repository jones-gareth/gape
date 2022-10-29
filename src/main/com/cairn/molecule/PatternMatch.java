package com.cairn.molecule;

import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.math3.util.Precision;
import org.apache.log4j.Logger;

import com.cairn.molecule.Ullman.MatchType;

/**
 * A class for matching a MolPattern object against a Molecule class.
 * 
 * @author Gareth Jones
 * @see com.cairn.molecule.MolPattern
 * @see com.cairn.molecule.Molecule
 * @see com.cairn.molecule.Ullman
 */
public class PatternMatch implements Ullman.UllmanCallback {
	static final int MAX_MATCHES = 100;

	// private volatile MolPattern pattern;

	protected volatile Molecule query, target;

	protected volatile int queryMatches[][];

	protected volatile int nMatches = 0;

	private volatile Ullman ull;

	private static final Logger logger = Logger.getLogger(PatternMatch.class);

	/**
	 * Empty constructor- to allow sub-class developmen
	 */
	protected PatternMatch() {
		;
	}

	/**
	 * Constructor. Creates a Query for matching a pattern against a molecule.
	 * 
	 * @param p
	 */
	public PatternMatch(MolPattern p) {
		setPattern(p);
	}

	/**
	 * Constructor from sln query
	 * 
	 * @param sln
	 */
	public PatternMatch(String sln) {
		MolPattern p = MolPattern.generateMolPattern(sln);
		setPattern(p);
	}

	/**
	 * Prepares the query pattern.
	 * 
	 * @param p
	 */
	private void setPattern(MolPattern p) {
		// pattern = p;
		if (p.getPatternMol() == null)
			p.patternToMolecule();
		query = p.getPatternMol();
		prepareQuery();
		queryMatches = new int[MAX_MATCHES][];
	}

	/**
	 * Sets up the query molecule. In particular we need to know cyclic atoms.
	 */
	private void prepareQuery() {
		query.assignAtomNeighbours();
		RingFinder rf = new RingFinder(query);
		rf.findRings();

		if (logger.isDebugEnabled()) {
			for (Atom atom : query.getAtoms()) {
				if (atom.isInRing()) {
					System.out.println("Cyclic Atom " + atom.info());
				}
			}
		}

	}

	/**
	 * Ullamn does basic atom-to-atom matching, but we need to check detailed
	 * bond environment, no of hydrogens, inclusion and exclusion atom lists and
	 * X and Y atoms. If all this passes then call the process routine
	 * 
	 * @param m
	 * @return
	 * 
	 * @see #process()
	 */
	private boolean checkMatch(Map<Atom, Atom> atomMapping) {

		logger.debug("checking match " + nMatches);

		if (queryMatches[nMatches] == null) {
			queryMatches[nMatches] = new int[query.getnAtoms()];
		}

		int match[] = queryMatches[nMatches];
		for (Entry<Atom, Atom> entry : atomMapping.entrySet()) {
			int queryNo = entry.getKey().getNo();
			int targetNo = entry.getValue().getNo();
			match[queryNo] = targetNo;
		}

		AtomType zType = null;
		AtomType yType = null;

		for (int i = 0; i < query.getnAtoms(); i++) {
			// Check hydrogen count matches
			QueryAtom queryAtom = (QueryAtom) query.getAtom(i);
			Atom targetAtom = target.getAtom(match[i]);

			if (queryAtom.getAtomType() == AtomType.Type.Z) {
				if (zType != null) {
					if (!zType.matchElementalType(targetAtom.getType()))
						return false;
				} else
					zType = targetAtom.getType();
			}

			if (queryAtom.getAtomType() == AtomType.Type.Y) {
				if (yType != null) {
					if (!yType.matchElementalType(targetAtom.getType()))
						return false;
				} else
					yType = targetAtom.getType();
				if (zType != null && zType.matchElementalType(yType))
					return false;
			}

			if (queryAtom.getnSlnHydrogens() != null) {
				if (targetAtom.getHydrogenCount() < queryAtom.getnSlnHydrogens()) {
					return false;
				}
			}

			// check special flags: fully connected and no hydrogens
			if (queryAtom.isNoHydrogens() && targetAtom.getHydrogenCount() > 0)
				return false;
			if (queryAtom.isFullyConnected() || queryAtom.getnConnected() > 0) {
				int nQuery = 0;
				if (queryAtom.isFullyConnected()) {
					nQuery = queryAtom.getnNotDummyNeighbours();
					if (queryAtom.getnSlnHydrogens() != null) {
						nQuery += queryAtom.getnSlnHydrogens();
					}
				} else
					nQuery = queryAtom.getnConnected();
				logger.debug("Connection test nQuery " + nQuery + " nTarget "
						+ targetAtom.getnNotDummyNeighbours());
				if (nQuery != targetAtom.getnNotDummyNeighbours())
					return false;
			}

			if (CollectionUtils.isNotEmpty(queryAtom.getIs())) {
				boolean ok = false;
				for (AtomType isType : queryAtom.getIs()) {
					logger.debug("Checking is types: ");
					if (isType.matchType(targetAtom.getType())) {
						logger.debug("matching type " + isType.getName());
						ok = true;
						break;
					}
				}
				if (!ok) {
					logger.debug("match failed");
					return false;
				}
			}
			if (CollectionUtils.isNotEmpty(queryAtom.getNot())) {
				for (AtomType notType : queryAtom.getNot()) {
					logger.debug("Checking not types: ");
					if (notType.matchType(targetAtom.getType())) {
						logger.debug("matching type " + notType.getName() + " returning");
						return false;
					}
				}
			}

			if (queryAtom.getQueryInRing() != null) {
				if (!queryAtom.getQueryInRing() && targetAtom.isInRing()) {
					return false;
				}
				if (queryAtom.getQueryInRing() && !targetAtom.isInRing()) {
					return false;
				}
			}

			if (queryAtom.getFormalCharge() != null) {
				if (targetAtom.getFormalCharge() == null) {
					return false;
				}
				if (targetAtom.getFormalCharge() != queryAtom.getFormalCharge()) {
					return false;
				}
			}

			if (queryAtom.getPartialCharge() != 0) {
				if (!Precision.equals(queryAtom.getPartialCharge(),
						targetAtom.getPartialCharge(), 0.1)) {
					return false;
				}
			}
		}

		// check bond environment
		for (Bond bond : query.getBonds()) {
			QueryBond queryBond = (QueryBond) bond;
			Atom targetAtom1 = target.getAtom(match[queryBond.getAtom1().getNo()]);
			Atom targetAtom2 = target.getAtom(match[queryBond.getAtom2().getNo()]);
			Bond targetBond = target.getBond(targetAtom1, targetAtom2);
			if (targetBond == null)
				return false;
			if (!queryBond.matchBond(targetBond))
				return false;
		}

		logger.debug("got match");
		nMatches++;
		process();
		return true;
	}

	/**
	 * This is called when we get a match- subclass to do something useful. It
	 * just prints out the match.
	 */
	public void process() {
		int match[] = queryMatches[nMatches - 1];
		logger.info("Match " + nMatches);
		for (int i = 0; i < query.getnAtoms(); i++) {
			Atom qAtom = query.getAtom(i);
			Atom tAtom = target.getAtom(match[i]);
			System.out.println(qAtom.info() + " --> " + tAtom.info());
		}
	}

	/*
	 * Called when the Ullman algorithm finds a match. Call checkMatch, which
	 * will then call process if the match is ok.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.molecule.Ullman.UllmanCallback#callback(boolean[][])
	 */
	@Override
	public void callback(Map<Atom, Atom> atomMapping) {
		checkMatch(atomMapping);
	}

	/**
	 * Call this to match the pattern to the molecule.
	 * 
	 * @param t
	 */
	public synchronized void match(Molecule t) {
		nMatches = 0;
		target = t;
		doUllman();
	}

	/**
	 * Call by match to do the work
	 * 
	 */
	private synchronized void doUllman() {
		ull = new Ullman(query, target, MatchType.PATTERN);

		ull.doUllman(this);
	}

	/**
	 * Test routine- mathces an sln to a mol2 file.
	 * 
	 * @param args
	 */
	public static void main(String args[]) {
		if (args.length != 2) {
			System.err.println("Usage: PatternMatch <sln_query> <mol2_file>");
			System.exit(0);
		}

		MolPattern p = MolPattern.generateMolPattern(args[0]);
		List<Molecule> molecules = Molecule.loadFiles(new String[] { args[1] });

		PatternMatch pM = new PatternMatch(p);
		for (Molecule mol : molecules) {
			System.out.println("Searching " + mol.getName());
			pM.match(mol);
		}
	}
}
