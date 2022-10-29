package com.cairn.molecule;

import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.math3.util.FastMath;
import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;
import com.cairn.molecule.Ullman.MatchType;
import com.cairn.molecule.Ullman.UllmanCallback;

/**
 * Determines the RMS between two molecules. Accounts for all isomorphisms (all
 * ways of mapping one structure to another) and returns the best rms.
 * 
 * @author Gareth Jones
 *
 */
public class Rms implements UllmanCallback {
	private static final Logger logger = Logger.getLogger(Rms.class);

	static {
		// logger.setLevel(Level.DEBUG);
	}

	private final boolean subgraph, matchElementalTypes;

	private Molecule moleculeA, moleculeB;
	private double rms;
	private boolean matched;

	public Rms(boolean subgraph, boolean matchElementalTypes) {
		super();
		this.subgraph = subgraph;
		this.matchElementalTypes = matchElementalTypes;
	}

	/**
	 * Determine the rms between two molecules
	 * 
	 * @param moleculeA
	 * @param moleculeB
	 * @return
	 */
	public synchronized double determineRms(Molecule moleculeA, Molecule moleculeB) {
		this.moleculeA = moleculeA;
		this.moleculeB = moleculeB;

		Ullman ullman = new Ullman(moleculeA, moleculeB,
				subgraph ? MatchType.FRAG_HEAVY : MatchType.MATCH);
		ullman.setMatchElementalTypes(matchElementalTypes);
		rms = Double.MAX_VALUE;
		ullman.doUllman(this);

		if (!matched) {
			throw new RuntimeException("Unable to map molecule " + moleculeA.getName()
					+ " to " + moleculeB.getName());
		}
		return rms;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.molecule.Ullman.UllmanCallback#callback(java.util.Map)
	 */
	@Override
	public void callback(Map<Atom, Atom> atomMapping) {

		// got an isomporphism, determine rms.

		matched = true;
		double sqrDistanceSum = .0;
		for (Entry<Atom, Atom> entry : atomMapping.entrySet()) {
			double[] coordA = moleculeA.getCoord(entry.getKey().getNo());
			double[] coordB = moleculeB.getCoord(entry.getValue().getNo());
			double sqrDistance = Coord.sqrDistance(coordA, coordB);
			sqrDistanceSum += sqrDistance;
		}

		double isoRms = FastMath.sqrt(sqrDistanceSum / atomMapping.size());
		logger.debug("Evaluated isomorphism, rms = " + isoRms);
		if (isoRms < rms) {
			rms = isoRms;
		}

	}

	/**
	 * Main method: takes molecule files and determines rms between first two
	 * structures.
	 * 
	 * @param args
	 */
	public static void main(String[] args) {

		boolean init = false, subgraph = false, matchElementalTypes = true;
		if (args.length < 1) {
			System.err.println("Usage: " + Rms.class.getName() + " <molecule files..>");
			System.exit(0);
		}

		List<Molecule> molecules = Molecule.loadFiles(args, init);
		if (molecules.size() > 2) {
			logger.warn(molecules.size()
					+ " molecules present.  Determining rms between the first two.");
		}
		if (molecules.size() < 2) {
			logger.error("Only " + molecules.size()
					+ " molecules present.  Unable to determine rms between 2 molecules");
			System.exit(0);
		}

		Rms rmsFitter = new Rms(subgraph, matchElementalTypes);
		double rms = rmsFitter.determineRms(molecules.get(0), molecules.get(1));
		System.out.println("Rms is " + rms);
	}
}
