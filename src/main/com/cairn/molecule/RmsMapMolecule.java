package com.cairn.molecule;

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;
import com.cairn.molecule.Ullman.MatchType;

/**
 * This inner class does the work of finding all isomorphisms of a molecule to
 * itself and selecting the one mapping with the best rmsd
 * 
 */
class RmsMapMolecule implements Ullman.UllmanCallback {

	private static final Logger logger = Logger.getLogger(RmsMapMolecule.class);

	private double xVals[][], yVals[][], trans[][];

	private Molecule molecule, template;

	private HashMap<Integer, Integer> mapping, bestMapping;

	private double bestRMS;

	private int nIsomorphisms;

	private final boolean matchElementalTypes, subgraph;

	public RmsMapMolecule(boolean matchElementalTypes, boolean subgraph) {
		super();
		this.matchElementalTypes = matchElementalTypes;
		this.subgraph = subgraph;
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
	public void callback(Map<Atom, Atom> mapping) {
		searchIsomorphism(mapping);
	}

	/**
	 * Test an isomorphism
	 * 
	 * @param m
	 */
	private void searchIsomorphism(Map<Atom, Atom> atomMapping) {

		int length = atomMapping.size();

		if (mapping == null) {
			mapping = new HashMap<Integer, Integer>();
			bestRMS = Double.MAX_VALUE;
			xVals = new double[length][4];
			yVals = new double[length][];
			trans = new double[4][4];
		}

		nIsomorphisms++;
		if (nIsomorphisms % 100 == 0)
			logger.info("Searching isomorphism " + nIsomorphisms);

		mapping.clear();

		// template is query in Ullamn atomsA, mapping key - molecule is
		// target in atomsB, mapping value
		atomMapping.entrySet().stream().forEach(entry -> {
			int no2 = entry.getKey().getNo();
			int no1 = entry.getValue().getNo();
			logger.debug("Atom " + no1 + "-->" + no2);
			mapping.put(no1, no2);
		});

		int i = 0;
		for (Entry<Integer, Integer> entry : mapping.entrySet()) {
			int no1 = entry.getKey();
			int no2 = entry.getValue();
			Coord.copy(molecule.getCoord(no1), xVals[i]);
			yVals[i] = template.getCoord(no2);
			i++;
		}

		Coord.leastSquaresFit(xVals, yVals, trans);

		double msSum = 0;
		for (i = 0; i < length; i++) {
			Coord.transPointInPlace(trans, xVals[i]);
			msSum += Coord.sqrDistance(xVals[i], yVals[i]);
		}
		double rms = Math.sqrt(msSum / (length));
		if (rms < bestRMS) {
			bestRMS = rms;
			bestMapping = mapping;
			mapping = new HashMap<Integer, Integer>();
		}
		logger.debug("Isomorphism " + nIsomorphisms + " rms " + rms);

	}

	public int getNIsomorphisms() {
		return nIsomorphisms;
	}

	public double getBestRms() {
		return bestRMS;
	}

	/**
	 * Search through all isomorphisms and returns the one with the smallest
	 * RMSD
	 * 
	 */
	public HashMap<Integer, Integer> searchIsomorphisms(Molecule molecule,
			Molecule template) {
		this.molecule = molecule;
		this.template = template;

		// checkSame(molecule, template);
		Ullman ullman = new Ullman(template, molecule,
				subgraph ? MatchType.FRAG_HEAVY : MatchType.MATCH);
		// ignore aromatic / sp2 distinctions.
		ullman.setIgnoreAromatic(true);
		// ignore hydrogens
		// ullman.setIgnoreHydrogensInMatch(true);
		ullman.setMatchElementalTypes(matchElementalTypes);

		ullman.doUllman(this);

		return bestMapping;
	}

}