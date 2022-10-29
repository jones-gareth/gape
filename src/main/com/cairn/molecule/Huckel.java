package com.cairn.molecule;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.log4j.Logger;

/**
 * Checks to see if a ring system (set of connected sp2 rings) is aromatic.
 * 
 * according to Huckels rule
 * 
 * @author Gareth Jones
 * 
 */
public class Huckel {
	private static final Logger logger = Logger.getLogger(Huckel.class);
	static {
		// logger.setLevel(Level.DEBUG);
	}

	private volatile boolean aromatic;
	private final Set<Atom> atoms;
	private final Molecule molecule;

	/**
	 * Determines all sp ring systems in a molecule and determines aromaticity
	 * according to Huckels rule. setSp2() must have been called for all rings
	 * in the molecule previously.
	 * 
	 * @param molecule
	 * @return List of all sp2 ring systems
	 */
	public static List<Huckel> huckel(Molecule molecule) {
		List<Huckel> ringSystems = new ArrayList<Huckel>();

		boolean[] used = new boolean[molecule.getnRings()];
		for (int i = 0; i < molecule.getnRings(); i++) {
			if (used[i])
				continue;
			Ring ring = molecule.getRings().get(i);
			if (!ring.isAllSp2())
				continue;
			List<Ring> rings = new ArrayList<Ring>();
			rings.add(ring);
			used[i] = true;

			boolean changed = false;
			do {
				changed = false;
				for (int j = i + 1; j < molecule.getnRings(); j++) {
					if (used[j])
						continue;
					Ring other = molecule.getRings().get(j);
					if (!other.isAllSp2())
						continue;
					if (inRingSystem(other, rings)) {
						rings.add(other);
						used[j] = true;
						changed = true;
					}
				}
			} while (changed);

			Huckel ringSystem = new Huckel(molecule, rings);
			ringSystem.checkAromatic();
			ringSystems.add(ringSystem);

		}

		return ringSystems;
	}

	/**
	 * @param ring
	 * @param others
	 * @return true if ring is bonded to any ring in others.
	 */
	private static boolean inRingSystem(Ring ring, List<Ring> others) {

		for (Ring other : others) {

			for (Bond bond1 : ring.getBonds()) {
				for (Bond bond2 : other.getBonds()) {
					if (bond1 == bond2) {
						return true;
					}
				}
			}
		}

		return false;

	}

	/**
	 * Ring system constructor
	 * 
	 * @param molecule
	 * @param rings
	 */
	private Huckel(Molecule molecule, List<Ring> rings) {
		this.molecule = molecule;
		atoms = new HashSet<Atom>();

		for (Ring ring : rings) {
			for (Atom atom : ring.getAtoms()) {
				if (!atoms.contains(atom))
					atoms.add(atom);
			}
		}
	}

	/**
	 * Checks to see that the ring system is aromatic or not. Marks all rings as
	 * aromatic if it is.
	 * 
	 * @return true if this ring system is aromatic
	 */
	private boolean checkAromatic() {

		boolean aromaticOk = true;
		if (atoms.size() < 5)
			aromaticOk = false;
		logger.debug("checking ring system for aromaticity.");

		int nElectrons = 0;

		for (Atom atom : atoms) {
			if (!aromaticOk) {
				logger.debug("not aromatic");
				break;
			}

			if (atom.getType().isPhosphorousType())
				aromaticOk = false;
			else if (atom.getType().isOxygenType())
				nElectrons += 2;
			else if (atom.getType().isSulphurType()) {
				if (sp2Ok(atom))
					nElectrons += 2;
				else
					aromaticOk = false;
			}

			else if (atom.getType().isCarbonType()) {
				if (sp2Ok(atom))
					nElectrons += 1;
				// -1 carbon, with two single bonds can contribute 1 electron
				// (not sure about this one if hydrogens are added automatically

				// else if (atom.nNotDummyNeighbours == 2 &&
				// atom.nSingleNeighbours == 2)
				// nElectrons += 1;

				else {

					aromaticOk = false;
				}
			}

			else if (atom.getType().isNitrogenType()) {
				if (sp2Ok(atom))
					nElectrons += 1;
				else if (atom.getnNotDummyNeighbours() == 3
						&& atom.getnSingleNeighbours() == 3)
					nElectrons += 2;
				else {
					aromaticOk = false;
				}
			}

			else {
				logger.debug("atom " + atom.info() + " is not aromatic");
				aromaticOk = false;
			}
		}

		if (aromaticOk) {
			if (nElectrons == 6 || nElectrons == 10 || nElectrons == 14
					|| nElectrons == 18)
				aromaticOk = true;
			else
				aromaticOk = false;
		}

		if (aromaticOk) {
			// 4n + 2 rule
			int test = nElectrons - 2;
			if (test % 4 == 0)

				if (nElectrons == 6 || nElectrons == 10 || nElectrons == 14
						|| nElectrons == 18) {
					logger.debug("Ring system passes huckels test");
					aromaticOk = true;
				} else {
					logger.debug("Ring system fails huckels test");
					aromaticOk = false;
				}
		}

		if (aromaticOk) {
			for (Ring ring : molecule.getRings()) {
				if (containsRing(ring)) {
					logger.debug("setting ring " + ring.info() + " aromatic");
					ring.setAromatic(true);
				}
			}
		}
		aromatic = aromaticOk;
		return aromaticOk;
	}

	/**
	 * @param atom
	 * @return true if this is a typical sp2 ring atom
	 */
	private boolean sp2Ok(Atom atom) {
		int nTotal = 0, nDouble = 0, nSingle = 0, nAromatic = 0;

		for (Atom other : atom.getNotDummyNeighbours()) {
			if (atoms.contains(other)) {
				Bond bond = molecule.getBond(atom, other);
				logger.debug("Got bond " + bond.info());
				if (bond.getBondType() == BondType.Type.SINGLE)
					nSingle++;
				else if (bond.getBondType() == BondType.Type.DOUBLE)
					nDouble++;
				else if (bond.getBondType() == BondType.Type.AR)
					nAromatic++;
				nTotal++;
			}
		}

		if (nTotal == 2) {
			if (nSingle == 1 && nDouble == 1)
				return true;
			if (nAromatic == 2)
				return true;
		}
		if (nTotal == 3) {
			if (nSingle == 2 && nDouble == 1)
				return true;
			if (nAromatic == 3)
				return true;
		}
		logger.debug("Atom " + atom.info() + " is not aromatic");
		return false;

	}

	/**
	 * This could be simpler now the class is constructed from a list of rings.
	 * 
	 * @param ring
	 * @return true if this ring is part of the ring system.
	 */
	private boolean containsRing(Ring ring) {

		for (Atom atom : ring.getAtoms()) {
			if (!atoms.contains(atom))
				return false;
		}
		return true;
	}

	/**
	 * @return the aromatic
	 */
	public boolean isAromatic() {
		return aromatic;
	}

	/**
	 * @return the atoms
	 */
	public Set<Atom> getAtoms() {
		return atoms;
	}

	/**
	 * @return the molecule
	 */
	public Molecule getMolecule() {
		return molecule;
	}

}
