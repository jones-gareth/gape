package com.cairn.molecule;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * This class finds the set of unconnected structures in a molecule structure.
 * Contains methods to generate an array of fully connected molecules.
 * 
 * @author Gareth Jones
 * 
 */
public class SplitMolecule {
	private Molecule molecule;

	private int moleculeNos[], nMolecules;

	private boolean inMolecule[];

	HashMap<Integer, Integer> fragmentSizes;

	/**
	 * Simple constructor
	 */
	public SplitMolecule() {
		;
	}

	/**
	 * Set current molecule
	 * 
	 * @param m
	 */
	public void setMolecule(Molecule m) {
		molecule = m;
	}

	/**
	 * Finds all the fully connected fragments and their sizes in the current
	 * molecule.
	 */
	void findFragments() {
		int nAtoms = molecule.getnAtoms();
		inMolecule = new boolean[nAtoms];
		moleculeNos = new int[nAtoms];
		fragmentSizes = new HashMap<Integer, Integer>();

		boolean moleculeFound = true;
		int moleculeNo = 0;
		while (moleculeFound) {
			moleculeFound = false;
			for (int i = 0; i < nAtoms; i++) {
				if (!inMolecule[i]) {
					moleculeFound = true;
					findFragment(i, moleculeNo);
					moleculeNo++;
				}
			}
		}

		nMolecules = moleculeNo;
		for (int i = 0; i < nMolecules; i++) {
			int size = fragmentSize(i);
			System.out.println("Fragment " + String.valueOf(i + 1) + " size " + size);
		}
	}

	/**
	 * Determines the size of a particular fragment
	 * 
	 * @param moleculeNo
	 * @return
	 */
	private int fragmentSize(int moleculeNo) {
		int nAtoms = molecule.getnAtoms();
		int count = 0;
		for (int i = 0; i < nAtoms; i++)
			if (inMolecule[i] && moleculeNos[i] == moleculeNo)
				count++;
		fragmentSizes.put(moleculeNo, count);
		return count;
	}

	/**
	 * Recursive method to find fully connected fragments. Marks this atom as
	 * being in the current fragment then marks all unmarked neighbours.
	 * 
	 * @param atomNo
	 * @param moleculeNo
	 */
	private void findFragment(int atomNo, int moleculeNo) {
		Atom atom = molecule.getAtom(atomNo);
		moleculeNos[atomNo] = moleculeNo;
		inMolecule[atomNo] = true;

		for (Atom neighbour : atom.getNotDummyNeighbours()) {
			int neighbourNo = neighbour.getNo();
			if (!inMolecule[neighbourNo])
				findFragment(neighbourNo, moleculeNo);
		}

	}

	/**
	 * Converts all the fully connected fragments found to an array of
	 * Molecules.
	 * 
	 * @return
	 */
	public Molecule[] getMoleculeFragments() {
		Molecule molecules[] = new Molecule[nMolecules];
		for (int i = 0; i < nMolecules; i++)
			molecules[i] = fragmentToMolecule(i);
		return molecules;
	}

	/**
	 * Converts a fully connected fragment into a Molecule class.
	 * 
	 * @param moleculeNo
	 * @return
	 */
	private Molecule fragmentToMolecule(int moleculeNo) {
		int nFragmentAtoms = fragmentSizes.get(moleculeNo);
		List<Atom> newAtoms = new ArrayList<>();
		List<double[]> newCoords = new ArrayList<>();
		HashMap<Atom, Atom> atomMap = new HashMap<Atom, Atom>();

		int atomNo = 0;
		int nLonePairs = 0;
		int nMolAtoms = molecule.getnAtoms();
		for (int i = 0; i < nMolAtoms; i++) {
			if (inMolecule[i] && moleculeNos[i] == moleculeNo) {
				Atom oldAtom = molecule.getAtom(i);
				Atom newAtom = oldAtom.copy(null);
				newAtoms.add(newAtom);
				newCoords.add(molecule.getCoord(i));
				newAtom.setNo(atomNo);
				atomMap.put(newAtom, oldAtom);
				atomNo++;
				if (oldAtom.getAtomType() == AtomType.Type.LP)
					nLonePairs++;
			}
		}

		ArrayList<Bond> bonds = new ArrayList<Bond>();
		int bondNo = 0;
		for (int i = 0; i < nFragmentAtoms; i++) {
			Atom atomA = newAtoms.get(i);
			for (int j = i + 1; j < nFragmentAtoms; j++) {
				Atom atomB = newAtoms.get(j);
				Atom oldAtomA = atomMap.get(atomA);
				Atom oldAtomB = atomMap.get(atomB);
				Bond bond = molecule.getBond(oldAtomA, oldAtomB);
				if (bond == null)
					continue;
				Bond newBond = new Bond(bondNo, atomA, atomB, bond.getBondType());
				bondNo++;
				bonds.add(newBond);
			}
		}

		String name = molecule.getName() + "_fragment_" + String.valueOf(moleculeNo + 1);
		Molecule mol = new Molecule(name, newAtoms, bonds, newCoords);
		mol.setnLonePairs(nLonePairs);

		System.out.println("Created molecule for fragment "
				+ String.valueOf(moleculeNo + 1));

		return mol;
	}

}
