package com.cairn.gape.molecule;

import com.cairn.gape.feature.LonePairAddition;
import com.cairn.molecule.Atom;
import com.cairn.molecule.AtomType;
import com.cairn.molecule.Bond;
import com.cairn.molecule.BondType;
import com.cairn.molecule.Solvate;

/**
 * Application to protonate an atom in a Molecule. The Solvate class will be
 * able to do this for most obvious cases.
 * 
 * @author Gareth Jones
 * @see Solvate
 * 
 */
class Protonation {

	public static void main(String args[]) {
		GaMolecule mol = new GaMolecule();

		if (args.length != 3) {
			System.out.println("Usage: Protonation <mol2file> <atom_no> <outfile>");
			System.exit(0);
		}

		mol.loadFile(args[0]);
		int atomNo = Integer.valueOf(args[1]).intValue();
		atomNo--;
		Atom atom = mol.getAtom(atomNo);
		System.out.println("Trying to add proton to " + atom.info());

		// Hijack the Lone Pair addition code to add the proton.
		double hBondLength = 1.008;
		int nConn = atom.getnNotDummyNeighbours();
		AtomType.Geometry geometry = atom.getType().getGeometry();
		double protonCoord[] = new double[4];
		double coord[] = mol.getCoord(atom.getNo());

		if (geometry == AtomType.Geometry.TRI) {
			if (nConn == 2) {
				double[] atom2 = mol
						.getCoord(atom.getNotDummyNeighbours().get(0).getNo());
				double[] atom3 = mol
						.getCoord(atom.getNotDummyNeighbours().get(1).getNo());
				LonePairAddition.addOnePairToTrigonal(coord, atom2, atom3, protonCoord,
						hBondLength);
			} else
				throw new RuntimeException("Unable to protonate Trigonal Atom");
		} else if (geometry == AtomType.Geometry.TET) {
			if (nConn == 3) {
				double[] atom2 = mol
						.getCoord(atom.getNotDummyNeighbours().get(0).getNo());
				double[] atom3 = mol
						.getCoord(atom.getNotDummyNeighbours().get(1).getNo());
				double[] atom4 = mol
						.getCoord(atom.getNotDummyNeighbours().get(2).getNo());
				LonePairAddition.addOnePairToTetrahedral(coord, atom2, atom3, atom4,
						protonCoord, hBondLength);
			} else
				throw new RuntimeException("Unable to protonate Tetrahedral Atom");

		}

		Atom proton = new Atom(mol, mol.getnAtoms(), AtomType.Type.H);
		Bond bond = new Bond(mol.getnBonds(), atom, proton, BondType.Type.SINGLE);
		mol.addAtom(proton, protonCoord);
		mol.addBond(bond);
		mol.update();

		mol.writeSybylMol2File(args[2], mol.getName() + " protonated");
	}
}
