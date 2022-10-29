package com.cairn.gape.molecule;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import com.cairn.molecule.Molecule;
import com.cairn.molecule.Molecule.MolReadException;

/**
 * Relabels the atoms to C1, C2, H1, H2 etc. The ordering is the same as the
 * input file.
 * 
 * @author Gareth Jones
 * 
 * @see com.cairn.molecule.Molecule#relabelAtoms()
 */
class RelabelMol2 {

	public static void main(String args[]) {
		if (args.length != 2) {
			System.err.println("Usage: RelabelMol2 <in file.mol2> <out file.mol2>");
			System.exit(0);
		}

		BufferedReader in = Molecule.openReader(args[0], Molecule.Source.FILE);
		if (in == null)
			System.exit(0);

		ArrayList<GaMolecule> molecules = new ArrayList<GaMolecule>();
		while (true) {
			GaMolecule molecule = new GaMolecule();
			try {
				molecule.loadSybylMol2(in);
				molecule.build();
				molecules.add(molecule);
			} catch (IOException | MolReadException ex) {
				break;
			}
		}
		try {
			in.close();
		} catch (IOException ex) {
			System.err.println("IO error closing " + ex);
		}

		for (int i = 0; i < molecules.size(); i++) {
			GaMolecule molecule = molecules.get(i);
			molecule.relabelAtoms();
		}

		try {
			FileWriter out = new FileWriter(args[1]);
			for (int i = 0; i < molecules.size(); i++) {
				GaMolecule molecule = molecules.get(i);
				molecule.writeSybylMol2(out, null);
			}
			out.close();
		} catch (java.io.IOException ex) {
			System.out.println(ex);
		}
	}
}