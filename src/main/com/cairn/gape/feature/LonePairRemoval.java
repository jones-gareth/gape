package com.cairn.gape.feature;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import com.cairn.gape.molecule.GaMolecule;
import com.cairn.molecule.Molecule;

/**
 * Application to remove lone pairs from a MOL2 file.
 * 
 * @author Gareth Jones
 *
 */
class LonePairRemoval {

	/**
	 * Removes all lone pars from structures in a file.
	 * 
	 * @param args
	 */
	public static void main(String args[]) {

		if (args.length != 1 && args.length != 2) {
			System.err.println("Usage: LonePairRemoval <mol2file> [<outfile>]");
			System.err.println();
			System.err.println("Removes Lone Pairs from file <mol2file>");
			System.err.println("Writes to <outfile> or \"<mol2file> no lp.mol2\".");
			System.exit(0);
		}

		if (!args[0].toUpperCase().endsWith(".MOL2")) {
			System.err.println(args[0] + " is not a mol2 file!!");
			System.exit(0);
		}

		try {
			List<GaMolecule> molecules = GaMolecule.loadFiles(new String[] { args[0] },
					Molecule.FileType.MOL2, Molecule.Source.FILE);

			String outName;
			if (args.length == 2)
				outName = args[1];
			else {
				String b = (new File(args[0])).getName();
				b = b.substring(0, b.length() - 5);
				outName = b + " no lp.mol2";
			}

			System.out.println("Writing to " + outName);
			FileWriter out = new FileWriter(new File(outName));
			for (GaMolecule molecule : molecules) {
				molecule.removeLonePairs();
				molecule.writeSybylMol2(out, null);
			}
			out.close();
		} catch (IOException ex) {
			System.err.println(ex);
		}
	}

}
