package com.cairn.molecule;

import java.util.List;

import com.cairn.molecule.Molecule.FileType;

/**
 * Convertes a multimol sdf file to a number of mol2 files.
 * 
 * @author Gareth Jones
 * 
 */
public class ConvertToMol2Files {

	public static void main(String args[]) {
		try {
			Molecule.setAddHydrogensFlag(true);

			if (args.length != 1) {
				System.out.println("Usage: " + ConvertToMol2Files.class.getName()
						+ " <sdfile>");
				System.exit(0);
			}

			String fileName = args[0];
			if (Molecule.getType(fileName) != FileType.SDF)
				throw new IllegalArgumentException("not an sdf file");

			List<Molecule> molecules = Molecule.loadFiles(new String[] { fileName });
			String base = fileName.substring(0, fileName.toUpperCase().indexOf(".SDF"));
			int no = 0;
			for (Molecule mol : molecules) {
				no++;
				String outFile = base + "_" + no + ".mol2";
				System.out.println("writing " + outFile);
				mol.writeSybylMol2File(outFile, "converted from molecule " + no + " in "
						+ fileName);
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}