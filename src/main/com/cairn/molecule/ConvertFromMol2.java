package com.cairn.molecule;

import java.util.ArrayList;
import java.util.List;

import com.cairn.molecule.AtomType.Type;
import com.cairn.molecule.Molecule.FileType;

/**
 * Convertes a multimol mol2 file to a multimol sdf file, stripping out dummy
 * atoms.
 * 
 * @author gjones
 * 
 */
public class ConvertFromMol2 {
	public static void main(String args[]) {
		try {
			Molecule.setAddHydrogensFlag(true);

			if (args.length != 1) {
				System.out.println("Usage: " + ConvertToMol2Files.class.getName()
						+ " <mol2file>");
				System.exit(0);
			}

			String fileName = args[0];
			if (Molecule.getType(fileName) != FileType.MOL2)
				throw new IllegalArgumentException("not an mol2 file");

			List<Molecule> molecules = Molecule.loadFiles(new String[] { fileName });
			String base = fileName.substring(0, fileName.toUpperCase().indexOf(".MOL2"));

			// remove dummy atoms
			for (Molecule molecule : molecules) {
				List<Atom> dummyAtoms = new ArrayList<Atom>();
				for (Atom atom : molecule.getAtoms()) {
					if (atom.getAtomType() == Type.DU)
						dummyAtoms.add(atom);
				}
				for (Atom atom : dummyAtoms)
					molecule.deleteAtom(atom);
			}

			String outFile = base + ".sdf";
			System.out.println("writing " + outFile);
			Molecule.writeMols(molecules, outFile);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
