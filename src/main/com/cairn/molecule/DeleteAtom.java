package com.cairn.molecule;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.cairn.molecule.Molecule.FileType;

/**
 * Deletes atoms from a structure file
 * 
 * @author Gareth Jones
 * 
 */
class DeleteAtom {
	public static void main(String args[]) {

		try {
			Molecule.setAssignTypesFlag(false);
			Molecule.setAddHydrogensFlag(false);
			Molecule.setSolvateFlag(false);
			Molecule.setChargeFlag(false);

			if (args.length != 2 && args.length != 3) {
				System.err.println("Usage: DeleteAtom <mol2file> <atom_ids..>");
				System.exit(0);
			}

			FileType fileType = Molecule.getType(args[0]);
			Molecule molecule = new Molecule(args[0], fileType, Molecule.Source.FILE);

			String outName;
			String b = (new File(args[0])).getName();
			b = b.substring(0, b.lastIndexOf('.'));
			outName = b + "_edit." + fileType.getSuffix();

			System.out.println("Writing to " + outName);

			List<Integer> atomIds = new ArrayList<Integer>();

			for (int i = 1; i < args.length; i++) {
				int no = Integer.valueOf(args[i]);
				no--;
				atomIds.add(no);

			}

			// delete atoms from the largest index down so we don't have to
			// renumber.
			Collections.sort(atomIds, Collections.reverseOrder());
			for (Integer atomId : atomIds) {
				System.out.println("Deleting atom " + atomId);
				molecule.deleteAtom(molecule.getAtom(atomId));
			}

			molecule.write(outName, null);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

}
