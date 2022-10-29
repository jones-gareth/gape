package com.cairn.molecule;

import static java.lang.System.out;

import java.util.List;

/**
 * Application to extract a structure from a file of structures and save in a
 * new file.
 * 
 * @author Gareth Jones
 * 
 */
public class ExtractStructure {
	public static void main(String[] args) {

		try {
			if (args.length != 1 && args.length != 3) {
				String name = ExtractStructure.class.getClass().getName();
				out.println("Usage: " + name
						+ "<molecule_name> <molecules_file> <output_file>");
				out.println("      or");
				out.println("       " + name + "<molecules_file>");
				System.exit(0);
			}

			if (args.length == 1)
				listStructures(args[0]);
			else
				findStructures(args[0], args[1], args[2]);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void listStructures(String file) {
		List<Molecule> molecules = Molecule.loadFiles(new String[] { file });
		for (Molecule molecule : molecules) {
			out.println(molecule.getName());
		}
	}

	public static void findStructures(String name, String molFile, String saveFile) {

		List<Molecule> molecules = Molecule.loadFiles(new String[] { molFile });
		boolean saved = false;
		for (Molecule molecule : molecules) {
			if (molecule.getName().equalsIgnoreCase(name)) {
				molecule.write(saveFile, null);
				saved = true;
				break;
			}

		}
		if (!saved)
			out.println("failed find any molecule named " + name);
	}
}
