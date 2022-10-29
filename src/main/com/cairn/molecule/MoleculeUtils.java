package com.cairn.molecule;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import com.cairn.molecule.Molecule.MolReadException;

/**
 * A simple class with a main method for common molecular operations
 * 
 * @author Gareth Jones
 * 
 */
public class MoleculeUtils {

	/**
	 * Transforms molecules by initializing- setting atom types, adding
	 * hydrogens and solvating
	 * 
	 * @param inFile
	 * @param outFile
	 */
	private static void prepare(String inFile, String outFile) {
		MolTransform transform = new MolTransform(inFile, outFile) {
			@Override
			public boolean transform() {
				molecule.init();
				return true;
			}
		};
		transform.run();
	}

	/**
	 * Transforms molecules by solvating only
	 * 
	 * @param inFile
	 * @param outFile
	 */
	private static void solvate(String inFile, String outFile) {
		final Solvate solvate = new Solvate();

		MolTransform transform = new MolTransform(inFile, outFile) {
			@Override
			public boolean transform() {
				molecule.assignAtomTypes();
				solvate.setMolecule(molecule);
				solvate.solvate();
				return true;
			}
		};
		transform.run();
	}

	/**
	 * Transforms molecules by adding hydrogens only
	 * 
	 * @param inFile
	 * @param outFile
	 */
	private static void addHydrogens(String inFile, String outFile) {
		final AddHydrogens addHydrogens = new AddHydrogens();
		MolTransform transform = new MolTransform(inFile, outFile) {
			@Override
			public boolean transform() {
				addHydrogens.setMolecule(molecule);
				addHydrogens.addHydrogens();
				molecule.assignAtomTypes();
				return true;
			}
		};
		transform.run();
	}

	/**
	 * Transforms molecules by removing lone pairs- assumes they've been added
	 * by GAPE- i.e. they're at the end of the atom and bond tables.
	 * 
	 * @param inFile
	 * @param outFile
	 */
	private static void removeLonePairs(String inFile, String outFile) {
		MolTransform transform = new MolTransform(inFile, outFile) {
			@Override
			public boolean transform() {
				molecule.removeLonePairs();
				molecule.update();
				return true;
			}
		};
		transform.run();
	}

	/**
	 * Splits a large file of molecules into separate molecule molecules.
	 * 
	 * @param inFile
	 */
	private static void split(String inFile) {
		MolTransform transform = new MolTransform(inFile) {
			@Override
			public boolean transform() {
				Molecule.FileType fileFormat = inFileType;
				String outName = molecule.getName().replace(' ', '_');
				if (fileFormat == Molecule.FileType.MOL2)
					outName += ".mol2";
				else
					outName += ".sdf";
				System.out.println("Writing to " + outName);
				molecule.write(outName, null);
				return true;
			}
		};

		transform.run();
	}

	private static void taffEnergy(String inFile) {
		Molecule.setAddHydrogensFlag(false);
		Molecule.setChargeFlag(true);
		Molecule.setSolvateFlag(false);

		MolTransform transform = new MolTransform(inFile) {
			@Override
			public boolean transform() {
				molecule.init();
				Taff taff = new Taff(molecule);
				double energy = taff.molEnergy();
				System.out.println(
						"Molecule " + molecule.getName() + " Total energy: " + energy);
				System.out.println("eVdw  " + taff.geteVdw());
				System.out.println("eBond " + taff.geteBond());
				System.out.println("eAng  " + taff.geteAng());
				System.out.println("eOop  " + taff.geteOop());
				System.out.println("eTor  " + taff.geteTor());
				return true;
			}
		};

		transform.run();
	}

	/**
	 * Creates an inclusion sphere in each molecule. Atoms outside the sphere
	 * are deleted.
	 * 
	 * @param inFile
	 * @param outFile
	 * @param x
	 * @param y
	 * @param z
	 * @param distance
	 */
	private static void coordinateSphere(String inFile, String outFile, final double x,
			final double y, final double z, final double distance) {
		MolTransform transform = new MolTransform(inFile, outFile) {
			@Override
			public boolean transform() {
				System.out.println("Exclusion sphere size " + distance + " at [" + x + ","
						+ y + "," + z + "]");
				molecule.inclusionSphere(x, y, z, distance);
				return true;
			}
		};
		transform.run();
	}

	/**
	 * Merges all the separate molecules in the input file into a single
	 * structure.
	 * 
	 * @param inFile
	 * @param outFile
	 */
	private static void merge(String inFile, String outFile) {
		List<Molecule> mols = Molecule.loadFiles(new String[] { inFile });
		System.out.println("molecule 0 " + " nAtoms " + mols.get(0).getnAtoms()
				+ " nBonds " + mols.get(0).getnBonds());
		int no = 0;
		for (Molecule mol : mols) {
			no++;
			System.out.println("merging molecule " + no + " nAtoms " + mol.getnAtoms()
					+ " nBonds " + mol.getnBonds());
			mols.get(0).merge(mol, false);
			System.out.println("molecule 0 " + " nAtoms " + mols.get(0).getnAtoms()
					+ " nBonds " + mols.get(0).getnBonds());
		}
		mols.get(0).build();
		System.out.println("Writing merged to" + outFile);
		mols.get(0).write(outFile, "Merged Molecules");
	}

	/**
	 * Splits all the input structures into fully connected fragments- such that
	 * each fragment has a separate structure entry in the output file
	 * 
	 * @param inFile
	 * @param outFile
	 */
	private static void findFragments(String inFile, String outFile) {

		List<Molecule> mols = Molecule.loadFiles(new String[] { inFile });

		SplitMolecule splitMolecule = new SplitMolecule();
		ArrayList<Molecule> allMolecules = new ArrayList<Molecule>();
		for (Molecule mol : mols) {
			// find fully connected fragments
			splitMolecule.setMolecule(mol);
			splitMolecule.findFragments();
			Molecule fragments[] = splitMolecule.getMoleculeFragments();
			for (int j = 0; j < fragments.length; j++)
				allMolecules.add(fragments[j]);
		}

		Molecule.write(allMolecules, outFile, "All fully connected fragments");
	}

	/**
	 * Application to provide command line access to the molecular transforms
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length == 0) {
			usage();
			return;
		}

		String op = args[0];

		if (op.equals("coordinate_sphere")) {
			String inFile = args[1];
			String outFile = args[2];
			double x = Double.parseDouble(args[3]);
			double y = Double.parseDouble(args[4]);
			double z = Double.parseDouble(args[5]);
			double distance = Double.parseDouble(args[6]);
			coordinateSphere(inFile, outFile, x, y, z, distance);
		}

		else if (op.equals("prepare")) {
			String inFile = args[1];
			String outFile = args[2];
			prepare(inFile, outFile);
		}

		else if (op.equals("merge")) {
			String inFile = args[1];
			String outFile = args[2];
			merge(inFile, outFile);
		}

		else if (op.equals("add_hydrogens")) {
			String inFile = args[1];
			String outFile = args[2];
			addHydrogens(inFile, outFile);
		}

		else if (op.equals("solvate")) {
			String inFile = args[1];
			String outFile = args[2];
			solvate(inFile, outFile);
		}

		else if (op.equals("find_fragments")) {
			String inFile = args[1];
			String outFile = args[2];
			findFragments(inFile, outFile);
		}

		else if (op.equals("remove_lone_pairs")) {
			String inFile = args[1];
			String outFile = args[2];
			removeLonePairs(inFile, outFile);
		}

		else if (op.equals("split")) {
			String inFile = args[1];
			split(inFile);
		}

		else if (op.equals("taff_energy")) {
			String inFile = args[1];
			taffEnergy(inFile);
		}

		else {
			System.out.println("Unknown Option " + op);
			usage();
		}
	}

	/**
	 * Prints out usage information.
	 */
	static void usage() {
		System.out.println("Usage:");
		System.out.println(
				"MoleculeUtils coordinate_sphere <in_file> <out_file> <x> <y> <z> <distance>");
		System.out.println("MoleculeUtils merge <in_file> <out_file>");
		System.out.println("MoleculeUtils prepare <in_file> <out_file>");
		System.out.println("MoleculeUtils add_hydrogens <in_file> <out_file>");
		System.out.println("MoleculeUtils solvate <in_file> <out_file>");
		System.out.println("MoleculeUtils find_fragments <in_file> <out_file>");
		System.out.println("MoleculeUtils split <in_file>");
		System.out.println("MoleculeUtils taff_energy <in_file>");

	}
}

/**
 * An abstract class for handling molecular transforms. It opens input and
 * output structure files. Call run to loop through all the structures applying
 * a transfomr.
 * 
 * @author Gareth Jones
 * 
 */
abstract class MolTransform {
	String inFile, outFile;

	BufferedReader in;

	BufferedWriter out;

	Molecule molecule;

	Molecule.FileType inFileType, outFileType;

	int moleculeNo;

	MolTransform() {
		;
	}

	/**
	 * Creates a transform that both reads and writes molecules.
	 * 
	 * @param _inFile
	 * @param _outFile
	 */
	MolTransform(String _inFile, String _outFile) {
		inFile = _inFile;
		outFile = _outFile;
		openReader();
		openWriter();
	}

	/**
	 * Creates a transform that just reads molecule.
	 * 
	 * @param _inFile
	 */
	MolTransform(String _inFile) {
		inFile = _inFile;
		openReader();
	}

	/**
	 * Use this method to apply operations to a molecule.
	 */
	abstract boolean transform();

	/**
	 * Reads each input structure, applying the transform and outputting the
	 * transformed molecule
	 */
	void run() {
		while (getNextMolecule()) {
			String name = molecule.getName();
			if (transform()) {
				if (out != null)
					outputMolecule();
				moleculeNo++;
				if (moleculeNo % 100 == 0)
					System.out.println("Processed " + moleculeNo + " molecules");
			} else
				System.out.println("Failed to transform molecule " + name);
		}
	}

	/**
	 * Opens a structure file for reading. Works with gzipped files.
	 */
	private void openReader() {
		try {
			if (inFile.toUpperCase().endsWith(".GZ")) {
				GZIPInputStream zipStream = new GZIPInputStream(
						new FileInputStream(inFile));
				in = new BufferedReader(new InputStreamReader(zipStream));
			} else {
				in = new BufferedReader(new FileReader(inFile));
			}
		} catch (FileNotFoundException ex) {
			ex.printStackTrace();
			System.exit(0);
		} catch (IOException ex) {
			ex.printStackTrace();
			System.exit(0);

		}
		inFileType = Molecule.getType(inFile);
	}

	/**
	 * Opens a structure file for writing. Works with gzipped files.
	 */
	private void openWriter() {
		try {
			if (outFile.toUpperCase().endsWith(".GZ")) {
				GZIPOutputStream zipStream = new GZIPOutputStream(
						new FileOutputStream(outFile));
				out = new BufferedWriter(new OutputStreamWriter(zipStream));
			} else {
				out = new BufferedWriter(new FileWriter(outFile));
			}
		} catch (IOException ex) {
			ex.printStackTrace();
			System.exit(0);
		}
		outFileType = Molecule.getType(outFile);
	}

	/**
	 * Gets the next molecule.
	 * 
	 * @return false if there are no more molecules.
	 */
	private boolean getNextMolecule() {
		molecule = new Molecule();

		try {
			if (inFileType == Molecule.FileType.MOL2)
				molecule.loadSybylMol2(in);
			else {
				molecule.loadSdfMol(in);
			}
		} catch (IOException | MolReadException ex) {
			finish();
			return false;
		}
		return true;
	}

	/**
	 * Writes out the transformed molecule.
	 */
	private void outputMolecule() {
		try {
			molecule.write(outFileType, out, "", false);
		} catch (IOException ex) {
			ex.printStackTrace();
			System.exit(0);
		}
	}

	/**
	 * Closes files
	 */
	void finish() {
		try {
			in.close();
			if (out != null)
				out.close();
		} catch (IOException ex) {
			ex.printStackTrace();
			System.exit(0);
		}
	}

}