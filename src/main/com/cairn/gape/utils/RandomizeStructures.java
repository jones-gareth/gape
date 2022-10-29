package com.cairn.gape.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.cairn.gape.ga.BaseSupervisor;
import com.cairn.molecule.Molecule;
import com.cairn.molecule.RotatableBond;

/**
 * An application to randomize a file of structures.
 * 
 * @author Gareth Jones
 * 
 */
public class RandomizeStructures {
	private static ThreadLocal<BaseSupervisor> supervisor = new ThreadLocal<BaseSupervisor>();

	private static void randomizeMolecule(Molecule molecule, boolean rigid) {
		if (!rigid) {
			for (RotatableBond rotatableBond : molecule.getRotatableBonds()) {
				int bVal = supervisor.get().randomInt(0, 255);
				rotatableBond.rotateBond(bVal);
			}
		}

		double x = 10 * supervisor.get().normalRand();
		double y = 10 * supervisor.get().normalRand();
		double z = 10 * supervisor.get().normalRand();

		for (double[] coord : molecule.getCoords()) {
			coord[0] += x;
			coord[1] += y;
			coord[2] += z;
		}
	}

	public static void main(String args[]) {
		try {

			Molecule.setAssignTypesFlag(true);
			Molecule.setAddHydrogensFlag(true);
			Molecule.setSolvateFlag(false);
			Molecule.setChargeFlag(false);

			supervisor.set(new BaseSupervisor());
			supervisor.get().setupRandomGenerator();
			if (args.length != 2 && args.length != 3) {
				System.err.println("Usage: " + RandomizeStructures.class.getName()
						+ " [-rigid] <file> <randomized_file>");
				System.exit(0);
			}

			boolean rigid = false;
			List<String> argsList = new ArrayList<>(Arrays.asList(args));
			if (argsList.get(0).equals("-rigid")) {
				argsList.remove(0);
				rigid = true;
			}

			String infile = argsList.remove(0);
			String outfile = argsList.remove(0);

			List<Molecule> molecules = Molecule.loadFiles(new String[] { infile });
			for (Molecule molecule : molecules) {
				if (!rigid) {
					molecule.assignRotatableBonds();
				}
				randomizeMolecule(molecule, rigid);
			}

			Molecule.write(molecules, outfile, "Randomized structure");

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
