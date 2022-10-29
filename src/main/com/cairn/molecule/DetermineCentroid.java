package com.cairn.molecule;

import java.util.List;

/**
 * Application class for determining the heavy atom centroid of a set of
 * molecules
 * 
 * @author Gareth Jones
 *
 */
public class DetermineCentroid {
	public static void main(String args[]) throws Exception {
		List<Molecule> molecules = Molecule.loadFiles(args, false);
		for (Molecule molecule : molecules) {
			double[] centroid = molecule.calculateCentroid(true);
			System.out.println(molecule.getName() + ": " + centroid[0] + " "
					+ centroid[1] + " " + centroid[2]);
		}
	}
}
