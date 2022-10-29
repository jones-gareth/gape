package com.cairn.molecule;

import java.util.List;

/**
 * A utility class for finding the
 * 
 * @author Gareth Jones
 *
 */
public class FindClosestConformer {

	// private static final Logger logger =
	// Logger.getLogger(FindClosetConformer.class);

	public static void main(String[] args) {
		boolean init = false;

		if (args.length == 0) {
			System.err.println("Usage: " + FindClosestConformer.class.getName()
					+ " <multi-conformer files..>");
			System.exit(0);
		}

		List<Molecule> molecules = Molecule.loadFiles(args, init);
		Molecule query = molecules.remove(0);
		double minRms = molecules.stream().map(m -> rms(query, m))
				.peek(r -> System.out.println("Rms is " + r)).mapToDouble(r -> r).min()
				.getAsDouble();
		System.out.println("Minimum rms is " + minRms);

	}

	private static double rms(Molecule query, Molecule target) {
		boolean subgraph = false, matchElementalTypes = false;
		RmsMapMolecule rmsFitter = new RmsMapMolecule(matchElementalTypes, subgraph);
		rmsFitter.searchIsomorphisms(query, target);
		double rms = rmsFitter.getBestRms();
		return rms;
	}
}
