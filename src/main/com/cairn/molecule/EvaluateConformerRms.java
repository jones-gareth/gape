package com.cairn.molecule;

import java.util.List;

import org.apache.log4j.Logger;

/**
 * A utility class to read a file containing multiconformers of a single
 * molecule.
 * 
 * @author Gareth Jones
 *
 */
public class EvaluateConformerRms {
	private static final Logger logger = Logger.getLogger(EvaluateConformerRms.class);

	public static void main(String[] args) {
		boolean init = false, subgraph = false, matchElementalTypes = false;
		if (args.length == 0) {
			System.err.println("Usage: " + EvaluateConformerRms.class.getName()
					+ " <multi-conformer file>");
			System.exit(0);
		}

		List<Molecule> molecules = Molecule.loadFiles(args, init);
		for (int i = 0; i < molecules.size(); i++) {
			for (int j = i + 1; j < molecules.size(); j++) {
				RmsMapMolecule rmsFitter = new RmsMapMolecule(matchElementalTypes,
						subgraph);
				rmsFitter.searchIsomorphisms(molecules.get(j), molecules.get(i));
				double rms = rmsFitter.getBestRms();
				logger.info("Rms between conformer " + String.valueOf(i + 1)
						+ " and conformer " + String.valueOf(j + 1) + " is " + rms);
				// System.out.println("{make_pair(" + i + ", " + j + "), " + rms
				// + " },");
			}
			break;
		}

	}
}
