package com.cairn.gape;

import java.util.ArrayList;
import java.util.List;

import com.cairn.gape.chromosome.SuperpositionChromosome;
import com.cairn.gape.feature.Feature;
import com.cairn.gape.molecule.GaMolecule;

/**
 * Scores one or more overlays using the GAPE scoring function. This gives
 * pretty close results using GAPE output files, though there are differences
 * (presumably only rounding errors!).
 * 
 * @author Gareth Jones
 * 
 */
public class RankSolutions {

	public void rank(String configFile, String molFiles[]) {

		SuperpositionChromosome solutions[] = new SuperpositionChromosome[molFiles.length];
		SuperpositionChromosome creator = new SuperpositionChromosome();

		int no = 0;
		for (String molFile : molFiles) {
			Superposition superposition = new Superposition();
			superposition.init(configFile);
			Feature.setRadius(superposition.finishFittingRadius);

			List<? extends GaMolecule> molecules = superposition
					.loadMolecules(new String[] { molFile });
			List<GaMolecule> okMolecules = new ArrayList<GaMolecule>();
			for (GaMolecule m : molecules) {
				if (m.getName().matches("GA \\d+ Pharmacophore all features")
						|| m.getName().matches("GA \\d+ Pharmacophore"))
					continue;
				m.clearSdfFields();
				m.setRandomize(false);
				okMolecules.add(m);
			}

			superposition.setupMolecules(okMolecules);

			SuperpositionChromosome c = creator.createChromosome();
			c.createEmpty(superposition);

			// set rigid to fit
			for (GaMolecule m : molecules)
				m.setFixed(true);
			c.fitMolecules(false);

			// then turn off to score
			for (GaMolecule m : molecules)
				m.setFixed(false);
			c.getFitness();

			System.out.println("Fitness for " + molFile + " is " + c.fitnessInfo());
			for (GaMolecule m : molecules)
				m.setFixed(true);

			solutions[no] = c;

			no++;
		}

	}

	public static void main(String[] args) {
		if (args.length < 2) {
			System.out.print("Usage: " + RankSolutions.class.getName()
					+ " <configFile> <molFiles..>");
			System.exit(0);
		}
		try {
			String configFile = args[0];
			String[] molFiles = new String[args.length - 1];
			for (int i = 1; i < args.length; i++)
				molFiles[i - 1] = args[i];
			RankSolutions rank = new RankSolutions();
			rank.rank(configFile, molFiles);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
