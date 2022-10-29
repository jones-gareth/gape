package com.cairn.gape.utils;

import java.util.List;

import com.cairn.common.utils.Coord;
import com.cairn.gape.Superposition;
import com.cairn.gape.chromosome.SuperpositionChromosome;
import com.cairn.gape.feature.Feature;
import com.cairn.gape.molecule.GaMolecule;
import com.cairn.molecule.Molecule;

/**
 * Simple class to test GAPE scoring function
 * 
 * @author Gareth Jones
 * 
 */
class ScoreOverlay extends Superposition {

	ScoreOverlay() {
	}

	/*
	 * Application which takes a file of overlaid compounds and scores them
	 * using the GAPE scoring function
	 * 
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Superposition#main(java.lang.String[])
	 */
	public static void main(String args[]) {
		ScoreOverlay sp = new ScoreOverlay();
		if (args.length < 2) {
			System.out.println("Usage: ScoreOverlay <conf file> <molecule files..>");
			System.exit(0);
		}

		String molFiles[] = new String[args.length - 1];
		for (int i = 1; i < args.length; i++)
			molFiles[i - 1] = args[i];

		List<GaMolecule> molecules = GaMolecule.loadFiles(molFiles,
				Molecule.FileType.MOL2, Molecule.Source.FILE);

		sp.scoreMolecules(args[0], molecules);
	}

	/**
	 * Scores the moleucles using GAPE scoring function. Does this by creating a
	 * single chomsome containing coordinates from the molecules.
	 * 
	 * @param configFile
	 * @param m
	 */
	void scoreMolecules(String configFile, List<GaMolecule> m) {

		for (GaMolecule molecule : m) {
			molecule.setRandomize(false);
		}

		init(configFile);
		setupMolecules(m);
		SuperpositionChromosome chrom = new SuperpositionChromosome();
		chrom.createEmpty(this);
		SuperpositionChromosome.setNf(getNumberFormat());

		// fool chromosome into being fitted
		int length = molecules.size();
		for (int i = 0; i < length; i++) {
			GaMolecule molecule = molecules.get(i);
			int nFeatures = molecule.getNFeatures();
			for (int j = 0; j < nFeatures; j++) {
				double[] f = molecule.getAllFeature(j).calculateCoordinate();
				Coord.copy(f, chrom.getFeatureCoordinates(i, j));
			}
		}

		SuperpositionChromosome.setCurrentChromosome(chrom);
		chrom.setFitted(true);
		chrom.setOk(true);

		Feature.setRadius(finishFittingRadius);
		chrom.getFitness();
		System.out.println(chrom.fitnessInfo());
		outputSolution(chrom, "Score overlay", true);
		if (!useFeatureClustering)
			System.out.println(chrom.mappingInfo());

	}

}
