package com.cairn.gape;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import com.cairn.common.utils.Timer;
import com.cairn.gape.chromosome.SuperpositionChromosome;
import com.cairn.gape.feature.Feature.FeatureType;
import com.cairn.gape.feature.FeatureMapping;
import com.cairn.gape.molecule.GaMolecule;
import com.cairn.molecule.Molecule;

/**
 * Does GAPE search of a database. Given a query structure and molecule files
 * does pairwise GA overlap of the query against every target structure.
 * 
 * @author Gareth Jones
 * 
 */
class StructureSearch extends Superposition {

	private static final Logger logger = Logger.getLogger(StructureSearch.class);
	static {
		// logger.setLevel(Level.DEBUG);
	}
	private volatile GaMolecule query, target;
	private volatile List<GaMolecule> database;

	private final List<SuperpositionChromosome> matches = new ArrayList<>();

	private volatile FileWriter allSolutionsWriter;

	/**
	 * Structure search application.
	 * 
	 * @param args
	 */
	public static void main(String args[]) {
		try {
			StructureSearch ss = new StructureSearch();

			if (args.length < 2) {
				System.out
						.println("\nUsage: StructureSearch <conf file> <molecule files..>");
				System.out.println();
				System.out
						.println("First structure is query, remainder are database files");
				System.exit(0);
			}

			ss.setupInfoOutput();
			String molFiles[] = new String[args.length - 1];
			for (int i = 1; i < args.length; i++)
				molFiles[i - 1] = args[i];
			List<GaMolecule> molecules = GaMolecule.loadFiles(molFiles,
					Molecule.FileType.MOL2, Molecule.Source.FILE);

			GaMolecule query = molecules.get(0);
			List<GaMolecule> database = molecules.subList(1, molecules.size());
			ss.fileFormat = Molecule.getType(molFiles[0]);

			ss.searchMolecules(args[0], query, database);
		} catch (Exception e) {
			logger.error("Fatal Exception ", e);
		}
	}

	/**
	 * Searches a single molecule against an array of molecules.
	 * 
	 * @param configFile
	 * @param q
	 * @param d
	 */
	public void searchMolecules(String configFile, GaMolecule q, List<GaMolecule> d) {

		totalTimer = new Timer();
		runTimer = new Timer();

		query = q;
		database = d;

		try {
			init(configFile);
			matches.clear();
			solutions.clear();
			allSolutionsWriter = new FileWriter("All_Solutions.mol2");

			setupQuery();

			for (int i = 0; i < database.size(); i++) {
				runTimer.reset();

				target = database.get(i);
				setupTarget();
				matchTarget();

				int no = i + 1;
				runTimer.interval();
				infoMessageLogger.infoMessageln("Database " + no + " "
						+ database.get(i).getName() + " Execution Time "
						+ runTimer.info() + "\n");
				runTimer.reset();

			}

			allSolutionsWriter.close();
			summarize();

			totalTimer.interval();
			infoMessageln("Total  Execution Time " + totalTimer.info());

		} catch (IOException ex) {
			System.err.println("StructureSearch: IOException: " + ex);
		}
	}

	/**
	 * Searches a single target against the query.
	 * 
	 */
	private void matchTarget() {
		SuperpositionChromosome.setAllocated(false);

		solutions.clear();
		for (int i = 0; i < nRuns; i++) {
			run(i);
			best = saveBest();
			solutions.add(best);

			int no = i + 1;
			infoMessageln("\n\nSolution " + no + " " + best.fitnessInfo() + "\n");

			pop.freePop();
		}

		Map<SuperpositionChromosome, Integer> initialPositions = new HashMap<>();
		for (int i = 0; i < solutions.size(); i++) {
			initialPositions.put(solutions.get(i), i);
		}
		Collections.sort(solutions, chromosomeComparator());

		try {
			for (int i = 0; i < solutions.size(); i++) {
				int r = i + 1;
				SuperpositionChromosome c = solutions.get(i);
				int no = initialPositions.get(c) + 1;
				c.rebuild();

				infoMessageln("Rank " + r + " solution " + no + " " + c.fitnessInfo());

				String name = target.getName();
				target.setName("GA ranked " + r + " " + name);
				target.writeSybylMol2(allSolutionsWriter, null, incLonePairs);
				target.setName(name);

				if (i == 0)
					matches.add(c);

			}
			allSolutionsWriter.flush();
		}

		catch (IOException ex) {
			throw new RuntimeException("writeMol2File IO exception: " + ex);
		}
	}

	/**
	 * Writes out a ranked summary of hits.
	 * 
	 */
	private void summarize() {
		Map<SuperpositionChromosome, Integer> initialPositions = new HashMap<>();
		for (int i = 0; i < matches.size(); i++) {
			initialPositions.put(matches.get(i), i);
		}
		Collections.sort(matches, chromosomeComparator());
		try {
			FileWriter summaryWriter = new FileWriter("Matches.mol2");
			String name = query.getName();
			query.setName("Query " + name);
			query.writeSybylMol2(summaryWriter, null);
			query.setName(name);

			for (int i = 0; i < matches.size(); i++) {

				SuperpositionChromosome c = matches.get(i);
				target = database.get(initialPositions.get(c));
				setupTarget();
				c.rebuild();

				int r = i + 1;
				int no = initialPositions.get(c) + 1;

				infoMessageln("Rank " + r + " Target " + no + " name " + target.getName()
						+ c.fitnessInfo());

				name = target.getName();
				target.setName(name + " Target " + no + " Rank " + r);
				target.writeSybylMol2(summaryWriter, null, incLonePairs);
				target.setName(name);

			}
			summaryWriter.close();
		} catch (IOException ex) {
			throw new RuntimeException("IOException " + ex);
		}
	}

	/**
	 * Initialises the query.
	 * 
	 */
	private void setupQuery() {
		commonMolSetup();

		query.setRigid(true);
		baseMolecule = query;
		fittingMolecule = query;

		baseMoleculeNo = 0;

		query.setProblem(this);
		query.setup();
		if (tordist != null)
			tordist.moleculeMatchTordist(query);

		molecules = new ArrayList<>();
		molecules.add(query);
		// create space for target
		molecules.add(null);
		binaryEntryPoints = new int[] { 0, 0 };
		integerEntryPoints = new int[] { 0 };
		integerStringLength = query.getNFeatures();
		integerStringRanges = new int[integerStringLength];
	}

	/**
	 * Initializes a target for matching against the query.
	 * 
	 */
	private void setupTarget() {
		if (target.getProblem() == null) {
			target.setProblem(this);
			target.setup();
			if (tordist != null)
				tordist.moleculeMatchTordist(target);
		}

		int nBonds = target.getnRotatableBonds();
		int nCorners = target.getNFreeCorners();
		int nBits = nBonds * 8 + nCorners;
		molecules.set(1, target);
		binaryStringLength = nBits;

		setupFeatureClustering();

		Map<FeatureType, FeatureMapping> baseFeatures = query.getFeatureMappings();
		int pos = 0;
		for (FeatureMapping otherFeatures : baseFeatures.values()) {
			for (int k = 0; k < baseFeatures.get(otherFeatures.getFeatureType())
					.getNFeatures(); k++) {
				integerStringRanges[pos] = otherFeatures.getNFeatures();
				pos++;
			}
		}

	}

}
