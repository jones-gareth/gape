package com.cairn.molecule;

import static java.lang.System.out;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import com.cairn.gape.molecule.GaMolecule;
import com.cairn.gape.molecule.MultiGaMolecule;
import com.cairn.gape.utils.InfoMessageLogger;

/**
 * Searches a conformer library for conformations of a given query and finds the
 * closest.
 * 
 * @author Gareth Jones
 * 
 */
public class SearchConformations {

	private static Logger logger;
	static {
		logger = Logger.getLogger(SearchConformations.class);
		// logger.setLevel(Level.DEBUG);
	}

	private final List<GaMolecule> queries;
	private final List<MultiGaMolecule> library;
	private final Map<String, Double> minRmsMap = new HashMap<String, Double>();

	public SearchConformations(List<GaMolecule> queries, List<MultiGaMolecule> library) {
		super();
		this.queries = queries;
		this.library = library;
	}

	/**
	 * Does the search
	 */
	public void search() {

		// override process rms in Superimpose to keep track of the closest
		// rms to each query.
		Superimpose superimpose = new Superimpose() {

			/*
			 * (non-Javadoc)
			 * 
			 * @see com.cairn.molecule.Superimpose#processRms(double, int)
			 */
			@Override
			public void processRms(double rms, int no) {
				String queryName = confX.getName();
				String targetName = confY.getName();

				logger.debug("Isomorphism " + no + " query " + queryName + " target "
						+ targetName + " rms " + rms);
				if (minRmsMap.containsKey(targetName)) {
					double minRms = minRmsMap.get(targetName);
					if (minRms > rms)
						minRmsMap.put(targetName, rms);
				} else {
					minRmsMap.put(targetName, rms);
				}
			}
		};

		// make sure we add hydrogens etc
		// for (Molecule query : queries) {
		// query.logLevel = 0;
		// query.init();
		// }
		// for (Molecule target : library) {
		// target.logLevel = 0;
		// target.init();
		// }

		// search
		for (GaMolecule query : queries) {
			for (MultiGaMolecule target : library) {
				logger.debug("searching " + query.getName() + " against "
						+ target.getName());
				for (int confNo = 0; confNo < target.nConformers(); confNo++) {
					target.setConformer(confNo);
					logger.debug("searching conformer number " + confNo);
					superimpose.search(query, target);
				}
			}
		}

		// report
		out.println();
		double rmsSum = 0, rmsCnt = 0, maxRms = 0, minRms = Double.MAX_VALUE;
		for (Molecule query : queries) {
			String name = query.getName();

			if (minRmsMap.containsKey(name)) {
				double rms = minRmsMap.get(name);
				out.println(name + " " + rms);
				rmsCnt++;
				rmsSum += rms;
				if (rms > maxRms)
					maxRms = rms;
				if (rms < minRms)
					minRms = rms;
			} else
				out.println("No matching conformers for " + name);

		}

		out.println();
		out.println("count " + rmsCnt);
		out.println("avg " + rmsSum / rmsCnt);
		out.println("min " + minRms);
		out.println("max " + maxRms);

	}

	public static void main(String[] args) {
		try {
			if (args.length != 2) {
				out.println("Usage: " + SearchConformations.class.getName()
						+ " <query file> <library file>");
				System.exit(0);
			}
			InfoMessageLogger infoMessageLogger = new InfoMessageLogger();
			List<GaMolecule> queries = GaMolecule.loadFiles(new String[] { args[0] },
					Molecule.FileType.UNKNOWN, Molecule.Source.FILE, infoMessageLogger,
					true);
			List<MultiGaMolecule> library = MultiGaMolecule.loadMultiMolFiles(
					new String[] { args[1] }, Molecule.FileType.UNKNOWN,
					Molecule.Source.FILE, infoMessageLogger, true);

			SearchConformations search = new SearchConformations(queries, library);
			search.search();

		} catch (Exception e) {
			e.printStackTrace();
		}

	}
}
