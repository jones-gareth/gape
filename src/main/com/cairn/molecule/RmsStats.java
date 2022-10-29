package com.cairn.molecule;

import java.util.DoubleSummaryStatistics;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;

/**
 * A application which takes a list of molecules. The first is a template
 * compound. All other molecules are compared to it and average, max and min rms
 * values determined.
 * 
 * @author Gareth Jones
 *
 */
public class RmsStats {
	private static final Logger logger = Logger.getLogger(RmsStats.class);

	public static void main(String[] args) {
		boolean init = false, subgraph = false, matchElementalTypes = true;
		if (args.length < 1) {
			System.err.println("Usage: " + RmsStats.class.getName()
					+ " <molecule files..>");
			System.exit(0);
		}

		List<Molecule> molecules = Molecule.loadFiles(args, init);

		if (molecules.size() < 2) {
			logger.error("Only " + molecules.size()
					+ " molecules present.  Unable to determine rms between 2 molecules");
			System.exit(0);
		}

		Rms rmsFitter = new Rms(subgraph, matchElementalTypes);
		Molecule template = molecules.get(0);
		molecules = molecules.subList(1, molecules.size());

		List<Double> rmsValues = molecules.stream()
				.map(m -> rmsFitter.determineRms(template, m))
				.peek(r -> System.out.println("Rms is " + r))
				.collect(Collectors.toList());
		DoubleSummaryStatistics stats = rmsValues.stream().collect(
				Collectors.summarizingDouble(r -> r));
		String summary = rmsValues.stream().map(r -> String.valueOf(r))
				.collect(Collectors.joining(" "));

		System.out.println("RMS values " + summary);
		System.out.println("Min " + stats.getMin() + " Avg " + stats.getAverage()
				+ " Max " + stats.getMax());
	}
}
