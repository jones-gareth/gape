package com.cairn.molecule;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Logger;

import com.cairn.gape.utils.CombinatorialChooser;
import com.cairn.gape.utils.CombinatorialChooser.CombinatorialChooserCallback;

public class ExhaustiveSuperpositionComparison extends IncrementalSuperpositionComparison
		implements CombinatorialChooserCallback<Integer> {

	private static final Logger logger = Logger
			.getLogger(ExhaustiveSuperpositionComparison.class);
	static {
		// logger.setLevel(Level.DEBUG);
	}
	private Solution bestSolution = null;

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.molecule.IncrementalSuperpositionComparison#buildOverlays
	 * ()
	 */
	@Override
	public boolean buildOverlays() {
		int nRequired = nTemplateMolecules / 2 + nTemplateMolecules % 2;
		int nMolecules = nMolecules();
		bestSolution = null;

		if (nRequired > nMolecules) {
			logger.info("Fail because insufficient number of molecules in overlay!");
		} else {
			IntStream.rangeClosed(nRequired, nMolecules).forEach(
					nChoose -> {
						logger.info("Enumerating overlays of size " + nChoose);
						CombinatorialChooser<Integer> chooser = CombinatorialChooser
								.combinatorialIndexChooser(nChoose, nMolecules, this);
						chooser.enumerate();
					});
		}

		// remove any old fitting file
		File fittedMoleculeFile = new File(fittedMoleculeFilename());
		if (fittedMoleculeFile.exists()) {
			fittedMoleculeFile.delete();
		}

		if (bestSolution != null) {
			logger.info("Best solution: " + StringUtils.join(bestSolution.values, ", ")
					+ " mean rms " + bestSolution.meanRms);
			System.out.print("**Best** **pass** " + bestSolution.values.size() + "/"
					+ nMolecules + "/" + nTemplateMolecules + " " + bestSolution.meanRms
					+ " ");
			String str = bestSolution.values.stream()
					.map((i) -> getMoleculeOverlay().get(i).getName())
					.collect(Collectors.joining(",", "[", "]"));
			System.out.println(str);

			// rebuild the best result and save it
			IntStream.range(0, nMolecules).forEach(no -> setFittingMolecule(no, false));
			bestSolution.values.forEach(value -> setFittingMolecule(value, true));
			alignMolecules();
			writeOverlay();
			overlaySummary();
		} else {
			System.out.println("**Best** **fail**");
		}

		return bestSolution != null;
	}

	public static void main(String args[]) {
		ExhaustiveSuperpositionComparison superpositonComparison = new ExhaustiveSuperpositionComparison();
		if (!superpositonComparison.processArguments(args)) {
			return;
		}
		superpositonComparison.buildOverlays();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.utils.CombinatorialChooser.CombinatorialChooserCallback
	 * #callback(java.util.List)
	 */
	@Override
	public void callback(List<Integer> values) {
		int nMolecules = nMolecules();
		IntStream.range(0, nMolecules).forEach(no -> setFittingMolecule(no, false));
		values.forEach(value -> setFittingMolecule(value, true));
		alignMolecules();

		boolean ok = true;
		double allFittingMoleculeRms = getFittingMoleculeRms();

		if (allFittingMoleculeRmsThreshold > .0
				&& allFittingMoleculeRms > allFittingMoleculeRmsThreshold) {
			ok = false;
		}

		double meanRms = -1;
		if (ok) {
			List<Double> testRmss = values.stream().map(no -> getMoleculeRmss(no))
					.collect(Collectors.toList());
			ok = testRmss.stream().allMatch(
					rms -> rms <= everyFittingMoleculeRmsThreshold);
			meanRms = testRmss.stream().mapToDouble(Double::doubleValue).average()
					.getAsDouble();
		}

		String label = StringUtils.join(values, ", ");
		logger.debug(label + ": ok " + ok + " all fitting molecule rms "
				+ allFittingMoleculeRms + " mean rms " + meanRms);
		if (ok) {
			logger.debug(label + " passes mean rms " + meanRms);
			Solution currentSolution = new Solution(values, meanRms);
			if (bestSolution == null || currentSolution.compareTo(bestSolution) == 1) {
				bestSolution = currentSolution;
				logger.info("New best solution: "
						+ StringUtils.join(bestSolution.values, ", ") + " mean rms "
						+ bestSolution.meanRms);
			}
		}
	}

	private static class Solution implements Comparable<Solution> {
		private final List<Integer> values;
		private final double meanRms;

		public Solution(List<Integer> values, double meanRms) {
			super();
			this.values = new ArrayList<>(values);
			this.meanRms = meanRms;
		}

		@Override
		public int compareTo(Solution o) {
			if (values.size() - o.values.size() == 0) {
				return Double.compare(o.meanRms, meanRms);
			}
			return Integer.compare(values.size(), o.values.size());
		}

	}
}
