package com.cairn.molecule;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.Parser;
import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.lang3.StringUtils;

/**
 * Attempts to identify the subset of compounds in the overlay which match the
 * template.
 * 
 * Start with a single molecule for fitting. Create an alignment and find the
 * closest structure (which is not currently used for fitting) and add it to the
 * set of fitting compounds, if the overall rms for the fitting molecules is
 * less than the threshold. Continue adding fitting molecules while the
 * threshold is not exceeded.
 * 
 * Repeat for all molecules to create a set of overlays the best is the largest
 * (if there is a tie resolve by rms).
 * 
 * @author Gareth Jones
 * 
 */
public class IncrementalSuperpositionComparison extends SuperpositionComparison {

	protected double allFittingMoleculeRmsThreshold = 2.0;
	protected double everyFittingMoleculeRmsThreshold = 2.0;

	private volatile boolean writeOverlay, printSummary;

	// private static final Logger logger = Logger
	// .getLogger(IncrementalSuperpositionComparison.class);

	public IncrementalSuperpositionComparison(boolean init, boolean subgraph,
			boolean mapNames) {
		super(init, subgraph, mapNames);
	}

	public IncrementalSuperpositionComparison() {
		super();
	}

	/**
	 * Determines if the incremental overlay is OK. The criteria is that the
	 * average rms for all molecules to the template is less or equal to
	 * allFittingMoleculeRmsThreshold and that each fitting molecule in the
	 * overlay is within everyFittingMoleculeRmsThreshold of the template
	 * molecule
	 * 
	 * @return
	 */
	private boolean overlayOk() {
		int nMolecules = nMolecules();

		if (allFittingMoleculeRmsThreshold > .0
				&& getFittingMoleculeRms() > allFittingMoleculeRmsThreshold) {
			return false;
		}

		if (everyFittingMoleculeRmsThreshold > .0) {
			for (int i = 0; i < nMolecules; i++) {
				if (isFittingMolecule(i)) {
					double testRms = getMoleculeRmss(i);
					if (testRms > everyFittingMoleculeRmsThreshold) {
						return false;
					}
				}
			}
		}

		return true;
	}

	/**
	 * Add the closest compound which is not used as a fitting molecule. Return
	 * true if the overlay is OK
	 * {@link IncrementalSuperpositionComparison#overlayOk()}
	 * 
	 * @return
	 */
	private boolean extendOverlay() {
		int nMolecules = nMolecules();

		double minRms = Double.MAX_VALUE;
		int closestNo = -1;
		for (int i = 0; i < nMolecules; i++) {
			if (isFittingMolecule(i))
				continue;
			double testRms = getMoleculeRmss(i);
			if (testRms < minRms) {
				minRms = testRms;
				closestNo = i;
			}
		}

		assert closestNo != -1 : "failed to find free molecule";

		setFittingMolecule(closestNo, true);
		alignMolecules();

		return overlayOk();
	}

	/**
	 * Build an overlay starting with a particular molecule.
	 * 
	 * @param start
	 *            the starting structure.
	 * 
	 * @return
	 */
	private OverlayResult buildOverlay(int start) {
		int nMolecules = nMolecules();
		for (int i = 0; i < nMolecules; i++) {
			setFittingMolecule(i, false);
		}
		setFittingMolecule(start, true);
		alignMolecules();

		OverlayResult overlayResult = new OverlayResult();
		Set<Integer> fittingMolecules = overlayResult.overlay;

		if (overlayOk()) {

			while (extendOverlay()) {
				for (int i = 0; i < nMolecules; i++) {
					if (isFittingMolecule(i))
						fittingMolecules.add(i);
				}
				overlayResult.rms = getFittingMoleculeRms();

				if (fittingMolecules.size() == nMolecules)
					break;

			}

		}

		return overlayResult;

	}

	/**
	 * Class to hold result for each overlay
	 * 
	 */
	private static class OverlayResult {
		double rms = .0;
		Set<Integer> overlay = new HashSet<Integer>();
	}

	/**
	 * Build all the overlays and identify the best one.
	 * 
	 */
	public boolean buildOverlays() {
		List<OverlayResult> overlays = new ArrayList<OverlayResult>();

		// for each molecule build an overlay
		int nMolecules = nMolecules();
		for (int start = 0; start < nMolecules; start++) {
			overlays.add(buildOverlay(start));
		}

		// sort overlays by size (then rms). The last overlay in the list is the
		// best.
		Collections.sort(overlays, new Comparator<OverlayResult>() {

			@Override
			public int compare(OverlayResult arg0, OverlayResult arg1) {
				Integer int1 = arg0.overlay.size();
				Integer int2 = arg1.overlay.size();
				int test = int1.compareTo(int2);
				if (test != 0)
					return test;
				Double d1 = arg0.rms;
				Double d2 = arg1.rms;
				return d1.compareTo(d2);
			}

		});

		List<Molecule> molecules = getMoleculeOverlay();

		boolean pass = false;

		// print out a summary
		for (int start = 0; start < nMolecules; start++) {

			OverlayResult overlayResult = overlays.get(start);
			List<Integer> overlay = new ArrayList<Integer>(overlayResult.overlay);
			Collections.sort(overlay);
			List<String> names = new ArrayList<String>();
			for (int moleculeNo : overlay) {
				names.add(molecules.get(moleculeNo).getName());
			}

			String info = "";
			if (start == nMolecules - 1) {
				info += "**Best** ";
				if (overlay.size() * 2 >= nTemplateMolecules) {
					info += "**pass** ";
					pass = true;
				} else {
					pass = false;
					info += "**fail** ";
				}
			}
			info += overlay.size() + "/" + nMolecules + "/" + nTemplateMolecules + ", "
					+ overlayResult.rms + " [" + StringUtils.join(names, ',') + "]";

			if (printSummary)
				System.out.println(info);
		}

		if (writeOverlay) {
			// rebuild the best result and save it
			Set<Integer> best = overlays.get(nMolecules - 1).overlay;
			if (CollectionUtils.isNotEmpty(best)) {
				for (int i = 0; i < nMolecules; i++) {
					setFittingMolecule(i, best.contains(i));
				}
				alignMolecules();
				writeOverlay();
				overlaySummary();
			}
		}

		return pass;
	}

	/**
	 * @return the writeOverlay
	 */
	public boolean isWriteOverlay() {
		return writeOverlay;
	}

	/**
	 * @param writeOverlay
	 *            the writeOverlay to set
	 */
	public void setWriteOverlay(boolean writeOverlay) {
		this.writeOverlay = writeOverlay;
	}

	/**
	 * @return the printSummary
	 */
	public boolean isPrintSummary() {
		return printSummary;
	}

	/**
	 * @return the allFittingMoleculeRmsThreshold
	 */
	public double getAllFittingMoleculeRmsThreshold() {
		return allFittingMoleculeRmsThreshold;
	}

	/**
	 * @return the everyFittingMoleculeRmsThreshold
	 */
	public double getEveryFittingMoleculeRmsThreshold() {
		return everyFittingMoleculeRmsThreshold;
	}

	/**
	 * @param allFittingMoleculeRmsThreshold
	 *            the allFittingMoleculeRmsThreshold to set
	 */
	public void setAllFittingMoleculeRmsThreshold(double allFittingMoleculeRmsThreshold) {
		this.allFittingMoleculeRmsThreshold = allFittingMoleculeRmsThreshold;
	}

	/**
	 * @param everyFittingMoleculeRmsThreshold
	 *            the everyFittingMoleculeRmsThreshold to set
	 */
	public void setEveryFittingMoleculeRmsThreshold(
			double everyFittingMoleculeRmsThreshold) {
		this.everyFittingMoleculeRmsThreshold = everyFittingMoleculeRmsThreshold;
	}

	/**
	 * @param printSummary
	 *            the printSummary to set
	 */
	public void setPrintSummary(boolean printSummary) {
		this.printSummary = printSummary;
	}

	protected boolean processArguments(String[] args) {
		boolean init = false, subGraph = false, matchByOrder = false;

		Option initOption = new Option("i", "init", false,
				"initialize structures using GAPE");
		Option subGraphOption = new Option("s", "subGraph", false,
				"consider template a substructure");
		Option matchByOrderOption = new Option("m", "matchByOrder", false,
				"Match in structure order (rather than attempting  to map names)");

		Option allFittingMoleculeRmsThresholdOption = new Option("a",
				"allFittingMoleculeRmsThreshold", true,
				"The rms threshold for all molecules in the overlay");
		Option everyFittingMoleculeRmsThresholdOption = new Option("e",
				"everyFittingMoleculeRmsThreshold", true,
				"The rms threshold for each molecule in the overlay");
		Option helpOption = new Option("h", "help", false, "Print help message");

		Options options = new Options();
		options.addOption(initOption);
		options.addOption(subGraphOption);
		options.addOption(matchByOrderOption);
		options.addOption(helpOption);
		options.addOption(allFittingMoleculeRmsThresholdOption);
		options.addOption(everyFittingMoleculeRmsThresholdOption);

		Parser parser = new BasicParser();
		CommandLine line = null;
		try {
			line = parser.parse(options, args);
		} catch (ParseException ex) {
			throw new RuntimeException(ex);
		}
		String[] extraArgs = line.getArgs();
		if (line.hasOption("help") || extraArgs.length != 2) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(getClass().getName()
					+ " <options> <molecules> <template>", options);
			return false;
		}

		if (line.hasOption("init")) {
			init = true;
		}
		if (line.hasOption("subGraph")) {
			subGraph = true;
		}
		if (line.hasOption("matchByOrder")) {
			matchByOrder = true;
		}

		this.setInit(init);
		this.setSubgraph(subGraph);
		setMapNames(!matchByOrder);
		setPrintSummary(true);
		setWriteOverlay(true);
		setMatchElementalTypes(true);

		if (line.hasOption("allFittingMoleculeRmsThreshold")) {
			double allFittingMoleculeRmsThreshold = Double.parseDouble(line
					.getOptionValue("allFittingMoleculeRmsThreshold"));
			setAllFittingMoleculeRmsThreshold(allFittingMoleculeRmsThreshold);
		}
		if (line.hasOption("everyFittingMoleculeRmsThreshold")) {
			double everyFittingMoleculeRmsThreshold = Double.parseDouble(line
					.getOptionValue("everyFittingMoleculeRmsThreshold"));
			setEveryFittingMoleculeRmsThreshold(everyFittingMoleculeRmsThreshold);
		}
		loadMolecules(extraArgs[0], extraArgs[1]);

		return true;
	}

	/**
	 * Application attempts to identify the subset of compounds in the overlay
	 * which match the template
	 * 
	 * @param args
	 */
	public static void main(String args[]) {
		IncrementalSuperpositionComparison superpositonComparison = new IncrementalSuperpositionComparison();
		if (!superpositonComparison.processArguments(args)) {
			return;
		}
		superpositonComparison.buildOverlays();
	}
}
