package com.cairn.molecule;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;

/**
 * A class to superimpose two overlays.
 * 
 * @author Gareth Jones
 * 
 */
public class SuperpositionComparison {
	private static final Logger logger = Logger.getLogger(SuperpositionComparison.class);

	static {
		// logger.setLevel(Level.DEBUG);
	}

	private List<Molecule> moleculeOverlay, templateOverlay, pharmMolecules, allMolecules;

	private List<HashMap<Integer, Integer>> mappings;
	private boolean[] fittingMolecules;
	private final List<Double> moleculeRmss = new ArrayList<Double>();
	private double fittingMoleculeRms, allMoleculeRms;
	private int nFittingAtoms;
	private String moleculeFile;

	private int nAllAtoms;
	private boolean init = false, subgraph = false, mapNames = false;
	protected boolean matchElementalTypes = false;
	protected int nTemplateMolecules;

	public SuperpositionComparison(boolean init, boolean subgraph, boolean mapNames) {
		super();
		this.init = init;
		this.subgraph = subgraph;
		this.mapNames = mapNames;
	}

	public SuperpositionComparison() {
		;
	}

	/**
	 * Application. Attempts to fit one overlay to another
	 * 
	 * @param args
	 */
	public static void main(String args[]) {
		if (args.length < 2) {
			System.out.println(
					"Usage: SuperpositionComparison [-p] [-s] [-m] [-matchElementalTypes]"
							+ "[-i <names>] [-x <names>] <molecules> <template>");
			System.out.println("-p: prepare input structures using GAPE");
			System.out.println("-m: attempt to map names");
			System.out.println(
					"-matchElementalTypes: match on elemental rather than hybridised types");
			System.out.println("-i <names>: include structures that contain these "
					+ "names in the overlay (comma separated names)");
			System.out.println("-x <names>: exclude structures that contain these "
					+ "names in the overlay (comma separated names)");
			System.out.println("-s: consider template as a substructure of molecules");
			System.exit(0);
		}

		int no = 0;
		boolean init = false, subgraph = false, mapNames = false,
				matchElementalTypes = false;
		String[] includeNames = null, excludeNames = null;
		while (args[no].startsWith("-")) {
			String arg = args[no].toLowerCase();

			if (arg.equals("-p"))
				init = true;
			else if (arg.equals("-s"))
				subgraph = true;
			else if (arg.equals("-m"))
				mapNames = true;
			else if (arg.equals("-matchelementaltypes"))
				matchElementalTypes = true;
			else if (arg.equals("-i")) {
				includeNames = StringUtils.split(args[++no], ",");
			} else if (arg.equals("-x")) {
				excludeNames = StringUtils.split(args[++no], ",");
			} else {
				System.out.println("unknown flag: " + arg);
				System.exit(0);
			}
			no++;
		}

		String moleculeFile = args[no++];
		String templateFile = args[no++];

		SuperpositionComparison superpositionComparison = new SuperpositionComparison(
				init, subgraph, mapNames);
		superpositionComparison.setMatchElementalTypes(matchElementalTypes);
		superpositionComparison.alignOverlays(moleculeFile, templateFile, includeNames,
				excludeNames);
	}

	/**
	 * Loads the overlay and template
	 * 
	 * @param moleculeFile
	 * @param templateFile
	 */
	public void loadMolecules(String moleculeFile, String templateFile) {

		this.moleculeFile = moleculeFile;

		Molecule.setChargeFlag(init);
		Molecule.setAddHydrogensFlag(init);
		Molecule.setSolvateFlag(init);

		moleculeOverlay = loadFile(moleculeFile, false);
		templateOverlay = loadFile(templateFile, true);
		nTemplateMolecules = templateOverlay.size();

		allMolecules = new ArrayList<Molecule>(moleculeOverlay);
		allMolecules.addAll(pharmMolecules);

		if (mapNames) {
			ArrayList<Molecule> mapped = new ArrayList<Molecule>();
			for (Molecule molecule : moleculeOverlay) {
				boolean matched = false;
				for (Molecule template : templateOverlay) {
					if (matchNames(molecule.getName(), template.getName())) {
						matched = true;
						logger.info(molecule.getName() + " -> " + template.getName());
						mapped.add(template);
						break;
					}
				}
				if (!matched)
					throw new RuntimeException(
							"Unable to find template for " + molecule.getName());
			}
			templateOverlay = mapped;
		}

		nAllAtoms = 0;
		// For each molecule find the optimal mapping
		if (moleculeOverlay.size() != templateOverlay.size())
			throw new RuntimeException("alignOverlays: molecule count mismatch");
		mappings = new ArrayList<HashMap<Integer, Integer>>();
		int no = 1;
		for (int i = 0; i < moleculeOverlay.size(); i++) {
			RmsMapMolecule mapMolecule = new RmsMapMolecule(matchElementalTypes,
					subgraph);
			Molecule molecule = moleculeOverlay.get(i);
			Molecule template = templateOverlay.get(i);
			HashMap<Integer, Integer> mapping = mapMolecule.searchIsomorphisms(molecule,
					template);
			int nIsomorphisms = mapMolecule.getNIsomorphisms();
			if (nIsomorphisms < 1)
				throw new RuntimeException("unable to map " + template.getName() + " to "
						+ molecule.getName());
			nAllAtoms += mapping.size();
			mappings.add(mapping);
			no *= nIsomorphisms;
			if (logger.isDebugEnabled())
				logger.debug("Mapping " + molecule.getName() + " to " + template.getName()
						+ " n isomorphisms " + nIsomorphisms + " best rms "
						+ mapMolecule.getBestRms());
		}
		logger.debug("total no of mappings " + no);
		logger.debug("total no of atom mappings " + nAllAtoms);

		fittingMolecules = new boolean[moleculeOverlay.size()];

	}

	/**
	 * @return the matchElementalTypes
	 */
	public boolean isMatchElementalTypes() {
		return matchElementalTypes;
	}

	/**
	 * @param matchElementalTypes
	 *            the matchElementalTypes to set
	 */
	public void setMatchElementalTypes(boolean matchElementalTypes) {
		this.matchElementalTypes = matchElementalTypes;
	}

	/**
	 * 
	 * This routine does the work of aligning the two overlays and determining
	 * RMS distances. We check all isomorphisms within each molecule.
	 * 
	 * @param init
	 * @param moleculeFile
	 * @param templateFile
	 */
	private void alignOverlays(String moleculeFile, String templateFile,
			String[] includeNames, String[] excludeNames) {
		loadMolecules(moleculeFile, templateFile);

		if (moleculeOverlay.size() == 0)
			throw new RuntimeException("no target molecules");
		if (templateOverlay.size() == 0)
			throw new RuntimeException("no template molecules");

		for (int i = 0; i < moleculeOverlay.size(); i++) {
			Molecule molecule = moleculeOverlay.get(i);
			Molecule template = templateOverlay.get(i);

			boolean fit = true;
			if (includeNames != null) {
				fit = matchName(molecule, includeNames);
			}
			if (fit && excludeNames != null) {
				fit = !matchName(molecule, excludeNames);
			}
			fittingMolecules[i] = fit;
			if (fit)
				System.out.println("Using " + molecule.getName() + " to "
						+ template.getName() + " in mapping");
			else {
				System.out.println("Not Using " + molecule.getName() + " to "
						+ template.getName() + " in mapping");
			}
		}

		alignMolecules();
		overlaySummary();
		writeOverlay();

	}

	/**
	 * Print summary information about the overlay.
	 */
	public void overlaySummary() {
		System.out.println();
		System.out.println("RMS for fitting molecules " + getFittingMoleculeRms());
		System.out.println("RMS for all molecules " + getAllMoleculeRms());
		System.out.println("Used " + getnFittingAtoms() + " atom mappings in fitting");
		System.out.println("Total number of atom mappings " + getnAllAtoms());
		System.out.println();

		System.out.println("Molecules used in the fitting");
		for (int i = 0; i < moleculeOverlay.size(); i++) {
			if (!isFittingMolecule(i))
				continue;
			Molecule molecule = moleculeOverlay.get(i);
			System.out
					.println("RMS for " + molecule.getName() + " " + getMoleculeRmss(i));
		}

		System.out.println("Molecules not used in the fitting");
		for (int i = 0; i < moleculeOverlay.size(); i++) {
			if (isFittingMolecule(i))
				continue;
			Molecule molecule = moleculeOverlay.get(i);
			System.out
					.println("RMS for " + molecule.getName() + " " + getMoleculeRmss(i));
		}
	}

	/**
	 * Aligns molecules to the template
	 * 
	 */
	public void alignMolecules() {

		logger.debug("Aligning molecules");
		// Now fit the first overlay to the second
		double xVals[][] = new double[nAllAtoms][];
		double yVals[][] = new double[nAllAtoms][];

		int no = 0;
		for (int i = 0; i < moleculeOverlay.size(); i++) {
			if (!fittingMolecules[i])
				continue;
			Molecule molecule = moleculeOverlay.get(i);
			Molecule template = templateOverlay.get(i);
			HashMap<Integer, Integer> mapping = mappings.get(i);
			for (int no1 : mapping.keySet()) {
				int no2 = mapping.get(no1);
				xVals[no] = molecule.getCoord(no1);
				yVals[no] = template.getCoord(no2);
				no++;
			}
		}
		nFittingAtoms = no;

		double trans[][] = new double[4][4];
		Coord.leastSquaresFit(xVals, yVals, trans);

		// transform all molecule coordinates
		for (Molecule molecule : allMolecules) {
			for (double[] coord : molecule.getCoords()) {
				Coord.transPointInPlace(trans, coord);
			}
		}

		int allMoleculeNo = 0, fittingMoleculeNo = 0;

		// determine per molecule rms
		moleculeRmss.clear();
		List<Double> fittingMoleculeRmss = new ArrayList<>();

		for (int i = 0; i < moleculeOverlay.size(); i++) {
			Molecule molecule = moleculeOverlay.get(i);
			Molecule template = templateOverlay.get(i);

			HashMap<Integer, Integer> mapping = mappings.get(i);
			double moleculeSum = 0;
			for (int no1 : mapping.keySet()) {
				int no2 = mapping.get(no1);
				double sqrDistance = Coord.sqrDistance(molecule.getCoord(no1),
						template.getCoord(no2));
				moleculeSum += sqrDistance;

			}
			no = mapping.size();
			double moleculeRms = Math.sqrt(moleculeSum / (no));
			logger.debug("RMS for molecule " + molecule.getName() + " is " + moleculeRms);
			moleculeRmss.add(moleculeRms);

			allMoleculeNo += no;
			if (fittingMolecules[i]) {
				fittingMoleculeRmss.add(moleculeRms);
				fittingMoleculeNo += no;
			}

		}

		assert nFittingAtoms == fittingMoleculeNo : "fitting atom count mismatch";
		assert nAllAtoms == allMoleculeNo : "all atom count mismatch";

		// determine total RMS
		allMoleculeRms = moleculeRmss.stream().mapToDouble(Double::doubleValue).average()
				.getAsDouble();
		fittingMoleculeRms = fittingMoleculeRmss.stream().mapToDouble(Double::doubleValue)
				.average().getAsDouble();
		logger.debug("RMS for all molecules is " + allMoleculeRms
				+ ", fitting molecules olnly is " + fittingMoleculeRms);
	}

	/**
	 * Create the filename for the fitted overlay
	 * 
	 * @return
	 */
	String fittedMoleculeFilename() {
		// Write out the fitted overlay
		File file = new File(moleculeFile);
		String base = file.getName();
		String fittedFile = "Fitted_" + base;
		return fittedFile;
	}

	/**
	 * Writes out the molecules
	 * 
	 */
	public void writeOverlay() {
		// Write out the fitted overlay
		String fittedFile = fittedMoleculeFilename();
		System.out.println("Writing fitting to " + fittedFile);
		try {
			Molecule.writeMols(allMolecules, fittedFile);
		} catch (IOException e) {
			throw new RuntimeException(
					"IOException on writing fitted molecules to " + fittedFile);
		}
	}

	/**
	 * @param molecule
	 * @param names
	 * @return true if the name of the molecule contains any of the strings in
	 *         the names.
	 */
	private boolean matchName(Molecule molecule, String[] names) {
		String name = molecule.getName();
		for (String test : names) {
			if (name.contains(test)) {
				return true;
			}
		}
		return false;
	}

	private static Pattern conformerMatchPattern = Pattern
			.compile("^(cons_|avg_)?(.*?)((_\\d+)+)?(_mol)?\\z");
	// private static Pattern conformerMatchPattern2 =
	// Pattern.compile("^(cons_)?(.*)\\z");
	// private static Pattern conformerMatchPattern3 = Pattern
	// .compile("^(cons_)?(.*_\\d+)_mol\\z");

	private static Pattern gapeMatchPattern = Pattern
			.compile("^ga (ranked )?\\d+ (.*?)((_\\d+)+)?( base molecule)?\\z");

	String baseName(String molName) {
		molName = molName.toLowerCase();

		Matcher m = conformerMatchPattern.matcher(molName);
		if (m.find())
			molName = m.group(2);
		// m = conformerMatchPattern2.matcher(molName);
		// if (m.find())
		// molName = m.group(2);
		// m = conformerMatchPattern3.matcher(molName);
		// if (m.find())
		// molName = m.group(2);
		m = gapeMatchPattern.matcher(molName);
		if (m.find())
			molName = m.group(2);

		return molName;
	}

	/**
	 * @param molName
	 * @param templateName
	 * @return true if template and molecule names are equivalent
	 */
	boolean matchNames(String molName, String templateName) {

		logger.debug("matching " + molName + " to template " + templateName);

		molName = baseName(molName);
		templateName = baseName(templateName);

		return molName.equals(templateName);
	}

	/**
	 * Load in a file stripping out Pharmacophore molecules. For the overlay to
	 * be fitted save the pharmacophore molecules and rotate them with the
	 * overlay.
	 * 
	 * @param file
	 * @return
	 */
	private List<Molecule> loadFile(String file, boolean template) {
		if (!template)
			pharmMolecules = new ArrayList<Molecule>();

		List<Molecule> allMols = Molecule.loadFiles(new String[] { file });
		List<Molecule> molecules = new ArrayList<Molecule>();
		for (Molecule molecule : allMols) {
			if (molecule.getName().toLowerCase().contains("pharm")) {
				if (!template) {
					pharmMolecules.add(molecule);
				}
				continue;
			}
			if (init)
				molecule.init();
			molecules.add(molecule);
		}
		return molecules;
	}

	/**
	 * @return the number of molecules in the overlay
	 */
	public int nMolecules() {
		return moleculeOverlay.size();
	}

	/**
	 * @return the init
	 */
	public boolean isInit() {
		return init;
	}

	/**
	 * @return the subgraph
	 */
	public boolean isSubgraph() {
		return subgraph;
	}

	/**
	 * @return the mapNames
	 */
	public boolean isMapNames() {
		return mapNames;
	}

	/**
	 * @param init
	 *            the init to set
	 */
	public void setInit(boolean init) {
		this.init = init;
	}

	/**
	 * @param subgraph
	 *            the subgraph to set
	 */
	public void setSubgraph(boolean subgraph) {
		this.subgraph = subgraph;
	}

	/**
	 * @param mapNames
	 *            the mapNames to set
	 */
	public void setMapNames(boolean mapNames) {
		this.mapNames = mapNames;
	}

	/**
	 * @return the moleculeOverlay
	 */
	public List<Molecule> getMoleculeOverlay() {
		return moleculeOverlay;
	}

	/**
	 * @return the templateOverlay
	 */
	public List<Molecule> getTemplateOverlay() {
		return templateOverlay;
	}

	/**
	 * @param moleculeNo
	 * @return true if this molecule is used to create the alignment
	 */
	public boolean isFittingMolecule(int moleculeNo) {
		return fittingMolecules[moleculeNo];
	}

	/**
	 * Set whether to use a fitting molecule in the alignment
	 * 
	 * @param moleculeNo
	 * @param fittingMolecule
	 */
	public void setFittingMolecule(int moleculeNo, boolean fittingMolecule) {
		fittingMolecules[moleculeNo] = fittingMolecule;
	}

	/**
	 * @return the moleculeRmss
	 */
	public double getMoleculeRmss(int moleculeNo) {
		return moleculeRmss.get(moleculeNo);
	}

	/**
	 * @return the fittingMoleculeRms
	 */
	public double getFittingMoleculeRms() {
		return fittingMoleculeRms;
	}

	/**
	 * @return the allMoleculeRms
	 */
	public double getAllMoleculeRms() {
		return allMoleculeRms;
	}

	/**
	 * @return the nFittingAtoms
	 */
	public int getnFittingAtoms() {
		return nFittingAtoms;
	}

	/**
	 * @return the nAllAtoms
	 */
	public int getnAllAtoms() {
		return nAllAtoms;
	}

}
