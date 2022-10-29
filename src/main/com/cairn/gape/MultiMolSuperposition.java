package com.cairn.gape;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.log4j.Logger;

import com.cairn.gape.molecule.GaMolecule;
import com.cairn.gape.molecule.MultiGaMolecule;
import com.cairn.molecule.Molecule;

/**
 * An extension to the superposition application where a conformer library is
 * used instead of fully flexible search.
 * 
 * @author Gareth Jones
 * 
 */
public class MultiMolSuperposition extends Superposition {
	private static Logger logger;
	static {
		logger = Logger.getLogger(MultiMolSuperposition.class);
		// logger.setLevel(Level.DEBUG);
	}

	// the length of the conformer integer string
	private volatile int conformerStringLength;
	private volatile int[] conformerStringRanges;

	// a lookup table for molecule number to position on integer string where
	// conformer number is stored
	private final Map<Integer, Integer> conformerEntryPoints = new HashMap<Integer, Integer>();
	private static final String VERSION = "0.0.5";

	/**
	 * Runs GAPE with multi-molecule conformer input. Requires a configuration
	 * file as an argument. Additional optional files contain input molecular
	 * structures. Just calls {@link #fitMolecules}.
	 * 
	 * @param args
	 *            filename arguments.
	 */
	public static void main(String args[]) {
		try {
			if (args.length < 1) {
				System.out.println("Usage: " + MultiMolSuperposition.class.getName()
						+ " -version | <conf file> [molecule files..]");
				System.exit(0);
			}
			if (args.length == 1 && args[0].equals("-version")) {
				System.out.println("version = " + VERSION);
				System.exit(1);
			}

			MultiMolSuperposition sp = new MultiMolSuperposition();
			sp.setupInfoOutput();
			System.out
					.print("\n\nGAPE [Multi Conformer Version]");
			System.out.print("Gareth Jones jones_gareth@outlook.com\n\n");

			String molFiles[] = ArrayUtils.remove(args, 0);

			sp.fitMolecules(args[0], molFiles);
			sp.infoMessageLogger.finish();
		} catch (Exception e) {
			logger.error("Fatal Exception ", e);
		}

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Superposition#fitMolecules(java.lang.String,
	 * java.lang.String[])
	 */
	@Override
	public void fitMolecules(String configFile, String[] molFiles) {
		super.fitMolecules(configFile, molFiles);
	}

	/**
	 * Loads in the multi-conformer molecules from a file;
	 * 
	 * @param molFiles
	 * @return
	 * 
	 */
	protected List<MultiGaMolecule> loadMultiGaMolecules(String[] molFiles) {
		List<MultiGaMolecule> m = MultiGaMolecule.loadMultiMolFiles(molFiles,
				Molecule.FileType.UNKNOWN, Molecule.Source.FILE, infoMessageLogger, true);
		fileFormat = Molecule.getType(molFiles[0]);
		return m;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Superposition#loadMolecules(java.lang.String[])
	 */
	@Override
	protected List<GaMolecule> loadMolecules(String[] molFiles) {
		List<GaMolecule> molecules = new ArrayList<>(loadMultiGaMolecules(molFiles));
		return molecules;
	}

	/**
	 * Initializes molecules for superpositon. Removes lone pairs. Checks to see
	 * if any molecules are rigid or fixed. Assigns any activities. Matches
	 * torsional distributions. Calls {@link #findBaseMolecule()} and
	 * {@link GaMolecule#setup()}. Finds binary string length and entry points.
	 */
	@Override
	public void setupMolecules(List<GaMolecule> m) {
		commonMolSetup();

		molecules = m;

		if (molecules.size() > MAX_MOLS) {
			infoMessageln("Maximum number of molecules that can be fitted is " + MAX_MOLS);
			infoMessageln("Fitting the first " + MAX_MOLS + " molecules");
			molecules = new ArrayList<>(molecules.subList(0, MAX_MOLS));
		}
		if (molecules.size() < 2) {
			throw new RuntimeException("GAPE needs at least 2 molecules: only "
					+ molecules.size() + " present");
		}

		for (GaMolecule molecule : molecules) {

			// Remove pharmacophore information
			infoMessageln("\n\nProcessing " + molecule.getName());
			molecule.removeLonePairs();

			molecule.setProblem(this);
			molecule.setup();

			if (isUseActivities()) {
				molecule.setActivity(getActivity(molecule.getName()));
				infoMessageln("pKa activity is " + molecule.getActivity());
			}

			if (tordist != null)
				tordist.moleculeMatchTordist(molecule);
		}

		findBaseMolecule();

		int nBits = 0;
		binaryEntryPoints = new int[molecules.size()];
		int i = 0;
		for (GaMolecule molecule : molecules) {
			binaryEntryPoints[i] = nBits;

			// note the molecule is not completely rigid as we need to rotate
			// hydrogens/lone pairs to fitting points.
			nBits += molecule.getConformationalBitLength();
			if (molecule.isRelaxMolecule() && molecule != fittingMolecule)
				nBits += 6 * 8;
			i++;
		}
		binaryStringLength = nBits;

		setupFeatureClustering();
	}

	/**
	 * Determine the length of the integer string portion. Add in space for
	 * conformers in a second integer string.
	 */
	@Override
	protected void determineIntegerStringLength() {
		super.determineIntegerStringLength();

		int no = 0;
		for (int i = 0; i < molecules.size(); i++) {
			MultiGaMolecule molecule = (MultiGaMolecule) molecules.get(i);
			if (molecule.nConformers() > 1) {
				int position = no;
				logger.debug("Conformer entry point for molecule " + i + " is position "
						+ position);
				conformerEntryPoints.put(i, position);
				no++;
			}
		}
		conformerStringLength = no;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Superposition#chromosomeClassName()
	 */
	@Override
	protected String chromosomeClassName() {
		return "com.cairn.gape.chromosome.MultiMolSuperpositionChromosome";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Superposition#findBaseMolecule()
	 */
	@Override
	public void findBaseMolecule() {
		super.findBaseMolecule();

		// augment information for integer string information to include
		// conformational data.
		int no = 0;
		conformerStringRanges = new int[conformerStringLength];
		for (GaMolecule molecule : molecules) {
			int nConformers = ((MultiGaMolecule) molecule).nConformers();
			if (nConformers > 1) {
				conformerStringRanges[no++] = nConformers;
			}
		}

	}

	/**
	 * @param moleculeNo
	 * @return the position on the integer string that stores the conformer
	 *         number for the molecule.
	 */
	public Integer getConformerEntryPoint(int moleculeNo) {
		return conformerEntryPoints.get(moleculeNo);
	}

	/**
	 * @return the conformerStringLength
	 */
	public int getConformerStringLength() {
		return conformerStringLength;
	}

	/**
	 * @return the conformerStringRanges
	 */
	public int[] getConformerStringRanges() {
		return conformerStringRanges;
	}

}