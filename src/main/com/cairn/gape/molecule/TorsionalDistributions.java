package com.cairn.gape.molecule;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.StringTokenizer;

import com.cairn.gape.utils.InfoMessageLogger;
import com.cairn.molecule.Fragments;
import com.cairn.molecule.Molecule;
import com.cairn.molecule.RotatableBond;

/**
 * Torsional distributions as used by GOLD/GASP etc
 * 
 * Torsional distribution routines
 * 
 * @author Gareth Jones
 * 
 */
public class TorsionalDistributions {
	private static final boolean DEBUG = false;

	private volatile int nDistributions = 0, nBins, maxPeakSize;

	// Table of torsional Distributions
	private volatile TorsionalDistribution distributions[];

	private volatile double granularity, deltaE = 2.5;

	// Tordist weights
	private volatile double elTypeWt, sybTypeWt, fragWt, hCntWt, linkageWt;

	private volatile boolean removeHighEnergy = false, useProbs = false,
			energyProbs = false;

	private volatile InfoMessageLogger infoMessageLogger;

	public TorsionalDistributions() {
		;
	}

	public TorsionalDistributions(InfoMessageLogger infoMessageLogger) {
		super();
		this.infoMessageLogger = infoMessageLogger;
	}

	// readTordistFile

	// Open the torsional distribution file and read it in.

	public void readTordistFile(String file) {
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(file)));
			infoMessageLogger.infoMessageln(2, "\nReading torsional distributions from "
					+ file);
			readTordistFile(in);
		} catch (IOException ex) {
			throw new RuntimeException("Failed to read " + file + " " + ex);
		}
	}

	public void readTordistFile() {
		try {
			InputStream stream = getClass().getResourceAsStream("TORDIST.txt");
			BufferedReader in = new BufferedReader(new InputStreamReader(stream));
			infoMessageLogger.infoMessageln(2,
					"\nReading default torsional distributions " + "from TORDIST.txt");
			readTordistFile(in);
		} catch (IOException ex) {
			throw new RuntimeException("Failed to read TORDIST.txt" + ex);
		}
	}

	private void readTordistFile(BufferedReader in) throws IOException {
		readTordistHeader(in);
		readTordistEntries(in);
	}

	// readTordistHeader

	// Read in the header of the torsional distributional file,
	// processing each line.
	private void readTordistHeader(BufferedReader in) throws IOException {
		boolean inHeader = false;

		String line = null;
		while ((line = in.readLine()) != null) {
			line = line.trim();
			if (line.startsWith("#") || line.equals(""))
				continue;
			if (line.startsWith("START_HEADER")) {
				inHeader = true;
				break;
			} else
				throw new RuntimeException("readTordistHeader: can't find header start");
		}
		if (!inHeader)
			throw new RuntimeException("readTordistHeader: can't find header start");

		while ((line = in.readLine()) != null) {
			line = line.trim();
			if (line.startsWith("#") || line.equals(""))
				continue;
			if (line.endsWith("END_HEADER"))
				return;
			processTordistHeaderLine(line);
		}
		throw new RuntimeException("readTordistHeader: can't find header end");
	}

	// processTordistHeaderLine

	// Process a header line in the torsional distribution file
	private void processTordistHeaderLine(String line) {

		StringTokenizer st = new StringTokenizer(line);
		String key = st.nextToken();
		String check = st.nextToken();
		if (!check.equals("="))
			throw new RuntimeException("processTordistHeaderLine: can't parse line: "
					+ line);
		String value = st.nextToken();

		// Maximum peak size, distributions with a maximum peak size
		// less than this will be discarded
		if (key.equals("MAX_PEAK_SIZE")) {
			maxPeakSize = Integer.valueOf(value).intValue();
			infoMessageLogger.infoMessageln(2,
					"Smallest size of largest peak in torsional " + "distributions "
							+ maxPeakSize);
			return;
		}

		// No of bins in the distribution
		if (key.equals("N_BINS")) {
			nBins = Integer.valueOf(value).intValue();
			if (nBins > TorsionalDistribution.MAX_BINS)
				throw new RuntimeException(
						"processTordistHeaderLine: maximum no of bins is "
								+ TorsionalDistribution.MAX_BINS);
			granularity = 360.0 / nBins;
			infoMessageLogger.infoMessageln(2,
					"No of bins in torsional distributions is " + nBins);
			return;
		}

		// Set to one to remove small high energy bars in the
		// distributions
		if (key.equals("REMOVE_HIGH_ENERGY")) {
			int rem = Integer.valueOf(value).intValue();
			removeHighEnergy = (rem == 1) ? true : false;
			if (removeHighEnergy)
				infoMessageLogger.infoMessageln(2, "Removing high energy torsions");
			return;
		}

		// Weight for sybyl atom type
		if (key.equals("SYB_TYPE_WT")) {
			sybTypeWt = Double.valueOf(value).doubleValue();
			return;
		}

		// This keyword isn't implemented yet: probabistic sampling
		// based on histograms.
		if (key.equals("USE_PROBS")) {
			int usep = Integer.valueOf(value).intValue();
			useProbs = (usep == 1) ? true : false;
			if (useProbs)
				infoMessageLogger.infoMessageln(2,
						"Using torsional histogram probabilities");
			return;
		}

		// Weight for fragment definition
		if (key.equals("FRAG_WT")) {
			fragWt = Double.valueOf(value).doubleValue();
			return;
		}

		// Delta E for energy-based bar removal
		if (key.equals("DELTA_E")) {
			deltaE = Double.valueOf(value).doubleValue();
			infoMessageLogger.infoMessageln(2, "Delta E set to " + deltaE);
			return;
		}

		// Weight for specifying a hydrogen count
		if (key.equals("H_CNT_WT")) {
			hCntWt = Double.valueOf(value).doubleValue();
			return;
		}

		// Weight for specifying a bond type
		if (key.equals("LINKAGE_WT")) {
			linkageWt = Double.valueOf(value).doubleValue();
			return;
		}

		// Weight for specifying an element rather than a specific
		// atom type.
		if (key.equals("EL_TYPE_WT")) {
			elTypeWt = Double.valueOf(value).doubleValue();
			return;
		}

		throw new RuntimeException(
				"processTordistHeaderLine: can't process header line:\n" + line);
	}

	/**
	 * Gets a new torsional distribution table entry.
	 * 
	 * @return
	 */
	private TorsionalDistribution getTordistTableEntry() {
		TorsionalDistribution td = new TorsionalDistribution(this);
		td.setEntryNo(nDistributions);
		nDistributions++;
		return td;
	}

	// readTordistEntries

	// Reads in all the torsional distributions from the file.
	void readTordistEntries(BufferedReader in) throws IOException {
		boolean inEntries = false;
		boolean endEntries = false;
		ArrayList<TorsionalDistribution> v = new ArrayList<TorsionalDistribution>();
		String line = null;

		while ((line = in.readLine()) != null) {
			if (line.startsWith("#") || line.equals(""))
				continue;
			if (!line.startsWith("START_TORSIONS"))
				throw new RuntimeException("readTordistHeader: can't find torsional "
						+ "distribution entries start");
			else {
				inEntries = true;
				break;
			}
		}
		if (!inEntries)
			throw new RuntimeException(
					"readTordistEntries: can't find torsional distribution "
							+ "entries start");

		while ((line = in.readLine()) != null) {
			if (line.startsWith("#") || line.equals(""))
				continue;
			if (line.startsWith("END_TORSIONS")) {
				endEntries = true;
				break;
			}

			// get a new torsional distribution data structure and fill
			// it up from the file

			TorsionalDistribution currentTordist = getTordistTableEntry();

			// The first line is the torsion name
			currentTordist.setTorsion(line);

			// The next line is the torsion definition
			line = in.readLine();
			currentTordist.readTordistFragmentEntry(line);

			// The third line is the torsional histogram
			while ((line = in.readLine()) != null) {
				if (line.startsWith("#") || line.equals(""))
					continue;
				else
					break;
			}
			currentTordist.readTordistHistogram(line);

			if (DEBUG) {
				System.out.print(currentTordist.info());
				System.out.print(currentTordist.histogramInfo());
			}
			v.add(currentTordist);
		}

		if (!endEntries)
			throw new RuntimeException(
					"readTordistHeader: can't find torsional distribution "
							+ "entries end");

		distributions = new TorsionalDistribution[nDistributions];
		for (int i = 0; i < nDistributions; i++)
			distributions[i] = v.get(i);
	}

	/**
	 * Check to see if any or the rotatable bond in molecule can be matched by
	 * any of the torsional distributions.
	 * 
	 * @param molecule
	 * @return
	 * 
	 */
	public int moleculeMatchTordist(Molecule molecule) {

		InfoMessageLogger infoMessageLogger = molecule.getInfoMessageLogger();

		Fragments.searchMolecule(molecule);

		int nMatched = 0;

		// loop over all rotatable bonds
		infoMessageLogger.infoMessageln(2,
				"Finding torsional distributions that occur in " + molecule.getName());
		for (RotatableBond rbond : molecule.getRotatableBonds()) {

			if (rbond.isFlipBond())
				continue;
			rbond.setTordistInfo(null);
			TorsionalDistribution matchedTordist = null;

			// TODO look again at this

			// We loop over all torsional distributions keeping the one
			// with the highest weight. On reflection I think that
			// this is wrong. Matching torsions should be "AND"ed
			// together. This requires a bit of work: we would now
			// have to loop over all torsions with the rotatable bond
			// in the center and develop code to AND two distributions
			// (taking account of the offset). This problem was
			// pointed out by the Zeneca GOLD Beta testers. Torsional
			// distributions with a fixed period would confuse things
			// here. I'd liked to have done this, but unfortunately
			// there isn't time.

			// loop over all torsional distributions
			for (int j = 0; j < distributions.length; j++) {
				TorsionalDistribution tordist = distributions[j];
				int match = 0;
				if ((match = tordist.matchTorsion(rbond, molecule)) != TorsionalDistribution.NO_MATCH) {
					// if we have a match then retain the distribution
					// with the highest weight
					if (rbond.getTordistInfo() == null
							|| rbond.getTordistInfo().getWeight() < tordist.getWeight()) {
						if (tordist.isMaxFailFlg()) {
							infoMessageLogger.infoMessageln(3,
									"WARNING: can't use tordist " + tordist.getTorsion()
											+ ":\nmaximum peak " + " size is less than "
											+ maxPeakSize);
							continue;
						}
						matchedTordist = tordist;
						rbond.setTordistInfo(new TordistInfo(rbond, tordist, molecule));
						if (match == TorsionalDistribution.REVERSE_MATCH)
							rbond.getTordistInfo().setReverse(true);
						else
							rbond.getTordistInfo().setReverse(false);
						rbond.getTordistInfo().setupTordistInfo();
					}
				}

			}

			if (rbond.getTordistInfo() != null) {

				if (infoMessageLogger.getLogLevel() > 3) {
					int no1 = rbond.getTordistInfo().getAtom1().getNo() + 1;
					int no2 = rbond.getTordistInfo().getAtom2().getNo() + 1;
					int no3 = rbond.getTordistInfo().getAtom3().getNo() + 1;
					int no4 = rbond.getTordistInfo().getAtom4().getNo() + 1;
					infoMessageLogger.infoMessageln("Rotatable bond [" + no1 + "," + no2
							+ "," + no3 + "," + no4 + "] matches torsion:");
				}
				infoMessageLogger.infoMessage(3, matchedTordist.info());
				nMatched++;
			}

		}

		infoMessageLogger.infoMessageln(2, nMatched + " torsional distributions matched");
		return nMatched;
	}

	public static void main(String[] args) {
		String torFile = args[0];
		String molFile = args[1];

		Molecule mol = new Molecule(molFile, Molecule.FileType.MOL2, Molecule.Source.FILE);
		mol.assignRotatableBonds();
		TorsionalDistributions td = new TorsionalDistributions();
		td.readTordistFile(torFile);
		td.moleculeMatchTordist(mol);

	}

	public TorsionalDistribution[] getDistributions() {
		return distributions;
	}

	public int getNBins() {
		return nBins;
	}

	public int getNDistributions() {
		return nDistributions;
	}

	protected double getDeltaE() {
		return deltaE;
	}

	protected int getMaxPeakSize() {
		return maxPeakSize;
	}

	protected double getHCntWt() {
		return hCntWt;
	}

	protected double getLinkageWt() {
		return linkageWt;
	}

	protected double getSybTypeWt() {
		return sybTypeWt;
	}

	protected double getElTypeWt() {
		return elTypeWt;
	}

	protected double getFragWt() {
		return fragWt;
	}

	protected double getGranularity() {
		return granularity;
	}

	/**
	 * @return the removeHighEnergy
	 */
	protected boolean isRemoveHighEnergy() {
		return removeHighEnergy;
	}

	/**
	 * @return the energyProbs
	 */
	protected boolean isEnergyProbs() {
		return energyProbs;
	}

	/**
	 * @param removeHighEnergy
	 *            the removeHighEnergy to set
	 */
	protected void setRemoveHighEnergy(boolean removeHighEnergy) {
		this.removeHighEnergy = removeHighEnergy;
	}

	/**
	 * @param energyProbs
	 *            the energyProbs to set
	 */
	protected void setEnergyProbs(boolean energyProbs) {
		this.energyProbs = energyProbs;
	}

}
