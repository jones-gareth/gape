package com.cairn.gape.molecule;

import java.util.StringTokenizer;

import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Logger;

import com.cairn.gape.molecule.TordistNode.ParentLinkage;
import com.cairn.molecule.Atom;
import com.cairn.molecule.Molecule;
import com.cairn.molecule.RotatableBond;

/**
 * A Torsional distribution.
 * 
 * Torsional distribution routines
 * 
 * @author Gareth Jones
 * 
 */
public class TorsionalDistribution {
	private final TorsionalDistributions distributions;

	private volatile String torsion;

	static final int MAX_BINS = 36;

	private volatile boolean maxFailFlg;

	private final TordistNode atom1, atom2, atom3, atom4;

	private final int count[];

	private volatile int entryNo;

	private volatile boolean expandFlg = false, periodFlg = false;

	private volatile double expandStart, expandEnd, periodStart, periodEnd, weight;

	private static final Logger logger = Logger.getLogger(TorsionalDistribution.class);

	public TorsionalDistribution(TorsionalDistributions td) {
		distributions = td;
		count = new int[distributions.getNBins()];
		atom1 = new TordistNode(distributions);
		atom2 = new TordistNode(distributions);
		atom3 = new TordistNode(distributions);
		atom4 = new TordistNode(distributions);
	}

	/**
	 * Reads in a torsional histogram line.
	 * 
	 * @param line
	 * 
	 */
	public void readTordistHistogram(String line) {
		int max = 0;
		maxFailFlg = false;
		StringTokenizer tk = new StringTokenizer(line);
		for (int i = 0; i < count.length; i++) {
			count[i] = Integer.valueOf(tk.nextToken()).intValue();
			if (count[i] > max)
				max = count[i];
		}

		if (max < distributions.getMaxPeakSize())
			maxFailFlg = true;

		// If necessary expand the histogram
		expandHistogram();
	}

	//
	/**
	 * Read in the torsional distribution definition.
	 * 
	 * @param line
	 * 
	 */
	protected void readTordistFragmentEntry(String line) {

		weight = .0;

		// process each of the 5 portions of the definition. The
		// definition is of the form node1 | node2 | node3 | node4 [|
		// directives

		StringTokenizer tk = new StringTokenizer(line, "|");

		/* process the 4 nodes */
		for (int entry_no = 1; entry_no < 5; entry_no++) {
			String entry = tk.nextToken();
			entry = entry.trim();
			TordistNode atomNode = null;
			switch (entry_no) {
			case 1:
				atomNode = atom1;
				break;
			case 2:
				atomNode = atom2;
				break;
			case 3:
				atomNode = atom3;
				break;
			case 4:
				atomNode = atom4;
				break;
			}

			atomNode.processNode(entry);
			// get the node weight
			weight += atomNode.getWeight();
		}

		/* process any directives */
		while (tk.hasMoreTokens()) {
			String entry = tk.nextToken();
			entry = entry.trim();
			processTorDirective(entry);
		}
	}

	/**
	 * Process torsional directives, these are expand and period
	 * 
	 * @param directive
	 * 
	 */
	private void processTorDirective(String directive) {
		if (directive.startsWith("expand"))
			processExpandDirective(directive);
		else if (directive.startsWith("period"))
			processPeriodDirective(directive);
		else
			throw new RuntimeException("process_tor_directive: unknown directive "
					+ directive);
	}

	/**
	 * Read an expand directive.
	 * 
	 * @param directive
	 */
	private void processExpandDirective(String directive) {
		expandFlg = true;
		StringTokenizer st = new StringTokenizer(directive);
		st.nextToken();
		expandStart = Double.valueOf(st.nextToken()).doubleValue();
		expandEnd = Double.valueOf(st.nextToken()).doubleValue();
	}

	/**
	 * Read a period directive.
	 * 
	 * @param directive
	 */
	private void processPeriodDirective(String directive) {
		periodFlg = true;
		StringTokenizer st = new StringTokenizer(directive);
		st.nextToken();
		periodStart = Double.valueOf(st.nextToken()).doubleValue();
		periodEnd = Double.valueOf(st.nextToken()).doubleValue();
	}

	private int angleToBinNo(double degrees) {
		double r = (degrees + 180) / distributions.getGranularity();
		return (int) r;
	}

	/**
	 * If the expand directive is set, part of the histogram should be expanded.
	 * In this case the distribution query is symmetric (so only a portion of
	 * the histogram was calculated), but the matching portion of the molecule
	 * may not be, so the distribution need to be expaned over 360deg.
	 * 
	 * 
	 */
	private void expandHistogram() {

		if (!expandFlg)
			return;

		if (logger.isDebugEnabled()) {
			logger.debug("expand start");
			logger.debug(histogramInfo());
		}

		double start = expandStart;
		double end = expandEnd;
		double diff = end - start;
		double step = distributions.getGranularity();

		// Here we expand by 180 degress.
		if (diff > 179.999 && diff < 180.001) {
			double sp = start;
			for (double angle = sp - step / 2; angle >= -180.0; angle -= step) {
				int b1 = angleToBinNo(angle);
				int b2 = angleToBinNo(2.0 * sp - angle);
				if (count[b1] == 0)
					count[b1] = count[b2];
			}
			sp = end;
			for (double angle = sp + step / 2; angle <= 180.0; angle += step) {
				int b1 = angleToBinNo(angle);
				int b2 = angleToBinNo(2.0 * sp - angle);
				if (count[b1] == 0)
					count[b1] = count[b2];
			}
		}

		// we can also expand over other periods, however the existing
		// histogram must start at 0deg
		else if (start > -.00001 && start < .00001) {
			double s = end;
			for (double angle = 0; angle < diff; angle += step) {
				double angle2 = angle + s;
				if (angle2 >= 180.0)
					angle2 = angle2 - 360.0;
				double angle3 = end - angle;
				int b1 = angleToBinNo(angle2);
				int b2 = angleToBinNo(angle3);
				if (count[b1] == 0)
					count[b1] = count[b2];
			}
			if (diff < 179.0) {
				if (!periodFlg)
					throw new RuntimeException("expandHistogram: symmetry is %f, "
							+ "but no period is set\n");
				if (periodStart < -.0001 || periodStart > .0001)
					throw new RuntimeException("expandHistogram: period starts at "
							+ periodStart + ", period must start at 0\n");
				double period = periodEnd;
				double d2 = 2.0 * diff;
				if (period < (d2 - .0001) || period > (d2 + .0001))
					throw new RuntimeException("expandHistogram: period sould be "
							+ period + " but is " + d2);
				double ip = period + 0.001;
				int nSections = 360 / ((int) ip);
				nSections--;
				for (int i = 0; i < nSections; i++) {
					s = period + i * period;
					for (int angle = 0; angle < period; angle += step) {
						double angle2 = angle + s;
						if (angle2 >= 180.0)
							angle2 = angle2 - 360.0;
						int b1 = angleToBinNo(angle2);
						int b2 = angleToBinNo(angle);
						if (count[b1] == 0)
							count[b1] = count[b2];
					}

				}
			}
		}

		else
			throw new RuntimeException("expand_histogram: expand range set to " + diff
					+ " expand start set to " + start
					+ " expand start must be 0 unless expand range is 180");

		if (logger.isDebugEnabled()) {
			logger.debug("expand end  ");
			logger.debug(histogramInfo());
		}
	}

	/**
	 * describe a tordist definition
	 * 
	 * @return
	 * 
	 */
	public String info() {
		// Print out the four toplevel nodes
		String rtn = "";
		if (StringUtils.isNotEmpty(torsion))
			rtn += torsion + "\n";
		rtn += atom1.info() + "|" + atom2.info() + "|" + atom3.info() + "|"
				+ atom4.info();

		// include any directives
		if (expandFlg)
			rtn += "| expand " + expandStart + " " + expandEnd;

		if (periodFlg)
			rtn += "| expand " + periodStart + " " + periodEnd;

		return rtn + "\n";
	}

	/**
	 * describe torsional distribution histogram
	 * 
	 * @return
	 */
	public String histogramInfo() {
		return histogramInfo(count);
	}

	public static String histogramInfo(int count[]) {
		String rtn = "";
		for (int i = 0; i < count.length; i++)
			rtn += count[i] + " ";
		return rtn + " ";
	}

	public static final int NO_MATCH = 0, MATCH = 1, REVERSE_MATCH = 2;

	private boolean atomNodesUsed[];

	private Atom no2Atom, no3Atom;

	/**
	 * Matches a rotatabale bond to a torsional distribution
	 * 
	 * @param rbond
	 * @param molecule
	 * @return
	 * 
	 */
	protected int matchTorsion(RotatableBond rbond, Molecule molecule) {
		// TODO- use pattern matching instead

		// use a depth first search for this so keep account of what
		// atoms are already matched in atom_nodes_used.
		atomNodesUsed = new boolean[molecule.getnAtoms()];

		Atom no2 = rbond.getAtom1();
		Atom no3 = rbond.getAtom2();
		no2Atom = no2;
		no3Atom = no3;
		atomNodesUsed[no2.getNo()] = atomNodesUsed[no3.getNo()] = true;
		int match = NO_MATCH;

		// Convert the 4 torsion nodes into two nodes by adding the
		// first to the neighbour list of the second and the fourth to
		// the third.
		atom2.addNeighbour(atom1);
		atom1.setParentLinkage(atom2.getParentLinkage());
		atom3.addNeighbour(atom4);

		// Now try to match the middle two nodes in the distribution
		// to the two atoms in the rotatable bond: there are two ways
		// of doing this.
		if (atom2.matchNode(null, no2, null, this)
				&& atom3.matchNode(null, no3, null, this))
			match = MATCH;
		else if (atom2.matchNode(null, no3, null, this)
				&& atom3.matchNode(null, no2, null, this))
			match = REVERSE_MATCH;

		// Unlink nodes 1 and 4 from 2 and 3
		atom2.removeLastNeighbour();
		atom3.removeLastNeighbour();
		atom1.setParentLinkage(ParentLinkage.PL_ANY);
		return match;
	}

	public TordistNode getAtom1() {
		return atom1;
	}

	public TordistNode getAtom2() {
		return atom2;
	}

	public TordistNode getAtom3() {
		return atom3;
	}

	public TordistNode getAtom4() {
		return atom4;
	}

	public double getWeight() {
		return weight;
	}

	public TorsionalDistributions getDistributions() {
		return distributions;
	}

	public int[] getCount() {
		return count;
	}

	public int getCount(int i) {
		return count[i];
	}

	protected boolean[] getAtomNodesUsed() {
		return atomNodesUsed;
	}

	protected boolean isAtomNodesUsed(int i) {
		return atomNodesUsed[i];
	}

	protected void setAtomNodesUsed(boolean[] atomNodesUsed) {
		this.atomNodesUsed = atomNodesUsed;
	}

	protected void setAtomNodesUsed(int i, boolean atomNodeUsed) {
		this.atomNodesUsed[i] = atomNodeUsed;
	}

	protected Atom getNo2Atom() {
		return no2Atom;
	}

	protected Atom getNo3Atom() {
		return no3Atom;
	}

	protected int getEntryNo() {
		return entryNo;
	}

	protected void setEntryNo(int entryNo) {
		this.entryNo = entryNo;
	}

	protected boolean isMaxFailFlg() {
		return maxFailFlg;
	}

	protected void setMaxFailFlg(boolean maxFailFlg) {
		this.maxFailFlg = maxFailFlg;
	}

	protected String getTorsion() {
		return torsion;
	}

	protected void setTorsion(String torsion) {
		this.torsion = torsion;
	}

}
