package com.cairn.gape.molecule;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;
import com.cairn.gape.utils.InfoMessageLogger;
import com.cairn.molecule.Atom;
import com.cairn.molecule.AtomType;
import com.cairn.molecule.Bond;
import com.cairn.molecule.BondType;
import com.cairn.molecule.Molecule;
import com.cairn.molecule.Ring;
import com.cairn.molecule.RotatableBond;

/**
 * A class that handles conformational analysis of large cycles. The cycle has a
 * bond broken, so as to open it and allow it to rotate. Constraints are used to
 * enforce ring closure.
 * 
 * @author Gareth Jones
 * 
 */
public class LargeCycle {

	private static Logger logger;
	static {
		logger = Logger.getLogger(LargeCycle.class);
		// logger.setLevel(Level.DEBUG);
	}
	private final GaMolecule molecule;

	private final Ring ring;

	private volatile Atom atom1, atom2;

	private volatile Atom dummy1, dummy2;

	private volatile Bond dummyBond1, dummyBond2;

	private static volatile boolean hideBreakOnWriting = true;

	private final InfoMessageLogger infoMessageLogger;

	/**
	 * Creates a Large cycle for a given ring.
	 * 
	 * @param m
	 * @param r
	 */
	LargeCycle(GaMolecule m, Ring r) {
		molecule = m;
		ring = r;
		infoMessageLogger = m.getInfoMessageLogger();
	}

	/**
	 * Check a ring and retrurns true if it's considered a large cycle. A large
	 * cycle has at least 9 atoms and three single bonds.
	 * 
	 * @param m
	 * @param ring
	 * @return
	 */
	static boolean isLargeCycle(GaMolecule m, Ring ring) {
		if (ring.getAtoms().size() <= 8)
			return false;
		int nSingle = 0;
		for (Bond b : ring.getBonds()) {
			if (b.getBondType() == BondType.Type.SINGLE && b.getRings().size() == 1)
				nSingle++;
		}

		if (logger.isDebugEnabled())
			logger.debug("Ring " + ring.info() + " n single bonds " + nSingle);

		if (nSingle <= 3)
			return false;

		return true;
	}

	/**
	 * Breaks the bond and replace with two dummy atoms that are used for ring
	 * closure constraints.
	 * 
	 */
	private void breakBond() {

		Bond best = null;
		int minNeighbours = Integer.MAX_VALUE;
		for (Bond bond : ring.getBonds()) {
			if (bond.getBondType() != BondType.Type.SINGLE)
				continue;
			if (bond.getRings().size() != 1)
				continue;
			int nNeigbours = bond.getAtom1().getnNotDummyNeighbours()
					+ bond.getAtom2().getnNotDummyNeighbours();
			if (nNeigbours < minNeighbours) {
				minNeighbours = nNeigbours;
				best = bond;
			}
		}

		atom1 = best.getAtom1();
		atom2 = best.getAtom2();

		infoMessageLogger.infoMessageln(
				2,
				"Breaking large cycle bond between " + atom1.info() + " and "
						+ atom2.info());

		dummy1 = new Atom(molecule, molecule.getnAtoms(), AtomType.Type.DU);
		dummy2 = new Atom(molecule, molecule.getnAtoms(), AtomType.Type.DU);
		double dummy1Coords[] = new double[4];
		double dummy2Coords[] = new double[4];

		Coord.copy(molecule.getCoord(atom2.getNo()), dummy1Coords);
		Coord.copy(molecule.getCoord(atom1.getNo()), dummy2Coords);

		molecule.deleteBond(best);
		molecule.addAtom(dummy1, dummy1Coords);
		molecule.addAtom(dummy2, dummy2Coords);

		dummyBond1 = new Bond(molecule.getnBonds(), atom1, dummy1, BondType.Type.SINGLE);
		dummyBond2 = new Bond(molecule.getnBonds() + 1, atom2, dummy2,
				BondType.Type.SINGLE);
		molecule.addBond(dummyBond1);
		molecule.addBond(dummyBond2);

	}

	/**
	 * Reforms the bond to allow writing of molecule files with bond reformed.
	 * 
	 */
	public void reformBond() {
		if (!hideBreakOnWriting)
			return;

		// Don't output dummy atoms and one dummy bond
		dummy1.setOutput(false);
		dummy2.setOutput(false);
		dummyBond2.setOutput(false);

		// Convert dummy bond to original bond
		dummyBond1.setAtom2(atom2);
	}

	/**
	 * Restores the broken bond once the molecule file is written.
	 * 
	 */
	public void restoreBond() {
		if (!hideBreakOnWriting)
			return;

		// Restore dummy bond
		dummyBond1.setAtom2(dummy2);

		// Restore output flags- though I don't think we really need to do this
		dummy1.setOutput(true);
		dummy2.setOutput(true);
		dummyBond2.setOutput(true);
	}

	/**
	 * Returns the sqr distance between the tow dummy atoms and the two atoms in
	 * the broken bond. To be used with a ring closure constraint.
	 * 
	 * @return
	 */
	double penalty() {
		double sqrDist1 = Coord.sqrDistance(molecule.getCoord(atom1.getNo()),
				molecule.getCoord(dummy2.getNo()));
		double sqrDist2 = Coord.sqrDistance(molecule.getCoord(atom2.getNo()),
				molecule.getCoord(dummy1.getNo()));

		return sqrDist1 + sqrDist2;
	}

	/**
	 * Identifies all large cycles in a molecule.
	 * 
	 * @param molecule
	 */
	static void findlargeCycles(GaMolecule molecule) {
		ArrayList<LargeCycle> largeCycles = new ArrayList<LargeCycle>();
		for (Ring ring : molecule.getRings()) {
			if (isLargeCycle(molecule, ring))
				largeCycles.add(new LargeCycle(molecule, ring));
		}

		for (LargeCycle cycle : largeCycles)
			cycle.breakBond();

		molecule.setLargeCycles(largeCycles);
	}

	/**
	 * Test routine to break cycles and randomize molecule.
	 * 
	 * @param args
	 */
	public static void main(String args[]) {
		if (args == null || args.length != 1) {
			System.out.println("Usage: LargeCycle <structure file>");
			System.exit(0);
		}
		LargeCycle.hideBreakOnWriting = false;
		String molFile = args[0];

		List<GaMolecule> mols = GaMolecule.loadFiles(new String[] { molFile },
				Molecule.FileType.UNKNOWN, Molecule.Source.FILE);
		for (GaMolecule mol : mols) {
			LargeCycle.findlargeCycles(mol);
			mol.update();
			mol.assignRotatableBonds(false);

			Random randomGenerator = new java.util.Random();
			for (RotatableBond rotatableBond : mol.getRotatableBonds()) {
				int bVal = randomGenerator.nextInt(256);
				System.out.println("bVal " + bVal);
				rotatableBond.rotateBond(bVal);
			}
		}
		GaMolecule.write(mols, "large_cycle_" + molFile,
				"Randomized with large cycles open");

	}

	public Atom getAtom1() {
		return atom1;
	}

	public Atom getAtom2() {
		return atom2;
	}

}
