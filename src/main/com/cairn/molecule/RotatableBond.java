package com.cairn.molecule;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Supplier;

import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;
import com.cairn.gape.molecule.TordistInfo;

/**
 * A class to model a rotatable bond.
 * 
 * @author Gareth Jones
 *
 */
public class RotatableBond {
	// the bond
	private volatile Bond bond;
	// two atoms on bond
	private volatile Atom atom1, atom2;
	// torsions around bond
	private volatile List<Torsion> torsions;
	// List of all atoms on each side of bond
	private volatile List<Atom> atom1List, atom2List;
	private volatile Molecule molecule;
	private volatile boolean flipBond;
	// Any torsional distribution for the bond
	private volatile TordistInfo tordistInfo;

	private static final ThreadLocal<double[][]> trans = ThreadLocal
			.withInitial(() -> new double[4][4]);
	private static final Logger logger = Logger.getLogger(RotatableBond.class);

	protected RotatableBond() {
		;
	}

	private RotatableBond(Bond b, Molecule m) {
		this(b, m, new ArrayList<Torsion>(), false);
	}

	RotatableBond(Bond b, Molecule m, List<Torsion> t, boolean flip) {
		bond = b;
		atom1 = b.getAtom1();
		atom2 = b.getAtom2();
		molecule = m;
		atom1List = molecule.findNodes(atom1, bond);
		atom2List = molecule.findNodes(atom2, bond);
		// set things up so that atom1list has the smallest number of
		// atoms- on rotation we'll move these atoms
		if (atom2List.size() < atom1List.size()) {
			Atom tmpAtom = atom1;
			atom1 = atom2;
			atom2 = tmpAtom;
			List<Atom> tmpList = atom1List;
			atom1List = atom2List;
			atom2List = tmpList;
		}
		flipBond = flip;
		torsions = t;
	}

	public String info() {
		String rtn = null;
		if (flipBond)
			rtn = "Flip ";
		else
			rtn = "Free ";
		rtn += "Rotatable Bond (" + bond.info() + ") N torsions " + torsions.size();
		return rtn;
	}

	// bVal should be an 8bit binary number i.e 0 to 255
	public double rotateBond(int bVal) {
		double angle = .0;

		if (tordistInfo != null) {
			double r = bVal / 256.0;
			angle = tordistInfo.getTordistRotationAngle(r);
			debug(() -> "Tordist ");
		} else {
			angle = bVal * Math.PI / 128.0;
		}

		if (flipBond)
			angle = bVal > 127 ? Math.PI : 0.0;

		rotateBond(angle);
		return angle;
	}

	public void rotateBond(double rAngle) {
		// rAngle = Torsion.setAngle(rAngle);
		debug(() -> "Angle " + rAngle);
		List<double[]> coords = molecule.getCoords();

		Coord.determineRotation(coords.get(atom2.getNo()), coords.get(atom1.getNo()),
				rAngle, trans.get());
		debug(() -> "Rotation Matrix" + Coord.info(trans.get()));

		for (Atom atom : atom1List) {
			int no = atom.getNo();
			debug(() -> "Moving Point from " + Coord.info(coords.get(no)));
			Coord.transPointInPlace(trans.get(), coords.get(no));
			debug(() -> "Moved Point to " + Coord.info(coords.get(no)));
		}
	}

	public double torsionalEnergy(Taff taff) {
		double energy = torsions.stream().map(t -> taff.torsionEnergy(t))
				.reduce(0.0, (sum, tEnergy) -> sum + tEnergy);
		debug(() -> "Rotable bond torsion energy " + energy);
		return energy;
	}

	static void flattenBond(Molecule m, Atom a2, Atom a3) {
		Atom a1 = a2.getNotDummyNeighbours().get(0);
		if (a1 == a3) {
			a1 = a2.getNotDummyNeighbours().get(1);
		}
		Atom a4 = a3.getNotDummyNeighbours().get(0);
		if (a4 == a2) {
			a4 = a3.getNotDummyNeighbours().get(1);
		}
		flattenBond(m, a1, a2, a3, a4);
	}

	static void flattenBond(Molecule m, Atom a1, Atom a2, Atom a3, Atom a4) {
		Bond b = m.findBond(a2, a3);
		debug(() -> "Flattening bond " + a1.info() + b.info() + a4.info());

		RotatableBond rb = new RotatableBond(b, m);
		double angle = Coord.torsion(m.getCoord(a1.getNo()), m.getCoord(a2.getNo()),
				m.getCoord(a3.getNo()), m.getCoord(a4.getNo()));
		angle = Torsion.setAngle(angle);
		logger.debug("Torsion Angle is " + angle);

		// Keep trans/cis geometry
		if (angle > Math.PI - Math.PI / 3.0)
			angle -= Math.PI;
		if (angle < -Math.PI + Math.PI / 3.0)
			angle -= Math.PI;

		angle *= -1;
		rb.rotateBond(angle);
		if (logger.isDebugEnabled()) {
			logger.debug("Correction Angle is " + Torsion.setAngle(angle));
			angle = Coord.torsion(m.getCoord(a1.getNo()), m.getCoord(a2.getNo()),
					m.getCoord(a3.getNo()), m.getCoord(a4.getNo()));
			logger.debug("Torsion Angle after transform is " + Torsion.setAngle(angle));
		}
	}

	private static void debug(Supplier<String> message) {
		if (logger.isDebugEnabled()) {
			logger.debug(message.get());
		}
	}

	/**
	 * @return the atom1
	 */
	public Atom getAtom1() {
		return atom1;
	}

	/**
	 * @return the atom2
	 */
	public Atom getAtom2() {
		return atom2;
	}

	/**
	 * @param atom1
	 *            the atom1 to set
	 */
	protected void setAtom1(Atom atom1) {
		this.atom1 = atom1;
	}

	/**
	 * @param atom2
	 *            the atom2 to set
	 */
	protected void setAtom2(Atom atom2) {
		this.atom2 = atom2;
	}

	/**
	 * @return the bond
	 */
	public Bond getBond() {
		return bond;
	}

	/**
	 * @return the atom1List
	 */
	public List<Atom> getAtom1List() {
		return atom1List;
	}

	/**
	 * @return the atom2List
	 */
	public List<Atom> getAtom2List() {
		return atom2List;
	}

	/**
	 * @param bond
	 *            the bond to set
	 */
	protected void setBond(Bond bond) {
		this.bond = bond;
	}

	/**
	 * @param atom1List
	 *            the atom1List to set
	 */
	protected void setAtom1List(List<Atom> atom1List) {
		this.atom1List = atom1List;
	}

	/**
	 * @param atom2List
	 *            the atom2List to set
	 */
	protected void setAtom2List(List<Atom> atom2List) {
		this.atom2List = atom2List;
	}

	/**
	 * @param molecule
	 *            the molecule to set
	 */
	protected void setMolecule(Molecule molecule) {
		this.molecule = molecule;
	}

	/**
	 * @return the trans
	 */
	protected static ThreadLocal<double[][]> getTrans() {
		return trans;
	}

	/**
	 * @return the tordistInfo
	 */
	public TordistInfo getTordistInfo() {
		return tordistInfo;
	}

	/**
	 * @return the flipBond
	 */
	public boolean isFlipBond() {
		return flipBond;
	}

	/**
	 * @param tordistInfo
	 *            the tordistInfo to set
	 */
	public void setTordistInfo(TordistInfo tordistInfo) {
		this.tordistInfo = tordistInfo;
	}

}
