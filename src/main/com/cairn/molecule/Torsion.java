package com.cairn.molecule;

import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;

/**
 * Represents a 4-atom bonded torsion
 * 
 * @author Gareth Jones
 *
 */
public class Torsion {
	private static final Logger logger = Logger.getLogger(Torsion.class);

	static {
		// logger.setLevel(Level.DEBUG);
	}

	private final Molecule molecule;
	private final Atom atom1;
	private final Atom atom2;
	private final Atom atom3;
	private final Atom atom4;
	private final Bond bond;
	private volatile TaffTorsion taff;
	private volatile boolean reverse = false;
	private volatile double energy;

	Torsion(Molecule m, Bond b, Atom a1, Atom a2, Atom a3, Atom a4) {
		molecule = m;
		bond = b;
		atom1 = a1;
		atom2 = a2;
		atom3 = a3;
		atom4 = a4;
	}

	double determineAngle() {
		double angle = Coord.torsion(molecule.getCoord(atom1.getNo()),
				molecule.getCoord(atom2.getNo()), molecule.getCoord(atom3.getNo()),
				molecule.getCoord(atom4.getNo()));
		return angle;
	}

	static double setAngle(double angle) {
		double pi = Math.PI;
		while (angle > pi)
			angle -= 2.0 * pi;
		while (angle <= -pi)
			angle += 2.0 * pi;
		return angle;
	}

	// public double energy() {
	// energy = energy(determineAngle());
	// return energy;
	// }

	public double energy(double angle) {
		energy = taff.energy(angle);
		logger.debug("Torsion " + info() + " angle " + angle + " energy " + energy);
		return energy;
	}

	public String info() {
		int no1 = atom1.getNo() + 1;
		int no2 = atom2.getNo() + 1;
		int no3 = atom3.getNo() + 1;
		int no4 = atom4.getNo() + 1;
		int bno = bond.getNo() + 1;

		String type1 = atom1.getType().getName();
		String type2 = atom2.getType().getName();
		String type3 = atom3.getType().getName();
		String type4 = atom4.getType().getName();
		String btype = bond.getType().getName();

		return "[ (" + no1 + "," + type1 + ") " + "(" + no2 + "," + type2 + ") " + "("
				+ no3 + "," + type3 + ") " + "(" + no4 + "," + type4 + ") ] " + bno + " "
				+ btype;
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
	 * @return the atom3
	 */
	public Atom getAtom3() {
		return atom3;
	}

	/**
	 * @return the atom4
	 */
	public Atom getAtom4() {
		return atom4;
	}

	/**
	 * @return the bond
	 */
	public Bond getBond() {
		return bond;
	}

	/**
	 * @return the energy
	 */
	public double getEnergy() {
		return energy;
	}

	/**
	 * @return the taff
	 */
	public TaffTorsion getTaff() {
		return taff;
	}

	/**
	 * @return the reverse
	 */
	public boolean isReverse() {
		return reverse;
	}

	/**
	 * @return the molecule
	 */
	public Molecule getMolecule() {
		return molecule;
	}

	/**
	 * @param taff
	 *            the taff to set
	 */
	void setTaff(TaffTorsion taff) {
		this.taff = taff;
	}

	/**
	 * @param reverse
	 *            the reverse to set
	 */
	void setReverse(boolean reverse) {
		this.reverse = reverse;
	}

}
