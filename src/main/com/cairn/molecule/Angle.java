package com.cairn.molecule;

/**
 * Represents an angle in a compound
 * 
 * @author Gareth Jones
 *
 */
class Angle {
	private final Atom atom1;
	private final Atom mid;
	private final Atom atom2;
	private TaffAngle taff;

	Angle(Atom a1, Atom m, Atom a2) {
		atom1 = a1;
		atom2 = a2;
		mid = m;
	}

	public String info() {
		return "Angle " + atom1.getNo() + " " + mid.getNo() + " " + atom2.getNo();
	}

	/**
	 * @return the atom1
	 */
	public Atom getAtom1() {
		return atom1;
	}

	/**
	 * @return the mid
	 */
	public Atom getMid() {
		return mid;
	}

	/**
	 * @return the atom2
	 */
	public Atom getAtom2() {
		return atom2;
	}

	/**
	 * @return the taff
	 */
	public TaffAngle getTaff() {
		return taff;
	}

	/**
	 * @param taff
	 *            the taff to set
	 */
	void setTaff(TaffAngle taff) {
		this.taff = taff;
	}

}
