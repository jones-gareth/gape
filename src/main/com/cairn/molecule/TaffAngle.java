package com.cairn.molecule;

/**
 * Forcfield parameters for an angle.
 * 
 * @author Gareth Jones
 *
 */
class TaffAngle {
	private final AtomType.Type atom1;
	private final AtomType.Type mid;
	private final AtomType.Type atom2;
	private volatile int weight;
	private final double angle;
	private final double k;

	TaffAngle(AtomType.Type a1, AtomType.Type a2, AtomType.Type a3, double a, double k) {
		atom1 = a1;
		mid = a2;
		atom2 = a3;
		angle = a;
		this.k = k;
		weight = 3;
		if (atom1 == AtomType.Type.WILD)
			weight--;
		if (mid == AtomType.Type.WILD)
			weight--;
		if (atom2 == AtomType.Type.WILD)
			weight--;
	}

	/**
	 * @return the weight
	 */
	public int getWeight() {
		return weight;
	}

	/**
	 * @return the angle
	 */
	public double getAngle() {
		return angle;
	}

	/**
	 * @return the k
	 */
	public double getK() {
		return k;
	}

	/**
	 * @return the atom1
	 */
	public AtomType.Type getAtom1() {
		return atom1;
	}

	/**
	 * @return the mid
	 */
	public AtomType.Type getMid() {
		return mid;
	}

	/**
	 * @return the atom2
	 */
	public AtomType.Type getAtom2() {
		return atom2;
	}

}
