package com.cairn.molecule;

class TaffBond {
	private final AtomType.Type atom1;
	private final AtomType.Type atom2;
	private final BondType.Type bond;
	private final double len;
	private final double k;

	TaffBond(AtomType.Type a1, AtomType.Type a2, BondType.Type b, double l, double k) {
		atom1 = a1;
		atom2 = a2;
		bond = b;
		len = l;
		this.k = k;
	}

	/**
	 * @return the len
	 */
	public double getLen() {
		return len;
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
	 * @return the atom2
	 */
	public AtomType.Type getAtom2() {
		return atom2;
	}

	/**
	 * @return the bond
	 */
	public BondType.Type getBond() {
		return bond;
	}

}
