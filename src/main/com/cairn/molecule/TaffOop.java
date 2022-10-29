package com.cairn.molecule;

/**
 * Forcefield parameters for an out of plane atom.
 * 
 * @author Gareth Jones
 *
 */
class TaffOop {
	private final AtomType.Type atom;
	private final double k;

	TaffOop(AtomType.Type a, double k) {
		atom = a;
		this.k = k;
	}

	/**
	 * @return the atom
	 */
	public AtomType.Type getAtom() {
		return atom;
	}

	/**
	 * @return the k
	 */
	public double getK() {
		return k;
	}

}
