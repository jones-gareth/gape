package com.cairn.molecule;

/**
 * Forcefield parameters for an Atom
 * 
 * @author Gareth Jones
 *
 */
public class TaffAtom {
	private final AtomType.Type atom;
	private final double r;
	private final double k;

	TaffAtom(AtomType.Type a, double r, double k) {
		atom = a;
		this.k = k;
		this.r = r;
	}

	/**
	 * @return the atom
	 */
	public AtomType.Type getAtom() {
		return atom;
	}

	/**
	 * @return the r
	 */
	public double getR() {
		return r;
	}

	/**
	 * @return the k
	 */
	public double getK() {
		return k;
	}

}
