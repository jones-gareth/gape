package com.cairn.molecule;

import org.apache.log4j.Logger;

/**
 * A class for torsion parameters in the Tripos force field.
 * 
 * @author Gareth Jones
 *
 */
class TaffTorsion {
	private final AtomType.Type atom1;
	private final AtomType.Type atom2;
	private final AtomType.Type atom3;
	private final AtomType.Type atom4;
	private final BondType.Type bond;
	private volatile int weight;
	private final double p;
	private final double k;
	private static final Logger logger = Logger.getLogger(TaffTorsion.class);

	static {
		// logger.setLevel(Level.DEBUG);
	}

	TaffTorsion(AtomType.Type a1, AtomType.Type a2, AtomType.Type a3, AtomType.Type a4,
			BondType.Type b, double k, double p) {
		atom1 = a1;
		atom2 = a2;
		atom3 = a3;
		atom4 = a4;
		this.p = p;
		this.k = k;
		bond = b;
		weight = 4;
		if (atom1 == AtomType.Type.WILD)
			weight--;
		if (atom2 == AtomType.Type.WILD)
			weight--;
		if (atom3 == AtomType.Type.WILD)
			weight--;
		if (atom4 == AtomType.Type.WILD)
			weight--;
	}

	double energy(double angle) {
		double sign = (p > 0) ? 1.0 : -1.0;
		double e = k * (1 + sign * Math.cos(p * angle));
		logger.debug("Taff torsion " + angle + " P " + p + " k " + k + " energy " + e);
		return e;
	}

	String info() {
		String type1 = AtomType.sybName(atom1);
		String type2 = AtomType.sybName(atom2);
		String type3 = AtomType.sybName(atom3);
		String type4 = AtomType.sybName(atom4);
		return "[" + type1 + "," + type2 + "," + type3 + "," + type4 + "] "
				+ BondType.sybName(bond) + " (" + p + "," + k + ")";
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
	 * @return the atom3
	 */
	public AtomType.Type getAtom3() {
		return atom3;
	}

	/**
	 * @return the atom4
	 */
	public AtomType.Type getAtom4() {
		return atom4;
	}

	/**
	 * @return the bond
	 */
	public BondType.Type getBond() {
		return bond;
	}

	/**
	 * @return the weight
	 */
	public int getWeight() {
		return weight;
	}

	/**
	 * @return the p
	 */
	public double getP() {
		return p;
	}

	/**
	 * @return the k
	 */
	public double getK() {
		return k;
	}

}
