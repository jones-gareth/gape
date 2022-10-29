package com.cairn.molecule;

import java.awt.Color;

import org.apache.commons.math3.util.FastMath;

import com.cairn.common.utils.UtilColor;

/**
 * Class for representing atom types and their properties. Uses TAFF atom types
 * as a base.
 * 
 * @author Gareth Jones
 * 
 */
public class AtomType {
	private final Type type;

	private final String name;

	private final Geometry geometry;

	private final double radius, weight, gaussianWidth, neutralBondOrder,
			neutralBondOrder2;

	private final Color color;

	/**
	 * Enumerated list of atom type identifiers
	 */
	public enum Type {
		ATM_NONE,

		C3, C2, C1, CAR, CCAT,

		N3, N2, N1, NAR, NAM, NPL3, N4,

		O3, O2, OCO2, OAR,

		S3, S2, SO, SO2,

		P3, H, CL, BR,

		I, SI, LP, DU, NA, K, F, CA, LI, AL, ANY, MG, ZN, FE, MN,

		ATOM_USER1, ATOM_USER2, ATOM_USER3, ATOM_USER4, ATOM_USER5,

		ATOM_USER6, ATOM_USER7, ATOM_USER8, ATOM_USER9, ATOM_USER10,

		HET, HAL, HEV, WILD, DUC,

		C, N, P, S, OANY,

		X, Y, Z
	};

	/**
	 * Geometry types
	 */
	public enum Geometry {
		LIN, TRI, TET, CO, NONE
	};

	/**
	 * Static list of available atom types. Atom types are not contructed
	 * outside of this array
	 */
	static public final AtomType[] types = {
			new AtomType(Type.C3, "C.3", 1.7, Geometry.TET, 12.0, UtilColor.PART_WHITE,
					4, -1),
			new AtomType(Type.C2, "C.2", 1.7, Geometry.TRI, 12.0, UtilColor.PART_WHITE,
					4, -1),
			new AtomType(Type.C1, "C.1", 1.7, Geometry.LIN, 12.0, UtilColor.PART_WHITE,
					4, -1),
			new AtomType(Type.CAR, "C.ar", 1.7, Geometry.TRI, 12.0, UtilColor.PART_WHITE,
					4, -1),
			new AtomType(Type.CCAT, "C.cat", 1.7, Geometry.TRI, 12.0,
					UtilColor.PART_WHITE, 4, -1),
			new AtomType(Type.C, "C", 1.7, Geometry.NONE, 12.0, UtilColor.PART_WHITE, 4,
					-1),

			new AtomType(Type.N3, "N.3", 1.55, Geometry.TET, 14.0, UtilColor.BLUE, 3, -1),
			new AtomType(Type.N2, "N.2", 1.55, Geometry.TRI, 14.0, UtilColor.BLUE, 3, -1),
			new AtomType(Type.N1, "N.1", 1.55, Geometry.LIN, 14.0, UtilColor.BLUE, 3, -1),
			new AtomType(Type.NAR, "N.ar", 1.55, Geometry.TRI, 14.0, UtilColor.BLUE, 3,
					-1),
			new AtomType(Type.NAM, "N.am", 1.55, Geometry.TRI, 14.0, UtilColor.BLUE, 3,
					-1),
			new AtomType(Type.NPL3, "N.pl3", 1.55, Geometry.TRI, 14.0, UtilColor.BLUE, 3,
					-1),
			new AtomType(Type.N4, "N.4", 1.55, Geometry.TET, 14.0, UtilColor.BLUE, 3, -1),
			new AtomType(Type.N, "N", 1.55, Geometry.NONE, 14.0, UtilColor.BLUE, 3, -1),

			new AtomType(Type.O3, "O.3", 1.52, Geometry.TET, 16.0, Color.red, 2, -1),
			new AtomType(Type.O2, "O.2", 1.52, Geometry.TRI, 16.0, Color.red, 2, -1),
			new AtomType(Type.OCO2, "O.co2", 1.52, Geometry.TRI, 16.0, Color.red, 2, -1),
			new AtomType(Type.OAR, "O.ar", 1.52, Geometry.TRI, 16.0, Color.red, 2, -1),
			new AtomType(Type.OANY, "O", 1.52, Geometry.NONE, 16.0, Color.red, 2, -1),

			new AtomType(Type.S3, "S.3", 1.8, Geometry.TET, 32.0, UtilColor.YELLOW, 2, 4),
			new AtomType(Type.S2, "S.2", 1.8, Geometry.TRI, 32.0, UtilColor.YELLOW, 2, 4),
			new AtomType(Type.SO, "S.o", 1.7, Geometry.TET, 32.0, UtilColor.YELLOW, 2, 4),
			new AtomType(Type.SO2, "S.o2", 1.7, Geometry.TET, 32.0, UtilColor.YELLOW, 2,
					4),
			new AtomType(Type.S, "S", 1.7, Geometry.NONE, 32.0, UtilColor.YELLOW, 2, 4),

			new AtomType(Type.P3, "P.3", 1.8, Geometry.TET, 32.0, UtilColor.BROWN, 3, 5),
			new AtomType(Type.P, "P", 1.8, Geometry.NONE, 32.0, UtilColor.BROWN, 3, 5),

			new AtomType(Type.H, "H", 1.5, Geometry.LIN, 1.0, UtilColor.CYAN, 1, -1),
			new AtomType(Type.CL, "Cl", 1.75, Geometry.TET, 34.4, UtilColor.GREEN, 1, -1),
			new AtomType(Type.BR, "Br", 1.85, Geometry.TET, 79.9, UtilColor.BROWN, 1, -1),
			new AtomType(Type.I, "I", 1.98, Geometry.TET, 126.9, UtilColor.BROWN, 1, -1),
			new AtomType(Type.SI, "Si", 1.2, Geometry.TET, 28.1, UtilColor.BROWN, 1, -1),

			new AtomType(Type.LP, "Lp", 1.0, Geometry.LIN, 0, UtilColor.PURPLE, -1, -1),
			new AtomType(Type.DU, "Du", 1.0, Geometry.NONE, 0, UtilColor.PURPLE, -1, -1),
			new AtomType(Type.NA, "Na", 1.2, Geometry.CO, 23.0, UtilColor.PURPLE, -1, -1),
			new AtomType(Type.K, "K", 1.2, Geometry.CO, 39.1, UtilColor.PURPLE, -1, -1),
			new AtomType(Type.F, "F", 1.47, Geometry.TET, 19.0, UtilColor.ORANGE, -1, -1),

			new AtomType(Type.CA, "Ca", 1.11, Geometry.CO, 40.1, UtilColor.PURPLE, -1, -1),
			new AtomType(Type.LI, "Li", 1.2, Geometry.CO, 6.94, UtilColor.PURPLE, -1, -1),
			new AtomType(Type.AL, "Al", 1.2, Geometry.CO, 27.0, UtilColor.PURPLE, -1, -1),

			new AtomType(Type.MG, "Mg", 0.74, Geometry.CO, 24.3, UtilColor.PURPLE, -1, -1),
			new AtomType(Type.ZN, "Zn", 0.83, Geometry.CO, 65.4, UtilColor.PURPLE, -1, -1),
			new AtomType(Type.FE, "Fe", 0.83, Geometry.CO, 55.8, UtilColor.BROWN, -1, -1),
			new AtomType(Type.MN, "Mn", 0.89, Geometry.CO, 54.9, UtilColor.PURPLE, -1, -1),

			new AtomType(Type.HET, "Het", 0.0, Geometry.NONE, 0.0, UtilColor.PURPLE, -1,
					-1),
			new AtomType(Type.HAL, "Hal", 0.0, Geometry.NONE, 0.0, UtilColor.PURPLE, -1,
					-1),
			new AtomType(Type.HEV, "Hev", 0.0, Geometry.NONE, 0.0, UtilColor.PURPLE, -1,
					-1),
			new AtomType(Type.DUC, "Du.C", 0.0, Geometry.NONE, 0.0, UtilColor.PURPLE, -1,
					-1),

			// Query atom types
			new AtomType(Type.WILD, "*", 0.0, Geometry.NONE, 0.0, UtilColor.PURPLE, -1,
					-1),
			new AtomType(Type.ANY, "Any", 0.0, Geometry.NONE, 0.0, UtilColor.PURPLE, -1,
					-1),
			new AtomType(Type.X, "X", 0.0, Geometry.NONE, 0.0, UtilColor.PURPLE, -1, -1),
			new AtomType(Type.Y, "Y", 0.0, Geometry.NONE, 0.0, UtilColor.PURPLE, -1, -1),
			new AtomType(Type.Z, "Z", 0.0, Geometry.NONE, 0.0, UtilColor.PURPLE, -1, -1) };

	/**
	 * Returns an atom type given a Tripos atom type name
	 * 
	 * @param name
	 * @return
	 */
	public static AtomType sybType(String name) {
		for (int i = 0; i < types.length; i++) {
			if (types[i].name.toUpperCase().equals(name.toUpperCase()))
				return types[i];
		}
		if (!name.equals("Q"))
			System.err.println("No type for syb atom " + name);
		return sybType("Du");
	}

	/**
	 * Returns the Sybyl name of this type identifier.
	 * 
	 * @param type
	 * @return
	 */
	public static String sybName(Type type) {
		for (int i = 0; i < types.length; i++) {
			if (types[i].type == type)
				return types[i].name;
		}
		System.err.println("No type for type ID " + type);
		return sybType("Du").name;
	}

	/**
	 * Returns the AtomType class given a type identifier
	 * 
	 * @param type
	 * @return
	 */
	public static AtomType sybType(Type type) {
		for (int i = 0; i < types.length; i++) {
			if (types[i].type == type)
				return types[i];
		}
		System.err.println("No type for type ID " + type);
		return sybType("Du");
	}

	/**
	 * Returns the AtomType class given a type identifier. For sd files we have
	 * elemental type identifiers N, C, OANY, P and S for the complex hybridised
	 * Tripos types.
	 * 
	 * @param name
	 * @return
	 */
	public static AtomType sdfType(String name) {
		return sybType("Name");
	}

	/**
	 * Private contructor.
	 * 
	 * @param t
	 * @param n
	 * @param r
	 * @param g
	 * @param w
	 * @param c
	 * @param nbo
	 * @param nbo2
	 */
	private AtomType(Type t, String n, double r, Geometry g, double w, Color c, int nbo,
			int nbo2) {
		type = t;
		name = n;
		radius = r;
		geometry = g;
		weight = w;
		color = c;
		neutralBondOrder = nbo;
		neutralBondOrder2 = nbo2;
		gaussianWidth = gaussianWidth();
	}

	/**
	 * Tests to see if two atom types match. Checks wild card, X, Y Hev etc.
	 * Also N for N.2, N.3 etc and C, OANY, P, S. Aromatic types are not
	 * equivalent to sp2 types.
	 * 
	 * @param a2
	 * @return
	 */
	public boolean matchType(AtomType a2) {
		return matchType(a2, false);
	}

	public boolean matchElementalType(AtomType a2, boolean ignoreAromatic) {

		if (isOxygenType() && a2.isOxygenType())
			return true;
		if (isCarbonType() && a2.isCarbonType())
			return true;
		if (isNitrogenType() && a2.isNitrogenType())
			return true;
		if (isPhosphorousType() && a2.isPhosphorousType())
			return true;
		if (isSulphurType() && a2.isSulphurType())
			return true;

		return matchType(a2, ignoreAromatic);

	}

	/**
	 * Tests to see if two atom types match. Checks wild card, X, Y Hev etc.
	 * Also N for N.2, N.3 etc and C, OANY, P, S.
	 * 
	 * @param ignoreAromatic
	 *            make aromatic and sp2 atoms equivalent
	 * @return
	 */
	public boolean matchType(AtomType a2, boolean ignoreAromatic) {
		Type type2 = a2.type;

		if (type == Type.ANY || type2 == Type.ANY)
			return true;
		if (type == Type.WILD || type2 == Type.WILD)
			return true;

		if (type == Type.X && a2.isHeavy())
			return true;
		if (type == Type.Y && a2.isHeavy())
			return true;
		if (type == Type.Z && a2.isHeavy())
			return true;
		if (type2 == Type.X && isHeavy())
			return true;
		if (type2 == Type.Y && isHeavy())
			return true;
		if (type2 == Type.Z && isHeavy())
			return true;

		if (type == Type.OANY && a2.isOxygenType())
			return true;
		if (type2 == Type.OANY && isOxygenType())
			return true;
		if (type == Type.C && a2.isCarbonType())
			return true;
		if (type2 == Type.C && isCarbonType())
			return true;
		if (type == Type.N && a2.isNitrogenType())
			return true;
		if (type2 == Type.N && isNitrogenType())
			return true;
		if (type == Type.P && a2.isPhosphorousType())
			return true;
		if (type2 == Type.P && isPhosphorousType())
			return true;
		if (type == Type.S && a2.isSulphurType())
			return true;
		if (type2 == Type.S && isSulphurType())
			return true;

		if (type == Type.HEV && a2.isHeavy())
			return true;
		if (type2 == Type.HEV && isHeavy())
			return true;

		if (type == Type.HAL && (type2 == Type.CL || type2 == Type.BR || type2 == Type.F))
			return true;
		if (type2 == Type.HAL && (type == Type.CL || type == Type.BR || type == Type.F))
			return true;

		if (ignoreAromatic) {
			if (type == Type.O3 && type2 == Type.OAR)
				return true;
			if (type == Type.CAR && type2 == Type.C2)
				return true;
			if (type == Type.NAR && (type2 == Type.N2 || type2 == Type.NPL3))
				return true;
			if (type2 == Type.O3 && type == Type.OAR)
				return true;
			if (type2 == Type.CAR && type == Type.C2)
				return true;
			if (type2 == Type.NAR && (type == Type.N2 || type2 == Type.NPL3))
				return true;
		}
		if (type == type2)
			return true;

		return false;
	}

	/**
	 * Returns the elemental type. That is the non-hybridized type e.g. C for
	 * C.3.
	 * 
	 * @return
	 */
	public AtomType getElementalType() {
		AtomType type = this;

		if (type.isOxygenType())
			type = sybType(Type.OANY);
		if (type.isCarbonType())
			type = sybType(Type.C);
		if (type.isNitrogenType())
			type = sybType(Type.N);
		if (type.isPhosphorousType())
			type = sybType(Type.S);
		if (type.isSulphurType())
			type = sybType(Type.P);

		return type;
	}

	/**
	 * Matches atom types based on the elemental type.
	 * 
	 * @param a2
	 * @return
	 */
	public boolean matchElementalType(AtomType a2) {
		AtomType type1 = this.getElementalType();
		AtomType type2 = a2.getElementalType();

		return type1.matchType(type2);
	}

	/**
	 * Returns true if type is an oxygen
	 * 
	 * @return
	 */
	public boolean isOxygenType() {
		if (type == Type.O2 || type == Type.O3 || type == Type.OCO2 || type == Type.OANY
				|| type == Type.OAR)
			return true;
		return false;
	}

	/**
	 * @return true is this is a hetero atom type
	 */
	public boolean isHeteroType() {
		if (isOxygenType() || isNitrogenType() || isPhosphorousType() || isSulphurType())
			return true;
		if (type == Type.CL || type == Type.BR || type == Type.F)
			return true;
		return false;
	}

	/**
	 * @return true if this type is a halogen.
	 */
	public boolean isHalogen() {
		if (type == Type.CL || type == Type.BR || type == Type.F)
			return true;
		return false;
	}

	/**
	 * Returns true if type is a carbon
	 * 
	 * @return
	 */
	public boolean isCarbonType() {
		if (type == Type.C1 || type == Type.C2 || type == Type.C3 || type == Type.CAR
				|| type == Type.CCAT || type == Type.C)
			return true;
		return false;
	}

	/**
	 * Returns true if type is a sulphur
	 * 
	 * @return
	 */
	public boolean isSulphurType() {
		if (type == Type.S2 || type == Type.S3 || type == Type.SO || type == Type.SO2
				|| type == Type.S)
			return true;
		return false;
	}

	/**
	 * Returns true if type is a phosphorous
	 * 
	 * @return
	 */
	public boolean isPhosphorousType() {
		if (type == Type.P3 || type == Type.P)
			return true;
		return false;
	}

	/**
	 * Returns true if type is a nitrogen
	 * 
	 * @return
	 */
	public boolean isNitrogenType() {
		if (type == Type.N1 || type == Type.N2 || type == Type.N3 || type == Type.NAR
				|| type == Type.N4 || type == Type.NPL3 || type == Type.NAM
				|| type == Type.N)
			return true;
		return false;
	}

	/**
	 * Returns true if this is a heavy atom.
	 * 
	 * @return
	 */
	public boolean isHeavy() {
		if (type == Type.H || type == Type.DU || type == Type.LP)
			return false;
		return true;
	}

	/**
	 * Returns true is this is not a dummy atom or metal ion.
	 * 
	 * @return
	 */
	public boolean isNotDummy() {
		if (type == Type.DU || type == Type.NA || type == Type.K || type == Type.CA
				|| type == Type.LI || type == Type.MG || type == Type.ZN
				|| type == Type.FE || type == Type.MN || type == Type.LP)
			return false;
		return true;
	}

	public static final double ATOMIC_N = 2.7;

	/**
	 * Determines the alpha parameter (or Gaussian width paramterer) for an
	 * atomic Gaussian.
	 * 
	 * @return
	 */
	private double gaussianWidth() {
		double r = radius;
		double alpha = FastMath.pow((3.0 * ATOMIC_N) / (4 * Math.PI * r * r * r),
				2.0 / 3.0) * Math.PI;
		return alpha;
	}

	/**
	 * Return the SDF type string for this type.
	 * 
	 * @return
	 */
	public String sdType() {
		if (isCarbonType())
			return "C";
		if (isNitrogenType())
			return "N";
		if (isOxygenType())
			return "O";
		if (isPhosphorousType())
			return "P";
		if (isSulphurType())
			return "S";
		if (type == Type.ANY)
			return "A";
		if (type == Type.DU)
			return "Q";
		return name;
	}

	/**
	 * Returns true if this is an elemental type for which a hybridised type is
	 * available.
	 * 
	 * @return
	 */
	public boolean isElementalType() {
		switch (type) {
		case C:
		case OANY:
		case N:
		case P:
		case S:
			return true;
		default:
			return false;

		}
	}

	/**
	 * @return the type
	 */
	public Type getType() {
		return type;
	}

	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	/**
	 * @return the geometry
	 */
	public Geometry getGeometry() {
		return geometry;
	}

	/**
	 * @return the radius
	 */
	public double getRadius() {
		return radius;
	}

	/**
	 * @return the weight
	 */
	public double getWeight() {
		return weight;
	}

	/**
	 * @return the gaussianWidth
	 */
	public double getGaussianWidth() {
		return gaussianWidth;
	}

	/**
	 * @return the neutralBondOrder
	 */
	public double getNeutralBondOrder() {
		return neutralBondOrder;
	}

	/**
	 * @return the neutralBondOrder2
	 */
	public double getNeutralBondOrder2() {
		return neutralBondOrder2;
	}

	/**
	 * @return the color
	 */
	public Color getColor() {
		return color;
	}

}
