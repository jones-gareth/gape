package com.cairn.molecule;

/**
 * Class to represent bond types
 * 
 * @author Gareth Jones
 * 
 */
public class BondType {
	private final Type type;

	private final String name;

	public enum Type {
		UNKNOWN, SINGLE, DOUBLE, TRIPLE, AR, AM, UNK, DU, NC
	}

	// Array of available bond types- bond types should always come from here
	// and not via the constructor
	public static final BondType types[] = { new BondType(Type.UNKNOWN, "unknown"),
			new BondType(Type.SINGLE, "1"), new BondType(Type.DOUBLE, "2"),
			new BondType(Type.TRIPLE, "3"), new BondType(Type.AR, "ar"),
			new BondType(Type.AM, "am"), new BondType(Type.UNK, "un"),
			new BondType(Type.DU, "du"), new BondType(Type.NC, "nc") };

	/**
	 * Constructor, from bond type and name. This should never be called -
	 * BondTypes should always come from the types array
	 * 
	 * @param t
	 * @param n
	 */
	private BondType(Type t, String n) {
		type = t;
		name = n;
	}

	/**
	 * Returns bond type given SDF integer type.
	 * 
	 * @param type
	 * @return
	 */
	public static BondType sdfType(int type) {
		switch (type) {
		case 1:
			return sybType(Type.SINGLE);
		case 2:
			return sybType(Type.DOUBLE);
		case 3:
			return sybType(Type.TRIPLE);
		case 4:
			return sybType(Type.AR);
		}
		return sybType(Type.DU);
	}

	/**
	 * Returns bond type given enumerated type for bond
	 * 
	 * @param type
	 * @return
	 */
	public static BondType sybType(Type type) {
		for (int i = 0; i < types.length; i++) {
			if (types[i].type == type)
				return types[i];
		}
		System.err.println("No type for type ID " + type);
		return sybType(Type.DU);
	}

	/**
	 * Returns bond type given Tripos MOL2 name
	 * 
	 * @param n
	 * @return
	 */
	public static BondType sybType(String n) {
		for (int i = 0; i < types.length; i++) {
			if (types[i].name.equals(n))
				return types[i];
		}
		System.out.println("No syb type for bond " + n);
		return sybType("du");
	}

	/**
	 * Returns the Tipos MOL2 name for this bond.
	 * 
	 * @param type
	 * @return
	 */
	static String sybName(Type type) {
		return sybType(type).name;
	}

	/**
	 * Find the sdf integer type for this bond type.
	 * 
	 * @return
	 */
	public int sdType() {
		if (type == BondType.Type.SINGLE)
			return 1;
		if (type == BondType.Type.DOUBLE)
			return 2;
		if (type == BondType.Type.TRIPLE)
			return 3;
		if (type == BondType.Type.AR)
			return 4;
		if (type == BondType.Type.AM)
			return 1;
		if (type == BondType.Type.UNK)
			return 8;
		if (type == BondType.Type.DU)
			return 1;
		if (type == BondType.Type.NC)
			return 1;

		return 1;
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

}
