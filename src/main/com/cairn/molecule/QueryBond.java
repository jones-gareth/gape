package com.cairn.molecule;

import com.cairn.molecule.BondType.Type;

public class QueryBond extends Bond {

	private Boolean queryInRing;

	// Hacks so that we can use this class for a query: doesn't fit
	// well with the rest of the bond representation

	private volatile boolean hasMultipleBondTypes;

	// these fields are set for queries with multiple bond types only
	private volatile boolean isSingle, isDouble, isTriple, isAromatic;

	public QueryBond(int no, Atom a1, Atom a2, Type t) {
		super(no, a1, a2, t);
	}

	/**
	 * Matches two bonds. Handles multiple bond types and treats amide bonds as
	 * single bonds.
	 * 
	 * @param bond2
	 * @return
	 */
	public boolean matchBond(Bond bond2) {

		if (queryInRing != null && queryInRing && !bond2.isInRing()) {
			return false;
		}
		if (queryInRing != null && !queryInRing && bond2.isInRing()) {
			return false;
		}

		BondType.Type type1 = getType().getType();
		BondType.Type type2 = bond2.getType().getType();
		if (type1 == BondType.Type.AM)
			type1 = BondType.Type.SINGLE;
		if (type2 == BondType.Type.AM)
			type2 = BondType.Type.SINGLE;

		if (type1 == type2)
			return true;

		if (hasMultipleBondTypes) {

			if (isSingle && type2 == BondType.Type.SINGLE)
				return true;

			if (isDouble && type2 == BondType.Type.DOUBLE)
				return true;

			if (isTriple && type2 == BondType.Type.TRIPLE)
				return true;

			if (isAromatic && type2 == BondType.Type.AR)
				return true;

		}

		return false;
	}

	/**
	 * @return the queryInRing
	 */
	public Boolean getQueryInRing() {
		return queryInRing;
	}

	/**
	 * @param queryInRing
	 *            the queryInRing to set
	 */
	public void setQueryInRing(Boolean queryInRing) {
		this.queryInRing = queryInRing;
	}

	/**
	 * @return the hasMultipleBondTypes
	 */
	public boolean isHasMultipleBondTypes() {
		return hasMultipleBondTypes;
	}

	/**
	 * @return the isSingle
	 */
	public boolean isSingle() {
		return isSingle;
	}

	/**
	 * @return the isDouble
	 */
	public boolean isDouble() {
		return isDouble;
	}

	/**
	 * @return the isTriple
	 */
	public boolean isTriple() {
		return isTriple;
	}

	/**
	 * @return the isAromatic
	 */
	public boolean isAromatic() {
		return isAromatic;
	}

	/**
	 * @param hasMultipleBondTypes
	 *            the hasMultipleBondTypes to set
	 */
	public void setHasMultipleBondTypes(boolean hasMultipleBondTypes) {
		this.hasMultipleBondTypes = hasMultipleBondTypes;
	}

	/**
	 * @param isSingle
	 *            the isSingle to set
	 */
	public void setSingle(boolean isSingle) {
		this.isSingle = isSingle;
	}

	/**
	 * @param isDouble
	 *            the isDouble to set
	 */
	public void setDouble(boolean isDouble) {
		this.isDouble = isDouble;
	}

	/**
	 * @param isTriple
	 *            the isTriple to set
	 */
	public void setTriple(boolean isTriple) {
		this.isTriple = isTriple;
	}

	/**
	 * @param isAromatic
	 *            the isAromatic to set
	 */
	public void setAromatic(boolean isAromatic) {
		this.isAromatic = isAromatic;
	}

}
