package com.cairn.molecule;

import java.util.List;

import com.cairn.molecule.AtomType.Type;

/**
 * An extension to Atom for use in queries.
 * 
 * @author Gareth Jones
 *
 */
public class QueryAtom extends Atom {
	// These are for query slns
	private volatile boolean fullyConnected, noHydrogens;
	private Boolean queryInRing;
	private volatile int nConnected;
	private volatile List<AtomType> is, not;
	private volatile Integer nSlnHydrogens;

	public QueryAtom(Molecule _molecule, int no, Type t) {
		super(_molecule, no, t);
	}

	/**
	 * @return the fullyConnected
	 */
	public boolean isFullyConnected() {
		return fullyConnected;
	}

	/**
	 * @return the noHydrogens
	 */
	public boolean isNoHydrogens() {
		return noHydrogens;
	}

	/**
	 * @return the nConnected
	 */
	public int getnConnected() {
		return nConnected;
	}

	/**
	 * @return the is
	 */
	public List<AtomType> getIs() {
		return is;
	}

	/**
	 * @return the not
	 */
	public List<AtomType> getNot() {
		return not;
	}

	/**
	 * @param fullyConnected
	 *            the fullyConnected to set
	 */
	void setFullyConnected(boolean fullyConnected) {
		this.fullyConnected = fullyConnected;
	}

	/**
	 * @param noHydrogens
	 *            the noHydrogens to set
	 */
	void setNoHydrogens(boolean noHydrogens) {
		this.noHydrogens = noHydrogens;
	}

	/**
	 * @param is
	 *            the is to set
	 */
	void setIs(List<AtomType> is) {
		this.is = is;
	}

	/**
	 * @param not
	 *            the not to set
	 */
	void setNot(List<AtomType> not) {
		this.not = not;
	}

	/**
	 * @param nConnected
	 *            the nConnected to set
	 */
	void setnConnected(int nConnected) {
		this.nConnected = nConnected;
	}

	/**
	 * @param nSlnHydrogens
	 *            the nSlnHydrogens to set
	 */
	void setnSlnHydrogens(Integer nSlnHydrogens) {
		this.nSlnHydrogens = nSlnHydrogens;
	}

	/**
	 * @return the nSlnHydrogens
	 */
	public Integer getnSlnHydrogens() {
		return nSlnHydrogens;
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

}
