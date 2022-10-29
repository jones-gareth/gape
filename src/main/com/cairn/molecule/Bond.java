package com.cairn.molecule;

import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

import com.cairn.molecule.BondType.Type;

/**
 * Class for representing atomic bonds
 * 
 * @author Gareth Jones
 * 
 */
public class Bond {
	private static final Logger logger = Logger.getLogger(Bond.class);
	// TODO subclass query bond

	private volatile Atom atom1, atom2;

	// set output true to include this bond in output structures
	private volatile boolean output = true;

	private volatile BondType type;

	private volatile String name;

	private volatile List<Ring> rings = new ArrayList<>();

	// this is used by queries and the ring finder algorithm
	private volatile boolean inRing;

	private volatile int no;

	private volatile TaffBond taff;

	public enum Stereo {
		NONE, UP, DOWN, EITHER, CIS_TRANS
	};

	private Stereo stereo = Stereo.NONE;

	/**
	 * Bonds can be fully rotatable, planar (0/180) or fixed.
	 */
	public enum Rotatable {
		NO, FULL, FLIP
	};

	private volatile Rotatable rotatable = Rotatable.NO;

	private volatile boolean canFlatten;

	/**
	 * Constructor
	 */
	public Bond() {
		output = true;
	}

	/**
	 * Bond constructor. The bond is between two atoms and the name is the
	 * Tripos bond type.
	 * 
	 * @param no
	 * @param a1
	 * @param a2
	 * @param n
	 */
	public Bond(int no, Atom a1, Atom a2, String n) {
		this();
		atom1 = a1;
		atom2 = a2;
		this.no = no;
		name = n;
		type = BondType.sybType(n);
	}

	/**
	 * Bond constructor. The bond is between two atoms with bond type t.
	 * 
	 * @param no
	 * @param a1
	 * @param a2
	 * @param t
	 */
	public Bond(int no, Atom a1, Atom a2, BondType.Type t) {
		this();
		atom1 = a1;
		atom2 = a2;
		this.no = no;
		type = BondType.sybType(t);
	}

	/**
	 * Bond constructor. The bond is between two atoms and t is the integer sdf
	 * bond type.
	 * 
	 * @param no
	 * @param a1
	 * @param a2
	 * @param t
	 */
	public Bond(int no, Atom a1, Atom a2, int t) {
		atom1 = a1;
		atom2 = a2;
		this.no = no;
		type = BondType.sdfType(t);
	}

	/**
	 * For use in creating a new molecule. Note that the bond atoms are in the
	 * old molecule so we need to get the new atoms from the new structure.
	 * 
	 * @param molecule
	 *            New molecule.
	 */
	public Bond copy(Molecule m) {
		Bond bond = new Bond();
		if (m != null) {
			bond.atom1 = m.getAtom(atom1.getNo());
			bond.atom2 = m.getAtom(atom2.getNo());
		}
		bond.name = name;
		bond.no = no;
		bond.type = type;
		bond.stereo = stereo;
		return bond;
	}

	/**
	 * Sets bond stereo chemistry. Bond can be up, down (or either). Used for SD
	 * files.
	 * 
	 * @param s
	 */
	public void setStereo(Stereo s) {
		stereo = s;
		// System.out.println("Setting stereo to "+stereo);
	}

	/**
	 * @return stereo setting
	 */
	public Stereo getStereo() {
		return stereo;
	}

	/**
	 * Generates a informative string for this bond.
	 * 
	 * @return
	 */
	public String info() {
		return atom1.info() + " Bond [" + type.getName() + "] " + String.valueOf(no + 1)
				+ " " + atom2.info();
	}

	/**
	 * Determines if this bond is rotatable or not. Set flipAmideBonds to allow
	 * amide bonds to flip (0/180) as rotatable bonds
	 * 
	 * @param flipAmideBonds
	 * @return
	 */
	Rotatable isRotatable(boolean flipAmideBonds) {
		AtomType.Type type1 = atom1.getType().getType();
		AtomType.Type type2 = atom2.getType().getType();

		rotatable = Rotatable.NO;

		if (isInRing())
			return rotatable;

		if (type.getType() != BondType.Type.SINGLE) {
			if (type.getType() == BondType.Type.AM)
				canFlatten = true;
			if (type.getType() == BondType.Type.AM && flipAmideBonds)
				rotatable = Rotatable.FLIP;

			return rotatable;
		}

		if (isTerminal())
			return rotatable;

		// bonds between C.2 and N.pl3 don't normally rotate
		// but they can flip by 180 degrees
		// However, N.pl3 linking rings is unlikely to be able to stay planar
		if (type1 == AtomType.Type.NPL3 && !atom1.isInRing()
				&& (type2 == AtomType.Type.C2 || type2 == AtomType.Type.CAR)) {
			// if (atom1.Npl3Link) {
			if (false) {
				rotatable = Rotatable.FULL;
			} else {
				canFlatten = true;
				if (atom1.isPlanarNH2group())
					rotatable = Rotatable.NO;
				else
					rotatable = Rotatable.FLIP;
			}
		} else if (type2 == AtomType.Type.NPL3 && !atom2.isInRing()
				&& (type1 == AtomType.Type.C2 || type1 == AtomType.Type.CAR)) {
			// if (atom2.Npl3Link) {
			if (false) {
				rotatable = Rotatable.FULL;
			} else {
				canFlatten = true;
				if (atom2.isPlanarNH2group())
					rotatable = Rotatable.NO;
				else
					rotatable = Rotatable.FLIP;
			}
		}

		// =N-N= is not rotatable
		else if (type1 == AtomType.Type.N2 && !atom1.isInRing()
				&& type2 == AtomType.Type.N2 && !atom2.isInRing()) {
			canFlatten = true;
			return rotatable;
		}

		// COOH is planar
		else if (type1 == AtomType.Type.C2 && type2 == AtomType.Type.O3
				&& atom1.isCOOHgroup() && atom2.isCOOHgroup()) {
			canFlatten = true;
			rotatable = Rotatable.FLIP;
		} else if (type2 == AtomType.Type.C2 && type1 == AtomType.Type.O3
				&& atom2.isCOOHgroup() && atom1.isCOOHgroup()) {
			canFlatten = true;
			rotatable = Rotatable.FLIP;
		}

		else
			rotatable = Rotatable.FULL;

		return rotatable;
	}

	/**
	 * Determines if this bond is a terminal structure bond.
	 * 
	 * @return
	 */
	public boolean isTerminal() {
		if (atom1.getnNeighbours() == 2 && atom1.getnLonePairs() == 1)
			return true;
		if (atom2.getnNeighbours() == 2 && atom2.getnLonePairs() == 1)
			return true;
		if (atom1.getnNeighbours() == 1)
			return true;
		if (atom2.getnNeighbours() == 1)
			return true;
		return false;
	}

	/**
	 * Find the correct bond type for this bond. The atom types must be properly
	 * set before calling this.
	 * 
	 * @return
	 */
	BondType.Type checkBondType() {

		if (atom1.getAtomType() == AtomType.Type.CCAT
				|| atom1.getAtomType() == AtomType.Type.CCAT) {
			logger.debug("Checking CCAT");
		}
		if (atom1.getType().getType() == AtomType.Type.LP
				|| atom2.getType().getType() == AtomType.Type.LP)
			return BondType.Type.SINGLE;

		// Acid groups- assign one double and one single bond
		if (atom1.getType().getType() == AtomType.Type.OCO2)
			return oco2Bond(atom1, atom2);
		if (atom2.getType().getType() == AtomType.Type.OCO2)
			return oco2Bond(atom2, atom1);

		if (atom1.getType().getType() == AtomType.Type.NPL3
				&& atom2.getType().getType() == AtomType.Type.CCAT)
			return ccatBond(atom2, atom1);
		if (atom2.getType().getType() == AtomType.Type.NPL3
				&& atom1.getType().getType() == AtomType.Type.CCAT)
			return ccatBond(atom1, atom2);

		if (atom1.getType().getType() == AtomType.Type.O2)
			return BondType.Type.DOUBLE;
		if (atom2.getType().getType() == AtomType.Type.O2)
			return BondType.Type.DOUBLE;

		if (atom1.getType().getType() == AtomType.Type.NAM && atom2.isAmideCarbon())
			return BondType.Type.AM;
		if (atom2.getType().getType() == AtomType.Type.NAM && atom1.isAmideCarbon())
			return BondType.Type.AM;
		if (atom1.getType().getType() == AtomType.Type.NAM
				&& atom2.isSulphonamideSulphur())
			return BondType.Type.SINGLE;
		if (atom2.getType().getType() == AtomType.Type.NAM
				&& atom1.isSulphonamideSulphur())
			return BondType.Type.SINGLE;

		return getType().getType();
	}

	/**
	 * Finds bondtype for charged bonds in O.co2 acids. We use the sd bond model
	 * where one oxygen has a double bond and the other has a single. Works OK
	 * in most cases, but the routine setChargedAcidGroups in Atom is best for
	 * this.
	 * 
	 * @param oco2Atom
	 * @param linkAtom
	 * @return
	 * 
	 * @see Atom#setChargedAcidGroups()
	 */
	private BondType.Type oco2Bond(Atom oco2Atom, Atom linkAtom) {
		Atom otherOco2Atom = null;
		for (Atom neighbour : linkAtom.getNotDummyNeighbours()) {
			if (neighbour != oco2Atom && neighbour.getType().isOxygenType()
					&& neighbour.getnNotDummyNeighbours() == 1)
				otherOco2Atom = neighbour;
		}

		if (otherOco2Atom == null) {
			return BondType.Type.DOUBLE;
		}

		// if one atom has a formal charge use that
		if (oco2Atom.getFormalCharge() != null && oco2Atom.getFormalCharge() == -1)
			return BondType.Type.SINGLE;
		if (otherOco2Atom.getFormalCharge() != null
				&& otherOco2Atom.getFormalCharge() == -1)
			return BondType.Type.DOUBLE;

		// Otherwise the first atom has the double bond to it
		if (oco2Atom.getNo() < otherOco2Atom.getNo())
			return BondType.Type.DOUBLE;
		else
			return BondType.Type.SINGLE;
	}

	/**
	 * Finds bondtype for charged bonds in guanidinium bases. We use the sd bond
	 * model where one nitrogen has a double bond and the others are single.
	 * 
	 * @param ccatAtom
	 * @param nAtom
	 * @return
	 */
	private BondType.Type ccatBond(Atom ccatAtom, Atom nAtom) {
		Atom firstNAtom = null, formalChargeAtom = null;

		// The atom with the double bond is the one with the formal charge or
		// the first nitrogen with two hydrogens.
		for (Atom neighbour : ccatAtom.getNotDummyNeighbours()) {
			if (neighbour.getType().isNitrogenType()) {
				if (neighbour.getFormalCharge() != null
						&& neighbour.getFormalCharge() == 1) {
					formalChargeAtom = neighbour;
					break;
				}
				int no = neighbour.getHydrogenCount();
				if (no == 2) {
					if (firstNAtom == null)
						firstNAtom = neighbour;
					else if (neighbour.getNo() < firstNAtom.getNo())
						firstNAtom = neighbour;
				}
			}
		}

		if (formalChargeAtom == null && firstNAtom == null)
			return null;
		else if (formalChargeAtom == nAtom) {
			logger.debug("Setting CCAT bond from " + ccatAtom.info() + " to "
					+ nAtom.info() + " to double ");
			return BondType.Type.DOUBLE;
		} else if (formalChargeAtom == null && firstNAtom == nAtom) {
			logger.debug("Setting CCAT bond from " + ccatAtom.info() + " to "
					+ nAtom.info() + " to double ");
			return BondType.Type.DOUBLE;
		} else {
			logger.debug("Setting CCAT bond from " + ccatAtom.info() + " to "
					+ nAtom.info() + " to single ");
			return BondType.Type.SINGLE;
		}

	}

	/**
	 * Returns the equivalent SYBYL bond type
	 * 
	 * @param molecule
	 * @return
	 */
	public Type sybylBondType(Molecule molecule) {
		if (atom1.getType().getType() == AtomType.Type.OCO2
				&& atom2.getType().isNotDummy()) {
			return Type.AR;
		}
		if (atom2.getType().getType() == AtomType.Type.OCO2
				&& atom1.getType().isNotDummy()) {
			return Type.AR;
		}
		if (atom1.getType().getType() == AtomType.Type.CCAT
				&& atom2.getType().isNitrogenType()) {
			return Type.AR;
		}
		if (atom2.getType().getType() == AtomType.Type.CCAT
				&& atom1.getType().isNitrogenType()) {
			return Type.AR;
		}

		return type.getType();
	}

	/**
	 * Find the sdf integer type for this bond.
	 * 
	 * @return
	 */
	public int sdType() {
		return type.sdType();
	}

	/**
	 * Adds a ring to this bond
	 * 
	 * @param r
	 */
	public void addRing(Ring r) {
		// System.out.println("Adding ring "+r.info()+" to bond "+info()+" no
		// "+String.valueOf(nRings+1));
		rings.add(r);
	}

	/**
	 * A convenience method to get the underlying bond type.
	 * 
	 * @return
	 */
	public Type getBondType() {
		return type.getType();
	}

	/**
	 * Returns true if this bond is in a ring.
	 * 
	 * @return
	 */
	public boolean isInRing() {
		return inRing || rings.size() > 0;
	}

	/**
	 * @param inRing
	 *            the inRing to set
	 */
	void setInRing(boolean inRing) {
		this.inRing = inRing;
	}

	/**
	 * Remove all ring information
	 */
	void setNotInRing() {
		rings.clear();
		inRing = false;
	}

	/**
	 * @return the atom1
	 */
	public Atom getAtom1() {
		return atom1;
	}

	/**
	 * @return the atom2
	 */
	public Atom getAtom2() {
		return atom2;
	}

	/**
	 * @return the output
	 */
	public boolean isOutput() {
		return output;
	}

	/**
	 * @return the type
	 */
	public BondType getType() {
		return type;
	}

	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	/**
	 * @return the rings
	 */
	public List<Ring> getRings() {
		return rings;
	}

	/**
	 * @return the no
	 */
	public int getNo() {
		return no;
	}

	/**
	 * @return the taff
	 */
	public TaffBond getTaff() {
		return taff;
	}

	/**
	 * @return the rotatable
	 */
	public Rotatable getRotatable() {
		return rotatable;
	}

	/**
	 * @return the canFlatten
	 */
	public boolean isCanFlatten() {
		return canFlatten;
	}

	/**
	 * @param type
	 *            the type to set
	 */
	public void setType(BondType type) {
		this.type = type;
	}

	/**
	 * @param output
	 *            the output to set
	 */
	public void setOutput(boolean output) {
		this.output = output;
	}

	/**
	 * @param atom2
	 *            the atom2 to set
	 */
	public void setAtom2(Atom atom2) {
		this.atom2 = atom2;
	}

	/**
	 * @param no
	 *            the no to set
	 */
	void setNo(int no) {
		this.no = no;
	}

	/**
	 * @param taff
	 *            the taff to set
	 */
	void setTaff(TaffBond taff) {
		this.taff = taff;
	}

}
