package com.cairn.molecule;

import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

/**
 * Class to represent an Atom
 * 
 * @author Gareth Jones
 * 
 * @see com.cairn.intranet.molecule.AtomType
 * 
 */
public class Atom {
	private static final Logger logger = Logger.getLogger(Atom.class);

	public static final boolean DEBUG = false;

	public enum AtomFragmentType {
		ATM_NONE, ATM_RIBOSE, ATM_ADENINE, ATM_URACIL, ATM_BENZENE, ATM_CYTOSINE
	};

	private volatile AtomType type;

	private volatile String label, name;

	private volatile List<Ring> rings = new ArrayList<>();

	private volatile int no;

	// Set output false if you don't want an atom to appear in output files
	private volatile boolean output = true;

	// SD fields
	private volatile Integer formalCharge, sdMassDifference, sdValence;

	private volatile double partialCharge;

	enum SdStereo {
		NONE, ODD, EVEN, EITHER
	};

	private volatile SdStereo sdStereo;

	private volatile boolean donorHydrogen = false;

	/**
	 * Mills and Dean Donor type.
	 */
	private volatile com.cairn.gape.feature.HydrogenBondingType donorType = null;

	/**
	 * Mills and Dean Acceptor type.
	 */
	private volatile com.cairn.gape.feature.HydrogenBondingType acceptorType = null;

	/**
	 * User defined feature.
	 */
	public volatile com.cairn.gape.feature.UserFeature userFeature = null;

	// These guys come from getNeighbours()
	private volatile Integer nNeighbours, nSingleNeighbours, nDoubleNeighbours,
			nTripleNeighbours, nAromaticNeighbours;
	// Call this hydrogenCount as nHydrogens already set for SLN
	private volatile int hydrogenCount;

	private final List<Atom> notDummyNeighbours = new ArrayList<>();
	private volatile List<Atom> lonePairs;

	private final List<AtomFragmentType> fragments = new ArrayList<>();

	// The count of lone pairs is only used while adding lone pairs- once added
	// it should be the same as the number of lone pairs in the lone pair list.
	private volatile int nLonePairs;

	private volatile double bondOrder;

	// These come from Taff.setup
	private volatile TaffAtom taff;

	private volatile TaffOop taffOop;

	private volatile List<Bond> bonds;

	private volatile List<Angle> angles;

	private volatile List<Torsion> torsions;

	private volatile Double energy, lastEnergy;

	private volatile Boolean converged = false;

	// flags for rotatable bond determination
	private volatile boolean coohGroup, planarNH2group, npl3Link;

	// this is for SLNs and also for the ring finder
	private volatile boolean inRing;

	// Tripos MOL2 substructure number and name
	private volatile Integer subNo;
	private volatile String subName;

	private volatile Molecule molecule;

	private Atom(Molecule _molecule) {
		output = true;
		molecule = _molecule;
	}

	/**
	 * Constructor.
	 * 
	 * @param _molecule
	 * @param no
	 * @param l
	 * @param n
	 */
	public Atom(Molecule _molecule, int no, String l, String n) {
		this(_molecule);
		assert (_molecule.getnAtoms() == no);
		this.no = no;
		label = l;
		name = n;
		type = AtomType.sybType(n);
	}

	/**
	 * @param _molecule
	 * @param l
	 * @param n
	 */
	public Atom(Molecule _molecule, String l, String n) {
		this(_molecule);
		no = _molecule.getnAtoms();
		label = l;
		name = n;
		type = AtomType.sybType(n);
	}

	/**
	 * Constructor. Molecule is set when constructing molecule
	 * 
	 * @param no
	 *            Atom number
	 * @param l
	 *            Atom label of name
	 * @param n
	 *            Sybyl atom type
	 */
	public Atom(int no, String l, String n) {
		this.no = no;
		label = l;
		name = n;
		type = AtomType.sybType(n);
	}

	/**
	 * Constructor. Atom label generated from atom number.
	 * 
	 * @param no
	 *            Atom number
	 * @param t
	 *            Atom type
	 */
	public Atom(Molecule _molecule, int no, AtomType.Type t) {
		this(_molecule);
		assert (_molecule.getnAtoms() == no);
		this.no = no;
		int n = no + 1;
		name = "Atom " + n;
		type = AtomType.sybType(t);
		label = type.getName();
	}

	/**
	 * Constructor. Atom label generated from atom number. Molecule is set when
	 * constructing molecule.
	 * 
	 * @param no
	 *            Atom number
	 * @param t
	 *            Atom type
	 */
	public Atom(int no, AtomType.Type t) {
		this.no = no;
		int n = no + 1;
		name = "Atom " + n;
		type = AtomType.sybType(t);
		label = type.getName();
	}

	/**
	 * Copies this Atom. For use in creating a new molecule.
	 */
	public Atom copy(Molecule newMolecule) {
		Atom newAtom = new Atom(newMolecule);
		newAtom.no = no;
		newAtom.label = label;
		newAtom.name = name;
		newAtom.type = type;
		newAtom.subName = subName;
		newAtom.subNo = subNo;

		newAtom.formalCharge = formalCharge;
		newAtom.sdMassDifference = sdMassDifference;
		newAtom.sdStereo = sdStereo;
		newAtom.sdValence = sdValence;

		return newAtom;
	}

	/**
	 * Returns true if Atom matches type t. Used by TAFF classes. Needs to be
	 * integrated with AtomType.matchType
	 * 
	 * @param t
	 * @return
	 * 
	 * @see AtomType
	 */
	boolean matchType(AtomType.Type t) {
		if (type.getType() == t)
			return true;
		if (t == AtomType.Type.WILD)
			return true;
		return false;
	}

	/**
	 * @return true is this is a heavy atom
	 */
	public boolean isHeavy() {
		return type.isHeavy();
	}

	/**
	 * @return true is this is a real atom
	 */
	public boolean isNotDummy() {
		return type.isNotDummy();
	}

	/**
	 * @return the number of heavy atoms bonded to this one
	 */
	public int getNHeavyNeighbours() {
		return notDummyNeighbours.size() - hydrogenCount;
	}

	/**
	 * Sets up TAFF forcefield data structures.
	 */
	void taffSetup() {
		angles = new ArrayList<Angle>();
		bonds = new ArrayList<Bond>();
		torsions = new ArrayList<Torsion>();

		for (Bond b : molecule.getBonds()) {
			if (b.getAtom1() == this || b.getAtom2() == this)
				bonds.add(b);
		}
		for (Angle a : molecule.getAngles()) {
			if (a.getAtom1() == this || a.getAtom2() == this || a.getMid() == this) {
				angles.add(a);
			}
		}
		for (Torsion t : molecule.getTorsions()) {
			if (t.getAtom1() == this || t.getAtom2() == this || t.getAtom3() == this
					|| t.getAtom4() == this) {
				torsions.add(t);
			}
		}

		lastEnergy = Double.MAX_VALUE;
		converged = false;
	}

	/**
	 * Frees up TAFF arrays
	 */
	public void taffFree() {
		taff = null;
		bonds = null;
		angles = null;
		torsions = null;
	}

	/**
	 * @param test
	 * @return true if test is bonded to this atom (and is not a dummy atom)
	 */
	public boolean isNotDummyNeighbour(Atom test) {
		return notDummyNeighbours.contains(test);
	}

	/**
	 * @param test
	 * @return true if test is bonded to this atom with a single bond
	 */
	public boolean isSingleNeighbour(Atom test) {
		Bond bond = molecule.getBond(this, test);
		if (bond != null && bond.getBondType() == BondType.Type.SINGLE)
			return true;
		else
			return false;
	}

	/**
	 * @param test
	 * @return true if test is bonded to this atom with a double bond
	 */
	public boolean isDoubleNeighbour(Atom test) {
		Bond bond = molecule.getBond(this, test);
		if (bond != null && bond.getBondType() == BondType.Type.DOUBLE)
			return true;
		else
			return false;
	}

	/**
	 * @param test
	 * @return true if test is bonded to this atom with a aromatic bond
	 */
	public boolean isAromaticNeighbour(Atom test) {
		Bond bond = molecule.getBond(this, test);
		if (bond != null && bond.getBondType() == BondType.Type.AR)
			return true;
		else
			return false;
	}

	/**
	 * @param fragment
	 * @return true if this atom is contained within the fragment
	 */
	public boolean inFragment(AtomFragmentType fragment) {
		return fragments.contains(fragment);
	}

	/**
	 * Add this fragment to the atom fragment list
	 * 
	 * @param frag
	 */
	public void addFragment(AtomFragmentType frag) {
		if (fragments.contains(frag)) {
			return;
		}
		fragments.add(frag);
	}

	/**
	 * @return a description of the fragments this atom belongs to.
	 */
	public String fragmentsInfo() {
		String rtn = "";
		if (fragments.size() == 0)
			return rtn;
		int id = no + 1;
		rtn = id + " " + name;
		for (AtomFragmentType fragment : fragments) {
			String name = Fragments.idToFragName(fragment);
			if (name != null)
				rtn += " " + name;
			else
				rtn += " unknown fragment " + fragment;
			rtn += "\n";
		}
		return rtn;
	}

	/**
	 * Creates neighbour lookup list for an atom. Also counts number of single,
	 * double, triple neighbours and number of hydrogens.
	 * 
	 * @param atom
	 */
	public void getNeighbours() {
		double bondOrder = .0;

		notDummyNeighbours.clear();

		int nSingle = 0, nDouble = 0, nTriple = 0, nNeighbours = 0, nAromatic = 0,
				nHydrogens = 0;

		for (Bond b : molecule.getBonds()) {
			Atom otherAtom = null;
			if (this == b.getAtom2())
				otherAtom = b.getAtom1();
			else if (this == b.getAtom1())
				otherAtom = b.getAtom2();
			else
				continue;
			BondType.Type type = b.getType().getType();

			nNeighbours++;

			if (otherAtom.isNotDummy()) {
				notDummyNeighbours.add(otherAtom);

				if (otherAtom.getAtomType() == AtomType.Type.H) {
					nHydrogens++;
				}

				switch (type) {
				case SINGLE:
					nSingle++;
					bondOrder += 1;
					break;
				case DOUBLE:
					nDouble++;
					bondOrder += 2;
					break;
				case TRIPLE:
					nTriple++;
					bondOrder += 3;
					break;
				case AM:
					nSingle++;
					bondOrder += 1;
					break;
				case AR:
					nAromatic++;
					// This is flawed !!
					bondOrder += 1.5;
					break;
				default:
					break;
				}
			}

		}

		this.nNeighbours = nNeighbours;
		nSingleNeighbours = nSingle;
		nDoubleNeighbours = nDouble;
		nTripleNeighbours = nTriple;
		nAromaticNeighbours = nAromatic;
		hydrogenCount = nHydrogens;
		this.bondOrder = bondOrder;
	}

	/**
	 * @return a descriptive label.
	 */
	public String info() {
		int id = no + 1;
		return new String(id + ":" + type.getName() + " [" + name + "]");
	}

	/**
	 * Finds and sets common group flags
	 */
	public void findGroups() {
		checkNpl3Link();
		checkPlanarNH2group();
		checkCOOHgroup();
	}

	/**
	 * Checks and sets the Npl3 link flags. Npl3 link is an Npl3 atom that links
	 * two rings
	 */
	private void checkNpl3Link() {
		if (type.getType() != AtomType.Type.NPL3)
			return;
		if (isInRing())
			return;
		if (getnNotDummyNeighbours() != 3)
			return;
		int cnt = 0;
		Atom ringAtoms[] = new Atom[3];
		for (int i = 0; i < 3; i++)
			if (notDummyNeighbours.get(i).isInRing()) {
				ringAtoms[cnt] = notDummyNeighbours.get(i);
				cnt++;
			}
		if (cnt < 2)
			return;
		npl3Link = true;
		for (int i = 0; i < cnt; i++)
			ringAtoms[i].npl3Link = true;

		if (DEBUG)
			System.out.println("Npl3 ring link group " + info());
	}

	/**
	 * Checks and sets flags for planar NH2 groups.
	 */
	private void checkPlanarNH2group() {
		if (type.getType() != AtomType.Type.NPL3)
			return;
		if (isInRing())
			return;
		if (hydrogenCount != 2)
			return;
		if (getnNotDummyNeighbours() != 3)
			return;
		planarNH2group = true;
		for (Atom neighbour : notDummyNeighbours) {
			if (neighbour.getAtomType() == AtomType.Type.H) {
				neighbour.planarNH2group = true;
			}
		}

		if (DEBUG)
			System.out.println("Planar NH2 group " + info());
	}

	/**
	 * Checks and sets flags for planar carboxylic acid groups.
	 */
	private void checkCOOHgroup() {
		if (type.getType() != AtomType.Type.C2)
			return;
		if (isInRing())
			return;
		if (getnNotDummyNeighbours() != 3)
			return;

		Atom o2atom = null, o3atom = null;
		for (Atom check : notDummyNeighbours) {
			if (check.getAtomType() == AtomType.Type.O3)
				o3atom = check;
			if (check.getAtomType() == AtomType.Type.O2)
				o2atom = check;
		}
		if (o3atom == null || o2atom == null)
			return;
		if (o3atom.hydrogenCount != 1)
			return;

		coohGroup = true;
		o2atom.coohGroup = true;
		o3atom.coohGroup = true;

		for (Atom o3AtomNeighbour : o3atom.getNotDummyNeighbours()) {
			if (o3AtomNeighbour.getAtomType() == AtomType.Type.H) {
				o3AtomNeighbour.coohGroup = true;
			}
		}

		if (DEBUG)
			System.out.println("COOH group " + info());
	}

	/**
	 * @param a2
	 * @return the bond betwween this atom and a2
	 */
	Bond getBond(Atom a2) {
		return molecule.getBond(this, a2);
	}

	/**
	 * @return true if this is an aromatic oxygen.
	 */
	public boolean isAromaticOxygen() {
		if (!type.isOxygenType())
			return false;
		if (getnNotDummyNeighbours() != 2)
			return false;
		if (nAromaticNeighbours != 2)
			return false;
		return true;
	}

	/**
	 * Checks the current atom type.
	 * 
	 * @return the predicted atom type.
	 */
	public AtomType.Type checkAtomType() {
		if (type.isCarbonType())
			return getCarbonType();
		if (type.isNitrogenType())
			return getNitrogenType();
		if (type.isOxygenType())
			return getOxygenType();
		if (type.isSulphurType())
			return getSulphurType();
		if (type.isPhosphorousType())
			return getPhosphorousType();
		else
			return getAtomType();
	}

	/**
	 * @return the type for a carbon
	 */
	private AtomType.Type getCarbonType() {
		if (!type.isCarbonType())
			return null;

		int no = getnNotDummyNeighbours();

		if (no == 4)
			return AtomType.Type.C3;

		if (no == 3) {
			if (bondOrder >= 5) {
				return AtomType.Type.CCAT;
			}

			// Look for arginine
			if (isArginineCarbon())
				return AtomType.Type.CCAT;

			if (isCarboxylateCarbon())
				return AtomType.Type.C2;

			if (nAromaticNeighbours == 2 && nSingleNeighbours == 1)
				return AtomType.Type.CAR;

			if (nDoubleNeighbours == 1 && nSingleNeighbours == 2)
				return AtomType.Type.C2;

			if (nAromaticNeighbours == 3)
				return AtomType.Type.CAR;
		}

		if (no == 2)
			return AtomType.Type.C1;

		return null;
	}

	/**
	 * Checks to see if this carbon is the center of an arginine group
	 * 
	 * @return
	 */
	boolean isArginineCarbon() {
		if (getnNotDummyNeighbours() != 3) {
			return false;
		}

		if (notDummyNeighbours.get(0).type.isNitrogenType()
				&& notDummyNeighbours.get(1).type.isNitrogenType()
				&& notDummyNeighbours.get(2).type.isNitrogenType()) {
			int nSingle = 0, nDouble = 0, nAromatic = 0;

			Bond b1 = getBond(notDummyNeighbours.get(0));
			if (b1.getBondType() == BondType.Type.SINGLE)
				nSingle++;
			else if (b1.getBondType() == BondType.Type.DOUBLE)
				nDouble++;
			else if (b1.getBondType() == BondType.Type.AR)
				nAromatic++;

			Bond b2 = getBond(notDummyNeighbours.get(1));
			if (b2.getBondType() == BondType.Type.SINGLE)
				nSingle++;
			else if (b2.getBondType() == BondType.Type.DOUBLE)
				nDouble++;
			else if (b2.getBondType() == BondType.Type.AR)
				nAromatic++;

			Bond b3 = getBond(notDummyNeighbours.get(2));
			if (b3.getBondType() == BondType.Type.SINGLE)
				nSingle++;
			else if (b3.getBondType() == BondType.Type.DOUBLE)
				nDouble++;
			else if (b3.getBondType() == BondType.Type.AR)
				nAromatic++;

			int n1 = notDummyNeighbours.get(0).getnNotDummyNeighbours();
			int n2 = notDummyNeighbours.get(1).getnNotDummyNeighbours();
			int n3 = notDummyNeighbours.get(2).getnNotDummyNeighbours();

			if ((nAromatic == 3 || (nSingle == 2 && nDouble == 1)) && n1 == 3 && n2 == 3
					&& n3 == 3)
				return true;
		}
		return false;
	}

	/**
	 * Checks to see if this carbon is in an amide group
	 * 
	 * @return
	 */
	boolean isAmideCarbon() {
		if (!type.isCarbonType())
			return false;

		if (getnNotDummyNeighbours() != 3)
			return false;
		if (nDoubleNeighbours != 1)
			return false;

		if (nSingleNeighbours != 2)
			return false;

		boolean hasNitrogen = false, hasOxygen = false;
		for (Atom test : notDummyNeighbours) {
			if (test.type.isNitrogenType())
				hasNitrogen = true;
			if (test.type.isOxygenType() && test.getnNotDummyNeighbours() == 1
					&& test.nDoubleNeighbours == 1)
				hasOxygen = true;
		}

		if (hasOxygen && hasNitrogen)
			return true;
		else
			return false;
	}

	/**
	 * Check to see if this carbon is in a nitro group.
	 * 
	 * @return
	 */
	boolean isNitroNitrogen() {
		if (!type.isNitrogenType())
			return false;

		if (getnNotDummyNeighbours() != 3)
			return false;

		int no = 0;
		for (Atom check : notDummyNeighbours) {
			if (check.type.isOxygenType() && check.getnNotDummyNeighbours() == 1)
				no++;
		}

		if (no == 2)
			return true;
		return false;
	}

	/**
	 * @return true if this oxygen is in a nitro group.
	 */
	public boolean isNitroOxygen() {
		if (!type.isOxygenType())
			return false;

		if (getnNotDummyNeighbours() != 1)
			return false;

		return notDummyNeighbours.get(0).isNitroNitrogen();
	}

	/**
	 * @return the predicted atom type for a nitrogen
	 */
	private AtomType.Type getNitrogenType() {
		if (!type.isNitrogenType()) {
			return null;
		}

		int no = getnNotDummyNeighbours();

		if (no == 1)
			return AtomType.Type.N1;

		if (no == 4)
			return AtomType.Type.N4;

		// NPL3 and CCAT
		if (nNeighbours == 3 && nSingleNeighbours == 2) {
			for (Atom checkAtom : notDummyNeighbours) {
				AtomType.Type check = checkAtom.getCarbonType();
				if (check == AtomType.Type.CCAT) {
					return AtomType.Type.NPL3;
				}
			}
		}

		if (no >= 2 && nDoubleNeighbours == 0) {

			// NAM
			for (int i = 0; i < no; i++) {
				if (notDummyNeighbours.get(i).isAmideCarbon()) {
					if (no == 2) {
						logger.warn("Negatively charged amide nitrogen!");
					}
					return AtomType.Type.NAM;
				}
				if (notDummyNeighbours.get(i).isSulphonamideSulphur()) {
					if (no == 2) {
						logger.warn("Negatively charged sulfonamide nitrogen!");
					}
					// never too sure if these guys should be NPL3 or NAM
					return AtomType.Type.NPL3;
				}
			}

			// NPL3
			if (nSingleNeighbours == 3) {
				for (Atom checkAtom : notDummyNeighbours) {
					AtomType.Type check = checkAtom.getCarbonType();
					if (check == AtomType.Type.C2 || check == AtomType.Type.CAR
							|| check == AtomType.Type.CCAT)
						return AtomType.Type.NPL3;
					// Can't call getNitrogenType here- possible infinite
					// recursion
					if (checkAtom.type.isNitrogenType()
							&& checkAtom.nDoubleNeighbours == 1)
						return AtomType.Type.NPL3;
				}
			}

			if (nSingleNeighbours == 3)
				return AtomType.Type.N3;

		}

		if (no == 3 && isNitroNitrogen())
			return AtomType.Type.N2;

		// Charged Planar Nitrogens
		if (no == 3 && nDoubleNeighbours == 1)
			return AtomType.Type.N2;
		if (no == 3 && nAromaticNeighbours == 2)
			return AtomType.Type.NAR;

		if (no == 2 && nDoubleNeighbours == 1)
			return AtomType.Type.N2;

		if (no == 2 && nAromaticNeighbours == 2)
			return AtomType.Type.NAR;

		if (no == 2 && nSingleNeighbours == 2)
			return AtomType.Type.NPL3;

		// Charged Linear Nitrogen
		if (no == 2 && nTripleNeighbours == 1 && nSingleNeighbours == 1)
			return AtomType.Type.N1;

		// this is presumably from a planar nitrogen in a Kekulized ring system.
		// This is take at face value, but additional checking should still be
		// done. Ring perception etc comes after atom typing
		if (no == 3 && nAromaticNeighbours == 3 && isInRing()) {
			return AtomType.Type.NAR;
		}
		return null;
	}

	/**
	 * Checks to see if this carbon is a carboxylate
	 * 
	 * @return
	 */
	private boolean isCarboxylateCarbon() {
		if (!type.isCarbonType())
			return false;

		if (getnNotDummyNeighbours() != 3)
			return false;

		int no = 0;
		for (Atom check : notDummyNeighbours) {
			if (check.type.isOxygenType() && check.getnNotDummyNeighbours() == 1)
				no++;
		}
		if (no == 2)
			return true;
		else
			return false;
	}

	/**
	 * Checks to see if this oxygen is in a carboxylate group
	 * 
	 * @return
	 */
	public boolean isCarboxylateOxygen() {
		if (getnNotDummyNeighbours() != 1)
			return false;
		if (notDummyNeighbours.get(0).isCarboxylateCarbon())
			return true;
		return false;
	}

	/**
	 * Checks to see if this sulphur is in a sulponamide group
	 * 
	 * @return
	 */
	boolean isSulphonamideSulphur() {
		if (!type.isSulphurType())
			return false;

		int nOxygens = nPSOxygens();
		if (nOxygens < 2)
			return false;
		for (Atom check : notDummyNeighbours) {
			if (check.type.isNitrogenType() && check.getnNotDummyNeighbours() == 3
					&& check.nSingleNeighbours == 3) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Finds the number of oxygens bound to this atom which are not bound to any
	 * other neighbours. Used for P and S base groups.
	 * 
	 * @return
	 */
	int nPSOxygens() {
		if (!type.isPhosphorousType() && !type.isSulphurType())
			return -1;
		int no = 0;
		for (Atom check : notDummyNeighbours) {
			if (check.type.isOxygenType() && check.getnNotDummyNeighbours() == 1) {
				no++;
			}
		}
		return no;
	}

	/**
	 * @return the predicted atom type for an oxygen atom
	 */
	private AtomType.Type getOxygenType() {
		if (!type.isOxygenType())
			return null;

		int no = getnNotDummyNeighbours();

		if (no == 2) {
			if (nSingleNeighbours == 2)
				return AtomType.Type.O3;

			// 0 Minus =O- or aromatic oxygen in Ring.
			if (isInRing() && nDoubleNeighbours == 1 && nSingleNeighbours == 1)
				return AtomType.Type.O2;
			if (isInRing() && nAromaticNeighbours == 2)
				return AtomType.Type.OAR;
		}

		if (no == 1) {
			Atom neighbour = notDummyNeighbours.get(0);

			if (neighbour.type.isCarbonType() && neighbour.isCarboxylateCarbon())
				return AtomType.Type.OCO2;

			else if (neighbour.type.isPhosphorousType())
				return AtomType.Type.OCO2;

			else if (neighbour.type.isSulphurType()) {
				no = neighbour.nPSOxygens();
				if (no == 1 || no == 2)
					return AtomType.Type.O2;

				if (no == 3) {
					// Arbitrary assignment of first neighbour to O.co2
					for (Atom check : notDummyNeighbours) {
						if (check.type.isOxygenType()
								&& check.getnNotDummyNeighbours() == 1) {
							if (check == this)
								return AtomType.Type.OCO2;
							else
								return AtomType.Type.O2;
						}
					}
				}

			}

			else if (nDoubleNeighbours == 1)
				return AtomType.Type.O2;

			// O Minus
			else if (nSingleNeighbours == 1)
				return AtomType.Type.O3;

		}

		return null;
	}

	/**
	 * @return the predicted atom type for a Phosphoros atom
	 */
	private AtomType.Type getPhosphorousType() {
		if (!type.isPhosphorousType())
			return null;
		return AtomType.Type.P3;
	}

	/**
	 * @return the predicted atom type for a sulphur atom
	 */
	private AtomType.Type getSulphurType() {
		if (!type.isSulphurType())
			return null;

		int no = getnNotDummyNeighbours();
		int nOxygens = nPSOxygens();
		if (nOxygens >= 2)
			return AtomType.Type.SO2;
		else if (nOxygens == 1)
			return AtomType.Type.SO;

		if (no == 4)
			return AtomType.Type.S3;
		if ((no == 1 || no == 2 || no == 3) && nDoubleNeighbours > 0)
			return AtomType.Type.S2;
		// if (no == 2) return AtomType.Type.S2;
		if (no == 2 || no == 3)
			return AtomType.Type.S3;

		return null;
	}

	/**
	 * Tries to predict formal charge from bond order. This method is flawed as
	 * aromatic nitrogens with three connections my be charged or neutral.
	 * 
	 * @return
	 */
	public double getChargeFromBondOrder() {
		if (type.getNeutralBondOrder() == -1)
			return .0;
		double charge1 = bondOrder - type.getNeutralBondOrder();
		if (type.getNeutralBondOrder2() == -1)
			return charge1;
		double charge2 = bondOrder - type.getNeutralBondOrder2();
		if (Math.abs(charge1) < Math.abs(charge2))
			return charge1;
		return charge2;
	}

	/**
	 * Adds a ring to this atom.
	 * 
	 * @param r
	 */
	public void addRing(Ring r) {
		rings.add(r);
	}

	/**
	 * Sets atom types, bond types and charges for charged acid groups.
	 */
	public void setChargedAcidGroups() {
		if (!type.isCarbonType() && !type.isPhosphorousType() && !type.isSulphurType())
			return;

		int nOxygens = 0;
		Atom oxygens[] = new Atom[4];
		Bond bonds[] = new Bond[4];
		for (Atom check : notDummyNeighbours) {
			if (check.type.isOxygenType() && check.getnNotDummyNeighbours() == 1) {
				oxygens[nOxygens] = check;
				bonds[nOxygens] = molecule.getBond(check, this);
				nOxygens++;
			}
		}

		if (nOxygens < 2)
			return;

		if (nOxygens == 2 && ((type.isCarbonType() && getnNotDummyNeighbours() == 3)
				|| (type.isPhosphorousType() && getnNotDummyNeighbours() == 4))) {
			oxygens[0].formalCharge = -1;
			oxygens[0].partialCharge = -0.5;
			oxygens[0].type = AtomType.sybType(AtomType.Type.OCO2);
			bonds[0].setType(BondType.sybType(BondType.Type.SINGLE));

			oxygens[1].formalCharge = null;
			oxygens[1].partialCharge = -0.5;
			oxygens[1].type = AtomType.sybType(AtomType.Type.OCO2);
			bonds[1].setType(BondType.sybType(BondType.Type.DOUBLE));

			oxygens[0].getNeighbours();
			oxygens[1].getNeighbours();
			getNeighbours();
		}

		else if (nOxygens == 2 && type.isSulphurType() && getnNotDummyNeighbours() == 4) {
			oxygens[0].formalCharge = null;
			oxygens[0].partialCharge = .0;
			oxygens[0].type = AtomType.sybType(AtomType.Type.O2);
			bonds[0].setType(BondType.sybType(BondType.Type.DOUBLE));

			oxygens[1].formalCharge = null;
			oxygens[1].partialCharge = .0;
			oxygens[1].type = AtomType.sybType(AtomType.Type.O2);
			bonds[1].setType(BondType.sybType(BondType.Type.DOUBLE));

			oxygens[0].getNeighbours();
			oxygens[1].getNeighbours();
			getNeighbours();
		}

		else if (nOxygens == 3 && type.isSulphurType()) {
			oxygens[0].formalCharge = null;
			oxygens[0].partialCharge = -0.333;
			oxygens[0].type = AtomType.sybType(AtomType.Type.OCO2);
			bonds[0].setType(BondType.sybType(BondType.Type.DOUBLE));

			oxygens[1].formalCharge = null;
			oxygens[1].partialCharge = -0.333;
			oxygens[1].type = AtomType.sybType(AtomType.Type.OCO2);
			bonds[1].setType(BondType.sybType(BondType.Type.DOUBLE));

			oxygens[2].formalCharge = -1;
			oxygens[2].partialCharge = -0.333;
			oxygens[2].type = AtomType.sybType(AtomType.Type.OCO2);
			bonds[2].setType(BondType.sybType(BondType.Type.SINGLE));

			oxygens[0].getNeighbours();
			oxygens[1].getNeighbours();
			oxygens[2].getNeighbours();
			getNeighbours();
		}

		else if (nOxygens == 3 && type.isPhosphorousType()) {
			oxygens[0].formalCharge = -1;
			oxygens[0].partialCharge = -0.667;
			oxygens[0].type = AtomType.sybType(AtomType.Type.OCO2);
			bonds[0].setType(BondType.sybType(BondType.Type.SINGLE));

			oxygens[1].formalCharge = null;
			oxygens[1].partialCharge = -0.667;
			oxygens[1].type = AtomType.sybType(AtomType.Type.OCO2);
			bonds[1].setType(BondType.sybType(BondType.Type.DOUBLE));

			oxygens[2].formalCharge = -1;
			oxygens[2].partialCharge = -0.667;
			oxygens[2].type = AtomType.sybType(AtomType.Type.OCO2);
			bonds[2].setType(BondType.sybType(BondType.Type.SINGLE));

			oxygens[0].getNeighbours();
			oxygens[1].getNeighbours();
			oxygens[2].getNeighbours();
			getNeighbours();
		}

	}

	/**
	 * Convenience method to return the base atom type
	 * 
	 * @return
	 */
	public AtomType.Type getAtomType() {
		return getType().getType();
	}

	/**
	 * Return true if this atom is in a ring.
	 * 
	 * @return
	 */
	public boolean isInRing() {
		return inRing || rings.size() > 0;
	}

	/**
	 * Remove all ring information
	 */
	void setNotInRing() {
		rings.clear();
		inRing = false;
	}

	/**
	 * @return the debug
	 */
	public static boolean isDebug() {
		return DEBUG;
	}

	/**
	 * @return the type
	 */
	public AtomType getType() {
		return type;
	}

	/**
	 * @return the label
	 */
	public String getLabel() {
		return label;
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
	 * @param inRing
	 *            the inRing to set
	 */
	void setInRing(boolean inRing) {
		this.inRing = inRing;
		if (inRing == false) {
			rings.clear();
		}
	}

	/**
	 * @return the output
	 */
	public boolean isOutput() {
		return output;
	}

	/**
	 * @return the formalCharge
	 */
	public Integer getFormalCharge() {
		return formalCharge;
	}

	/**
	 * @return the sdMassDifference
	 */
	public Integer getSdMassDifference() {
		return sdMassDifference;
	}

	/**
	 * @return the sdValence
	 */
	public Integer getSdValence() {
		return sdValence;
	}

	/**
	 * @return the partialCharge
	 */
	public double getPartialCharge() {
		return partialCharge;
	}

	/**
	 * @return the sdStereo
	 */
	public SdStereo getSdStereo() {
		return sdStereo;
	}

	/**
	 * @return the donorHydrogen
	 */
	public boolean isDonorHydrogen() {
		return donorHydrogen;
	}

	/**
	 * @return the donorType
	 */
	public com.cairn.gape.feature.HydrogenBondingType getDonorType() {
		return donorType;
	}

	/**
	 * @return the acceptorType
	 */
	public com.cairn.gape.feature.HydrogenBondingType getAcceptorType() {
		return acceptorType;
	}

	/**
	 * @return the userFeature
	 */
	public com.cairn.gape.feature.UserFeature getUserFeature() {
		return userFeature;
	}

	/**
	 * @return the nNeighbours
	 */
	public Integer getnNeighbours() {
		return nNeighbours;
	}

	/**
	 * @return the nSingleNeighbours
	 */
	public Integer getnSingleNeighbours() {
		return nSingleNeighbours;
	}

	/**
	 * @return the nDoubleNeighbours
	 */
	public Integer getnDoubleNeighbours() {
		return nDoubleNeighbours;
	}

	/**
	 * @return the nTripleNeighbours
	 */
	public Integer getnTripleNeighbours() {
		return nTripleNeighbours;
	}

	/**
	 * @return the nAromaticNeighbours
	 */
	public Integer getnAromaticNeighbours() {
		return nAromaticNeighbours;
	}

	/**
	 * @return the hydrogenCount
	 */
	public int getHydrogenCount() {
		return hydrogenCount;
	}

	/**
	 * @return the notDummyNeighbours
	 */
	public List<Atom> getNotDummyNeighbours() {
		return notDummyNeighbours;
	}

	/**
	 * @return the notDummyNeighbours
	 */
	public int getnNotDummyNeighbours() {
		return notDummyNeighbours.size();
	}

	/**
	 * @return the lonePairs
	 */
	public List<Atom> getLonePairs() {
		return lonePairs;
	}

	/**
	 * @return the fragments
	 */
	public List<AtomFragmentType> getFragments() {
		return fragments;
	}

	/**
	 * @return the nLonePairs
	 */
	public int getnLonePairs() {
		return nLonePairs;
	}

	/**
	 * @return the bondOrder
	 */
	public double getBondOrder() {
		return bondOrder;
	}

	/**
	 * @return the taff
	 */
	public TaffAtom getTaff() {
		return taff;
	}

	/**
	 * @return the taffOop
	 */
	public TaffOop getTaffOop() {
		return taffOop;
	}

	/**
	 * @return the bonds
	 */
	public List<Bond> getBonds() {
		return bonds;
	}

	/**
	 * @return the angles
	 */
	public List<Angle> getAngles() {
		return angles;
	}

	/**
	 * @return the torsions
	 */
	public List<Torsion> getTorsions() {
		return torsions;
	}

	/**
	 * @return the energy
	 */
	public Double getEnergy() {
		return energy;
	}

	/**
	 * @return the lastEnergy
	 */
	public Double getLastEnergy() {
		return lastEnergy;
	}

	/**
	 * @return the converged
	 */
	public Boolean getConverged() {
		return converged;
	}

	/**
	 * @return the cOOHgroup
	 */
	public boolean isCOOHgroup() {
		return coohGroup;
	}

	/**
	 * @return the planarNH2group
	 */
	public boolean isPlanarNH2group() {
		return planarNH2group;
	}

	/**
	 * @return the npl3Link
	 */
	public boolean isNpl3Link() {
		return npl3Link;
	}

	/**
	 * @return the subNo
	 */
	public Integer getSubNo() {
		return subNo;
	}

	/**
	 * @return the subName
	 */
	public String getSubName() {
		return subName;
	}

	/**
	 * @param formalCharge
	 *            the formalCharge to set
	 */
	public void setFormalCharge(Integer formalCharge) {
		this.formalCharge = formalCharge;
	}

	/**
	 * @param partialCharge
	 *            the partialCharge to set
	 */
	public void setPartialCharge(double partialCharge) {
		this.partialCharge = partialCharge;
	}

	/**
	 * @param type
	 *            the type to set
	 */
	public void setType(AtomType type) {
		this.type = type;
	}

	/**
	 * @param molecule
	 *            the molecule to set
	 */
	void setMolecule(Molecule molecule) {
		this.molecule = molecule;
	}

	/**
	 * @param donorType
	 *            the donorType to set
	 */
	public void setDonorType(com.cairn.gape.feature.HydrogenBondingType donorType) {
		this.donorType = donorType;
	}

	/**
	 * @param acceptorType
	 *            the acceptorType to set
	 */
	public void setAcceptorType(
			com.cairn.gape.feature.HydrogenBondingType acceptorType) {
		this.acceptorType = acceptorType;
	}

	/**
	 * @param donorHydrogen
	 *            the donorHydrogen to set
	 */
	public void setDonorHydrogen(boolean donorHydrogen) {
		this.donorHydrogen = donorHydrogen;
	}

	/**
	 * @param output
	 *            the output to set
	 */
	public void setOutput(boolean output) {
		this.output = output;
	}

	/**
	 * @param nLonePairs
	 *            the nLonePairs to set
	 */
	public void setnLonePairs(int nLonePairs) {
		this.nLonePairs = nLonePairs;
	}

	/**
	 * @param lonePairs
	 *            the lonePairs to set
	 */
	public void setLonePairs(List<Atom> lonePairs) {
		if (lonePairs != null) {
			assert lonePairs.size() == nLonePairs;
		}
		this.lonePairs = lonePairs;
	}

	/**
	 * @param sdMassDifference
	 *            the sdMassDifference to set
	 */
	void setSdMassDifference(Integer sdMassDifference) {
		this.sdMassDifference = sdMassDifference;
	}

	/**
	 * @param sdValence
	 *            the sdValence to set
	 */
	void setSdValence(Integer sdValence) {
		this.sdValence = sdValence;
	}

	/**
	 * @param sdStereo
	 *            the sdStereo to set
	 */
	void setSdStereo(SdStereo sdStereo) {
		this.sdStereo = sdStereo;
	}

	/**
	 * @param subNo
	 *            the subNo to set
	 */
	void setSubNo(Integer subNo) {
		this.subNo = subNo;
	}

	/**
	 * @param subName
	 *            the subName to set
	 */
	void setSubName(String subName) {
		this.subName = subName;
	}

	/**
	 * @param no
	 *            the no to set
	 */
	void setNo(int no) {
		this.no = no;
	}

	/**
	 * @param label
	 *            the label to set
	 */
	void setLabel(String label) {
		this.label = label;
	}

	/**
	 * @param userFeature
	 *            the userFeature to set
	 */
	public void setUserFeature(com.cairn.gape.feature.UserFeature userFeature) {
		this.userFeature = userFeature;
	}

	/**
	 * @param energy
	 *            the energy to set
	 */
	public void setEnergy(Double energy) {
		this.energy = energy;
	}

	/**
	 * @param taff
	 *            the taff to set
	 */
	void setTaff(TaffAtom taff) {
		this.taff = taff;
	}

	/**
	 * @param taffOop
	 *            the taffOop to set
	 */
	void setTaffOop(TaffOop taffOop) {
		this.taffOop = taffOop;
	}

}
