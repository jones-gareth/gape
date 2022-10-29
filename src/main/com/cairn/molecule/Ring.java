package com.cairn.molecule;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;
import com.google.common.primitives.Ints;

/**
 * Represents a ring in a molecule
 * 
 * @author gjones
 * 
 */
public class Ring {

	private static Logger logger;
	static {
		logger = Logger.getLogger(Ring.class);
		// logger.setLevel(Level.DEBUG);
	}
	private final List<Atom> atoms;

	private final List<Integer> orderedAtomIds;

	private final List<Bond> bonds;

	private volatile double center[] = new double[4];

	private volatile double display[] = new double[4];

	private volatile double normals[][] = new double[2][4];

	// Strict aromatic ring
	private boolean aromatic;

	// Planar Ring
	private boolean allSp2;

	Molecule molecule;

	Ring(int nAtoms, int ids[], Molecule m) {
		molecule = m;
		atoms = new ArrayList<>();
		bonds = new ArrayList<>();
		for (int i = 0; i < nAtoms; i++) {
			Atom atom = m.getAtom(ids[i]);
			atoms.add(atom);
			atom.addRing(this);
		}

		orderedAtomIds = atoms.stream().map(atom -> atom.getNo()).sorted()
				.collect(Collectors.toList());

		addBonds();
		if (m.getCoords() != null)
			findCenter();
	}

	boolean sameRing(Ring other) {
		if (molecule != other.molecule) {
			return false;
		}
		return orderedAtomIds.equals(other.orderedAtomIds);
	}

	/**
	 * @param other
	 * @return true if two rings are fused (contain a common bond).
	 */
	boolean connected(Ring other) {
		if (molecule != other.molecule) {
			return false;
		}
		for (Bond bond1 : bonds) {
			for (Bond bond2 : other.bonds) {
				if (bond1.getNo() == bond2.getNo()) {
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * Find all the bonds in the ring
	 */
	private void addBonds() {

		Molecule m = molecule;

		for (int i = 0; i < atoms.size(); i++) {
			Atom atom2 = i == atoms.size() - 1 ? atoms.get(0) : atoms.get(i + 1);
			Bond b = m.getBond(atoms.get(i), atom2);
			bonds.add(b);
			b.addRing(this);
			// System.out.println("Bond Type "+b.type.name);
		}

	}

	/**
	 * Set atom types and bond types for a ring previously determined to be
	 * aromatic.
	 */
	public void setAromaticTypes() {
		if (!aromatic)
			return;

		// set aromatic bond types if Huckels rules passed
		for (int i = 0; i < atoms.size(); i++) {

			Atom atom = atoms.get(i);
			bonds.get(i).setType(BondType.sybType(BondType.Type.AR));

			AtomType.Type newType = null;
			if (atom.getType().isOxygenType())
				newType = AtomType.Type.OAR;
			else if (atom.getType().isSulphurType())
				;
			else if (atom.getType().isCarbonType())
				newType = AtomType.Type.CAR;
			else if (atom.getType().isNitrogenType()) {
				newType = AtomType.Type.NAR;
			}

			if (newType != null)
				atom.setType(AtomType.sybType(newType));
		}
	}

	/**
	 * Check to see if the ring is aromatic or contains all sp2 groups(i.e is
	 * planar).
	 */
	public void percieveSp2() {
		// if (aromatic)
		// allSp2 = true;

		boolean notallSp2 = false;

		allSp2 = false;

		for (Atom atom : atoms) {
			AtomType type = atom.getType();
			if (type.getGeometry() != AtomType.Geometry.TRI) {
				// S.2 is often mistyped as S.3!
				if (type.getType() == AtomType.Type.S3)
					;
				else if (type.isOxygenType())
					;
				else
					notallSp2 = true;
			}
		}

		if (!notallSp2) {
			logger.debug("Ring " + info() + " is sp2 planar");
			allSp2 = true;
		}
	}

	/**
	 * Check to see if these atom ids equate to this ring.
	 * 
	 * @param atomIds
	 *            List of atom numbers.
	 */
	public boolean matchesList(int[] atomIds) {
		if (atomIds.length != atoms.size()) {
			return false;
		}
		List<Integer> idList = Ints.asList(atomIds);
		return orderedAtomIds.stream().allMatch(x -> idList.contains(x));
	}

	/**
	 * Find the center of the ring
	 * 
	 * @return
	 */
	public double[] findCenter() {
		Molecule m = molecule;
		List<double[]> points = new ArrayList<>();
		for (Atom atom : atoms) {
			points.add(m.getCoord(atom.getNo()));
		}
		Coord.center(points, center);
		center[3] = 1.0;
		Coord.copy(center, display);
		return center;
	}

	/**
	 * Determine ring center for a separate coordinate system.
	 * 
	 * @param coords
	 * @return
	 */
	public double[] findCenter(List<double[]> coords) {
		List<double[]> points = new ArrayList<>();
		for (Atom atom : atoms) {
			points.add(coords.get(atom.getNo()));
		}
		Coord.center(points, display);
		return display;
	}

	/**
	 * Return neighbours in this ring for the atom
	 * 
	 * @param a
	 * @return
	 */
	public Atom[] ringNeighbours(Atom a) {
		Atom neighbours[] = new Atom[2];
		int no = 0;
		for (Bond bond : bonds) {
			if (bond.getAtom1() == a)
				neighbours[no++] = bond.getAtom2();
			if (bond.getAtom2() == a)
				neighbours[no++] = bond.getAtom1();
		}
		return neighbours;
	}

	/**
	 * Return true if this atom is in a ring.
	 * 
	 * @param a
	 * @return
	 */
	boolean inRing(Atom a) {
		for (Atom atom : atoms) {
			if (a == atom) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Return the bond linking two ring atoms.
	 * 
	 * @param a1
	 * @param a2
	 * @return
	 */
	public Bond findBond(Atom a1, Atom a2) {
		assert (atoms.contains(a1));
		assert (atoms.contains(a2));
		Bond bond = molecule.findBond(a1, a2);
		assert (bonds.contains(bond));
		return bond;
	}

	/**
	 * Return an informational string about the ring
	 * 
	 * @return
	 */
	public String info() {
		String rtn = "Ring [";
		rtn += StringUtils.join(orderedAtomIds, ",");
		rtn += "]";
		return rtn;
	}

	static private final ThreadLocal<double[]> point1 = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[4];
		}
	};
	static private final ThreadLocal<double[]> point2 = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[4];
		}
	};
	static private final ThreadLocal<double[]> n = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[4];
		}
	};
	static private final ThreadLocal<double[]> normal = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[4];
		}
	};

	/**
	 * Create ring normals for the ring
	 * 
	 * @param normalLength
	 * @return
	 */
	public double[][] ringNormals(double normalLength) {
		Coord.zero(normal.get());
		int nAtoms = atoms.size();
		for (int i = 0; i < nAtoms; i++) {
			double a1[] = molecule.getCoord(atoms.get(i).getNo());
			double a2[] = null;
			if (i == (nAtoms - 1)) {
				a2 = molecule.getCoord(atoms.get(0).getNo());
			} else {
				a2 = molecule.getCoord(atoms.get(i + 1).getNo());
			}

			Coord.subtract(a1, center, point1.get());
			Coord.subtract(a2, center, point2.get());
			Coord.vectorProduct(point1.get(), point2.get(), n.get());
			Coord.addInPlace(normal.get(), n.get());
		}

		Coord.setLength(normal.get(), normalLength);
		Coord.add(center, normal.get(), normals[0]);
		Coord.subtract(center, normal.get(), normals[1]);
		normals[0][3] = normals[1][3] = 1.0;
		return normals;
	}

	/**
	 * @return the aromatic
	 */
	public boolean isAromatic() {
		return aromatic;
	}

	/**
	 * @return the allSp2
	 */
	public boolean isAllSp2() {
		return allSp2;
	}

	/**
	 * @return the display
	 */
	public double[] getDisplay() {
		return display;
	}

	/**
	 * @return the center
	 */
	public double[] getCenter() {
		return center;
	}

	/**
	 * @return the atoms
	 */
	public List<Atom> getAtoms() {
		return atoms;
	}

	/**
	 * @return the bonds
	 */
	public List<Bond> getBonds() {
		return bonds;
	}

	/**
	 * @param aromatic
	 *            the aromatic to set
	 */
	public void setAromatic(boolean aromatic) {
		this.aromatic = aromatic;
	}

}
