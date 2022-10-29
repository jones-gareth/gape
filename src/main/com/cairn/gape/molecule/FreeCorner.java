package com.cairn.gape.molecule;

import java.util.Arrays;
import java.util.List;

import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;
import com.cairn.molecule.Atom;
import com.cairn.molecule.AtomType;
import com.cairn.molecule.Bond;
import com.cairn.molecule.BondType;
import com.cairn.molecule.Molecule;
import com.cairn.molecule.Ring;
import com.cairn.molecule.RotatableBond;
import com.cairn.molecule.Torsion;

/**
 * @author Gareth Jones
 * 
 *         See J. Mol Graph. V11 N2 Payne & Glen pp74-91
 * 
 *         Also Computers Chem V15 N1 Senn pp 93-94 for intersection point of
 *         three spheres
 */
public class FreeCorner {
	private static final Logger logger = Logger.getLogger(FreeCorner.class);

	// Free corner is at X
	// Ring goes A-B-X-C-D
	private volatile Atom atomX, atomA, atomB, atomC, atomD;

	private volatile Bond bondAB, bondBX, bondXC, bondCD;

	private volatile Molecule molecule;

	private volatile Ring ring;

	// When flipping a free corner if an interatomic distance
	// distortion of greater than CORNER_FLAP_TOL occurs then the
	// corner is not flipped.

	static final double CORNER_FLAP_TOL = 0.15;

	private final boolean freeCorner;

	private volatile CornerRotation rotationCD, rotationXC, rotationAB;

	private final double altPoint[] = new double[4];

	public FreeCorner(Atom x, Molecule m) {
		freeCorner = checkFreeCorner(x, m);
	}

	private boolean checkFreeCorner(Atom x, Molecule m) {
		atomX = x;
		molecule = m;

		if (atomX.getRings().size() != 1)
			return false;
		ring = atomX.getRings().get(0);
		logger.debug(ring.info());
		if (ring.getAtoms().size() < 4)
			return false;

		logger.debug("Atom X " + atomX.info());
		Atom[] neighbours = ring.ringNeighbours(atomX);
		atomC = neighbours[0];
		if (atomC.getRings().size() > 1) {
			return false;
		}
		atomB = neighbours[1];
		if (atomB.getRings().size() > 1) {
			return false;
		}
		logger.debug("Atom C " + atomC.info());
		logger.debug("Atom B " + atomB.info());

		bondBX = ring.findBond(atomB, atomX);
		bondXC = ring.findBond(atomX, atomC);
		logger.debug("Bond BX " + bondBX.info());
		logger.debug("Bond XC " + bondXC.info());
		if (bondBX.getBondType() != BondType.Type.SINGLE)
			return false;
		if (bondXC.getBondType() != BondType.Type.SINGLE)
			return false;

		neighbours = ring.ringNeighbours(atomC);
		atomD = neighbours[0] == atomX ? neighbours[1] : neighbours[0];
		neighbours = ring.ringNeighbours(atomB);
		atomA = neighbours[0] == atomX ? neighbours[1] : neighbours[0];

		bondAB = ring.findBond(atomA, atomB);
		bondCD = ring.findBond(atomC, atomD);
		if (bondAB.getBondType() != BondType.Type.SINGLE)
			return false;
		if (bondCD.getBondType() != BondType.Type.SINGLE)
			return false;

		rotationCD = new CornerRotation(atomC, bondCD, atomB);
		rotationXC = new CornerRotation(atomX, bondXC, atomB);
		rotationAB = new CornerRotation(atomB, bondAB, atomX);

		return true;
	}

	public static void main(String args[]) {
		if (args.length != 1) {
			logger.error("Usage: FreeCorner <mol2file>");
			System.exit(1);
		}
		GaMolecule mol = new GaMolecule();
		mol.loadFile(args[0]);
		mol.findFreeCorners();
		for (int i = 0; i < mol.getNFreeCorners(); i++) {
			logger.info("Flipping free corner " + String.valueOf(i + 1));
			mol.getFreeCorner(i).flipCorner();
		}

		mol.writeSybylMol2File("Free Corner Flip " + mol.getName() + ".mol2", "");
	}

	private final double m[] = new double[4];

	private void findOtherCorner() {
		double a[] = molecule.getCoord(atomA.getNo());
		double b[] = molecule.getCoord(atomB.getNo());
		double x[] = molecule.getCoord(atomX.getNo());
		Coord.midPoint(molecule.getCoord(atomC.getNo()),
				molecule.getCoord(atomD.getNo()), m);

		int no = 0;
		boolean fitted = false;

		double r1sqr = 0, r2sqr = 0, r3sqr = 0, x12 = 0, x13 = 0, y12 = 0, y13 = 0, y23 = 0, z12 = 0, z13 = 0, l1Sqr = 0, l2Sqr = 0, l3Sqr = 0, p12Sqr = 0, p13Sqr = 0, epsilon = 0;

		while (!fitted) {
			r1sqr = Coord.sqrDistance(x, a);
			r2sqr = Coord.sqrDistance(x, b);
			r3sqr = Coord.sqrDistance(x, m);

			if (logger.isDebugEnabled()) {
				logger.debug("r1sqr " + r1sqr);
				logger.debug("r2sqr " + r2sqr);
				logger.debug("r3sqr " + r3sqr);
			}

			x12 = a[0] - b[0];
			x13 = a[0] - m[0];
			y12 = a[1] - b[1];
			y13 = a[1] - m[1];
			y23 = b[1] - m[1];
			z12 = a[2] - b[2];
			z13 = a[2] - m[2];

			if (logger.isDebugEnabled()) {
				logger.debug("x12 " + x12);
				logger.debug("x13 " + x13);
				logger.debug("y12 " + y12);
				logger.debug("y13 " + y13);
				logger.debug("y23 " + y23);
				logger.debug("z12 " + z12);
				logger.debug("z13 " + z13);
			}

			l1Sqr = a[0] * a[0] + a[1] * a[1] + a[2] * a[2] - r1sqr;
			l2Sqr = b[0] * b[0] + b[1] * b[1] + b[2] * b[2] - r2sqr;
			l3Sqr = m[0] * m[0] + m[1] * m[1] + m[2] * m[2] - r3sqr;

			if (logger.isDebugEnabled()) {
				logger.debug("l1Sqr " + l1Sqr);
				logger.debug("l2Sqr " + l2Sqr);
				logger.debug("l3Sqr " + l3Sqr);
			}

			p12Sqr = 0.5 * (l1Sqr - l2Sqr);
			p13Sqr = 0.5 * (l1Sqr - l3Sqr);

			if (logger.isDebugEnabled()) {
				logger.debug("p12Sqr " + p12Sqr);
				logger.debug("p13Sqr " + p13Sqr);
				logger.debug("a[0] " + a[0]);
				logger.debug("b[0] " + b[0]);
				logger.debug("m[0] " + m[0]);
			}
			double denom = -y23 * a[0] + y13 * b[0] - y12 * m[0];
			epsilon = 1.0 / denom;
			logger.debug("denom " + denom);
			logger.debug("epsilon " + epsilon);

			double test = (epsilon > 0.0) ? epsilon : -epsilon;
			if (test < 1.0e-5) {
				a = molecule.getCoord(atomB.getNo());
				b = molecule.getCoord(atomA.getNo());
			} else {
				fitted = true;
			}
			no++;
			if (no > 2)
				throw new RuntimeException("Epsilon Test failed");
		}

		double alpha = epsilon * (y12 * p13Sqr - y13 * p12Sqr) - a[0];
		double beta = epsilon * (x13 * p12Sqr - x12 * p13Sqr) - a[1];
		double gamma = epsilon * (y13 * z12 - y12 * z13);
		double delta = epsilon * (x12 * z13 - x13 * z12);

		double c1 = 1.0 + gamma * gamma + delta * delta;
		double c2 = -2.0 * (a[2] - alpha * gamma - beta * delta);
		double c3 = a[2] * a[2] - r1sqr + alpha * alpha + beta * beta;

		solveQuadratic(c1, c2, c3);

		if (nearlyEqual(quadraticRoots[0], x[2]))
			altPoint[2] = quadraticRoots[1];
		else if (nearlyEqual(quadraticRoots[1], x[2]))
			altPoint[2] = quadraticRoots[0];
		else
			throw new RuntimeException("find other corner: didn't find original point");

		altPoint[0] = alpha + a[0] + gamma * altPoint[2];
		altPoint[1] = beta + a[1] + delta * altPoint[2];

		logger.debug("Other corner " + Coord.info(altPoint));
	}

	double quadraticRoots[] = new double[2];

	/**
	 * quadraticRoots[] are solutions of ax^2 + bx + c. A check is made for
	 * roots close to zero.
	 * 
	 * @param a
	 * @param b
	 * @param c
	 * @throws GaException
	 */
	void solveQuadratic(double a, double b, double c) {
		double sqrTerm = b * b - 4.0 * a * c;
		if (logger.isDebugEnabled()) {
			logger.debug("Solve Qaud a " + a + " b " + b + " c " + c);
			logger.debug("Solve Qaud b**2-4ac " + sqrTerm);
		}
		if (sqrTerm < .0) {
			double check = b * b / (4 * a * c);
			if (check > 0.9999999 && check < 1.00000001)
				sqrTerm = .0;
			else
				throw new RuntimeException("solveQuadratic: complex root\n");
		}

		double x = Math.sqrt(sqrTerm);
		quadraticRoots[0] = (-b + x) / (2.0 * a);
		quadraticRoots[1] = (-b - x) / (2.0 * a);
	}

	/**
	 * Returns true if two distances are nearly the same. The parameter
	 * CORNER_FLAP_TOL is used to check for this. This routine is used to
	 * determine which of the two alternate points for the free corner
	 * corresponds to the current corner.
	 * 
	 * @param p1
	 * @param p2
	 * @return
	 */
	boolean nearlyEqual(double p1, double p2) {
		double diff = p1 - p2;
		if (diff > CORNER_FLAP_TOL || -diff > CORNER_FLAP_TOL)
			return false;

		return true;
	}

	private final double yVec[] = new double[4], zVec[] = new double[4];

	/**
	 * Flips the free corner. See Payne and Glen for all the gory details
	 * 
	 * @throws GaException
	 */
	public void flipCorner() {
		List<double[]> coords = molecule.getCoords();

		if (logger.isDebugEnabled()) {
			logger.debug("Fliping corner");
			logger.debug("Atom A " + atomA.info());
			logger.debug("Atom B " + atomB.info());
			logger.debug("Atom X " + atomX.info());
			logger.debug("Atom C " + atomC.info());
			logger.debug("Atom D " + atomD.info());
		}

		findOtherCorner();

		for (int i = 0; i < 4; i++) {
			Atom otherAtom = null;
			if (i == 0)
				otherAtom = atomA;
			else if (i == 1)
				otherAtom = atomB;
			else if (i == 2)
				otherAtom = atomC;
			else if (i == 3)
				otherAtom = atomD;

			double d1 = Math.sqrt(Coord.sqrDistance(altPoint,
					coords.get(otherAtom.getNo())));
			double d2 = Math.sqrt(Coord.sqrDistance(coords.get(atomX.getNo()),
					coords.get(otherAtom.getNo())));
			double diff = d1 - d2;
			logger.debug("Point distance check " + diff);
			if (diff > CORNER_FLAP_TOL || -diff > CORNER_FLAP_TOL) {
				logger.debug("corner not flipped (" + i + ")");
				return;
			}
		}

		Coord.copy(coords.get(atomX.getNo()), yVec);
		Coord.copy(coords.get(atomB.getNo()), zVec);

		Atom atomY = null, atomZ = null;
		if (logger.isDebugEnabled()) {
			atomY = new Atom(molecule, molecule.getnAtoms(), AtomType.Type.DU);
			atomZ = new Atom(molecule, molecule.getnAtoms() + 1, AtomType.Type.DU);
			molecule.addAtom(atomY, yVec);
			molecule.addAtom(atomZ, zVec);
			Bond b1 = new Bond(molecule.getnBonds(), atomX, atomZ, "1");
			Bond b2 = new Bond(molecule.getnBonds() + 1, atomB, atomY, "1");
			molecule.addBond(b1);
			molecule.addBond(b2);
			coords = molecule.getCoords();
		}

		double rot1 = Coord.torsion(altPoint, coords.get(atomA.getNo()),
				coords.get(atomB.getNo()), yVec);
		logger.debug("rotationAB " + rot1);
		rotationAB.rotateBond(-rot1);
		Coord.transPointInPlace(rotationAB.getTransMatrix(), yVec);

		if (logger.isDebugEnabled()) {
			double dist = Coord.distance(altPoint, yVec);
			if (dist > CORNER_FLAP_TOL)
				logger.debug("flipCorner: (1) point mismatch, distance " + dist);
			else
				logger.debug("flipCorner: (1) distance check " + dist);
		}

		double rot2 = Coord.torsion(altPoint, coords.get(atomC.getNo()),
				coords.get(atomD.getNo()), coords.get(atomX.getNo()));
		logger.debug("rotationCD " + rot2);
		rotationCD.rotateBond(rot2);
		Coord.transPointInPlace(rotationCD.getTransMatrix(), zVec);

		if (logger.isDebugEnabled()) {
			double dist = Coord.distance(altPoint, coords.get(atomX.getNo()));
			if (dist > CORNER_FLAP_TOL)
				logger.debug("flipCorner: (2) point mismatch, distance " + dist);
			else
				logger.debug("flipCorner: (2) distance check " + dist);
		}

		double rot3 = Coord.torsion(zVec, coords.get(atomX.getNo()),
				coords.get(atomC.getNo()), coords.get(atomB.getNo()));
		logger.debug("rotationXC " + rot3);
		rotationXC.rotateBond(-rot3);
		Coord.transPointInPlace(rotationXC.getTransMatrix(), zVec);

		if (logger.isDebugEnabled()) {
			double dist = Coord.distance(coords.get(atomB.getNo()), zVec);
			if (dist > CORNER_FLAP_TOL)
				logger.debug("flipCorner: (3) point mismatch, distance " + dist);
			else
				logger.debug("flipCorner: (3) distance check " + dist);
			Coord.copy(zVec, coords.get(atomZ.getNo()));
			Coord.copy(yVec, coords.get(atomY.getNo()));
		}
	}

	void addTorsions(List<Torsion> v, Molecule m) {
		for (Torsion torsion : m.getTorsions()) {
			Bond bond = null;
			for (int j = 0; j < 4; j++) {
				if (j == 0)
					bond = bondAB;
				else if (j == 1)
					bond = bondBX;
				else if (j == 2)
					bond = bondXC;
				else if (j == 3)
					bond = bondCD;
			}
			if (torsion.getBond() == bond && !v.contains(torsion)) {
				v.add(torsion);
			}
		}
	}

	private class CornerRotation extends RotatableBond {
		private final Atom checkAtom, rootAtom;

		CornerRotation(Atom root, Bond b, Atom check) {
			checkAtom = check;
			rootAtom = root;

			setBond(b);
			if (b.getAtom1() == root) {
				setAtom1(b.getAtom1());
				setAtom2(b.getAtom2());
			} else {
				setAtom1(b.getAtom2());
				setAtom2(b.getAtom1());
			}
			setMolecule(FreeCorner.this.molecule);
			setAtom1List(FreeCorner.this.molecule.findNodes(rootAtom, b,
					Arrays.asList(new Atom[] { checkAtom })));
			setAtom2List(FreeCorner.this.molecule.findNodes(getAtom2(), b,
					Arrays.asList(new Atom[] { checkAtom })));

		}

		private double[][] getTransMatrix() {
			return getTrans().get();
		}

	}

	public boolean isFreeCorner() {
		return freeCorner;
	}

	public CornerRotation getRotationAB() {
		return rotationAB;
	}

	public CornerRotation getRotationCD() {
		return rotationCD;
	}

	public CornerRotation getRotationXC() {
		return rotationXC;
	}

}
