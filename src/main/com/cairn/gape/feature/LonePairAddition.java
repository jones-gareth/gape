package com.cairn.gape.feature;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.mutable.MutableBoolean;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Precision;
import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;
import com.cairn.gape.molecule.GaMolecule;
import com.cairn.molecule.Atom;
import com.cairn.molecule.AtomType;
import com.cairn.molecule.AtomType.Type;
import com.cairn.molecule.Bond;
import com.cairn.molecule.Molecule;

/**
 * Class to add Lone Pairs to acceptors.
 * 
 * The routines use geometric transformations to position the acceptor at the
 * origin and a heavy atom (bonded to the acceptor) along the z-axis. Lone pairs
 * were than added to complete the desired geometry (sp1, sp2 or sp3).
 * 
 * @author Gareth Jones
 * 
 */
public class LonePairAddition {
	private static final Logger logger = Logger.getLogger(LonePairAddition.class);
	static {
		// logger.setLevel(Level.DEBUG);
	}
	private volatile double lonePairLength = 2.9;

	private volatile Molecule molecule;

	private Atom currentAcceptor;

	static private final ThreadLocal<double[]> dummy = ThreadLocal
			.withInitial(() -> new double[4]);

	static private final ThreadLocal<double[]> dummy1 = ThreadLocal
			.withInitial(() -> new double[4]);

	static private final ThreadLocal<double[]> dummy2 = ThreadLocal
			.withInitial(() -> new double[4]);

	static private final ThreadLocal<double[]> centroid = ThreadLocal
			.withInitial(() -> new double[4]);

	static private final ThreadLocal<double[]> lonePair1 = ThreadLocal
			.withInitial(() -> new double[4]);

	static private final ThreadLocal<double[]> lonePair2 = ThreadLocal
			.withInitial(() -> new double[4]);

	static private final ThreadLocal<double[]> lonePair3 = ThreadLocal
			.withInitial(() -> new double[4]);

	static private final ThreadLocal<double[]> lpVec = ThreadLocal
			.withInitial(() -> new double[4]);

	static private final ThreadLocal<double[]> rescaleAtom2 = ThreadLocal
			.withInitial(() -> new double[4]);

	static private final ThreadLocal<double[]> rescaleAtom3 = ThreadLocal
			.withInitial(() -> new double[4]);

	static private final ThreadLocal<double[]> rescaleAtom4 = ThreadLocal
			.withInitial(() -> new double[4]);

	static private final ThreadLocal<double[][]> in = ThreadLocal
			.withInitial(() -> new double[4][4]);

	static private final ThreadLocal<double[][]> out = ThreadLocal
			.withInitial(() -> new double[4][4]);

	static private final ThreadLocal<double[][]> rot = ThreadLocal
			.withInitial(() -> new double[4][4]);

	/**
	 * Should not be initialized outside of this class
	 */
	private LonePairAddition() {

	}

	/**
	 * Set the molecule for current operations
	 * 
	 * @param m
	 */
	private void setMolecule(Molecule m) {
		molecule = m;
	}

	// single instance to do the work
	static private ThreadLocal<LonePairAddition> lonePairAddidion = new ThreadLocal<LonePairAddition>() {
		@Override
		protected LonePairAddition initialValue() {
			return new LonePairAddition();
		}
	};

	public static void addLonePairs(Molecule m) {
		addLonePairs(m, 2.9);
	}

	/**
	 * Adds lone pairs to molecule.
	 * 
	 * @param m
	 * @param distance
	 */
	public static void addLonePairs(Molecule m, double distance) {
		logger.debug("Molecule " + m.getName() + " nAtoms " + m.getnAtoms()
				+ " nLonePairs " + m.getnLonePairs());
		lonePairAddidion.get().setMolecule(m);
		lonePairAddidion.get().lonePairLength = distance;
		lonePairAddidion.get().addLonePairs();
		logger.debug("Molecule " + m.getName() + " nAtoms " + m.getnAtoms()
				+ " nLonePairs " + m.getnLonePairs());
	}

	/**
	 * Regenerates lone pair co-ordinates
	 * 
	 * @param m
	 */
	public static void updateLonePairs(Molecule m) {
		lonePairAddidion.get().setMolecule(m);
		lonePairAddidion.get().updateLonePairs();
	}

	/**
	 * Returns true if this atom is an acceptor
	 * 
	 * @param a
	 * @return
	 */
	static public boolean isAcceptor(Atom a) {

		HydrogenBondingType hBondType = a.getAcceptorType();

		if (hBondType != null) {
			logger.debug(a.info() + " is an acceptor");
			return true;
		}

		return false;
	}

	/**
	 * @param a
	 * @return the geometry of an acceptor atom
	 */
	static public HydrogenBondingType.AcceptorGeometry getAcceptorGeometry(Atom a) {
		HydrogenBondingType hBondType = a.getAcceptorType();

		HydrogenBondingType.AcceptorGeometry geometry = null;
		if (hBondType != null)
			geometry = hBondType.getGeometry();

		logger.debug(a.info() + " has geometry " + geometry);
		return geometry;
	}

	/**
	 * Routine for test and debugging. Add lone pairs to a single molecule and
	 * saves structure.
	 * 
	 * @param args
	 */
	public static void main(String args[]) {

		if (args.length != 1) {
			System.out.println("Usage: LonePairAddition <molfile>");
			System.exit(0);
		}

		List<GaMolecule> mols = GaMolecule.loadFiles(args, Molecule.FileType.UNKNOWN,
				Molecule.Source.FILE);

		HydrogenBondingType.loadParameters();
		mols.forEach(mol -> {
			HydrogenBondingType.searchMolecule(mol);
			addLonePairs(mol);
			mol.writeSybylMol2File(mol.getName() + " lone pairs.mol2", "Lone Pairs Added");
		});
	}

	/**
	 * Adds one LP to atom1 on a line that runs from atom2 through atom1.
	 * 
	 * @param atom1
	 * @param atom2
	 * @param lp
	 * @param lonePairLength
	 */
	static public void addOnePairToLinear(double atom1[], double atom2[], double lp[],
			double lonePairLength) {
		Coord.subtract(atom1, atom2, lpVec.get());
		Coord.setLength(lpVec.get(), lonePairLength);
		Coord.add(atom1, lpVec.get(), lp);
		lp[3] = 1.0;
	}

	/**
	 * origin is the donor atom and atom1 is bonded to it. Three lone pairs are
	 * added to complete sp3 geometry.
	 * 
	 * @param origin
	 * @param atom1
	 * @param lp1
	 * @param lp2
	 * @param lp3
	 * @param lonePairLength
	 */
	static public void addThreePairsToTetrahedral(double origin[], double atom1[],
			double lp1[], double lp2[], double lp3[], double lonePairLength) {

		// Set up transformation that tranforms donor to origin and atom2
		// along z-axis
		setUpTrans(origin, atom1, in.get(), out.get());
		// Set a dummy atom in xz plane such that distance from origin to
		// dummy atom is the length of the lone pair and the angle from
		// atom1-origin-dummy atom is tetrahedral (109.47)

		dummy.get()[0] = lonePairLength * Math.cos((19.47 / 180 * Math.PI));
		dummy.get()[1] = 0.0;
		dummy.get()[2] = lonePairLength * Math.sin((19.47 / 180 * Math.PI));
		dummy.get()[3] = 1.0;

		// generate 1st lone pair at 0 degrees
		double angle = 0;
		setUpRotation(angle, rot.get());
		Coord.transPoint(rot.get(), dummy.get(), lonePair1.get());
		angle = 2 * Math.PI / 3;
		// rotate by 120 degrees to generate 2nd and 3rd lone pairs.
		setUpRotation(angle, rot.get());
		Coord.transPoint(rot.get(), lonePair1.get(), lonePair2.get());
		Coord.transPoint(rot.get(), lonePair2.get(), lonePair3.get());

		// Apply inverse transformations to generate the final
		// lone pair co-ordinates
		Coord.transPoint(out.get(), lonePair1.get(), lp1);
		Coord.transPoint(out.get(), lonePair2.get(), lp2);
		Coord.transPoint(out.get(), lonePair3.get(), lp3);
	}

	/**
	 * Acceptor is at origin and is bonded to atom1 and atom2. Two lone pairs
	 * are added to complete the sp3 geometry.
	 * 
	 * @param origin
	 * @param atom1
	 * @param atom2
	 * @param lp1
	 * @param lp2
	 * @param lonePairLength
	 */
	static public void addTwoPairsToTetrahedral(double origin[], double atom1[],
			double atom2[], double lp1[], double lp2[], double lonePairLength) {

		// set up the the transformation that moves the acceptor to the
		// origin and atom1 on the z-axis
		setUpTrans(origin, atom1, in.get(), out.get());

		// dummy is the co-ordinates of atom2 under this transformation
		Coord.transPoint(in.get(), atom2, dummy.get());
		/* find the angle dummy makes with the zy-axis */
		double angle = Math.asin(dummy.get()[0]
				/ (Math.sqrt((dummy.get()[0] * dummy.get()[0] + dummy.get()[1]
						* dummy.get()[1]))));
		if (dummy.get()[1] < 0.0)
			angle = Math.PI - angle;
		angle *= -1.0;
		// dummy1 is dummy rotated onto the zy plane
		setUpRotation(angle, rot.get());
		Coord.transPoint(rot.get(), dummy.get(), dummy1.get());
		// atom1 lies along -ve z-axis; atom2 should be in +ve quadrant
		// of zy plane

		// Alpha is the angle between atom1-acceptor-atom2
		double alpha = Math.PI
				/ 2
				+ FastMath.asin(dummy1.get()[2]
						/ Math.sqrt(dummy1.get()[1] * dummy1.get()[1] + dummy1.get()[2]
								* dummy1.get()[2]));
		// Now add the two lone pairs in a plane perpendicular to the
		// zy plane that passes through the origin and bisects
		// alpha. Beta is the angle out between this plane and the y
		// axis.
		double beta = 0.5 * (2.0 * Math.PI - alpha) - 0.5 * Math.PI;
		// The angle between the two lone pairs is 109.47 so we can determine
		// their projections in the x and zy planes, as one lone pair will be
		// added above the zy plane and one lone pair will be added an equal
		// distance below it
		double xLen = lonePairLength * Math.sin(54.735 * Math.PI / 180);
		double zyLen = lonePairLength * Math.cos(54.735 * Math.PI / 180);
		lonePair1.get()[0] = xLen;
		lonePair2.get()[0] = -xLen;
		// The projections of the lone pairs on the z and y axis are
		// determined from beta
		lonePair1.get()[1] = lonePair2.get()[1] = -1.0 * zyLen * Math.cos(beta);
		lonePair1.get()[2] = lonePair2.get()[2] = 1.0 * zyLen * Math.sin(beta);
		lonePair1.get()[3] = lonePair2.get()[3] = 1.0;

		/* Move lone pairs back to the original co-ordinate system. */
		setUpRotation(-1.0 * angle, rot.get());
		Coord.transPoint(rot.get(), lonePair1.get(), dummy.get());
		Coord.transPoint(rot.get(), lonePair2.get(), dummy1.get());

		Coord.transPoint(out.get(), dummy.get(), lp1);
		Coord.transPoint(out.get(), dummy1.get(), lp2);
	}

	/**
	 * Add one lone pair to an acceptor bonded to three atoms in sp3 geometry.
	 * The acceptor is at co-ordinates atom and three atoms (atom2, atom3 and
	 * atom4) are bonded to it. If C is the centroid of atoms 2, 3 and 4 and C'
	 * is reflection through the origin, O, then the lone pair is added on OC'
	 * 
	 * @param atom
	 * @param atom2
	 * @param atom3
	 * @param atom4
	 * @param lp
	 * @param lonePairLength
	 */
	static public void addOnePairToTetrahedral(double atom[], double atom2[],
			double atom3[], double atom4[], double lp[], double lonePairLength) {
		logger.debug("adding one pair to tetrahedral");

		// rescale bondLengths to 1
		rescaleBonds(atom, atom2, rescaleAtom2.get());
		rescaleBonds(atom, atom3, rescaleAtom3.get());
		rescaleBonds(atom, atom4, rescaleAtom4.get());

		if (logger.isDebugEnabled()) {
			logger.debug("Atom " + Coord.info(atom));
			logger.debug("Atom2 " + Coord.info(atom2));
			logger.debug("Atom3 " + Coord.info(atom3));
			logger.debug("Atom4 " + Coord.info(atom4));
			logger.debug("Rescale atom2 " + Coord.info(rescaleAtom2.get()));
			logger.debug("Rescale atom3 " + Coord.info(rescaleAtom3.get()));
			logger.debug("Rescale atom4 " + Coord.info(rescaleAtom4.get()));
		}

		centroid.get()[0] = (rescaleAtom2.get()[0] + rescaleAtom3.get()[0] + rescaleAtom4
				.get()[0]) / 3.0;
		centroid.get()[1] = (rescaleAtom2.get()[1] + rescaleAtom3.get()[1] + rescaleAtom4
				.get()[1]) / 3.0;
		centroid.get()[2] = (rescaleAtom2.get()[2] + rescaleAtom3.get()[2] + rescaleAtom4
				.get()[2]) / 3.0;
		centroid.get()[3] = 1.0;

		// find lone pair
		Coord.subtract(atom, centroid.get(), lpVec.get());
		Coord.setLength(lpVec.get(), lonePairLength);
		Coord.add(atom, lpVec.get(), lp);
		lp[3] = 1.0;
	}

	/**
	 * Add one lone pair to the sp2 acceptor with co-ordinates origin. The
	 * acceptor is bonded to atoms atom1 and atom2.
	 * 
	 * @param origin
	 * @param atom1
	 * @param atom2
	 * @param lp
	 * @param lonePairLength
	 * @return
	 */
	public static boolean addOnePairToTrigonal(double origin[], double atom1[],
			double atom2[], double lp[], double lonePairLength) {

		logger.debug("Adding one pair to trigonal");

		// Set up the transformation that puts <origin> at the origin
		// and atom1 along the z-axis.
		setUpTrans(origin, atom1, in.get(), out.get());

		// move atom2
		double[] atom2Moved = dummy.get();
		Coord.transPoint(in.get(), atom2, atom2Moved);
		// <angle> is the angle atom2 makes with the zy plane */
		double angle = FastMath.asin(atom2Moved[0]
				/ (Math.sqrt(atom2Moved[0] * atom2Moved[0] + atom2Moved[1]
						* atom2Moved[1])));
		if (atom2Moved[1] < 0.0)
			angle = Math.PI - angle;
		angle *= -1.0;
		// new atom2 is rotated to dummy1 which lies on the zy plane
		setUpRotation(angle, rot.get());
		atom2Moved = dummy1.get();
		Coord.transPoint(rot.get(), dummy.get(), atom2Moved);
		// atom1 lies along -ve z-axis; atom2 should be in +ve
		// quadrant of zy plane
		if (atom2Moved[1] < 0 || atom2Moved[2] < 0)
			return false;

		// alpha is the angle between atom1-origin-atom2
		double alpha = Math.PI
				/ 2
				+ Math.asin(atom2Moved[2]
						/ Math.sqrt(atom2Moved[1] * atom2Moved[1] + atom2Moved[2]
								* atom2Moved[2]));
		// beta is the angle we want the lone pair to make with the
		// z-axis
		double beta = 0.5 * (2.0 * Math.PI - alpha) - 0.5 * Math.PI;
		lonePair1.get()[0] = 0;
		lonePair1.get()[1] = -1.0 * lonePairLength * Math.cos(beta);
		lonePair1.get()[2] = lonePairLength * Math.sin(beta);
		lonePair1.get()[3] = 1.0;

		// move lone pair back
		setUpRotation(-1.0 * angle, rot.get());
		Coord.transPoint(rot.get(), lonePair1.get(), dummy.get());
		Coord.transPoint(out.get(), dummy.get(), lp);
		return true;
	}

	/**
	 * Adds two pair to the acceptor with sp2 geometry at co-ordinates <atom1>.
	 * The acceptor is bonded to <origin> and <origin> is connected to atom2.
	 * 
	 * @param atom1
	 * @param origin
	 * @param atom2
	 * @param lp1
	 * @param lp2
	 * @param lonePairLength
	 * @return
	 */
	public static boolean addTwoPairsToTrigonal(double atom1[], double origin[],
			double atom2[], double lp1[], double lp2[], double lonePairLength) {

		// Set up the transformation that moves origin to the origin
		// and puts atom1 on the z-axis
		setUpTrans(origin, atom1, in.get(), out.get());

		// move atom1 and atom2 to dummy2 and dummy */
		Coord.transPoint(in.get(), atom2, dummy.get());

		Coord.transPoint(in.get(), atom1, dummy2.get());

		if (logger.isDebugEnabled()) {
			logger.debug("addTwoPairsToTrigonal: atom1 moved to "
					+ Coord.info(dummy2.get()));
			logger.debug("addTwoPairsToTrigonal: atom2 moved to "
					+ Coord.info(dummy.get()));
		}

		// angle is the angle atom2 makes with the zy plane
		double angle = FastMath.asin(dummy.get()[0]
				/ (Math.sqrt(dummy.get()[0] * dummy.get()[0] + dummy.get()[1]
						* dummy.get()[1])));
		if (dummy.get()[1] < 0.0)
			angle = Math.PI - angle;
		angle *= -1.0;
		setUpRotation(angle, rot.get());
		Coord.transPoint(rot.get(), dummy.get(), dummy1.get());

		if (logger.isDebugEnabled())
			logger.debug("addTwoPairsToTrigonal: atom2 moved to "
					+ Coord.info(dummy1.get()));

		// atom1 lies along -ve z-axis; atom2 should be in +ve
		// quadrant of zy plane
		if (dummy1.get()[1] < 0 || dummy1.get()[2] < 0)
			return false;

		// can now add lone pairs onto atom1 in the zy plane so that
		// the angle origin-atom1-lp is 120 degrees
		double alpha = Math.PI / 3.0;
		double l = lonePairLength * Math.sin(alpha);
		lonePair1.get()[0] = lonePair2.get()[0] = 0;
		lonePair1.get()[1] = l;
		lonePair2.get()[1] = -l;
		lonePair1.get()[2] = lonePair2.get()[2] = dummy2.get()[2] - lonePairLength
				* Math.cos(alpha);
		lonePair1.get()[3] = lonePair2.get()[3] = 1.0;

		// move lone pairs back
		setUpRotation(-1.0 * angle, rot.get());
		Coord.transPoint(rot.get(), lonePair1.get(), dummy1.get());
		Coord.transPoint(rot.get(), lonePair2.get(), dummy2.get());
		Coord.transPoint(out.get(), dummy1.get(), lp1);
		Coord.transPoint(out.get(), dummy2.get(), lp2);
		return true;
	}

	/**
	 * Adds two lone pairs to the sp2 acceptor atom1 which is bonded to the atom
	 * at co-ordinates origin. The plane in which the lone pairs are added is
	 * undefined
	 * 
	 * @param atom1
	 * @param origin
	 * @param lp1
	 * @param lp2
	 * @param lonePairLength
	 */
	@SuppressWarnings("unused")
	private static void addTwoPairsRandomlyToTrigonal(double atom1[], double origin[],
			double lp1[], double lp2[], double lonePairLength) {
		double alpha = Math.PI / 3.0;
		addTwoPairsRandomlyToTrigonal(atom1, origin, lp1, lp2, alpha, lonePairLength);
	}

	/**
	 * Adds two lone pairs to the sp2 acceptor atom1 which is bonded to the atom
	 * at co-ordinates origin. The plane in which the lone pairs are added is
	 * undefined
	 * 
	 * alpha argument is double the angle between the lone pairs, so we can use
	 * this transform for the cone feature
	 * 
	 * @param atom1
	 * @param origin
	 * @param lp1
	 * @param lp2
	 * @param alpha
	 * @param lonePairLength
	 */
	static void addTwoPairsRandomlyToTrigonal(double atom1[], double origin[],
			double lp1[], double lp2[], double alpha, double lonePairLength) {
		// get the transformation that moves <origin> to the origin
		// and atom1 to the z-axis
		setUpTrans(origin, atom1, in.get(), out.get());

		// Dummy is atom1 transformed
		Coord.transPoint(in.get(), atom1, dummy.get());
		// Atom1 lies along -ve z-axis. Add the lone pairs in the zy
		// plane to atom1 such that the angle lone-pair-atom1-origin
		// is 120 degrees
		double l = lonePairLength * Math.sin(alpha);
		lonePair1.get()[0] = lonePair2.get()[0] = 0;
		lonePair1.get()[1] = l;
		lonePair2.get()[1] = -l;
		lonePair1.get()[2] = lonePair2.get()[2] = dummy.get()[2] - lonePairLength
				* Math.cos(alpha);
		lonePair1.get()[3] = lonePair2.get()[3] = 1.0;

		// Angle could be random- but we just make it 0
		double angle = 0;
		setUpRotation(angle, rot.get());
		Coord.transPoint(rot.get(), lonePair1.get(), dummy1.get());
		Coord.transPoint(rot.get(), lonePair2.get(), dummy2.get());
		// move lone pairs back
		Coord.transPoint(out.get(), dummy1.get(), lp1);
		Coord.transPoint(out.get(), dummy2.get(), lp2);
	}

	static private final ThreadLocal<double[][]> xRot = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};
	static private final ThreadLocal<double[][]> yRot = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};
	static private final ThreadLocal<double[][]> trans = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};
	static private final ThreadLocal<double[][]> invTrans = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};
	static private final ThreadLocal<double[][]> invXrot = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};
	static private final ThreadLocal<double[][]> invYrot = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};
	static private final ThreadLocal<double[][]> product1 = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};
	static private final ThreadLocal<double[]> diff = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[4];
		}
	};

	/**
	 * Sets up two transformation matrices in and out. In moves point p1 to
	 * origin and alignes points p1 and p2 on z-axis. Out is the reverse
	 * transformation.
	 * 
	 * @param p1
	 * @param p2
	 * @param in
	 * @param out
	 */
	private static void setUpTrans(double p1[], double p2[], double in[][],
			double out[][]) {

		// Pretty much a copy of Coord.determineRotation

		double a, b, c, v, x, y, z;
		x = p1[0];
		y = p1[1];
		z = p1[2];
		Coord.subtract(p1, p2, diff.get());
		Coord.unit(diff.get());
		a = diff.get()[0];
		b = diff.get()[1];
		c = diff.get()[2];

		// trans: translates point one to the origin. xRot and yRot:
		// rotation about x and y axis such that point2 is aligned on
		// the z-axis. invTrans, invXrot and invYrot: inverse
		// translations and rotations.

		// Ref: Newman & Sproull, Principles of Interactive
		// computer Graphics, McGraw Hill, 1981 pp 346-348.

		Coord.identity(trans.get());
		trans.get()[3][0] = -x;
		trans.get()[3][1] = -y;
		trans.get()[3][2] = -z;
		Coord.identity(invTrans.get());
		invTrans.get()[3][0] = x;
		invTrans.get()[3][1] = y;
		invTrans.get()[3][2] = z;
		v = Math.sqrt(b * b + c * c);
		Coord.identity(xRot.get());
		Coord.identity(invXrot.get());
		if (v > .0) {
			xRot.get()[1][1] = xRot.get()[2][2] = invXrot.get()[1][1] = invXrot.get()[2][2] = c
					/ v;
			xRot.get()[1][2] = invXrot.get()[2][1] = b / v;
			xRot.get()[2][1] = invXrot.get()[1][2] = -b / v;
		}
		Coord.identity(yRot.get());
		Coord.identity(invYrot.get());
		yRot.get()[0][0] = yRot.get()[2][2] = invYrot.get()[0][0] = invYrot.get()[2][2] = v;
		yRot.get()[0][2] = invYrot.get()[2][0] = a;
		yRot.get()[2][0] = invYrot.get()[0][2] = -a;

		// Form in and out from matrix products of component translation
		// and rotations.
		Coord.product(trans.get(), xRot.get(), product1.get());
		Coord.product(product1.get(), yRot.get(), in);
		Coord.product(invYrot.get(), invXrot.get(), product1.get());
		Coord.product(product1.get(), invTrans.get(), out);
	}

	/**
	 * Initialises the matrix ,rot, to perform a rotation of angle radians about
	 * the z-axis.
	 * 
	 * @param angle
	 * @param rot
	 */
	private static void setUpRotation(double angle, double rot[][]) {
		Coord.identity(rot);
		// Normal format of z-axis rotation matrix
		rot[2][2] = rot[3][3] = 1;
		rot[0][0] = rot[1][1] = Math.cos(angle);
		rot[0][1] = -1.0 * Math.sin(angle);
		rot[1][0] = Math.sin(angle);
	}

	/**
	 * Returns the true number of lone pairs that an acceptor atom has. Differes
	 * from countLonePairs in that an acceptor with no geometry or cone geometry
	 * will return the true number of lone pairs (probably 3).
	 * 
	 * @param a
	 * @return
	 */
	private static int countTrueLonePairs(Atom a) {

		if (!isAcceptor(a))
			return 0;
		AtomType type = a.getType();

		int nConnections = 0;
		for (Atom neighbour : a.getNotDummyNeighbours()) {
			if (neighbour.getAtomType() != AtomType.Type.LP) {
				nConnections++;
			}
		}

		int total = 0;
		// hardwire nitro oxygens
		if (a.isNitroOxygen())
			total = 3;
		else if (type.getGeometry() == AtomType.Geometry.LIN)
			total = 2;
		else if (type.getGeometry() == AtomType.Geometry.TRI)
			total = 3;
		else if (type.getGeometry() == AtomType.Geometry.TET)
			total = 4;

		return total - nConnections;
	}

	/**
	 * Returns a count of lone pairs that need to be added to an atom. Take
	 * account of a reduced representation whereby an acceptor that accepts in a
	 * cone is represented by one (rather that 3) lone pairs
	 * 
	 * @param a
	 * @return
	 */
	private static int countLonePairs(Atom a) {
		HydrogenBondingType.AcceptorGeometry geometry = getAcceptorGeometry(a);

		if (AcceptorAtomFeature.USE_ACCEPTOR_GEOMETRY) {
			if (geometry == null)
				return 0;

			switch (geometry) {
			case DIR:
				return countTrueLonePairs(a);
			case PLANE:
				int test = countTrueLonePairs(a);
				if (test != 2)
					throw new RuntimeException(
							"Can't create plane lone pair geometry on acceptor atom "
									+ a.info());
				return 2;
			case CONE:
			case AG_NONE:
				// For the cone and none the lone pair is just a
				// representation of forward geometry and not a real
				// lone pair.
				return 1;
			}
		} else {
			return countTrueLonePairs(a);
		}
		return 0;
	}

	/**
	 * Goes though all acceptors in molecule and adds all the lone pairs in the
	 * required geometry.
	 * 
	 */
	private void addLonePairs() {
		int nLonePairs = 0;
		for (Atom a : molecule.getAtoms()) {
			a.setnLonePairs(countLonePairs(a));
			logger.debug(a.info() + " has " + a.getnLonePairs() + " LPs");
			nLonePairs += a.getnLonePairs();
		}

		molecule.setnLonePairs(nLonePairs);

		List<Atom> testAtoms = new ArrayList<>(molecule.getAtoms());
		for (Atom test : testAtoms) {
			currentAcceptor = test;
			if (!isAcceptor(currentAcceptor))
				continue;
			HydrogenBondingType.AcceptorGeometry geometry = getAcceptorGeometry(currentAcceptor);

			logger.debug("Adding LPs to" + currentAcceptor.info());

			addLonePairsToMolecule();

			if (AcceptorAtomFeature.USE_ACCEPTOR_GEOMETRY) {
				switch (geometry) {
				case DIR:
					addLonePairsToAtom();
					break;
				case PLANE:
					addLonePairsToAtom();
					break;
				case CONE:
				case AG_NONE:
					addOnePairToAtom();
					break;
				}
			} else
				addLonePairsToAtom();

			logger.debug(currentAcceptor.info() + " n lone pairs "
					+ currentAcceptor.getnLonePairs());
		}

	}

	/**
	 * Goes through all acceptors in a molecule and puts the lone pairs in the
	 * right place.
	 */
	void updateLonePairs() {
		// Remove Lone Pairs
		int nLonePairs = 0;
		for (Atom atom : new ArrayList<>(molecule.getAtoms())) {
			if (atom.getAtomType() == AtomType.Type.LP) {
				molecule.deleteAtom(atom);
			}
		}

		int no = 0;
		// Re-add lone-pairs
		for (Atom test : new ArrayList<>(molecule.getAtoms())) {
			currentAcceptor = test;
			if (!isAcceptor(currentAcceptor))
				continue;
			HydrogenBondingType.AcceptorGeometry geometry = getAcceptorGeometry(currentAcceptor);

			logger.debug("Adding LPs to" + currentAcceptor.info());
			no += addLonePairsToMolecule();

			if (AcceptorAtomFeature.USE_ACCEPTOR_GEOMETRY) {
				switch (geometry) {
				case DIR:
				case PLANE:
					addLonePairsToAtom();
					break;
				case CONE:
				case AG_NONE:
					addOnePairToAtom();
					break;
				}
			} else
				addLonePairsToAtom();

			int nAcceptorLonePairs = currentAcceptor.getnLonePairs();
			logger.debug(currentAcceptor.info() + " n lone pairs " + nAcceptorLonePairs);
			nLonePairs += nAcceptorLonePairs;
		}

		assert no == nLonePairs : "lone pair count error";
		molecule.setnLonePairs(no);
	}

	void addOnePairToAtom() {
		int no = currentAcceptor.getnNotDummyNeighbours();
		if (no == 1)
			addLpToSp1();
		else if (no == 2)
			add1LpToSp2();
		else if (no == 3)
			add1LpToSp3();
		else
			throw new RuntimeException("Unable to add a single LP to an atom with " + no
					+ " neighbours");
	}

	void addLonePairsToAtom() {
		AtomType.Geometry geometry = currentAcceptor.getType().getGeometry();
		int no = currentAcceptor.getnLonePairs();

		// Hardwire carboxylic acid
		if (currentAcceptor.isCarboxylateOxygen())
			addLpToOco2();
		// hardwire nitro group
		else if (currentAcceptor.isNitroOxygen())
			addLpToOco2();
		else {
			if (geometry == AtomType.Geometry.TET) {
				if (no == 1)
					add1LpToSp3();
				else if (no == 2)
					add2LpToSp3();
				else if (no == 3)
					add3LpToSp3();
			} else if (geometry == AtomType.Geometry.TRI) {
				if (no == 1)
					add1LpToSp2();
				else if (no == 2)
					add2LpToSp2();
			} else if (geometry == AtomType.Geometry.LIN) {
				if (no == 1)
					addLpToSp1();
			} else {
				throw new RuntimeException("Unable to add lone pairs to atom "
						+ currentAcceptor.info());
			}
		}

		assert checkLonePairLengths(currentAcceptor) : "Lone pair addition error";
	}

	/**
	 * Checks that the right number of lone pairs have been added and they are
	 * the right distance from the acceptor
	 * 
	 * @param acceptor
	 * @return
	 */
	private boolean checkLonePairLengths(Atom acceptor) {
		MutableBoolean ok = new MutableBoolean(true);
		MutableInt count = new MutableInt();

		molecule.getBonds()
				.forEach(
						bond -> {
							if ((bond.getAtom1() == acceptor && bond.getAtom2()
									.getAtomType() == Type.LP)
									|| (bond.getAtom2() == acceptor && bond.getAtom1()
											.getAtomType() == Type.LP)) {
								double testLen = Coord.distance(
										molecule.getCoord(bond.getAtom1().getNo()),
										molecule.getCoord(bond.getAtom2().getNo()));
								if (!Precision.equals(testLen, lonePairLength, 0.01)) {
									logger.warn("Lone pair length error: expected "
											+ lonePairLength + " got " + testLen);
									ok.setValue(false);
								}
								count.increment();
							}
						});

		int no = currentAcceptor.getnLonePairs();
		if (count.intValue() != no) {
			logger.warn("Lone pair count error: expected " + no + " got "
					+ count.intValue());
			ok.setValue(false);
		}

		return ok.booleanValue();

	}

	/**
	 * Adds one lone pair to an sp3 acceptor.
	 * 
	 */
	void add1LpToSp3() {

		if (currentAcceptor.getnNotDummyNeighbours() != 3)
			throw new RuntimeException("add1LpToSp3: atom " + currentAcceptor.info()
					+ " has " + currentAcceptor.getnNotDummyNeighbours() + " connections");

		double coord[] = molecule.getCoord(currentAcceptor.getNo());
		double coord2[] = molecule.getCoord(currentAcceptor.getNotDummyNeighbours()
				.get(0).getNo());
		double coord3[] = molecule.getCoord(currentAcceptor.getNotDummyNeighbours()
				.get(1).getNo());
		double coord4[] = molecule.getCoord(currentAcceptor.getNotDummyNeighbours()
				.get(2).getNo());
		double lp[] = molecule.getCoord(molecule.getnAtoms() - 1);

		addOnePairToTetrahedral(coord, coord2, coord3, coord4, lp, lonePairLength);
	}

	/**
	 * Adds one lone pair to an sp2 atom.
	 * 
	 */
	void add1LpToSp2() {

		if (currentAcceptor.getnNotDummyNeighbours() != 2)
			throw new RuntimeException("add1LpToSp2: atom " + currentAcceptor.info()
					+ " has " + currentAcceptor.getnNotDummyNeighbours() + " connections");

		double coord[] = molecule.getCoord(currentAcceptor.getNo());
		double coord2[] = molecule.getCoord(currentAcceptor.getNotDummyNeighbours()
				.get(0).getNo());
		double coord3[] = molecule.getCoord(currentAcceptor.getNotDummyNeighbours()
				.get(1).getNo());
		double lp[] = molecule.getCoord(molecule.getnAtoms() - 1);

		addOnePairToTrigonal(coord, coord2, coord3, lp, lonePairLength);
	}

	/**
	 * Adds one lone pair to an sp1 atom.
	 * 
	 */
	void addLpToSp1() {
		if (currentAcceptor.getnNotDummyNeighbours() != 1)
			throw new RuntimeException("add1LpToSp1: atom " + currentAcceptor.info()
					+ " has " + currentAcceptor.getnNotDummyNeighbours() + " connections");

		double coord[] = molecule.getCoord(currentAcceptor.getNo());
		double coord2[] = molecule.getCoord(currentAcceptor.getNotDummyNeighbours()
				.get(0).getNo());
		double lp[] = molecule.getCoord(molecule.getnAtoms() - 1);

		addOnePairToLinear(coord, coord2, lp, lonePairLength);
	}

	/**
	 * Adds two lone pairs to an sp3 atom.
	 * 
	 */
	void add2LpToSp3() {
		if (currentAcceptor.getnNotDummyNeighbours() != 2)
			throw new RuntimeException("add2LpToSp3: atom " + currentAcceptor.info()
					+ " has " + currentAcceptor.getnNotDummyNeighbours() + " connections");

		double coord[] = molecule.getCoord(currentAcceptor.getNo());
		double coord2[] = molecule.getCoord(currentAcceptor.getNotDummyNeighbours()
				.get(0).getNo());
		double coord3[] = molecule.getCoord(currentAcceptor.getNotDummyNeighbours()
				.get(1).getNo());
		double lp1[] = molecule.getCoord(molecule.getnAtoms() - 1);
		double lp2[] = molecule.getCoord(molecule.getnAtoms() - 2);

		addTwoPairsToTetrahedral(coord, coord2, coord3, lp1, lp2, lonePairLength);
	}

	/**
	 * Adds two lone pairs to an sp2 atom.
	 * 
	 */
	void add2LpToSp2() {

		if (currentAcceptor.getnNotDummyNeighbours() != 1)
			throw new RuntimeException("add2LpToSp3: atom " + currentAcceptor.info()
					+ " has " + currentAcceptor.getnNotDummyNeighbours() + " connections");

		double coord[] = molecule.getCoord(currentAcceptor.getNo());
		Atom otherAtom = currentAcceptor.getNotDummyNeighbours().get(0);
		Atom thirdAtom = null;
		for (Atom test : otherAtom.getNotDummyNeighbours()) {
			if (test != currentAcceptor) {
				thirdAtom = test;
			}
		}
		double coord2[] = molecule.getCoord(otherAtom.getNo());
		double coord3[] = molecule.getCoord(thirdAtom.getNo());
		double lp1[] = molecule.getCoord(molecule.getnAtoms() - 1);
		double lp2[] = molecule.getCoord(molecule.getnAtoms() - 2);

		if (logger.isDebugEnabled()) {
			logger.debug("add2LpToSp2 other atom is " + otherAtom.info());
			logger.debug("add2LpToSp2 third atom is " + thirdAtom.info());
		}
		if (!addTwoPairsToTrigonal(coord, coord2, coord3, lp1, lp2, lonePairLength))
			throw new RuntimeException("failed to add two pairs to trigonal for atom "
					+ currentAcceptor.info());
	}

	/**
	 * Adds two lone pairs to the carboxylic oxygen O.co2 atom. Also used to add
	 * lone pairs to nitro oxygen.
	 * 
	 */
	void addLpToOco2() {

		if (currentAcceptor.getnNotDummyNeighbours() != 1)
			throw new RuntimeException("addLpToOco2: atom " + currentAcceptor.info()
					+ " has " + currentAcceptor.getnNotDummyNeighbours() + " connections");

		double coord[] = molecule.getCoord(currentAcceptor.getNo());
		Atom otherAtom = currentAcceptor.getNotDummyNeighbours().get(0);
		Atom thirdAtom = null;
		for (Atom test : otherAtom.getNotDummyNeighbours()) {
			if (test != currentAcceptor && test.getType().isOxygenType()
					&& test.getnNotDummyNeighbours() == 1) {
				thirdAtom = test;
			}
		}
		if (thirdAtom == null)
			throw new RuntimeException("can't find other oco2 paired to atom "
					+ currentAcceptor.info());

		double coord2[] = molecule.getCoord(otherAtom.getNo());
		double coord3[] = molecule.getCoord(thirdAtom.getNo());
		double lp1[] = molecule.getCoord(molecule.getnAtoms() - 1);
		double lp2[] = molecule.getCoord(molecule.getnAtoms() - 2);

		addTwoPairsToTrigonal(coord, coord2, coord3, lp1, lp2, lonePairLength);
	}

	/**
	 * Adds three lone pairs to an sp3 atom.
	 * 
	 */
	void add3LpToSp3() {

		if (currentAcceptor.getnNotDummyNeighbours() != 1)
			throw new RuntimeException("add3LpToSp3: atom " + currentAcceptor.info()
					+ " has " + currentAcceptor.getnNotDummyNeighbours() + " connections");

		double coord[] = molecule.getCoord(currentAcceptor.getNo());
		double coord2[] = molecule.getCoord(currentAcceptor.getNotDummyNeighbours()
				.get(0).getNo());
		double lp1[] = molecule.getCoord(molecule.getnAtoms() - 1);
		double lp2[] = molecule.getCoord(molecule.getnAtoms() - 2);
		double lp3[] = molecule.getCoord(molecule.getnAtoms() - 3);

		addThreePairsToTetrahedral(coord, coord2, lp1, lp2, lp3, lonePairLength);
	}

	/**
	 * Updates molecule to reflect the addition of a lone pair to
	 * currentAcceptor.
	 */
	int addLonePairsToMolecule() {
		int no = currentAcceptor.getnLonePairs();
		List<Atom> lonePairs = new ArrayList<>();
		for (int i = 0; i < no; i++) {
			Atom lonePair = new Atom(molecule, molecule.getnAtoms(), AtomType.Type.LP);
			Bond bond = new Bond(molecule.getnBonds(), currentAcceptor, lonePair, "1");
			lonePairs.add(lonePair);
			molecule.addAtom(lonePair, new double[4]);
			molecule.addBond(bond);
			lonePair.getNeighbours();
		}
		currentAcceptor.setLonePairs(lonePairs);
		currentAcceptor.getNeighbours();
		return no;
	}

	static private final ThreadLocal<double[]> vec = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[4];
		}
	};

	/**
	 * rescaleAtom2 is on the line between atom1 and atom2, unit distance from
	 * atom1.
	 * 
	 * @param atom1
	 * @param atom2
	 * @param rescaleAtom2
	 */
	static void rescaleBonds(double atom1[], double atom2[], double rescaleAtom2[]) {
		Coord.subtract(atom2, atom1, vec.get());
		Coord.unit(vec.get());
		Coord.add(atom1, vec.get(), rescaleAtom2);
		rescaleAtom2[3] = 1.0;
	}

}
