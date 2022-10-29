package com.cairn.gape.utils;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.util.FastMath;
import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;
import com.cairn.gape.molecule.GaMolecule;
import com.cairn.molecule.Atom;
import com.cairn.molecule.AtomType;
import com.cairn.molecule.Molecule;

/**
 * @author Gareth Jones
 * 
 *         Rewritten for efficiency so that the Gaussian object is no longer
 *         used.
 * 
 *         G = n * exp(-alpha(r-rN)
 * 
 */
public class GaussianList {
	// number of gaussians in the list
	private volatile int no;

	// n[] stores the gaussian normalization parameter
	// alpha[] the width exponent and rN[][] the gaussian centers
	// intersetVolume is continuously modified when intersecting two Gaussian
	// Lists.
	private final double n[], alpha[], rN[][];

	public static final boolean USE_EXP_LOOKUP = true;

	private static final double VOL_CUTOFF = 1.0e-2, EXP_LOOKUP_MIN = -20;

	private static final int EXP_LOOKUP_SIZE = 10000;

	private volatile static double expLookup[];

	private static final Logger logger = Logger.getLogger(GaussianList.class);

	// private double intersectVolume;

	static {

		// create exponential and sqrt lookup tables

		if (USE_EXP_LOOKUP) {
			expLookup = new double[EXP_LOOKUP_SIZE];
			for (int i = 0; i < EXP_LOOKUP_SIZE; i++) {
				double s = i * EXP_LOOKUP_MIN / EXP_LOOKUP_SIZE;
				expLookup[i] = FastMath.exp(s);
			}

			if (logger.isDebugEnabled()) {
				for (int i = 0; i < EXP_LOOKUP_SIZE; i++) {
					double s = i * EXP_LOOKUP_MIN / EXP_LOOKUP_SIZE;
					double check = FastMath.exp(s);
					double val = lookupExp(s);
					double diff = check - val;
					logger.debug("Real " + check + " lookup " + val + " diff " + diff);
				}
			}
		}

	}

	/**
	 * @param v
	 * @return exponent of v using lookup tables.
	 */
	static public double lookupExp(double v) {
		assert USE_EXP_LOOKUP;

		if (v > 0)
			return 1;
		if (v <= EXP_LOOKUP_MIN)
			return 0;
		double p = v * EXP_LOOKUP_SIZE / EXP_LOOKUP_MIN + 0.5;
		if (p >= EXP_LOOKUP_SIZE)
			return .0;
		return expLookup[(int) p];
	}

	/**
	 * Create storage for a new list.
	 * 
	 * @param size
	 */
	private GaussianList(int size) {
		no = 0;
		n = new double[size];
		alpha = new double[size];
		rN = new double[size][4];
	}

	/**
	 * @param n
	 *            normalization factor
	 * @param alpha
	 *            width exponent
	 * @return The volumn of a Gaussian
	 */
	static private double gaussianVolume(double n, double alpha) {
		double f = Math.PI / alpha;
		double vol = n * f * Math.sqrt(f);
		return vol;
	}

	/**
	 * @param n1
	 *            normalization factor of first Gaussian
	 * @param alpha1
	 *            width exponent of first Gaussian
	 * @param rN1
	 *            center of first Gaussian
	 * @param n2
	 *            normalization factor of second Gaussian
	 * @param alpha2
	 *            width exponent of second Gaussian
	 * @param rN2
	 *            ceneter of second Gaussian
	 * @return The intersection volumn of two Gaussians
	 */
	static private double intersectionVolume(double n1, double alpha1, double rN1[],
			double n2, double alpha2, double rN2[]) {
		double rSqr = Coord.sqrDistance(rN1, rN2);
		double val = -((alpha1 * alpha2) / (alpha1 + alpha1)) * rSqr;
		double n = .0;
		if (USE_EXP_LOOKUP)
			n = n1 * n2 * lookupExp(val);
		else
			n = n1 * n2 * FastMath.exp(val);
		double alpha = alpha1 + alpha2;
		return gaussianVolume(n, alpha);
	}

	/**
	 * @param n1
	 *            normalization factor
	 * @param alpha1
	 *            width exponent
	 * @param rN1
	 *            center
	 * @return The total overlap volume of all Gaussian in this list with a
	 *         single Gaussian.
	 */
	public double overlapVolume(double n1, double alpha1, double rN1[]) {
		double vol = 0;
		for (int i = 0; i < no; i++) {
			vol += intersectionVolume(n[i], alpha[i], rN[i], n1, alpha1, rN1);
		}
		logger.debug("volume " + vol);
		return vol;
	}

	private final double r1Scale[] = new double[4], r2Scale[] = new double[4];

	/**
	 * Intersects two Gausians to create ane Gaussian in this list. The Gausian
	 * is only created if there is a minimum overlap. Increments Gaussian count
	 * if a new Gaussian is created. Increments intersectVolume.
	 * 
	 * @param n1
	 *            normalization factor of first Gaussian
	 * @param alpha1
	 *            width exponent of first Gaussian
	 * @param rN1
	 *            center of first Gaussian
	 * @param n2
	 *            normalization factor of second Gaussian
	 * @param alpha2
	 *            width exponent of second Gaussian
	 * @param rN2
	 *            center of second Gaussian
	 */
	private void intersectGaussians(double n1, double alpha1, double rN1[], double n2,
			double alpha2, double rN2[]) {
		double rSqr = Coord.sqrDistance(rN1, rN2);

		double val = -((alpha1 * alpha2) / (alpha1 + alpha1)) * rSqr;
		double n3 = .0;
		if (USE_EXP_LOOKUP)
			n3 = n1 * n2 * lookupExp(val);
		else
			n3 = n1 * n2 * FastMath.exp(val);
		double alpha3 = alpha1 + alpha2;
		double vol = gaussianVolume(n3, alpha3);
		// intersectVolume += vol;

		if (vol < VOL_CUTOFF)
			return;

		// new Gaussian parameters
		n[no] = n3;
		alpha[no] = alpha3;

		Coord.copy(rN1, r1Scale);
		Coord.copy(rN2, r2Scale);
		Coord.multiply(r1Scale, alpha1);
		Coord.multiply(r2Scale, alpha2);
		Coord.add(r1Scale, r2Scale, rN[no]);
		// new gaussian center
		Coord.multiply(rN[no], 1 / (alpha1 + alpha2));
		// increment no of Gaussians
		no++;
	}

	/**
	 * Intersects two Gaussians and returns parameters of Gaussian created by
	 * intersection.
	 * 
	 * @param n1
	 *            normalization factor of first Gaussian
	 * @param alpha1
	 *            width exponent of first Gaussian
	 * @param rN1
	 *            center of first Gaussian
	 * @param n2
	 *            normalization factor of second Gaussian
	 * @param alpha2
	 *            width exponent of second Gaussian
	 * @param rN2
	 *            center of second Gaussian
	 * @param rN3
	 *            The center of the intersection gaussian will be put here.
	 * @return a array containing normalization factor, width exponent and
	 *         volume of intersection Gaussian
	 */
	@Deprecated
	public static double[] intersectGaussians(double n1, double alpha1, double rN1[],
			double n2, double alpha2, double rN2[], double rN3[]) {
		double rSqr = Coord.sqrDistance(rN1, rN2);
		double val = -((alpha1 * alpha2) / (alpha1 + alpha1)) * rSqr;
		double n3 = .0;
		if (USE_EXP_LOOKUP)
			n3 = n1 * n2 * lookupExp(val);
		else
			n3 = n1 * n2 * FastMath.exp(val);
		double alpha3 = alpha1 + alpha2;
		double vol = gaussianVolume(n3, alpha3);

		double rtn[] = new double[3];
		rtn[0] = n3;
		rtn[1] = alpha3;
		rtn[2] = vol;

		double r1Scale[] = new double[4], r2Scale[] = new double[4];
		Coord.copy(rN1, r1Scale);
		Coord.copy(rN2, r2Scale);
		Coord.multiply(r1Scale, alpha1);
		Coord.multiply(r2Scale, alpha2);
		Coord.add(r1Scale, r2Scale, rN3);
		Coord.multiply(rN3, 1 / (alpha1 + alpha2));

		return rtn;
	}

	/**
	 * @return a new Gaussian list by intersecting this list with itself. Each
	 *         Gaussian is intersected with every other Gaussian, but not with
	 *         itself. Use to determine corrections in volume calculations.
	 */
	public GaussianList intersection() {
		GaussianList overlay = new GaussianList(no * no / 2);

		logger.debug("Max No " + no * no / 2);
		for (int i = 0; i < no; i++) {
			for (int j = i + 1; j < no; j++) {
				overlay.intersectGaussians(n[i], alpha[i], rN[i], n[j], alpha[j], rN[j]);
			}
		}
		logger.debug("Actual no " + overlay.no);
		return overlay;
	}

	/**
	 * @param b
	 * @return a new Gaussian list by intersecting this list with another
	 *         Gaussian list. Each Gaussian in this list is intersected with
	 *         every Gaussian in the other list. Used a a first step in
	 *         determining intersection volumes.
	 */
	public GaussianList intersection(GaussianList b) {
		GaussianList overlay = new GaussianList(no * b.no);

		for (int i = 0; i < no; i++)
			for (int j = 0; j < b.no; j++) {
				overlay.intersectGaussians(n[i], alpha[i], rN[i], b.n[j], b.alpha[j],
						b.rN[j]);
			}
		logger.debug("no " + overlay.no);
		return overlay;
	}

	/**
	 * @param b
	 * @return The intersection volume (simple sum of volume of intersection
	 *         Gaussians) determined by intersecting each Gaussian in this list
	 *         is intersected with every Gaussian in the other list.
	 */
	public double overlapVolume(GaussianList b) {
		double val = .0;
		for (int i = 0; i < no; i++)
			for (int j = 0; j < b.no; j++) {
				val += intersectionVolume(n[i], alpha[i], rN[i], b.n[j], b.alpha[j],
						b.rN[j]);
			}
		logger.debug("no " + b.no);
		return val;
	}

	// GaussianList intersection(GaussianList a, GaussianList b) {
	// GaussianList overlay = new GaussianList(no*b.no+no+a.no);

	// int pos = 0;
	// for (int i=0; i<no; i++)
	// for (int j=0; j<a.no; j++) {

	// Gaussian other = a.gaussians[j];
	// if (gaussians[i].parent1 == other) continue;
	// if (gaussians[i].parent2 == other) continue;

	// Gaussian g = gaussians[i].overlap(other);
	// if (g.volume > VOL_CUTOFF) {
	// overlay.gaussians[pos] = g;
	// pos++;
	// }
	// }

	// for (int i=0; i<no; i++)
	// for (int j=0; j<b.no; j++) {

	// Gaussian other = b.gaussians[j];
	// if (gaussians[i].parent1 == other) continue;
	// if (gaussians[i].parent2 == other) continue;

	// Gaussian g = gaussians[i].overlap(other);
	// if (g.volume > VOL_CUTOFF) {
	// overlay.gaussians[pos] = g;
	// pos++;
	// }
	// }

	// if (DEBUG) System.out.println("no "+pos);
	// overlay.no = pos;
	// return overlay;
	// }

	/**
	 * @param mol
	 * @return A list of atomic gaussians for this molecule.
	 */
	public static GaussianList atomicGaussians(GaMolecule mol) {
		GaussianList list = new GaussianList(mol.getnAtoms());
		int no = 0;
		for (Atom atom : mol.getAtoms()) {
			if (!atom.isNotDummy())
				continue;
			list.n[no] = AtomType.ATOMIC_N;
			list.alpha[no] = atom.getType().getGaussianWidth();
			list.rN[no] = mol.getCoord(atom.getNo());
			no++;
		}
		list.no = no;
		logger.debug("Atomic list no " + list.no);
		return list;
	}

	/**
	 * @param mol
	 * @return A list of atomic gaussians for this molecule. Restricts gaussians
	 *         to carbon and halogens.
	 */
	public static GaussianList hydrophobicAtomicGaussians(GaMolecule mol) {
		GaussianList list = new GaussianList(mol.getnAtoms());
		int no = 0;
		for (Atom atom : mol.getAtoms()) {
			if (!atom.isNotDummy())
				continue;
			// restrict overlay to hydrophobic atoms-
			AtomType type = atom.getType();
			if (type.isCarbonType() || type.isHalogen()) {
				list.n[no] = AtomType.ATOMIC_N;
				list.alpha[no] = atom.getType().getGaussianWidth();
				list.rN[no] = mol.getCoord(atom.getNo());
				no++;
			}
		}
		list.no = no;
		logger.debug("Atomic list no " + list.no);
		return list;
	}

	/**
	 * @return sum of volumes of all Gaussians in this list.
	 */
	public double totalVolume() {
		double vol = .0;
		for (int i = 0; i < no; i++)
			vol += gaussianVolume(n[i], alpha[i]);
		return vol;
	}

	/**
	 * @return a molecular represtation of this list of Gaussians so we can
	 *         display in a Molecular viewer. Each gaussian center is a dummy
	 *         atom.
	 */
	private Molecule toMolecule() {
		List<Atom> atoms = new ArrayList<>();
		List<double[]> coords = new ArrayList<>();
		for (int i = 0; i < no; i++) {
			String l = "N_" + n[i] + "_A_" + alpha[i];
			atoms.add(new Atom(i, l, "Du"));
			coords.add(rN[i]);
		}
		Molecule mol = new Molecule("Gaussian List", atoms, coords);
		return mol;
	}

	/**
	 * Saves this list of Gaussians in MOL2 format to a file. Each Gassian
	 * center is represented by a dummy atom.
	 * 
	 * @param name
	 */
	private void saveAsMol2(String name) {
		Molecule mol = toMolecule();
		mol.writeSybylMol2File(name, "Gaussian file");
	}

	/**
	 * Test routine. Takes two molecules and determines common volume using
	 * atomic aussian overlap.
	 * 
	 * @param args
	 */
	public static void main(String args[]) {
		GaMolecule molA, molB;

		if (args.length != 2) {
			System.out.println("Usage: GaussianList <molA.mol2> <molB.mol2>");
			System.exit(0);
		}
		molA = new GaMolecule();
		molA.loadFile(args[0]);
		molB = new GaMolecule();
		molB.loadFile(args[1]);

		GaussianList atoms = GaussianList.atomicGaussians(molA);
		GaussianList otherAtoms = GaussianList.atomicGaussians(molB);

		double volumeIntegral = molA.volumeIntegral(molB);
		System.out.println("Volume Integral " + volumeIntegral);

		atoms.saveAsMol2("molA_gaussians.mol2");
		otherAtoms.saveAsMol2("molB_gaussians.mol2");

		if (false) {
			// run this to visualize overlap Gaussians.
			GaussianList overlay = atoms.intersection(otherAtoms);
			overlay.saveAsMol2("overlay_gaussians.mol2");
			double volume = overlay.totalVolume();
			System.out.println("1st order Gaussian Integral " + volume);
		}

		double volume = atoms.overlapVolume(otherAtoms);
		logger.debug("1st order Gaussian Integral " + volume);

		// haven't figured out how to get all the higher order contribs yet!!

		// GAPE gaussianIntegral method stop here

		GaussianList overlay = atoms.intersection(otherAtoms);
		volume = overlay.totalVolume();
		// Subtract 2 order intersections
		overlay = overlay.intersection();
		double diff = overlay.totalVolume();
		volume -= diff;
		logger.debug("2nd order Gaussian Integral " + diff + " new vol " + volume);

		// Add 3 order intersections
		overlay = overlay.intersection();
		diff = overlay.totalVolume();
		volume += diff;
		logger.debug("3rd order Gaussian Integral " + diff + " new vol " + volume);

		logger.debug("Gaussian Integral " + volume);

		// Subtract 2 order intersections- still haven't figured out
		// how to do all this!!

		// overlay = overlay.intersection(atoms, otherAtoms);
		// overlay.saveAsMol2("second_overlay_gaussians.mol2");
		// double diff = overlay.totalVolume();
		// volume -= diff;
		// System.out.println
		// ("2nd order Gaussian Integral "+diff+" new vol "+volume);

	}
}
