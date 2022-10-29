package com.cairn.molecule;

import org.apache.commons.math3.util.FastMath;
import org.apache.log4j.Logger;

import com.cairn.common.stats.Simplex;
import com.cairn.common.utils.Coord;
import com.cairn.gape.utils.InfoMessageLogger;

/**
 * A class to implement the Tripos force field
 * 
 * @author Gareth Jones
 *
 */
public class Taff {
	// TODO move all this to taff package and split up classes

	private static final Logger logger = Logger.getLogger(Taff.class);

	static {
		// logger.setLevel(Level.DEBUG);
	}

	private static ThreadLocal<TaffAtoms> tfAtoms = ThreadLocal
			.withInitial(() -> new TaffAtoms());

	private static ThreadLocal<TaffBonds> tfBonds = ThreadLocal
			.withInitial(() -> new TaffBonds());

	private static ThreadLocal<TaffAngles> tfAngles = ThreadLocal
			.withInitial(() -> new TaffAngles());

	private static ThreadLocal<TaffTorsions> tfTorsions = ThreadLocal
			.withInitial(() -> new TaffTorsions());

	private static ThreadLocal<TaffOops> tfOops = ThreadLocal
			.withInitial(() -> new TaffOops());

	private volatile boolean ignore14 = false;
	private volatile boolean ignoreVdwAttractive = false;
	private volatile double tolerance = 1e-5;
	private volatile double eVdw, eAng, eBond, eOop, eTor, energy;
	private volatile Molecule mol;
	private volatile Atom currentAtom;
	private volatile int nCycles = 100;
	protected volatile InfoMessageLogger infoMessageLogger;

	// May need to make these threadlocal.
	private final double trans[][] = new double[4][4];
	private final double xRot[][] = new double[4][4];
	private final double yRot[][] = new double[4][4];
	private final double product1[][] = new double[4][4];
	private final double product2[][] = new double[4][4];
	private final double product3[][] = new double[4][4];
	private final double rot[][] = new double[4][4];
	private final double diff[] = new double[4];
	private final double dummy[] = new double[4];
	private final double dummy1[] = new double[4];

	public Taff(Molecule mol) {
		this.mol = mol;
		infoMessageLogger = mol.getInfoMessageLogger();
		init();
	}

	/**
	 * Bind the forcefield to the molecule
	 */
	void init() {
		tfBonds.get().addParameters(mol);

		tfAtoms.get().addParameters(mol);

		if (mol.getnAngles() == 0)
			mol.assignAngles();
		tfAngles.get().addParameters(mol);

		if (mol.getnTorsions() == 0)
			mol.assignTorsions();
		tfTorsions.get().addParameters(mol);

		tfOops.get().addParameters(mol);
		for (Atom atom : mol.getAtoms()) {
			atom.taffSetup();
		}
	}

	public static void main(String args[]) {
		if (args.length == 0) {
			System.err.println("Usage: Taff <mol2file> [<max_cycles>]");
			System.exit(0);
		}

		Molecule molecule = new Molecule();
		molecule.loadFile(args[0]);
		Taff taff = new Taff(molecule);
		if (args.length == 2)
			taff.nCycles = Integer.valueOf(args[1]).intValue();
		taff.minimizeAtoms();
		String out = "Minimized " + molecule.baseName + ".mol2";
		System.out.println("Molecule is in " + out);
		molecule.writeSybylMol2File(out, "Minimized Molecule: " + molecule.getName());

	}

	void minimizeMolecule() {
		TaffSimplexMolecule simplex = new TaffSimplexMolecule();
		for (int i = 0; i < 1; i++) {
			simplex.initSimplex();
			simplex.optimize();
			simplex.result();
			infoMessageLogger.infoMessageln(3, "Energy is " + energy);
		}
		infoMessageLogger.infoMessageln(3, "VDW energy " + eVdw);
		infoMessageLogger.infoMessageln(3, "Bond Stretch " + eBond);
		infoMessageLogger.infoMessageln(3, "Angle Bend " + eAng);
		infoMessageLogger.infoMessageln(3, "Out of Plane " + eOop);
		infoMessageLogger.infoMessageln(3, "Energy is " + energy);
	}

	protected void redraw() {
		;
	}

	protected void minimizeAtoms() {
		molEnergy();
		int no = nCycles;
		TaffSimplexAtom simplex = new TaffSimplexAtom();
		double last = Double.MIN_VALUE;

		for (int i = 1; i <= no; i++) {
			mol.getAtoms().stream().forEach(a -> {
				currentAtom = a;
				simplex.initSimplex();
				simplex.optimize();
				simplex.result();
			});

			molEnergy();
			infoMessageLogger.infoMessageln(3, "Cycle " + i + " Energy is " + energy);
			if (Math.abs(last - energy) < tolerance && i > no / 2) {
				infoMessageLogger.infoMessageln(3, "Converged");
				break;
			}
			last = energy;
		}
		molEnergy();

		infoMessageLogger.infoMessageln(3, "VDW energy " + eVdw);
		infoMessageLogger.infoMessageln(3, "Bond Stretch " + eBond);
		infoMessageLogger.infoMessageln(3, "Angle Bend " + eAng);
		infoMessageLogger.infoMessageln(3, "Out of Plane " + eOop);
		infoMessageLogger.infoMessageln(3, "Energy is " + energy);
	}

	public double moleculeVdwEnergy() {
		double eVdw = 0;
		int nAtoms = mol.getnAtoms();
		for (int i = 0; i < nAtoms; i++)
			for (int j = i + 1; j < nAtoms; j++) {
				Atom a1 = mol.getAtom(i);
				Atom a2 = mol.getAtom(j);
				if (a1.getTaff() == null || a2.getTaff() == null)
					continue;
				if (mol.bonded(a1, a2))
					continue;
				if (mol.bonded13(a1, a2))
					continue;
				if (ignore14)
					if (mol.bonded14(a1, a2))
						continue;
				double sqrD = Coord.sqrDistance(mol.getCoord(a1.getNo()),
						mol.getCoord(a2.getNo()));
				double vdw = a1.getTaff().getR() + a2.getTaff().getR();
				sqrD = sqrD / (vdw * vdw);
				double pow6 = sqrD * sqrD * sqrD;
				double pow12 = pow6 * pow6;
				double kMean = Math.sqrt(a1.getTaff().getK() * a2.getTaff().getK());
				double atomVdw = kMean * (1.0 / pow12 - 2.0 / pow6);
				if (ignoreVdwAttractive && atomVdw < .0)
					atomVdw = 0;
				eVdw += atomVdw;
			}
		return eVdw;
	}

	double atomVdwEnergy(Atom atom) {
		double eVdw = 0;
		for (Atom a2 : mol.getAtoms()) {
			if (atom == a2)
				continue;
			if (atom.getTaff() == null || a2.getTaff() == null)
				continue;
			if (mol.bonded(atom, a2))
				continue;
			if (mol.bonded13(atom, a2))
				continue;
			if (ignore14)
				if (mol.bonded14(atom, a2))
					continue;
			double sqrD = Coord.sqrDistance(mol.getCoord(atom.getNo()),
					mol.getCoord(a2.getNo()));
			double vdw = atom.getTaff().getR() + a2.getTaff().getR();
			sqrD = sqrD / (vdw * vdw);
			double pow6 = sqrD * sqrD * sqrD;
			double pow12 = pow6 * pow6;
			double kMean = Math.sqrt(atom.getTaff().getK() * a2.getTaff().getK());
			eVdw += kMean * (1.0 / pow12 - 2.0 / pow6);
		}
		return eVdw;
	}

	double bondEnergy(Bond b) {
		if (b.getTaff() == null)
			return .0;
		double d = Coord.distance(mol.getCoord(b.getAtom1().getNo()),
				mol.getCoord(b.getAtom2().getNo()));
		double diff = d - b.getTaff().getLen();
		double e = b.getTaff().getK() * (diff * diff);
		logger.debug("Length " + d + " taff len " + b.getTaff().getLen() + " k "
				+ b.getTaff().getK());
		logger.debug("Bond energy between " + String.valueOf(b.getAtom1().getNo() + 1)
				+ " and " + String.valueOf(b.getAtom2().getNo() + 1) + " : " + e);
		return e;
	}

	private static double v1[] = new double[3];
	private static double v2[] = new double[3];

	double angleEnergy(Angle a) {
		if (a.getTaff() == null)
			return .0;
		double p1[] = mol.getCoord(a.getAtom1().getNo());
		double mid[] = mol.getCoord(a.getMid().getNo());
		double p2[] = mol.getCoord(a.getAtom2().getNo());
		Coord.subtract(p1, mid, v1);
		Coord.subtract(p2, mid, v2);
		double angle = Coord.angle(v1, v2);
		angle = (angle * 180.0) / Math.PI;
		double diff = a.getTaff().getAngle() - angle;
		double e = a.getTaff().getK() * diff * diff;
		logger.debug("Angle " + angle + " taff angle " + a.getTaff().getAngle() + " k "
				+ a.getTaff().getK());
		logger.debug("Angle energy between " + String.valueOf(a.getAtom1().getNo() + 1)
				+ ", " + String.valueOf(a.getMid().getNo() + 1) + " and "
				+ String.valueOf(a.getAtom2().getNo() + 1) + ": " + e);
		return e;
	}

	public double torsionEnergy(Torsion t) {
		if (t.getTaff() == null)
			return .0;
		double a[] = mol.getCoord(t.getAtom1().getNo());
		double b[] = mol.getCoord(t.getAtom2().getNo());
		double c[] = mol.getCoord(t.getAtom3().getNo());
		double d[] = mol.getCoord(t.getAtom4().getNo());
		double angle = .0;
		if (t.isReverse())
			angle = Coord.torsion(d, c, b, a);
		else
			angle = Coord.torsion(a, b, c, d);
		double energy = t.energy(angle);
		return energy;
	}

	double oopEnergy(Atom atom) {
		double d = 0;
		if (atom.getTaffOop() == null)
			return .0;
		if (atom.getnNotDummyNeighbours() != 3)
			return .0;
		double p1[] = mol.getCoord(atom.getNotDummyNeighbours().get(0).getNo());
		double p2[] = mol.getCoord(atom.getNotDummyNeighbours().get(1).getNo());
		double p3[] = mol.getCoord(atom.getNotDummyNeighbours().get(2).getNo());
		double p[] = mol.getCoord(atom.getNo());
		// p1 to origin p2 along z-axis
		double x = p1[0];
		double y = p1[1];
		double z = p1[2];
		Coord.identity(trans);
		trans[3][0] = -x;
		trans[3][1] = -y;
		trans[3][2] = -z;
		Coord.subtract(p1, p2, diff);
		Coord.unit(diff);
		double a = diff[0];
		double b = diff[1];
		double c = diff[2];
		double v = Math.sqrt(b * b + c * c);
		Coord.identity(xRot);
		if (v > .0) {
			xRot[1][1] = xRot[2][2] = c / v;
			xRot[1][2] = b / v;
			xRot[2][1] = -b / v;
		}
		Coord.identity(yRot);
		yRot[0][0] = yRot[2][2] = v;
		yRot[0][2] = a;
		yRot[2][0] = -a;
		Coord.product(trans, xRot, product1);
		Coord.product(product1, yRot, product2);
		Coord.transPoint(product2, p3, dummy);
		// <angle> is the angle atom3 makes with the zy plane
		double angle = FastMath
				.asin(dummy[0] / (Math.sqrt(dummy[0] * dummy[0] + dummy[1] * dummy[1])));
		if (dummy[1] < 0.0)
			angle = Math.PI - angle;
		angle *= -1.0;

		/* atom3 is rotated onto the zy plane */
		rot[2][2] = rot[3][3] = 1;
		rot[0][0] = rot[1][1] = Math.cos(angle);
		rot[0][1] = -1.0 * Math.sin(angle);
		rot[1][0] = Math.sin(angle);
		Coord.product(product2, rot, product3);
		Coord.transPoint(product3, p, dummy1);
		d = dummy1[0];
		logger.trace("D " + d);

		return atom.getTaffOop().getK() * d * d;
	}

	double atomEnergy(Atom atom) {
		double eVdw = atomVdwEnergy(atom);
		double eBond = atom.getBonds().stream().mapToDouble(b -> bondEnergy(b)).sum();
		double eAng = atom.getAngles().stream().mapToDouble(a -> angleEnergy(a)).sum();
		double eOop = oopEnergy(atom);
		double eTor = atom.getTorsions().stream().mapToDouble(t -> torsionEnergy(t))
				.sum();

		double energy = eVdw + eBond + eAng + eOop + eTor;
		logger.debug("Energy is " + energy);
		if (Double.isNaN(energy)) {
			logger.error("NaN energy");
			energy = Double.MAX_VALUE;
		}
		atom.setEnergy(energy);
		return energy;
	}

	public double molEnergy() {
		synchronized (mol) {
			double eVdw = moleculeVdwEnergy();
			double eBond = mol.getBonds().stream().mapToDouble(b -> bondEnergy(b)).sum();
			double eAng = mol.getAngles().stream().mapToDouble(a -> angleEnergy(a)).sum();
			double eOop = mol.getAtoms().stream().mapToDouble(a -> oopEnergy(a)).sum();
			double eTor = mol.getTorsions().stream().mapToDouble(t -> torsionEnergy(t))
					.sum();
			this.eVdw = eVdw;
			this.eBond = eBond;
			this.eAng = eAng;
			this.eOop = eOop;
			this.eTor = eTor;
			this.energy = eVdw + eBond + eAng + eOop + eTor;
			logger.debug("Energy is " + energy);
			if (Double.isNaN(energy)) {
				logger.error("NaN energy");
				energy = Double.MAX_VALUE;
			}
			return energy;
		}
	}

	static boolean matchType(AtomType.Type moleculeAtomType, AtomType.Type ffAtomType) {
		assert (moleculeAtomType != AtomType.Type.WILD);
		assert (ffAtomType != AtomType.Type.OAR);

		if (ffAtomType == AtomType.Type.WILD) {
			return true;
		}
		// O.ar is not in any of the TAFF forcefield files- use O.3 as a
		// surrogate
		else if (moleculeAtomType == AtomType.Type.OAR) {
			return AtomType.Type.O3 == ffAtomType;
		} else {
			return moleculeAtomType == ffAtomType;
		}
	}

	/**
	 * @return the ignore14
	 */
	public boolean isIgnore14() {
		return ignore14;
	}

	/**
	 * @return the ignoreVdwAttractive
	 */
	public boolean isIgnoreVdwAttractive() {
		return ignoreVdwAttractive;
	}

	/**
	 * @return the mol
	 */
	protected Molecule getMol() {
		return mol;
	}

	/**
	 * @return the eVdw
	 */
	protected double geteVdw() {
		return eVdw;
	}

	/**
	 * @return the eAng
	 */
	protected double geteAng() {
		return eAng;
	}

	/**
	 * @return the eBond
	 */
	protected double geteBond() {
		return eBond;
	}

	/**
	 * @return the eOop
	 */
	protected double geteOop() {
		return eOop;
	}

	/**
	 * @return the eTor
	 */
	protected double geteTor() {
		return eTor;
	}

	/**
	 * @return the energy
	 */
	protected double getEnergy() {
		return energy;
	}

	/**
	 * @return the nCycles
	 */
	protected int getnCycles() {
		return nCycles;
	}

	/**
	 * @return the currentAtom
	 */
	protected Atom getCurrentAtom() {
		return currentAtom;
	}

	/**
	 * @param currentAtom
	 *            the currentAtom to set
	 */
	protected void setCurrentAtom(Atom currentAtom) {
		this.currentAtom = currentAtom;
	}

	/**
	 * @return the tolerance
	 */
	protected double getTolerance() {
		return tolerance;
	}

	/**
	 * @param ignore14
	 *            the ignore14 to set
	 */
	public void setIgnore14(boolean ignore14) {
		this.ignore14 = ignore14;
	}

	/**
	 * @param ignoreVdwAttractive
	 *            the ignoreVdwAttractive to set
	 */
	public void setIgnoreVdwAttractive(boolean ignoreVdwAttractive) {
		this.ignoreVdwAttractive = ignoreVdwAttractive;
	}

	class TaffSimplexMolecule extends Simplex {
		private int nEvals = 0;

		@Override
		public double value(double[] test) {
			int i = 0;
			for (double[] coord : mol.getCoords()) {
				coord[0] = test[i++];
				coord[1] = test[i++];
				coord[2] = test[i++];
			}
			molEnergy();
			nEvals++;
			if (nEvals % 100 == 0) {
				logger.info("Iter " + nEvals + " energy " + energy);
				redraw();
			}
			return energy;
		}

		public void initSimplex() {
			setMaxEvals(5000);
			int ndim = mol.getnAtoms() * 3;
			double[] steps = new double[ndim];
			double[] initialGuess = new double[ndim];
			int j = 0;
			for (double[] coord : mol.getCoords()) {
				steps[j] = 0.25;
				initialGuess[j++] = coord[0];
				steps[j] = 0.25;
				initialGuess[j++] = coord[1];
				steps[j] = 0.25;
				initialGuess[j++] = coord[2];
			}

			setStartingPoint(initialGuess);
			setSteps(steps);
		}

		void result() {
			value(getMinimum());
		}
	}

	public class TaffSimplexAtom extends Simplex {
		public TaffSimplexAtom() {
			;
		}

		@Override
		public double value(double[] test) {
			double coord[] = mol.getCoord(currentAtom.getNo());
			coord[0] = test[0];
			coord[1] = test[1];
			coord[2] = test[2];
			double e = atomEnergy(currentAtom);
			return e;
		}

		public void initSimplex() {
			setSteps(new double[] { 0.25, 0.25, 0.25 });
			double coord[] = mol.getCoord(currentAtom.getNo());
			setStartingPoint(new double[] { coord[0], coord[1], coord[2] });
		}

		public void result() {
			value(getMinimum());
		}
	}

}
