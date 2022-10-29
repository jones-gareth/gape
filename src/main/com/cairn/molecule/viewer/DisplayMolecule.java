package com.cairn.molecule.viewer;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Optional;

import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;
import com.cairn.molecule.Atom;
import com.cairn.molecule.AtomType;
import com.cairn.molecule.Bond;
import com.cairn.molecule.BondType;
import com.cairn.molecule.Molecule;
import com.cairn.molecule.Ring;
import com.cairn.molecule.viewer.MolView.LabelType;

public class DisplayMolecule extends Molecule {

	private static Logger logger = Logger.getLogger(DisplayMolecule.class);
	static {
		// logger.setLevel(Level.DEBUG);
	}

	private Bounds bounds;
	private double coords2Display[][];
	private List<double[]> display;
	private boolean colorByType = true;
	private Color color = Color.white;
	private static final double BOND_WIDTH = .15;
	private static boolean displayLonePairs = false, prepare = false;
	static final int MAX_DISPLAY = 200;

	public DisplayMolecule() {
		infoMessageLogger.setLogLevel(4);
	}

	public DisplayMolecule(Molecule m) {
		// we don't copy!!
		this();
		createReference(m);
		createDisplayCoordinates();
		resetDisplay();
	}

	public DisplayMolecule(String file) {
		super(file, FileType.MOL2, Source.FILE);
	}

	public DisplayMolecule(String file, FileType type, Source source) {
		super(file, type, source);
	}

	public DisplayMolecule(BufferedReader in, FileType type) {
		super(in, type);
	}

	@Override
	public void init() {
		if (prepare)
			super.init();
		else
			super.build();
		createDisplayCoordinates();
		resetDisplay();
	}

	public void createDisplayCoordinates() {
		display = new ArrayList<>();
		for (double[] coord : getCoords()) {
			display.add(Arrays.copyOf(coord, 4));
		}
	}

	public static List<DisplayMolecule> loadFiles(String files[], FileType startType,
			Source source) {

		ArrayList<DisplayMolecule> molecules = new ArrayList<DisplayMolecule>();
		for (int i = 0; i < files.length; i++) {
			if (files[i].startsWith("-"))
				continue;
			logger.debug("loading " + files[i]);
			BufferedReader in = openReader(files[i], source);
			if (in == null)
				continue;
			FileType type = startType;
			if (type == FileType.UNKNOWN) {
				if (files[i].toUpperCase().endsWith(".MOL2")) {
					type = FileType.MOL2;
				} else {
					type = FileType.SDF;
				}
			}
			boolean loading = true;
			while (loading) {
				DisplayMolecule molecule = new DisplayMolecule();
				molecule.infoMessageLogger.setLogLevel(5);

				try {
					if (type == FileType.MOL2) {
						molecule.loadSybylMol2(in);
					} else if (type == FileType.SDF) {
						molecule.loadSdfMol(in);
					}
				} catch (IOException ex) {
					logger.error("IO error reading " + ex);
					break;
				} catch (MolReadException ex) {
					logger.debug("Format error or end of stream reading " + ex);
					break;
				}
				molecule.infoMessageLogger.infoMessageln("display molecule: loaded mol "
						+ molecule.getName());
				molecule.init();

				molecules.add(molecule);
				if (molecules.size() > MAX_DISPLAY) {
					loading = false;
					logger.warn("More than " + MAX_DISPLAY
							+ " molecules present - only displaying " + MAX_DISPLAY
							+ " molecules");
				}
			}
			try {
				in.close();
			} catch (IOException ex) {
				logger.error("IO error closing " + ex);
			}
		}

		return molecules;
	}

	public void resetDisplay() {
		coords2Display = new double[4][4];
		Coord.identity(coords2Display);
	}

	synchronized public void shake() {
		for (double[] v : getCoords()) {
			double r1 = Math.random() - .5;
			double r2 = Math.random() - .5;
			double r3 = Math.random() - .5;
			v[0] += r1;
			v[1] += r2;
			v[2] += r3;
		}
		recalculateRingCenters();
	}

	public void recalculateRingCenters() {
		for (Ring ring : getRings()) {
			ring.findCenter();
		}
	}

	void createMoleculeBounds() {
		bounds = new Bounds();
	}

	public class Bounds {
		private double xMax, xMin, yMax, yMin, zMax, zMin;

		Bounds() {
			logger.debug("Creating bounds");
			if (getnAtoms() == 0) {
				xMin = yMin = zMin = -1;
				xMax = yMax = zMax = 1;
				return;
			}
			xMax = xMin = getCoord(0)[0];
			yMax = yMin = getCoord(0)[1];
			zMax = zMin = getCoord(0)[2];
			for (double[] coord : getCoords()) {
				double x = coord[0];
				double y = coord[1];
				double z = coord[2];
				logger.debug("x " + x + " y " + y + " z " + z);
				if (x > xMax)
					xMax = x;
				if (y > yMax)
					yMax = y;
				if (z > zMax)
					zMax = z;
				if (x < xMin)
					xMin = x;
				if (y < yMin)
					yMin = y;
				if (z < zMin)
					zMin = z;
			}
			logger.debug("Bounds: " + xMin + " " + yMin + " " + zMin + " " + xMax + " "
					+ yMax + " " + zMax);
		}

		/**
		 * @return the xMax
		 */
		public double getxMax() {
			return xMax;
		}

		/**
		 * @return the xMin
		 */
		public double getxMin() {
			return xMin;
		}

		/**
		 * @return the yMax
		 */
		public double getyMax() {
			return yMax;
		}

		/**
		 * @return the yMin
		 */
		public double getyMin() {
			return yMin;
		}

		/**
		 * @return the zMax
		 */
		public double getzMax() {
			return zMax;
		}

		/**
		 * @return the zMin
		 */
		public double getzMin() {
			return zMin;
		}

	}

	void translateDisplay(double trans[][]) {
		if (trans == null)
			return;

		double[][] product = new double[4][4];
		Coord.product(coords2Display, trans, product);
		double temp[][] = coords2Display;
		coords2Display = product;
		product = temp;

		logger.debug("translateDisplay: translating Atoms");
		translateAtoms();
	}

	void translateAtoms() {
		logger.debug("translating Atoms");

		for (int i = 0; i < getnAtoms(); i++) {
			Coord.transPoint(coords2Display, getCoord(i), display.get(i));
		}
		for (Ring ring : getRings()) {
			Coord.transPoint(coords2Display, ring.getCenter(), ring.getDisplay());
		}
	}

	void setDisplay(double[][] trans) {
		if (trans == null)
			return;
		Coord.copy(trans, coords2Display);
		logger.debug("setDisplay: translating Atoms");
		translateAtoms();
	}

	private class BondCompare implements Comparator<Bond> {
		double mid1[] = new double[4];
		double mid2[] = new double[4];

		@Override
		public int compare(Bond b1, Bond b2) {
			Coord.midPoint(display.get(b1.getAtom1().getNo()),
					display.get(b1.getAtom2().getNo()), mid1);
			Coord.midPoint(display.get(b2.getAtom1().getNo()),
					display.get(b2.getAtom2().getNo()), mid2);
			if (mid1[2] > mid2[2])
				return 1;
			if (mid1[2] < mid2[2])
				return -1;
			return 0;
		}
	}

	private class AtomCompare implements Comparator<Atom> {
		@Override
		public int compare(Atom a1, Atom a2) {
			if (display.get(a1.getNo())[2] > display.get(a2.getNo())[2])
				return 1;
			if (display.get(a1.getNo())[2] < display.get(a2.getNo())[2])
				return -1;
			return 0;
		}
	}

	synchronized void display(MolView vp, double trans[][]) {
		HashMap<Atom, Boolean> displayedAtoms = new HashMap<Atom, Boolean>();
		translateDisplay(trans);
		ArrayList<Bond> bondArrayList = new ArrayList<>(getBonds());
		Collections.sort(bondArrayList, new BondCompare());

		double m[] = new double[4];
		for (Bond bond : bondArrayList) {
			Atom a1 = bond.getAtom1();
			Atom a2 = bond.getAtom2();
			if (!displayLonePairs && a1.getAtomType() == AtomType.Type.LP)
				continue;
			if (!displayLonePairs && a2.getAtomType() == AtomType.Type.LP)
				continue;
			if (bond.getStereo() == Bond.Stereo.NONE) {
				if (vp.isHeavy() && a1.getAtomType() != AtomType.Type.DU && !a1.isHeavy())
					continue;
				if (vp.isHeavy() && a1.getAtomType() != AtomType.Type.DU && !a2.isHeavy())
					continue;
			}
			if (vp.getBondLabelType() != LabelType.NONE) {
				String l = null;
				if (vp.getBondLabelType() == LabelType.TYPE)
					l = bond.getType().getName();
				else if (vp.getBondLabelType() == LabelType.ID)
					l = String.valueOf(bond.getNo() + 1);
				vp.setColor(Color.white);
				Coord.midPoint(display.get(a1.getNo()), display.get(a2.getNo()), m);
				vp.drawString(l, m);
			}
			if (vp.isBallAndStick())
				displayThickBond(vp, a1, a2, bond);
			else if (bond.isInRing() && bond.getBondType() == BondType.Type.AR)
				insideRingBond(vp, a1, a2, bond, true);
			else if (bond.isInRing() && bond.getBondType() == BondType.Type.DOUBLE)
				insideRingBond(vp, a1, a2, bond, false);
			else if (!bond.isInRing() && bond.getBondType() == BondType.Type.DOUBLE)
				acyclicDoubleBond(vp, a1, a2, bond);
			else if (!bond.isInRing() && bond.getBondType() == BondType.Type.TRIPLE)
				acyclicTripleBond(vp, a1, a2, bond);
			else if (bond.getStereo() == Bond.Stereo.UP)
				stereoBondUp(vp, a1, a2, bond);
			else if (bond.getStereo() == Bond.Stereo.DOWN)
				stereoBondDown(vp, a1, a2, bond);
			else
				displaySingleBond(vp, a1, a2, bond);
			if (display.get(a1.getNo())[2] < display.get(a2.getNo())[2]) {
				if (!displayedAtoms.containsKey(a1))
					atomLabel(vp, a1, bond);
				if (!displayedAtoms.containsKey(a2))
					atomLabel(vp, a2, bond);
			} else {
				if (!displayedAtoms.containsKey(a2))
					atomLabel(vp, a2, bond);
				if (!displayedAtoms.containsKey(a1))
					atomLabel(vp, a1, bond);
			}
			displayedAtoms.put(a1, true);
			displayedAtoms.put(a2, true);
		}
		List<Atom> atomArrayList = new ArrayList<>();
		for (Atom atom : getAtoms()) {
			if (vp.isHeavy()) {
				// display dummy atoms even if heavy is set -should really check
				// that the neighbours are not light atoms!
				if (atom.getAtomType() == AtomType.Type.DU && atom.getnNeighbours() > 0)
					continue;
				if (atom.getNHeavyNeighbours() > 0)
					continue;
			} else {
				if (atom.getnNeighbours() > 0)
					continue;
			}
			if (!displayLonePairs && atom.getAtomType() == AtomType.Type.LP)
				continue;
			atomArrayList.add(atom);
		}

		Collections.sort(atomArrayList, new AtomCompare());
		for (Atom a : atomArrayList) {
			if (!atomLabel(vp, a, null)) {
				setColor(vp, a);
				vp.fillOval(display.get(a.getNo()), 0.5, 0.5);
			}
		}
	}

	void setColor(MolView vp, Atom a) {
		setColor(vp, a, display.get(a.getNo()));
	}

	void setColor(MolView vp, Atom a, double v[]) {
		if (!colorByType && a.getType().isCarbonType())
			vp.setColor(color, v);
		else
			vp.setColor(a.getType().getColor(), v);
	}

	boolean atomLabel(MolView vp, Atom a, Bond bond) {
		// if (vp.ballAndStick) {
		// setColor(vp, a);
		// double r = a.type.radius/2;
		// vp.fillOval(display[a.no], r, r);
		// }
		if (vp.getAtomLabelType() != LabelType.NONE) {
			String l = null;
			if (vp.getAtomLabelType() == LabelType.TYPE) {
				l = a.getType().getName();
			} else if (vp.getAtomLabelType() == LabelType.LABEL) {
				l = a.getLabel();
			} else if (vp.getAtomLabelType() == LabelType.ID) {
				l = String.valueOf(a.getNo() + 1);
			} else if (vp.getAtomLabelType() == LabelType.HETERO) {
				l = a.getType().getName();
				int s = l.indexOf('.');
				if (s > 0)
					l = l.substring(0, s);
				if (l.equals("C"))
					l = "";
				if (l.equals("H")) {
					if (bond != null && bond.getStereo() == Bond.Stereo.NONE)
						l = "";
				}
			} else if (vp.getAtomLabelType() == LabelType.PARTIAL_CHARGE) {
				l = "";
				if (a.getPartialCharge() != .0) {
					l = String.valueOf(a.getPartialCharge());
					if (a.getPartialCharge() > 0)
						l = "+" + l;
				}
			} else if (vp.getAtomLabelType() == LabelType.FORMAL_CHARGE) {
				l = "";
				if (a.getFormalCharge() != null) {
					l = String.valueOf(a.getFormalCharge());
					if (a.getFormalCharge() > 0)
						l = "+" + l;
				}
			} else if (vp.getAtomLabelType() == LabelType.SUBSTR) {
				l = a.getSubName();
			}
			if (l == null || l.equals("")) {
				return false;
			}
			// Color c = colorByType ? a.type.color : color;
			Color c = color;
			vp.drawCenteredString(l, display.get(a.getNo()), c);
			return true;
		}
		return false;
	}

	void stereoBondUp(MolView vp, Atom a1, Atom a2, Bond b) {
		logger.debug("stereo bond up");

		double bvec[] = new double[4];
		double zvec[] = new double[4];
		double pvec[] = new double[4];
		double p1[] = new double[4];
		double p2[] = new double[4];
		double mid[] = new double[4];
		double mid1[] = new double[4];
		double mid2[] = new double[4];
		double points1[][] = new double[3][];
		double points2[][] = new double[4][];

		Coord.subtract(display.get(a2.getNo()), display.get(a1.getNo()), bvec);
		zvec[0] = 0;
		zvec[1] = 0;
		zvec[2] = 1;
		zvec[3] = 1;
		Coord.vectorProduct(bvec, zvec, pvec);
		Coord.setLength(pvec, BOND_WIDTH);
		Coord.add(display.get(a2.getNo()), pvec, p1);
		Coord.subtract(display.get(a2.getNo()), pvec, p2);
		Coord.midPoint(display.get(a1.getNo()), display.get(a2.getNo()), mid);
		Coord.setLength(pvec, BOND_WIDTH / 2.0);
		Coord.add(mid, pvec, mid1);
		Coord.subtract(mid, pvec, mid2);

		points1[0] = display.get(a1.getNo());
		points1[1] = mid1;
		points1[2] = mid2;
		points2[0] = mid2;
		points2[1] = mid1;
		points2[2] = p1;
		points2[3] = p2;
		setColor(vp, a1);
		vp.drawPolygon(points1, true);
		setColor(vp, a2);
		vp.drawPolygon(points2, true);
	}

	void stereoBondDown(MolView vp, Atom a1, Atom a2, Bond b) {
		logger.debug("stereo bond down");

		double bvec[] = new double[4];
		double zvec[] = new double[4];
		double pvec[] = new double[4];
		double mid[] = new double[4];
		double mid1[] = new double[4];
		double mid2[] = new double[4];

		Coord.subtract(display.get(a2.getNo()), display.get(a1.getNo()), bvec);
		zvec[0] = 0;
		zvec[1] = 0;
		zvec[2] = 1;
		zvec[3] = 1;
		Coord.vectorProduct(bvec, zvec, pvec);
		int nTicks = 6;
		double step = Coord.mag(bvec) / nTicks;
		double h = BOND_WIDTH / nTicks;
		for (int i = 1; i <= nTicks; i++) {
			Coord.setLength(bvec, step * i);
			Coord.add(display.get(a1.getNo()), bvec, mid);
			Coord.setLength(pvec, h * i);
			Coord.add(mid, pvec, mid1);
			Coord.subtract(mid, pvec, mid2);
			if (i <= nTicks / 2)
				setColor(vp, a1);
			else
				setColor(vp, a2);
			vp.drawLine(mid1, mid2);
		}
	}

	void acyclicDoubleBond(MolView vp, Atom a1, Atom a2, Bond b) {
		double bvec[] = new double[4];
		double p1[] = new double[4];
		double p2[] = new double[4];
		double p4[] = new double[4];
		double p3[] = new double[4];
		double t2[] = new double[4];
		double bt[] = new double[4];
		double v1[] = new double[4];
		double v2[] = new double[4];
		double v3[] = new double[4];
		double v4[] = new double[4];
		double up[] = new double[4];
		double down[] = new double[4];

		logger.debug("Entering acyclicDoubleBond");

		boolean p1_bond = false, p2_bond = false, p3_bond = false, p4_bond = false;

		Coord.subtract(display.get(a2.getNo()), display.get(a1.getNo()), bvec);
		for (Atom test : a1.getNotDummyNeighbours()) {
			if (test == a2)
				continue;
			Coord.subtract(display.get(test.getNo()), display.get(a1.getNo()), t2);
			double angle = Coord.angle(t2, bvec);
			if (logger.isDebugEnabled()) {
				double a = angle * 180.0 / Math.PI;
				logger.debug("Angle 1 " + a);
			}
			if (angle > .75 * Math.PI)
				continue;
			Coord.subtract(display.get(test.getNo()), display.get(a1.getNo()), bt);
			Coord.setLength(bt, Coord.parallelDistance(bvec, bt, BOND_WIDTH / 2));
			if (p1_bond) {
				p2_bond = true;
				Coord.copy(bt, p2);
			} else {
				p1_bond = true;
				Coord.copy(bt, p1);
			}
		}

		Coord.subtract(display.get(a1.getNo()), display.get(a2.getNo()), bvec);
		for (Atom test : a2.getNotDummyNeighbours()) {
			if (test == a1)
				continue;
			Coord.subtract(display.get(test.getNo()), display.get(a2.getNo()), t2);
			double angle = Coord.angle(t2, bvec);
			if (logger.isDebugEnabled()) {
				double a = angle * Math.PI / 180.0;
				logger.debug("Angle 2 " + a);
			}
			if (angle > .75 * Math.PI)
				continue;
			Coord.subtract(display.get(test.getNo()), display.get(a2.getNo()), bt);
			Coord.setLength(bt, Coord.parallelDistance(bvec, bt, BOND_WIDTH / 2));
			if (p3_bond) {
				p4_bond = true;
				Coord.copy(bt, p4);
			} else {
				p3_bond = true;
				Coord.copy(bt, p3);
			}
		}

		if (p1_bond) {
			if (p3_bond) {
				double angle = Coord.angle(p1, p3);
				if (angle > .8 * Math.PI) {
					double tmp[] = p3;
					p3 = p4;
					p4 = tmp;
					p3_bond = p4_bond;
					p4_bond = true;
				}
			}
			if (p4_bond) {
				double angle = Coord.angle(p1, p4);
				if (angle < .8 * Math.PI) {
					double tmp[] = p4;
					p4 = p3;
					p3 = tmp;
					p4_bond = p3_bond;
					p3_bond = true;
				}
			}
		}
		if (p2_bond) {
			if (p3_bond) {
				double angle = Coord.angle(p2, p3);
				if (angle < .8 * Math.PI) {
					double tmp[] = p3;
					p3 = p4;
					p4 = tmp;
					p3_bond = p4_bond;
					p4_bond = true;
				}
			}
			if (p4_bond) {
				double angle = Coord.angle(p2, p4);
				if (angle > .8 * Math.PI) {
					double tmp[] = p4;
					p4 = p3;
					p3 = tmp;
					p4_bond = p3_bond;
					p3_bond = true;
				}
			}
		}

		boolean v1_bond = false, v2_bond = false, v3_bond = false, v4_bond = false;

		if (p1_bond) {
			Coord.add(display.get(a1.getNo()), p1, v1);
			v1_bond = true;
		}
		if (p2_bond) {
			Coord.add(display.get(a1.getNo()), p2, v2);
			v2_bond = true;
		}
		if (p3_bond) {
			Coord.add(display.get(a2.getNo()), p3, v3);
			v3_bond = true;
		}
		if (p4_bond) {
			Coord.add(display.get(a2.getNo()), p4, v4);
			v4_bond = true;
		}

		if (!v2_bond && v1_bond) {
			genNewPoint(a1, a2, v1, v2);
			v2_bond = true;
		}
		if (!v1_bond && v2_bond) {
			genNewPoint(a1, a2, v2, v1);
			v1_bond = true;
		}
		if (!v4_bond && v3_bond) {
			genNewPoint(a2, a1, v3, v4);
			v4_bond = true;
		}
		if (!v3_bond && v4_bond) {
			genNewPoint(a2, a1, v4, v3);
			v3_bond = true;
		}

		if (!v1_bond && v3_bond) {
			genNewPoint2(a1, a2, v3, v1);
			v1_bond = true;
		}
		if (!v2_bond && v4_bond) {
			genNewPoint2(a1, a2, v4, v2);
			v2_bond = true;
		}
		if (!v3_bond && v1_bond) {
			genNewPoint2(a2, a1, v1, v3);
			v3_bond = true;
		}
		if (!v4_bond && v2_bond) {
			genNewPoint2(a2, a1, v2, v4);
			v4_bond = true;
		}

		if (!v1_bond) {
			Coord.subtract(display.get(a1.getNo()), display.get(a2.getNo()), bvec);
			up[0] = bvec[1];
			up[1] = -bvec[0];
			up[2] = 0;
			down[0] = -bvec[1];
			down[1] = bvec[0];
			down[2] = 0;
			Coord.setLength(up, BOND_WIDTH / 2);
			Coord.setLength(down, BOND_WIDTH / 2);
			Coord.add(display.get(a1.getNo()), up, v1);
			Coord.add(display.get(a1.getNo()), down, v2);
			Coord.add(display.get(a2.getNo()), up, v3);
			Coord.add(display.get(a2.getNo()), down, v4);
		}

		if (logger.isDebugEnabled()) {
			logger.debug("v1 " + Coord.info(v1));
			logger.debug("v3 " + Coord.info(v3));
			logger.debug("v2 " + Coord.info(v2));
			logger.debug("v4 " + Coord.info(v4));
		}

		displayBondLine(vp, a1, a2, v1, v3);
		displayBondLine(vp, a1, a2, v2, v4);
	}

	void acyclicTripleBond(MolView vp, Atom a1, Atom a2, Bond b) {
		double bvec[] = new double[4];
		double v1[] = new double[4];
		double v2[] = new double[4];
		double v3[] = new double[4];
		double v4[] = new double[4];
		double up[] = new double[4];
		double down[] = new double[4];

		Coord.subtract(display.get(a1.getNo()), display.get(a2.getNo()), bvec);
		up[0] = bvec[1];
		up[1] = -bvec[0];
		up[2] = 0;
		down[0] = -bvec[1];
		down[1] = bvec[0];
		down[3] = 0;
		Coord.setLength(up, BOND_WIDTH / 2);
		Coord.setLength(down, BOND_WIDTH / 2);
		Coord.add(display.get(a1.getNo()), up, v1);
		Coord.add(display.get(a1.getNo()), down, v2);
		Coord.add(display.get(a2.getNo()), up, v3);
		Coord.add(display.get(a2.getNo()), down, v4);
		displayBondLine(vp, a1, a2, v1, v3);
		displayBondLine(vp, a1, a2, display.get(a1.getNo()), display.get(a2.getNo()));
		displayBondLine(vp, a1, a2, v2, v4);
	}

	void genNewPoint(Atom a1, Atom a2, double v[], double np[]) {
		double v2[] = new double[4];
		double v3[] = new double[4];
		double bvec[] = new double[4];

		Coord.subtract(display.get(a1.getNo()), display.get(a2.getNo()), bvec);
		v[3] = 1.0;
		Coord.rotatePoint(display.get(a2.getNo()), display.get(a1.getNo()), v, Math.PI,
				v2);

		Coord.subtract(v2, display.get(a1.getNo()), v3);
		double dist = Coord.perpendicularDistance(bvec, v3, BOND_WIDTH / 2);
		logger.debug("genNewPoint: dist " + dist);
		Coord.setLength(bvec, dist);
		Coord.subtract(v2, bvec, np);
	}

	void genNewPoint2(Atom a1, Atom a2, double v[], double np[]) {
		double v1[] = new double[4];
		double bvec[] = new double[4];

		Coord.subtract(display.get(a1.getNo()), display.get(a2.getNo()), bvec);
		Coord.subtract(v, display.get(a1.getNo()), v1);
		double dist = Coord.mag(bvec)
				+ Coord.perpendicularDistance(bvec, v1, BOND_WIDTH / 2);
		logger.debug("genNewPoint2: dist " + dist);
		Coord.setLength(bvec, dist);
		Coord.add(v, bvec, np);
	}

	void insideRingBond(MolView vp, Atom a1, Atom a2, Bond b, boolean dashed) {
		double v1[] = new double[4];
		double v2[] = new double[4];
		double v3[] = new double[4];
		double v4[] = new double[4];
		double mid[] = new double[4];
		double bond1[] = new double[4];
		double bond2[] = new double[4];

		logger.debug("Entering insideRingBond " + a1.info() + " " + a2.info());
		Coord.midPoint(display.get(a1.getNo()), display.get(a2.getNo()), mid);
		displaySingleBond(vp, a1, a2, b);
		Optional<Ring> optional = b.getRings().stream().filter(r -> r.isAromatic())
				.findFirst();
		Ring ring = optional.isPresent() ? optional.get() : b.getRings().get(0);
		double center[] = ring.getDisplay();
		if (logger.isDebugEnabled()) {
			logger.debug("Center 1 " + Coord.info(center));
			center = ring.findCenter(this.display);
			logger.debug("Center 2 " + Coord.info(center));
		}
		Coord.subtract(center, display.get(a1.getNo()), v1);
		Coord.subtract(center, display.get(a2.getNo()), v2);
		Coord.subtract(display.get(a1.getNo()), display.get(a2.getNo()), bond1);
		Coord.subtract(display.get(a2.getNo()), display.get(a1.getNo()), bond2);
		double sin1 = Coord.sin(v1, bond2);
		double sin2 = Coord.sin(v2, bond1);
		logger.debug("sin1 " + sin1);
		logger.debug("sin2 " + sin2);
		if (sin1 < .0)
			sin1 = -sin1;
		if (sin2 < .0)
			sin2 = -sin2;
		if (sin1 < .01)
			return;
		if (sin2 < .01)
			return;
		double d1 = BOND_WIDTH / sin1;
		double d2 = BOND_WIDTH / sin2;
		logger.debug("d1 " + d2);
		logger.debug("d2 " + d1);
		Coord.setLength(v1, d1);
		Coord.setLength(v2, d2);
		Coord.add(display.get(a1.getNo()), v1, v3);
		Coord.add(display.get(a2.getNo()), v2, v4);
		if (logger.isDebugEnabled()) {
			logger.debug("a1 " + Coord.info(display.get(a1.getNo())));
			logger.debug("a2 " + Coord.info(display.get(a2.getNo())));
			logger.debug("v3 " + Coord.info(v3));
			logger.debug("v4 " + Coord.info(v4));
		}
		if (dashed)
			displayDashedBondLine(vp, a1, a2, v3, v4);
		else
			displayBondLine(vp, a1, a2, v3, v4);
	}

	void displayBondLine(MolView vp, Atom a1, Atom a2, double v1[], double v2[]) {
		double mid[] = new double[4];
		Coord.midPoint(v1, v2, mid);
		setColor(vp, a1, v1);
		vp.drawLine(v1, mid);
		setColor(vp, a2, v2);
		vp.drawLine(mid, v2);
	}

	void displayDashedBondLine(MolView vp, Atom a1, Atom a2, double v1[], double v2[]) {
		double mid[] = new double[4];
		Coord.midPoint(v1, v2, mid);
		setColor(vp, a1, v1);
		vp.drawDashedLine(v1, mid, 0.1);
		setColor(vp, a2, v2);
		vp.drawDashedLine(mid, v2, 0.1);
	}

	void displaySingleBond(MolView vp, Atom a1, Atom a2, Bond b) {
		double mid[] = new double[4];
		Coord.midPoint(display.get(a1.getNo()), display.get(a2.getNo()), mid);
		setColor(vp, a1);
		vp.drawLine(display.get(a1.getNo()), mid);
		setColor(vp, a2);
		vp.drawLine(mid, display.get(a2.getNo()));
	}

	void displayThickBond(MolView vp, Atom a1, Atom a2, Bond b) {
		double mid[] = new double[4];
		Coord.midPoint(display.get(a1.getNo()), display.get(a2.getNo()), mid);
		setColor(vp, a1);
		vp.drawThickLine(display.get(a1.getNo()), mid, 0.15);
		setColor(vp, a2);
		vp.drawThickLine(mid, display.get(a2.getNo()), 0.15);
	}

	/**
	 * @return the bounds
	 */
	public Bounds getBounds() {
		return bounds;
	}

	/**
	 * @return the display
	 */
	public List<double[]> getDisplay() {
		return display;
	}

	/**
	 * @return the colorByType
	 */
	public boolean isColorByType() {
		return colorByType;
	}

	/**
	 * @return the color
	 */
	public Color getColor() {
		return color;
	}

	/**
	 * @param colorByType
	 *            the colorByType to set
	 */
	public void setColorByType(boolean colorByType) {
		this.colorByType = colorByType;
	}

	/**
	 * @param color
	 *            the color to set
	 */
	public void setColor(Color color) {
		this.color = color;
	}

	/**
	 * @return the displayLonePairs
	 */
	public static boolean isDisplayLonePairs() {
		return displayLonePairs;
	}

	/**
	 * @return the prepare
	 */
	public static boolean isPrepare() {
		return prepare;
	}

	/**
	 * @param displayLonePairs
	 *            the displayLonePairs to set
	 */
	public static void setDisplayLonePairs(boolean displayLonePairs) {
		DisplayMolecule.displayLonePairs = displayLonePairs;
	}

	/**
	 * @param prepare
	 *            the prepare to set
	 */
	public static void setPrepare(boolean prepare) {
		DisplayMolecule.prepare = prepare;
	}

}
