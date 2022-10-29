package com.cairn.gape.feature;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;
import com.cairn.gape.molecule.GaMolecule;
import com.cairn.gape.utils.InfoMessageLogger;
import com.cairn.molecule.Atom;
import com.cairn.molecule.Molecule;
import com.cairn.molecule.Ring;

/**
 * Class to represent aromatic or other planar rings as pharmacophore features.
 * 
 * @author Gareth Jones
 * 
 */
public class AromaticRingFeature extends Feature {
	private static Logger logger;
	private static boolean logDebug;

	static {
		logger = Logger.getLogger(AromaticRingFeature.class);
		// logger.setLevel(Level.DEBUG);
		logDebug = logger.isDebugEnabled();
	}

	private volatile Ring ring;

	// list of atoms in the ring ordered by id
	private volatile List<Atom> atoms;

	// two ring normals
	private volatile double normals[][];

	private static ThreadLocal<Double> normalLen = new ThreadLocal<Double>() {
		@Override
		protected Double initialValue() {
			return 3.0;
		}
	};

	/**
	 * Constructor
	 */
	private AromaticRingFeature() {
		featureType = FeatureType.AROMATIC_RING;
		featureSetName = "AROMATIC_RING";
		atomFeature = false;
	}

	/**
	 * Creates a feature for the molecular ring
	 * 
	 * @param m
	 * @param r
	 * @param no
	 */
	public AromaticRingFeature(GaMolecule m, Ring r, int no) {
		this();
		molecule = m;
		ring = r;
		featureSetNo = no;

		atoms = new ArrayList<>(r.getAtoms());

		// sort atoms list by id
		Collections.sort(atoms, (a1, a2) -> Integer.compare(a1.getNo(), a2.getNo()));

	}

	/**
	 * Creates a ring feature from a pharmacophore description string.
	 * 
	 * @param options
	 * @param FeatureDefinion
	 * @param geometry
	 * @param _moleculeNo
	 */
	public AromaticRingFeature(HashMap<String, String> options, String FeatureDefinion,
			PharmFeatureGeometry geometry, int _moleculeNo) {
		this();

		virtual = true;
		moleculeNo = _moleculeNo;

		pharmFeatureGeometry = geometry;
		VectorPharmFeatureGeometry g = (VectorPharmFeatureGeometry) geometry;
		normals = new double[][] { g.points[0], g.points[1] };
		coordinate = new double[4];
		coordinate[3] = 1.0;
		Coord.midPoint(g.points[0], g.points[1], coordinate);
	}

	/**
	 * Searches a molecules for planar rings and creates features.
	 * 
	 * @param m
	 * @return a list of ring features
	 */
	public static List<Feature> findFeatures(GaMolecule m) {
		InfoMessageLogger infoMessageLogger = m.getInfoMessageLogger();
		ArrayList<Feature> v = new ArrayList<Feature>();
		infoMessageLogger.infoMessage(3, "Looking for Aromatic Rings: ");
		int no = 0;
		for (Ring ring : m.getRings()) {
			// if (ring.aromatic) {
			if (ring.isAllSp2()) {
				infoMessageLogger.infoMessage(3, ring.info() + " ");
				AromaticRingFeature feature = new AromaticRingFeature(m, ring, no);
				v.add(feature);

				no++;
			}
		}

		infoMessageLogger.infoMessage(3, "");
		return v;
	}

	/*
	 * Finds fitting point for ring feature. The fitting point is the ring
	 * center. Ring normals are also calculated. (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Feature#getCoordinate()
	 */
	@Override
	public double[] calculateCoordinate() {
		if (virtual)
			return coordinate;

		coordinate = ring.findCenter();
		normals = ring.ringNormals(normalLen.get());

		return coordinate;
	}

	/*
	 * A simple label for the ring feature (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Feature#info()
	 */
	@Override
	public String info() {
		if (virtual)
			return "ring [virtual]";
		return ring.info();
	}

	/*
	 * A descriptive label that can be used in output files and to store
	 * pharmacophre feature. (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Feature#pharmLabel()
	 */
	@Override
	public String pharmLabel() {
		if (virtual)
			return featureSetName + " VIRTUAL";

		List<String> atomIds = atoms.stream().map(a -> String.valueOf(a.getNo() + 1))
				.collect(Collectors.toList());
		String alist = StringUtils.join(atomIds, "_");
		return featureSetName + " ATOMS_" + alist;
	}

	/*
	 * Determines similarity between two ring features. (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Feature#score(com.cairn.gape.Feature)
	 */
	@Override
	public double score(Feature f) {
		AromaticRingFeature other = (AromaticRingFeature) f;

		// Center score
		double vol = Feature.score(Coord.sqrDistance(coordinate, other.coordinate));

		// Normal scores
		double nVol = .0;
		for (int i = 0; i < 2; i++) {
			double v = 0;
			for (int j = 0; j < 2; j++) {
				double t = Feature.score(Coord.sqrDistance(normals[i], other.normals[j]));
				if (t > v)
					v = t;
			}
			if (logDebug)
				logger.debug("+" + v);
			nVol += v;
		}
		nVol = nVol / 2.0;

		// the score is the lesser of center and normal scores
		double score = (vol < nVol) ? vol : nVol;

		if (logDebug)
			logger.debug(" Score: " + score + " vol " + vol + " nVol " + nVol);

		if (score > maximumGaussianScore.get())
			score = maximumGaussianScore.get();

		if (logDebug)
			logger.debug(info() + " " + other.info() + " score " + score);

		geometricScore = score;
		return score;
	}

	/**
	 * Writes out ring normal fitting points to a mol2 file. Used for debugging.
	 * 
	 * @param mol
	 * @param file
	 */
	private static void aromaticRingFeatures2Mol2(GaMolecule mol, String file) {
		List<Feature> features = mol.getFeatureMappings(FeatureType.AROMATIC_RING)
				.getFeatures();

		List<Atom> atoms = new ArrayList<>(features.size() * 3);
		List<double[]> coords = new ArrayList<>(features.size() * 3);
		int n = 0;
		for (int i = 0; i < features.size(); i++) {
			String l = String.valueOf(i + 1);
			AromaticRingFeature f = (AromaticRingFeature) features.get(i);
			atoms.add(new Atom(n, "Ring_" + l, "Du"));
			coords.add(f.calculateCoordinate());
			n++;

			atoms.add(new Atom(mol, n, "Ring_" + l + "_normal_1", "Du"));
			coords.add(f.normals[0]);
			n++;

			atoms.add(new Atom(mol, n, "Ring_" + l + "_normal_2", "Du"));
			coords.add(f.normals[1]);
			n++;
		}
		Molecule vpMol = new Molecule("Aromatic Ring Fitting Points", atoms, coords);

		logger.info(file);
		vpMol.writeSybylMol2File(file, "Aromatic Ring Points");
	}

	/**
	 * Used for testing and debugging
	 * 
	 * @param args
	 */
	public static void main(String args[]) {
		GaMolecule molA = new GaMolecule();
		GaMolecule molB = new GaMolecule();

		molA.loadFile(args[0]);
		molA.findFeatures();
		aromaticRingFeatures2Mol2(molA, "aromaticA.mol2");
		List<Feature> features = molA.getFeatureMappings(FeatureType.AROMATIC_RING)
				.getFeatures();
		features.stream().forEach(f -> f.calculateCoordinate());

		molB.loadFile(args[1]);
		molB.findFeatures();
		aromaticRingFeatures2Mol2(molB, "aromaticB.mol2");
		features = molB.getFeatureMappings(FeatureType.AROMATIC_RING).getFeatures();
		features.stream().forEach(f -> f.calculateCoordinate());

		double score = molA.aromaticRingScore(molB);
		logger.info("score " + score);

	}

	/*
	 * Returns the appropriate pharmacophore feature geometry definition for
	 * this ring feature. (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Feature#getPharmFeatureGeometry()
	 */
	@Override
	public PharmFeatureGeometry getPharmFeatureGeometry() {
		if (virtual)
			return pharmFeatureGeometry;

		pharmFeatureGeometry = new VectorPharmFeatureGeometry(normals[0], normals[1]);
		return pharmFeatureGeometry;
	}

	public static double getNormalLen() {
		return normalLen.get();
	}

	public static void setNormalLen(double normalLen) {
		AromaticRingFeature.normalLen.set(normalLen);
	}

}
