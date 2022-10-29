package com.cairn.gape;

import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;
import com.cairn.gape.chromosome.SuperpositionChromosome;
import com.cairn.gape.feature.Feature;
import com.cairn.gape.molecule.GaMolecule;
import com.cairn.molecule.Atom;
import com.cairn.molecule.Ullman;
import com.cairn.molecule.Ullman.MatchType;

/**
 * This class handles the overlay of two chromomsomes
 * 
 * @author Gareth Jones
 */
public class RigidFit implements Ullman.UllmanCallback {
	private final SuperpositionChromosome base;

	private volatile SuperpositionChromosome other;

	private final GaMolecule fittingMolecule;

	private double trans[][];

	private volatile double[][] refY, refX, newX;

	private final int nPoints;

	private volatile double minRMS;

	private volatile Ullman ull;

	private static final Logger logger = Logger.getLogger(RigidFit.class);

	/**
	 * Constructor. The base chromosome should be the currently fitted
	 * chromosome.
	 * 
	 * @param _base
	 *            Chromosome that fitting should be done to.
	 */
	public RigidFit(SuperpositionChromosome _base) {
		logger.debug("Creating rigid fitting");

		base = _base;
		trans = new double[4][4];
		fittingMolecule = base.getFittingMolecule();
		nPoints = fittingMolecule.getnAtoms();
		newX = new double[nPoints][4];
		refY = new double[nPoints][4];
		for (int i = 0; i < nPoints; i++) {
			Coord.copy(base.getFittingMolecule().getCoord(i), refY[i]);
		}
	}

	/**
	 * Here we fit molecules in another chromosome to the base coordinates. The
	 * other chomosome should be the current fitted chromosome.
	 * 
	 * @param _other
	 */
	public void fit(SuperpositionChromosome _other) {
		logger.debug("Fitting chromosome");

		other = _other;
		refX = new double[nPoints][];
		for (int i = 0; i < nPoints; i++) {
			refX[i] = other.getFittingMolecule().getCoord(i);
		}
		minRMS = Double.MAX_VALUE;

		ull = new Ullman(fittingMolecule, fittingMolecule, MatchType.MATCH);
		ull.doUllman(this);

		logger.debug("Updating coordinates to min rms " + minRMS);

		List<? extends GaMolecule> molecules = other.getMolecules();
		int nMolecules = other.getMolecules().size();

		// Atom Coordinates
		for (int i = 0; i < nMolecules; i++) {
			GaMolecule molecule = other.getMolecules(i);
			for (double[] coord : molecule.getCoords()) {
				Coord.transPointInPlace(trans, coord);
			}
		}

		// rebuild feature coordinates:
		// In feature list

		for (int i = 0; i < nMolecules; i++) {
			GaMolecule molecule = molecules.get(i);
			int nFeatures = molecule.getNFeatures();
			for (int j = 0; j < nFeatures; j++) {
				Feature f = molecule.getAllFeature(j);
				// recalculate feature coordinate
				f.calculateCoordinate();
				// copy to chromosome.
				Coord.copy(f.getSavedCoordinate(), other.getFeatureCoordinates(i, j));

			}
		}

	}

	/**
	 * Handles base molecule match to itself.
	 */
	@Override
	public void callback(Map<Atom, Atom> atomMapping) {
		logger.debug("Found isomorphism");

		double mapX[][] = new double[nPoints][];
		double mapY[][] = new double[nPoints][];
		int mapping[] = new int[nPoints];
		for (int i = 0; i < mapping.length; i++)
			mapping[i] = -1;

		int no = 0;
		for (Entry<Atom, Atom> entry : atomMapping.entrySet()) {
			int no1 = entry.getKey().getNo();
			int no2 = entry.getValue().getNo();
			mapX[no] = refX[no1];
			mapY[no] = refY[no2];
			mapping[no1] = no2;
			no++;
		}

		double trans[][] = new double[4][4];
		Coord.leastSquaresFit(mapX, mapY, trans);

		for (int i = 0; i < nPoints; i++)
			Coord.transPoint(trans, refX[i], newX[i]);

		double d = .0, cnt = 0;
		for (int i = 0; i < nPoints; i++) {
			if (mapping[i] == -1)
				continue;
			double p1[] = newX[i];
			double p2[] = refY[mapping[i]];
			for (int j = 0; j < 3; j++) {
				double diff = p1[j] - p2[j];
				d += diff * diff;
			}
			cnt++;
		}
		double rms = Math.sqrt(d / cnt);
		logger.debug("RMS " + rms);
		if (rms < minRMS) {
			minRMS = rms;
			this.trans = trans;
		}

	}

}