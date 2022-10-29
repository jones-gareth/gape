package com.cairn.gape.chromosome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.util.FastMath;
import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;
import com.cairn.gape.Superposition;
import com.cairn.gape.feature.Feature;
import com.cairn.gape.feature.Feature.FeatureType;
import com.cairn.gape.feature.FeatureOverlay;
import com.cairn.gape.feature.Pharmacophore;
import com.cairn.gape.feature.UserFeatureSet;
import com.cairn.gape.ga.BaseSupervisor;
import com.cairn.gape.ga.GaSupervisor;
import com.cairn.gape.molecule.GaMolecule;
import com.cairn.molecule.Atom;
import com.cairn.molecule.Molecule;

/**
 * This class represents the GA chromosome encoding for the GAPE algorithm
 * 
 * @author Gareth Jones
 * @see Chromosome
 * 
 */
public class SuperpositionChromosome extends BinaryAndIntegerChromosome {
	private static final Logger logger = Logger.getLogger(SuperpositionChromosome.class);
	static {
		// logger.setLevel(Level.DEBUG);
	}
	private volatile boolean hasFitness = false, fitted = false, ok;

	private volatile double fitness, volumeIntegral, conformationalEnergy,
			donorHydrogenScore, acceptorAtomScore, aromaticRingScore, constraintScore,
			userFeatureScores[], conformationalEnergies[], volumeIntegrals[][],
			donorHydrogenScores[][], acceptorAtomScores[][], aromaticRingScores[][];

	protected GaMolecule baseMolecule, fittingMolecule;

	protected volatile List<? extends GaMolecule> molecules;

	private volatile int baseMoleculeNo;

	private static final ThreadLocal<SuperpositionChromosome> currentChromosome = new ThreadLocal<SuperpositionChromosome>();

	private volatile double featureCoordinates[][][];

	private static final ThreadLocal<NumberFormat> nf = new ThreadLocal<NumberFormat>();

	private String _info;

	private Feature matchedFeatures[];

	/**
	 * Creates sufficient chromosomes for a GAPE run
	 * 
	 * @param problem
	 * @throws GaException
	 */
	@Override
	public void allocate(GaSupervisor gaSupervisor) {
		problem = gaSupervisor;
		Superposition superposition = (Superposition) problem;

		int no = superposition.getNChromosomes() + 10 + superposition.getNRuns();
		synchronized (chromosomeList) {
			chromosomeList.clear();
			for (int i = 0; i < no; i++) {
				logger.debug("Created chrom  " + i + " of " + no);
				SuperpositionChromosome c = createChromosome();
				c.createEmpty(superposition);
				chromosomeList.add(c);
				c.free = true;
			}
		}

		logger.debug("allocated " + no + " chromosomes");

		nf.set(problem.getNumberFormat());
	}

	/**
	 * @return a new superposition chromomosome
	 */
	public SuperpositionChromosome createChromosome() {
		return new SuperpositionChromosome();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.chromosome.BinaryAndIntegerChromosome#freeChromosome
	 * ()
	 */
	@Override
	public void freeChromosome() {
		super.freeChromosome();
		hasFitness = fitted = false;
	}

	/**
	 * Allocates space for chromosome data structures. Including binary and
	 * integer strings, feature co-ordinates and scores.
	 * 
	 * @param s
	 * @throws GaException
	 */
	public void createEmpty(Superposition s) {
		problem = s;
		Superposition superposition = (Superposition) problem;

		if (superposition.getBinaryStringLength() > 0)
			binaryString = new BinaryStringChromosome(superposition,
					superposition.getBinaryStringLength());

		if (superposition.getIntegerStringLength() > 0) {
			integerString = new IntegerStringChromosome(superposition,
					superposition.getIntegerStringLength(),
					superposition.getIntegerStringRanges());
			integerString.setAllowNullDefault(true);
		}

		molecules = superposition.getMolecules();
		baseMolecule = superposition.getBaseMolecule();
		baseMoleculeNo = superposition.getBaseMoleculeNo();
		fittingMolecule = superposition.getFittingMolecule();

		userFeatureScores = new double[superposition.getnUserFeatureTypes()];

		int nMolecules = molecules.size();
		conformationalEnergies = new double[nMolecules];
		volumeIntegrals = new double[nMolecules][nMolecules];
		donorHydrogenScores = new double[nMolecules][nMolecules];
		acceptorAtomScores = new double[nMolecules][nMolecules];
		aromaticRingScores = new double[nMolecules][nMolecules];
		featureCoordinates = new double[nMolecules][][];
		for (int i = 0; i < nMolecules; i++) {
			molecules.get(i).getAtomicGaussians();
			featureCoordinates[i] = new double[molecules.get(i).getNFeatures()][4];
		}
		nf.set(problem.getNumberFormat());
	}

	/**
	 * Resets the molecular coordinates from reference coordinates.
	 * 
	 * @param moleculeNo
	 */
	protected void copyReferenceCoordinates(int moleculeNo) {
		GaMolecule molecule = molecules.get(moleculeNo);

		molecule.copyReferenceCorordinates();
	}

	/**
	 * Decodes chromosome and overlays molecules and copies feature coordinates
	 * to chromosome.
	 * 
	 * @param remap
	 *            set if the chromosome has already been fitted. Fitting is done
	 *            in two passes with a remapping step after the first- if you
	 *            remap a fitted chromosome the overlay may change as you'd be
	 *            doing three passes of fitting.
	 * 
	 * @throws GaException
	 * 
	 * @see GaMolecule#fitMolecule(IntegerStringChromosome, int, GaMolecule,
	 *      boolean)
	 */
	public void fitMolecules(boolean remap) {
		fitted = false;
		Superposition superposition = (Superposition) problem;

		int nMolecules = molecules.size();
		for (int i = 0; i < nMolecules; i++) {
			GaMolecule molecule = molecules.get(i);
			copyReferenceCoordinates(i);
			molecule.generateConformation(binaryString,
					superposition.getBinaryEntryPoints()[i]);
		}

		logger.debug("Fitting chromosome " + getChromosomeId());

		int no = 0;
		for (int i = 0; i < nMolecules; i++) {
			GaMolecule molecule = molecules.get(i);
			if (molecule == fittingMolecule)
				continue;
			if (molecule.isFixed())
				continue;
			ok = molecule.fitMolecule(integerString,
					superposition.getIntegerEntryPoints()[no], fittingMolecule, remap);
			if (!ok)
				return;
			if (molecule.isRelaxMolecule()) {
				int start = superposition.getBinaryEntryPoints()[i]
						+ molecule.getConformationalBitLength();
				molecule.relax(binaryString, start);
			}
			no++;
		}

		// ok may not be set if all compounds are rigid
		ok = true;
		for (int i = 0; i < nMolecules; i++) {
			GaMolecule molecule = molecules.get(i);
			for (int j = 0; j < molecule.getNFeatures(); j++) {
				double[] f = molecule.getAllFeature(j).calculateCoordinate();
				Coord.copy(f, featureCoordinates[i][j]);
			}
		}

		logger.debug("Fitted chromosome " + getChromosomeId());

		currentChromosome.set(this);
		fitted = true;
	}

	/**
	 * Determines conformational energies for each molecule in the overlay.
	 * 
	 * @return
	 * 
	 */
	public double getConformationalEnergy() {
		conformationalEnergy = .0;
		int nMolecules = molecules.size();
		for (int i = 0; i < nMolecules; i++) {
			GaMolecule molecule = molecules.get(i);
			conformationalEnergies[i] = molecule.conformationalEnergy();
			conformationalEnergy += conformationalEnergies[i];
			logger.debug("conformational energy for molecule " + i + " "
					+ conformationalEnergies[i]);
		}
		conformationalEnergy = conformationalEnergy / nMolecules;
		return conformationalEnergy;
	}

	/**
	 * Determines pairwise volume integral scores. If compareAll is set we do
	 * all pairs of molecules, otherwise just do pairs with the base molecule.
	 * 
	 * @return
	 * 
	 * @see GaMolecule#gaussianIntegral(GaMolecule)
	 */
	public double getVolumeIntegral() {
		volumeIntegral = .0;
		double cnt = .0;
		Superposition superposition = (Superposition) problem;

		int nMolecules = molecules.size();
		for (int i = 0; i < nMolecules; i++) {
			GaMolecule moleculeA = molecules.get(i);
			for (int j = i + 1; j < nMolecules; j++) {
				GaMolecule moleculeB = molecules.get(j);
				// we get a cleaner overlay if we compare against everything- at
				// a CPU cost.
				if (false && !superposition.isCompareAll() && i != baseMoleculeNo
						&& j != baseMoleculeNo)
					continue;

				double vol = moleculeA.gaussianIntegral(moleculeB);
				// double vol = moleculeA.volumeIntegral(moleculeB);
				vol = (moleculeA.getWeight() + moleculeB.getWeight()) * vol / 2.0;

				volumeIntegrals[i][j] = vol;
				volumeIntegral += vol;
				logger.debug("Volume Integral for molecules " + i + " and " + j + " "
						+ vol);
				cnt++;
			}
		}

		volumeIntegral = volumeIntegral / cnt;
		return volumeIntegral * 2;
	}

	/**
	 * Determine the donor score.
	 * 
	 * The default is to do the 3D clustering scoring.
	 * 
	 * Otherwise, unless compareAll is set we just call featureOverlay.
	 * CompareAll hasn't been used for a while- it should still work, but
	 * there's no guarantee.
	 * 
	 * Can also use the 3D clustering model.
	 * 
	 * @return
	 * 
	 * @see #featureOverlay(int)
	 * @see FeatureOverlay#groupFeatures(int)
	 */
	public double donorHydrogenScore() {
		Superposition superposition = (Superposition) problem;

		// check for 3D clustering scoring
		if (superposition.isUseFeatureClustering()) {
			donorHydrogenScore = superposition.getFeatureOverlay().groupFeatures(
					FeatureType.DONOR_INTERACTION_POINT);
			return donorHydrogenScore;
		}

		if (!superposition.isCompareAll()) {
			donorHydrogenScore = featureOverlay(FeatureType.DONOR_INTERACTION_POINT);
			return donorHydrogenScore;
		}

		donorHydrogenScore = .0;
		double cnt = 0;
		int nMolecules = molecules.size();
		for (int i = 0; i < nMolecules; i++) {
			GaMolecule moleculeA = molecules.get(i);
			for (int j = i + 1; j < nMolecules; j++) {
				GaMolecule moleculeB = molecules.get(j);
				if (!superposition.isCompareAll() && i != baseMoleculeNo
						&& j != baseMoleculeNo)
					continue;

				double score = moleculeA.donorHydrogenScore(moleculeB);
				score = (moleculeA.getWeight() + moleculeB.getWeight()) * score / 2.0;

				donorHydrogenScores[i][j] = score;
				donorHydrogenScore += score;
				logger.debug("Donor Hydrogen for molecule " + i + " and " + j + " "
						+ score);
				cnt++;

			}
		}
		donorHydrogenScore = donorHydrogenScore / cnt;
		return donorHydrogenScore;
	}

	/**
	 * Determines constraint score for the overlay. The constraint score is not
	 * scaled by the number of molecules.
	 * 
	 * @return
	 */
	public double constraintScore() {
		double score = 0;
		for (GaMolecule molecule : molecules) {
			score += molecule.constraintDistance();
		}
		constraintScore = score;
		return score;
	}

	/**
	 * Determine the Acceptor score. The default is to do the 3D clustering
	 * scoring.
	 * 
	 * Otherwise, unless compareAll is set we just call featureOverlay.
	 * CompareAll hasn't been used for a while- it should still work, but
	 * there's no guarantee.
	 * 
	 * Can also use the 3D clustering model.
	 * 
	 * @return
	 * 
	 * @see #featureOverlay(int)
	 * @see FeatureOverlay#groupFeatures(int)
	 */
	public double acceptorAtomScore() {
		Superposition superposition = (Superposition) problem;

		// check for 3D clustering scoring
		if (superposition.isUseFeatureClustering()) {
			acceptorAtomScore = superposition.getFeatureOverlay().groupFeatures(
					FeatureType.ACCEPTOR_ATOM);
			return acceptorAtomScore;
		}

		if (superposition.isScalePharmacophore() && !superposition.isCompareAll()) {
			acceptorAtomScore = featureOverlay(FeatureType.ACCEPTOR_ATOM);
			return acceptorAtomScore;
		}

		acceptorAtomScore = .0;
		double cnt = 0;
		int nMolecules = molecules.size();
		for (int i = 0; i < nMolecules; i++) {
			GaMolecule moleculeA = molecules.get(i);
			for (int j = i + 1; j < nMolecules; j++) {
				GaMolecule moleculeB = molecules.get(j);

				if (!superposition.isCompareAll() && i != baseMoleculeNo
						&& j != baseMoleculeNo)
					continue;

				double score = moleculeA.acceptorAtomScore(moleculeB);
				score = (moleculeA.getWeight() + moleculeB.getWeight()) * score / 2.0;

				acceptorAtomScores[i][j] = score;
				acceptorAtomScore += score;
				logger.debug("Acceptor Atom for molecule " + i + " and " + j + " "
						+ score);
				cnt++;

			}
		}
		acceptorAtomScore = acceptorAtomScore / cnt;
		return acceptorAtomScore;
	}

	/**
	 * Determine the aromatic ring score. The default is to do the 3D clustering
	 * scoring.
	 * 
	 * Otherwise, unless compareAll is set we just call featureOverlay.
	 * CompareAll hasn't been used for a while- it should still work, but
	 * there's no guarantee.
	 * 
	 * Can also use the 3D clustering model.
	 * 
	 * @return
	 * 
	 * @see #featureOverlay(int)
	 * @see FeatureOverlay#groupFeatures(int)
	 */
	public double aromaticRingScore() {
		Superposition superposition = (Superposition) problem;

		// check for 3D clustering scoring
		if (superposition.isUseFeatureClustering()) {
			aromaticRingScore = superposition.getFeatureOverlay().groupFeatures(
					FeatureType.AROMATIC_RING);
			return aromaticRingScore;
		}

		if (superposition.isScalePharmacophore() && !superposition.isCompareAll()) {
			aromaticRingScore = featureOverlay(FeatureType.AROMATIC_RING);
			return aromaticRingScore;
		}

		aromaticRingScore = .0;
		double cnt = 0;
		int nMolecules = molecules.size();
		for (int i = 0; i < nMolecules; i++) {
			GaMolecule moleculeA = molecules.get(i);
			for (int j = i + 1; j < nMolecules; j++) {
				GaMolecule moleculeB = molecules.get(j);
				if (!superposition.isCompareAll() && i != baseMoleculeNo
						&& j != baseMoleculeNo)
					continue;

				double score = moleculeA.aromaticRingScore(moleculeB);
				score = (moleculeA.getWeight() + moleculeB.getWeight()) * score / 2.0;

				aromaticRingScores[i][j] = score;
				aromaticRingScore += score;
				logger.debug("Aromatic ring for molecule " + i + " and " + j + " "
						+ score);
				cnt++;
			}
		}
		aromaticRingScore = aromaticRingScore / cnt;
		return aromaticRingScore;
	}

	/**
	 * Determine the UserFeature score. The default is to do the 3D clustering
	 * scoring.
	 * 
	 * Otherwise, unless compareAll is set we just call featureOverlay.
	 * CompareAll hasn't been used for a while- it should still work, but
	 * there's no guarantee.
	 * 
	 * Can also use the 3D clustering model.
	 * 
	 * @return
	 * 
	 * @see #featureOverlay(int)
	 * @see FeatureOverlay#groupFeatures(int)
	 */
	public double userFeatureScore(FeatureType featureType) {
		int userFeatureNo = Feature.userFeatureNo(featureType);
		Superposition superposition = (Superposition) problem;

		// check for 3D clustering scoring
		if (superposition.isUseFeatureClustering()) {
			userFeatureScores[userFeatureNo] = superposition.getFeatureOverlay()
					.groupFeatures(featureType);
			return userFeatureScores[userFeatureNo];
		}

		if (superposition.isScalePharmacophore() && !superposition.isCompareAll()) {
			userFeatureScores[userFeatureNo] = featureOverlay(featureType);
			return userFeatureScores[userFeatureNo];
		}

		double featureScore = .0;
		double cnt = 0;
		int nMolecules = molecules.size();
		for (int i = 0; i < nMolecules; i++) {
			GaMolecule moleculeA = molecules.get(i);
			for (int j = i + 1; j < nMolecules; j++) {
				GaMolecule moleculeB = molecules.get(j);
				if (!superposition.isCompareAll() && i != baseMoleculeNo
						&& j != baseMoleculeNo)
					continue;

				double score = moleculeA.featureOverlay(featureType, moleculeB);
				score = (moleculeA.getWeight() + moleculeB.getWeight()) * score / 2.0;

				featureScore += score;
				logger.debug("Aromatic ring for molecule " + i + " and " + j + " "
						+ score);
				cnt++;

			}
		}
		featureScore = featureScore / cnt;
		userFeatureScores[userFeatureNo] = featureScore;
		return featureScore;
	}

	/**
	 * Feature overlay with Base molecule. Either one to many or the single
	 * match (one-to-one). Feature clustering is typically better.
	 * 
	 * @param featureType
	 *            integer to identify feature type
	 * @return
	 * 
	 * @see Feature
	 * @see #singleMatchFeatureOverlay(int)
	 */
	double featureOverlay(FeatureType featureType) {
		Superposition superposition = (Superposition) problem;

		if (superposition.isSingleMatchOnly())
			return singleMatchFeatureOverlay(featureType);
		else
			return multipleMatchFeatureOverlay(featureType);
	}

	/**
	 * Older feature match routine .This routine allows one-to-many mapping to
	 * the base molecule.
	 * 
	 * While the one-to-many match should work there's no guarantee and it may
	 * be incompatible with some of the pharmacophore generation features.
	 * 
	 * @param featureType
	 * @return
	 */
	double multipleMatchFeatureOverlay(FeatureType featureType) {
		List<Feature> baseFeatures = baseMolecule.getFeatureMappings(featureType)
				.getFeatures();

		Superposition superposition = (Superposition) problem;

		double score = .0;
		int nFeatures = baseFeatures.size();
		int nMolecules = molecules.size();

		for (int i = 0; i < nFeatures; i++) {

			Feature baseFeature = baseFeatures.get(i);
			int nMatched = 0;
			double fScore = 0;

			for (int j = 0; j < nMolecules; j++) {
				if (j == baseMoleculeNo)
					continue;
				GaMolecule otherMolecule = molecules.get(j);

				boolean matched = false;
				List<Feature> otherFeatures = otherMolecule.getFeatureMappings(
						featureType).getFeatures();

				int nOtherFeatures = otherFeatures.size();
				for (int k = 0; k < nOtherFeatures; k++) {

					double s = baseFeature.score(otherFeatures.get(k));
					s = (baseMolecule.getWeight() + otherMolecule.getWeight()) * s / 2.0;
					fScore += s;

					if (!matched
							&& baseFeature.getGeometricScore() > superposition
									.getGeometricWeightCriteria()) {
						matched = true;
						nMatched++;
						logger.debug("nMatched (1) " + nMatched + " g score "
								+ baseFeature.getGeometricScore() + " "
								+ baseFeature.info());
					}
				}
			}

			// Scale up if pharmacophore
			if (superposition.isScalePharmacophore()) {
				if (nMatched > 1) {
					logger.debug("nMatched (2) " + nMatched);
					fScore *= FastMath.pow(nMatched, superposition.getPharmFactor());
					baseFeature.setPharmPoint(true);
					baseFeature.setnMatched(nMatched);
				} else if (nMolecules == 2 && nMatched == 1) {
					baseFeature.setPharmPoint(true);
					baseFeature.setnMatched(nMatched);
				} else {
					baseFeature.setPharmPoint(false);
				}
			}
			score += fScore;
		}

		score = score / ((double) nMolecules - 1);
		return score;
	}

	/**
	 * Feature overlay with Base molecule. Here we find a one to one mapping
	 * between a base molecule feature and features in each other molecule. This
	 * is the method that seems to give the best results in conjunction with
	 * fitting to the base molecule (as opposed to compareAll).
	 * 
	 * @param featureType
	 *            The feature we're scoring
	 * @return
	 */
	double singleMatchFeatureOverlay(FeatureType featureType) {

		// Initializep
		for (GaMolecule molecule : molecules) {
			for (Feature feature : molecule.getFeatureMappings(featureType).getFeatures()) {
				feature.setBestScore(.0);
				feature.setBestGeometricScore(.0);
				feature.setBaseFeature(null);
			}
		}

		int nMolecules = molecules.size();

		// For each base molecule feature
		for (Feature baseFeature : baseMolecule.getFeatureMappings(featureType)
				.getFeatures()) {

			// Loop through all the other molecules
			for (int j = 0; j < nMolecules; j++) {
				if (j == baseMoleculeNo)
					continue;
				GaMolecule molecule = molecules.get(j);

				Feature bestFeature = null;
				double bestScore = .0;
				double bestGeometricScore = .0;

				// Look at all equivalent features in the molecule and
				// find the feature with the best geometric score
				for (Feature otherFeature : molecule.getFeatureMappings(featureType)
						.getFeatures()) {

					double score = baseFeature.score(otherFeature);
					double geometricScore = baseFeature.getGeometricScore();

					if (geometricScore > bestGeometricScore) {
						bestGeometricScore = geometricScore;
						bestScore = score;
						bestFeature = otherFeature;
					}
				}

				// mark the best feature found as mapping to the base
				// molecule.
				if (bestFeature != null
						&& bestGeometricScore > bestFeature.getBestGeometricScore()) {
					bestFeature.setBaseFeature(baseFeature);
					bestFeature.setBestGeometricScore(bestGeometricScore);
					bestFeature.setBestScore(bestScore);
				}
			}
		}

		// Seperate routine to sum up score
		return doMapping(false, featureType);
	}

	/**
	 * sums up the scores from the mapping found by singleMatchFeatureOverlay.
	 * This is a separate routine so that we can call it later and print out the
	 * mapping. Returns total score- if set a scale up factor for pharmacophore
	 * elucidation is applied. Also identifies all the pharmacophore points in
	 * the overlay- a pharmacophore point has to appear in at least half the the
	 * molecules with good geometry.
	 * 
	 * @param createInfo
	 *            if this is set a description about the mapping is stored in
	 *            the private String _info
	 * @param featureType
	 *            this is the feature we're scoring
	 * @return
	 */
	double doMapping(boolean createInfo, FeatureType featureType) {
		int nMolecules = molecules.size();

		if (matchedFeatures == null)
			matchedFeatures = new Feature[nMolecules];
		Superposition superposition = (Superposition) problem;

		assert !superposition.isUseFeatureClustering();

		for (GaMolecule molecule : molecules) {
			for (Feature feature : molecule.getFeatureMappings(featureType).getFeatures()) {
				feature.setPharmPoint(false);
				feature.setnMatched(0);
			}
		}

		logger.debug("Finding pharmacophore for feature" + featureType);

		double score = .0;
		for (Feature baseFeature : baseMolecule.getFeatureMappings(featureType)
				.getFeatures()) {

			int nMatched = 0;
			double fScore = .0;

			if (createInfo) {
				_info = baseMolecule.getName() + " Base Feature " + baseFeature.info()
						+ " :";
			}

			for (GaMolecule molecule : molecules) {
				if (molecule == baseMolecule) {
					continue;
				}

				for (Feature otherFeature : molecule.getFeatureMappings(featureType)
						.getFeatures()) {

					// Is this feature in the other molecule mapped to
					// the base molecule feature?
					if (otherFeature.getBaseFeature() != null
							&& otherFeature.getBaseFeature() == baseFeature) {

						double weight = (baseMolecule.getWeight() + molecule.getWeight()) / 2.0;
						fScore += otherFeature.getBestScore() * weight;

						if (createInfo)
							_info += molecule.getName() + " " + otherFeature.info() + " "
									+ nf.get().format(otherFeature.getBestScore())
									+ " [wt: " + nf.get().format(weight) + "]\n";

						// Criteria for good geometric match
						if (otherFeature.getBestGeometricScore() > superposition
								.getGeometricWeightCriteria()) {
							matchedFeatures[nMatched++] = otherFeature;
						}
					}

				}
			}

			// Scale up pharmacophore
			if (superposition.isScalePharmacophore()) {

				// Pharmacopore criteria- expect multiple matches
				if (nMatched > 1 || (nMolecules == 2 && nMatched == 1)) {
					logger.debug("nMatched (2) " + nMatched);
					if (nMolecules > 2)
						fScore *= FastMath.pow(nMatched, superposition.getPharmFactor());
					// Mark base feature as a pharmacophore point
					baseFeature.setPharmPoint(true);

					// Mark features in other molecules as pharmacophore points
					// don't do this anymore- it doesn't seem to make any sense.

					// for (int j = 0; j < nMatched; j++) {
					// matchedFeatures[j].setPharmPoint(true);
					// GaMolecule mol = matchedFeatures[j].getMolecule();
					// if (DEBUG)
					// System.err.println("Pharm Point Mol " + mol.name
					// + " " + mol.getNPharmPoints() + " "
					// + matchedFeatures[j].pharmLabel());
					// }
				}

				if (createInfo && fScore > .0)
					_info += "Number matched " + nMatched + " score "
							+ nf.get().format(fScore);
			}
			score += fScore;

		}

		// Normalize to number of pairs with base molecule
		score = score / ((double) nMolecules - 1);
		if (createInfo)
			_info += " Scaled score " + nf.get().format(score) + "\n";

		return score;
	}

	/**
	 * Returns a string containing mapping information for all features
	 * 
	 * @return
	 * 
	 * @see #doMapping(boolean, int)
	 */
	public String mappingInfo() {
		if (this != currentChromosome.get())
			throw new RuntimeException("getFitness: not currently fitted chromosome");
		Superposition superposition = (Superposition) problem;

		if (!superposition.isSingleMatchOnly())
			return "";

		String rtn = "\n";
		doMapping(true, FeatureType.DONOR_INTERACTION_POINT);
		rtn += _info + "\n";
		doMapping(true, FeatureType.ACCEPTOR_ATOM);
		rtn += _info + "\n";
		doMapping(true, FeatureType.AROMATIC_RING);
		rtn += _info;
		for (int i = 0; i < superposition.getnUserFeatureTypes(); i++) {
			UserFeatureSet featureSet = superposition.getUserFeatureSet(Feature
					.userFeatureType(i));
			rtn += "\nUser Feature " + featureSet.getFeatureSetName();
			doMapping(true, featureSet.getFeatureType());
			rtn += _info;
		}
		return rtn;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Chromosome#getFitness()
	 * 
	 * Returns the fitness of the chromosome. Chromosome must be decoded and
	 * overlay generated. Array of molecules must correspond to this chromosome.
	 * 
	 * @see #getConformationalEnergy()
	 * 
	 * @see #getVolumeIntegral()
	 * 
	 * @see #singleMatchFeatureOverlay()
	 */
	@Override
	public double getFitness() {
		if (hasFitness)
			return fitness;

		if (!ok)
			throw new RuntimeException("getFitness: chromosome is not OK "
					+ getChromosomeId());
		if (free)
			throw new RuntimeException("getFitness: chromosome is not current "
					+ getChromosomeId());
		if (!fitted)
			throw new RuntimeException("getFitness: chromosome is not fitted "
					+ getChromosomeId());

		if (this != currentChromosome.get())
			throw new RuntimeException("getFitness: not currently fitted chromosome");

		fitness = .0;
		Superposition superposition = (Superposition) problem;
		// Conformational Energies
		if (superposition.getConfWt() > 0)
			getConformationalEnergy();

		// Volume Overlay
		if (superposition.getVolumeWt() > 0)
			getVolumeIntegral();

		// All the features now:

		// Clear pharmacophore information
		for (GaMolecule molecule : molecules) {
			for (Feature feature : molecule.getAllFeatures()) {
				feature.setPharmPoint(false);
				feature.setnMatched(0);
			}
		}

		// For donors, acceptors and aromatic rings given the standard
		// configuration of compareAll not set and singleMatchOnly set we'll end
		// up calling singleMatchFeatureOverlay for each of these.

		// Donor Hydrogens
		if (superposition.getDonorHydrogenWt() > 0)
			donorHydrogenScore();

		// Acceptor Atom
		if (superposition.getAcceptorAtomWt() > 0)
			acceptorAtomScore();

		// Aromatic Ring
		if (superposition.getAromaticRingWt() > 0)
			aromaticRingScore();

		// Constraints (including breaking large cycles)
		if (superposition.getConstraintWt() > 0)
			constraintScore();

		fitness = superposition.getDonorHydrogenWt() * donorHydrogenScore
				+ superposition.getAcceptorAtomWt() * acceptorAtomScore
				+ superposition.getAromaticRingWt() * aromaticRingScore
				+ superposition.getVolumeWt() * volumeIntegral
				- superposition.getConfWt() * conformationalEnergy
				- superposition.getConstraintWt() * constraintScore;

		// UserFeatures
		for (int i = 0; i < superposition.getnUserFeatureTypes(); i++) {
			FeatureType featureType = Feature.userFeatureType(i);
			userFeatureScore(featureType);
			UserFeatureSet featureSet = superposition.getUserFeatureSet(featureType);
			fitness += featureSet.getFeatureWeight() * userFeatureScores[i];
		}

		logger.debug("Fitness " + fitness);

		hasFitness = true;
		return fitness;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Chromosome#fitnessInfo()
	 * 
	 * Returns a summary of this chromosomes fitness
	 */
	@Override
	public String fitnessInfo() {
		Superposition superposition = (Superposition) problem;

		String rtn = "";
		if (BaseSupervisor.USE_SPRINTF)
			rtn += String
					.format("Fit %7.1f [Donor %2.1f Acc %2.1f Ring %2.1f Conf %5.1f Vol %5.1f Const %5.1",
							fitness, donorHydrogenScore, acceptorAtomScore,
							aromaticRingScore, conformationalEnergy, volumeIntegral,
							constraintScore);
		else
			rtn += "Fit " + nf.get().format(fitness) + " [Donor "
					+ nf.get().format(donorHydrogenScore) + " Acc "
					+ nf.get().format(acceptorAtomScore) + " Ring "
					+ nf.get().format(aromaticRingScore) + " Conf "
					+ nf.get().format(conformationalEnergy) + " Vol "
					+ nf.get().format(volumeIntegral) + " Const "
					+ nf.get().format(constraintScore);

		for (int i = 0; i < superposition.getnUserFeatureTypes(); i++) {
			FeatureType featureType = Feature.userFeatureType(i);
			UserFeatureSet featureSet = superposition.getUserFeatureSet(featureType);
			if (BaseSupervisor.USE_SPRINTF)
				rtn += String.format(" %s $5.1f", featureSet.getFeatureSetName(),
						userFeatureScores[i]);
			else
				rtn += " " + featureSet.getFeatureSetName() + " "
						+ nf.get().format(userFeatureScores[i]);
		}

		return rtn + "]\n";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Chromosome#rebuild()
	 * 
	 * Rebuilds by fitting molecules and calculating fitness.
	 */
	@Override
	public double rebuild() {
		hasFitness = false;
		fitMolecules(false);
		return getFitness();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.Chromosome#distance(com.cairn.gape.Chromosome)
	 * 
	 * For GAPE the distance between two chromosomes is defined as the sum of
	 * square distances between all features.
	 */
	@Override
	public double distance(Chromosome c2) {
		// actually mean inter-features Square distance
		SuperpositionChromosome c = (SuperpositionChromosome) c2;

		int no = 0;
		double sqrSum = .0;
		int nMolecules = molecules.size();
		for (int i = 0; i < nMolecules; i++) {
			int nFeatures = molecules.get(i).getAllFeatures().size();
			for (int j = 0; j < nFeatures; j++) {
				sqrSum += Coord.sqrDistance(featureCoordinates[i][j],
						c.featureCoordinates[i][j]);
				no++;
			}
		}

		double ms = sqrSum / no;
		logger.debug("Inter chromosome sqr distance " + ms);

		return ms;
	}

	/**
	 * Outputs a Pharmacophore of Dummy atoms to a file
	 * 
	 * @param file
	 * @param name
	 * 
	 * @see #outputPharmMol(Writer, String)
	 */
	public void outputPharmMol(String file, String name) {
		try {
			FileWriter out = new FileWriter(new File(file));
			outputPharmMol(out, name);
			out.close();
		}

		catch (IOException ex) {
			throw new RuntimeException("outputPharm IO exception: " + ex);
		}
	}

	/**
	 * Outputs a Pharmacophore of Dummy atoms to a Writer
	 * 
	 * @param out
	 *            Writer to open structure file
	 * @param name
	 *            Molecule name
	 */
	public void outputPharmMol(Writer out, String name) throws IOException {
		Superposition superposition = (Superposition) problem;

		if (this != currentChromosome.get())
			throw new RuntimeException("outputPharm: not currently fitted chromosome");
		GaMolecule pharmMol = new GaMolecule();
		pharmMol.setName(name);

		// add scores to pharmacophore mol
		pharmMol.removeSdfField("FITNESS");
		pharmMol.addSdfField("FITNESS", String.valueOf(fitness));
		if (superposition.getVolumeWt() > .0) {
			pharmMol.removeSdfField("VOLUME_INTEGRAL");
			pharmMol.addSdfField("VOLUME_INTEGRAL", String.valueOf(volumeIntegral));
		}
		if (superposition.getConfWt() > .0) {
			pharmMol.removeSdfField("CONFORMATIONAL_ENERGY");
			pharmMol.addSdfField("CONFORMATIONAL_ENERGY",
					String.valueOf(conformationalEnergy));
		}
		if (superposition.getDonorHydrogenWt() > .0) {
			pharmMol.removeSdfField("DONOR_HYDROGEN_SCORE");
			pharmMol.addSdfField("DONOR_HYDROGEN_SCORE",
					String.valueOf(donorHydrogenScore));
		}
		if (superposition.getAcceptorAtomWt() > .0) {
			pharmMol.removeSdfField("ACCEPTOR_ATOM_SCORE");
			pharmMol.addSdfField("ACCEPTOR_ATOM_SCORE", String.valueOf(acceptorAtomScore));
		}
		if (superposition.getAromaticRingWt() > .0) {
			pharmMol.removeSdfField("AROMATIC_RING_SCORE");
			pharmMol.addSdfField("AROMATIC_RING_SCORE", String.valueOf(aromaticRingScore));
		}
		if (superposition.getConstraintWt() > .0) {
			pharmMol.removeSdfField("CONSTRAINT_SCORE");
			pharmMol.addSdfField("CONSTRAINT_SCORE", String.valueOf(constraintScore));
		}

		int nPharmPoints = getNPharmPoints();

		List<Atom> atoms = new ArrayList<>();
		List<double[]> coords = new ArrayList<>();

		for (GaMolecule molecule : molecules) {
			// if feature clustering is not used then the pharmacophore points
			// are just in the base molecule.
			if (!superposition.isUseFeatureClustering() && molecule != baseMolecule)
				continue;
			// otherwise they can be in any molecule.
			molecule.getFeatureMappings().values().stream()
					.flatMap(m -> m.getFeatures().stream()).forEach(f -> {
						if (f.isPharmPoint()) {
							atoms.add(new Atom(pharmMol, f.atomLabel(), "Du"));
							coords.add(f.getSavedCoordinate());
						}
					});
		}

		pharmMol.setCoords(coords);
		pharmMol.setAtoms(atoms);

		if (atoms.size() != nPharmPoints)
			throw new IllegalStateException(
					"outputPharm: incorrect number of pharm points");
		if (superposition.getFileFormat() == Molecule.FileType.SDF)
			pharmMol.writeSdfMol(out, "Pharmacophore Points");
		else
			pharmMol.writeSybylMol2(out, "Pharmacophore Points");

		// create second molecule which shows geometric features
		Pharmacophore pharmacophore = new Pharmacophore(this);
		Molecule geomMol = pharmacophore.toMolecule();
		geomMol.setName(name + " all features");
		if (superposition.getFileFormat() == Molecule.FileType.SDF)
			geomMol.writeSdfMol(out, "Pharmacophore Geometry");
		else
			geomMol.writeSybylMol2(out, "Pharmacophore Geometry");

	}

	/**
	 * Adds Pharmacophore Information to a structure file. The pharmacophore
	 * information put in the comments section (MOl2) or in SDF fields
	 */
	public void addPharmInfo() {
		if (this != currentChromosome.get())
			throw new RuntimeException("outputPharm: not currently fitted chromosome");

		int total = (int) molecules.stream()
				.flatMap(m -> m.getFeatureMappings().values().stream())
				.flatMap(m -> m.getFeatures().stream()).filter(Feature::isPharmPoint)
				.peek(f -> {
					f.getMolecule().addPharmDescription(f.featureLabel());
				}).count();

		assert total == getNPharmPoints();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.Chromosome#sameNiche(com.cairn.gape.Chromosome)
	 * 
	 * Uses distance function to determine if two chromosomes share a niche.
	 * 
	 * @see #distance
	 */
	@Override
	public boolean sameNiche(Chromosome c) {
		Superposition superposition = (Superposition) problem;

		if (superposition.isNichesOn()) {
			double ms = distance(c);
			if (ms < superposition.getNicheDistance() * superposition.getNicheDistance())
				return true;
		}
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Chromosome#ok()
	 * 
	 * Return true if we successfully built an overlay.
	 */
	@Override
	public boolean ok() {
		if (!fitted)
			fitMolecules(true);
		return ok;
	}

	public static SuperpositionChromosome getCurrentChromosome() {
		return currentChromosome.get();
	}

	public List<? extends GaMolecule> getMolecules() {
		return molecules;
	}

	public GaMolecule getMolecules(int i) {
		return molecules.get(i);
	}

	public static java.text.NumberFormat getNf() {
		return nf.get();
	}

	public static void setNf(java.text.NumberFormat nf) {
		SuperpositionChromosome.nf.set(nf);
	}

	public double[][][] getFeatureCoordinates() {
		return featureCoordinates;
	}

	public double[] getFeatureCoordinates(int i, int j) {
		return featureCoordinates[i][j];
	}

	public void setFeatureCoordinates(double[][][] featureCoordinates) {
		this.featureCoordinates = featureCoordinates;
	}

	public static void setCurrentChromosome(SuperpositionChromosome currentChromosome) {
		SuperpositionChromosome.currentChromosome.set(currentChromosome);
	}

	/**
	 * should only be used for debugging and testing.
	 * 
	 * @param ok
	 */
	public void setOk(boolean ok) {
		this.ok = ok;
	}

	/**
	 * should only be used for debugging and testing.
	 * 
	 * @param fitted
	 */
	public void setFitted(boolean fitted) {
		this.fitted = fitted;
	}

	public GaMolecule getFittingMolecule() {
		return fittingMolecule;
	}

	/**
	 * @return the number of pharmacophore points in the overlay.
	 * @throws GaException
	 */
	private int getNPharmPoints() {
		Superposition superposition = (Superposition) problem;

		if (this != currentChromosome.get())
			throw new RuntimeException("getNPharmPoints: not currently fitted chromosome");
		int nPharmPoints = 0;
		for (GaMolecule molecule : molecules)
			nPharmPoints += molecule.getNPharmPoints();

		if (superposition.isUseFeatureClustering())
			assert nPharmPoints == superposition.getFeatureOverlay().getNPharmPoints();

		return nPharmPoints;
	}
}
