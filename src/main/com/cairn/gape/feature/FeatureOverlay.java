package com.cairn.gape.feature;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.math3.util.FastMath;

import com.cairn.common.utils.Coord;
import com.cairn.gape.Superposition;
import com.cairn.gape.feature.Feature.FeatureType;
import com.cairn.gape.molecule.GaMolecule;

/**
 * A class to encode pharmacophore point scoring without reference to the base
 * molecule. Uses 3D relocation clustering of feature fitting points. When
 * scoring a cluster of fitting points each feature is tried in turn and the
 * best scoring feature is selected as a "base" feature.
 * 
 * @author Gareth Jones
 * 
 */
public class FeatureOverlay {

	private final List<? extends GaMolecule> molecules;

	private final int nMolecules;

	private final Superposition problem;

	// minDistance and maxDistance are used to join and split feature points.
	private volatile double score, minDistance, maxDistance, minSqrDistance,
			maxSqrDistance;

	// set this to true to require all features in a set to belong to different
	// molecules.
	private final boolean differentMoleculeSet = true;

	private final Map<FeatureType, FeaturePointSet> featurePointSets = new HashMap<>();

	private static final boolean DEBUG = false;

	private static final int MAX_RELOCATIONS = 100;

	/**
	 * Creates an object for doing the feature overlay. Initializes
	 * FeaturePointSets for built-in and user features. Sets the pharmacophore
	 * counts for molecules to zero. If you're re-scoring (such as in a GA run)
	 * you need to set the molecule pharmacophore counts to zero before calling
	 * scoring functions.
	 * 
	 * @param problem
	 */
	public FeatureOverlay(Superposition problem) {
		this.problem = problem;
		molecules = problem.getMolecules();
		nMolecules = problem.getNMolecules();

		assert nMolecules > 2;

		featurePointSets.clear();

		// Donor Hydrogen
		if (problem.getDonorHydrogenWt() > 0)
			setup(FeatureType.DONOR_INTERACTION_POINT);

		// Acceptor Atom
		if (problem.getAcceptorAtomWt() > 0)
			setup(FeatureType.ACCEPTOR_ATOM);

		// Aromatic Ring
		if (problem.getAromaticRingWt() > 0)
			setup(FeatureType.AROMATIC_RING);

		// UserFeatures
		IntStream.range(0, problem.getnUserFeatureTypes()).forEach(i -> {
			FeatureType type = FeatureType.valueOf("USER_FEATURES" + i);
			setup(type);
		});
	}

	/**
	 * Initializes a feature set. Creates a FeaturePointSet.
	 * 
	 * @param featureType
	 */
	private void setup(FeatureType featureType) {
		// find all features in this set
		List<Feature> features = molecules.stream()
				.flatMap(m -> m.getFeatureMappings(featureType).getFeatures().stream())
				.collect(Collectors.toList());

		// create FeaturePointSet
		FeaturePointSet featurePointSet = new FeaturePointSet(featureType, features);
		featurePointSets.put(featureType, featurePointSet);
	}

	/**
	 * Takes the current saved feature coordinates and determines GAPE score.
	 * More typically we'll call scoreFeature for each feature.
	 * 
	 * @return feature score.
	 */
	public double score() {
		score = .0;

		// Donor Hydrogen
		if (problem.getDonorHydrogenWt() > 0)
			score += problem.getDonorHydrogenWt()
					* groupFeatures(FeatureType.DONOR_INTERACTION_POINT);

		// Acceptor Atom
		if (problem.getAcceptorAtomWt() > 0)
			score += problem.getAcceptorAtomWt()
					* groupFeatures(FeatureType.ACCEPTOR_ATOM);

		// Aromatic Ring
		if (problem.getAromaticRingWt() > 0)
			score += problem.getAromaticRingWt()
					* groupFeatures(FeatureType.AROMATIC_RING);

		// UserFeatures
		IntStream.range(1, 10).forEach((no) -> {
			FeatureType featureType = FeatureType.valueOf("USER_FEATURES" + no);
			double userfeatureScore = groupFeatures(featureType);
			UserFeatureSet featureSet = problem.getUserFeatureSet(featureType);
			score += featureSet.getFeatureWeight() * userfeatureScore;
		});

		return score;
	}

	/**
	 * Determines GAPE score for a feature. First clusters in 3D space then
	 * scores clusters.
	 * 
	 * @param featureSetNo
	 * @return
	 */
	public double groupFeatures(FeatureType featureType) {

		FeaturePointSet featurePointSet = featurePointSets.get(featureType);
		assert featurePointSet.featureType == featureType;
		// group points
		featurePointSet.groupPoints();
		// score
		double featureScore = featurePointSet.score();
		return featureScore;
	}

	/**
	 * Assumes that all features have been grouped and scored.
	 * 
	 * @return total number of pharmacophore points for this overlay.
	 */
	public int getNPharmPoints() {
		int nPharmPoints = 0;
		for (FeaturePointSet featurePointSet : featurePointSets.values()) {
			if (featurePointSet == null)
				continue;
			int no = featurePointSet.nPharmPoints;
			assert no == featurePointSet.getPharmCount();
			nPharmPoints += no;
		}
		return nPharmPoints;
	}

	/**
	 * @param featureSetNo
	 * @return number of pharmacophore points for this feature set.
	 */
	public int getNPharmPoints(FeatureType featureType) {
		return featurePointSets.get(featureType).nPharmPoints;
	}

	/**
	 * This class is used to store a set of feature points for a given feature
	 * type
	 * 
	 */
	private class FeaturePointSet {
		// set of points for the feature and free list
		private final Set<FeaturePoint> featurePoints, freeFeaturePoints;

		private final FeatureType featureType;

		private int nPharmPoints;

		private double score;

		List<Feature> features;

		/**
		 * Creates a feature point set for this feature type and the given
		 * features.
		 * 
		 * @param featureSetNo
		 * @param features
		 */
		private FeaturePointSet(FeatureType featureSetNo, List<Feature> features) {
			this.featureType = featureSetNo;
			this.features = features;

			featurePoints = new HashSet<FeaturePoint>();
			freeFeaturePoints = new HashSet<FeaturePoint>();

			// In the worst case we need a number of feature points equal to the
			// total number of features across all molecules
			for (int i = 0; i < features.size(); i++) {
				FeaturePoint featurePoint = new FeaturePoint();
				freeFeaturePoints.add(featurePoint);
			}
			assert featurePoints.size() + freeFeaturePoints.size() == features.size();
		}

		/**
		 * Assumes that all point have been grouped.
		 * 
		 * @return pharmacophore score for the set.
		 */
		public double score() {
			assert featurePoints.size() > 0;
			assert getPharmCount() == 0;

			// remove pharm point information here so that we always know
			// molecule pharm counts and feature counts are the same. This is
			// also done in groupPoints
			nPharmPoints = 0;
			for (Feature feature : features) {
				feature.setPharmPoint(false);
				feature.setnMatched(0);
			}

			score = 0;
			for (FeaturePoint featurePoint : featurePoints) {
				score += featurePoint.score();
				if (featurePoint.pharmPoint)
					nPharmPoints++;
			}
			return score;
		}

		/**
		 * @return number of pharmacophore points from checking feature array
		 */
		private int getPharmCount() {
			int no = 0;
			for (Feature feature : features) {
				if (feature.isPharmPoint())
					no++;
			}
			return no;
		}

		/**
		 * Groups all the points. First performs initial assignment of feature
		 * poonts then does relocation until convergence.
		 * 
		 * @see #relocate()
		 * @see #addFeature(Feature)
		 */
		public void groupPoints() {

			// clear out any previous grouping and scoring
			for (FeaturePoint featurePoint : featurePoints) {
				freeFeaturePoints.add(featurePoint);
			}
			featurePoints.clear();
			assert featurePoints.size() + freeFeaturePoints.size() == features.size();

			// Features are cleared in score()
			// for (Feature f : features)
			// f.setPharmPoint(false);

			// feature coordinates need to be this close (maxDistance) to form a
			// cluster
			maxDistance = Feature.getRadius() * 2;

			// features this close (minDistance) from the same molecule will be
			// allowed in the same cluster- otherwise features from the same
			// molecule will always occupy their own cluster
			if (differentMoleculeSet)
				minDistance = 0;
			else
				minDistance = problem.getFinishFittingRadius();
			maxSqrDistance = maxDistance * maxDistance;
			minSqrDistance = minDistance * minDistance;

			// / initial assignment
			for (Feature feature : features)
				addFeature(feature);

			// relocations
			boolean relocate = true;
			int nRelocations = 0;
			while (relocate) {
				nRelocations++;
				relocate = relocate();
				if (nRelocations == MAX_RELOCATIONS) {
					// if (DEBUG)
					System.out.println("groupPoints no convergance after "
							+ MAX_RELOCATIONS + " relocations");
					break;
				}
			}

			if (differentMoleculeSet)
				for (FeaturePoint featurePoint : featurePoints)
					assert featurePoint.checkSingleMolecules();
		}

		/**
		 * Adds a feature to the initial grouping.
		 * 
		 * @param feature
		 */
		private void addFeature(Feature feature) {
			GaMolecule molecule = feature.getMolecule();
			boolean added = false;

			// loop though any extiting points
			for (FeaturePoint featurePoint : featurePoints) {
				double sqrDistance = featurePoint.sqrDistance(feature);
				boolean containsMolecule = featurePoint.containsMoleculeFeature(molecule);
				if (containsMolecule && differentMoleculeSet)
					continue;
				else if (containsMolecule) {
					// criteria if the current point already contains a feature
					// from this molecule.
					if (sqrDistance < minSqrDistance) {
						// add updating center
						featurePoint.addPoint(feature, true);
						added = true;
						break;
					}
				} else {
					// criteria if the current point does not contain any
					// feature from this molecule.
					if (sqrDistance < maxSqrDistance) {
						// add updating center
						featurePoint.addPoint(feature, true);
						added = true;
						break;
					}
				}
			}

			if (!added) {
				// no close points so create a new feature point with this
				// feature.
				FeaturePoint newFeaturePoint = freeFeaturePoints.iterator().next();
				freeFeaturePoints.remove(newFeaturePoint);
				newFeaturePoint.seed(feature);
				featurePoints.add(newFeaturePoint);
			}
		}

		/**
		 * Performs a single cycle of relocation. Moves each feature to a closer
		 * point (if it's not already assigned to the closest point).
		 * 
		 * @return true if any features are relocated.
		 */
		private boolean relocate() {
			boolean relocate = false;

			assert featurePoints.size() + freeFeaturePoints.size() == features.size();

			for (Feature feature : features) {
				// foreach feature find the closest featurePoint and the point
				// containing that feature.
				FeaturePoint closestPoint = null, currentPoint = null;
				double minDistance = Double.MAX_VALUE;
				for (FeaturePoint featurePoint : featurePoints) {
					double distance = featurePoint.sqrDistance(feature);
					if (distance < minDistance) {
						minDistance = distance;
						closestPoint = featurePoint;
					}
					if (featurePoint.hasFeature(feature))
						currentPoint = featurePoint;

				}

				assert closestPoint != null : "Unable to find closest point";
				assert currentPoint != null : "Unable to find current point";

				// check if feature is alredy in closest point.
				if (closestPoint == currentPoint)
					continue;

				// check constraint that each point only contains one feature
				// from a molecule
				if (differentMoleculeSet
						&& closestPoint.containsMoleculeFeature(feature.getMolecule()))
					continue;

				// otherwise move it to closest point.
				relocate = true;
				currentPoint.removePoint(feature, false);
				closestPoint.addPoint(feature, false);

				// if the relocation leaves an empty point remove it.
				if (currentPoint.getNFeatures() == 0) {
					featurePoints.remove(currentPoint);
					freeFeaturePoints.add(currentPoint);
					assert featurePoints.size() + freeFeaturePoints.size() == features
							.size();
				}

			}

			// update centers
			for (FeaturePoint featurePoint : featurePoints)
				featurePoint.recalculateCenter();

			return relocate;
		}
	}

	/**
	 * This class represents a cluster of pharmacophore features.
	 * 
	 */
	private class FeaturePoint {
		private final double center[] = new double[] { 0, 0, 0, 1 },
				sum[] = new double[] { 0, 0, 0, 1 };

		private double score = 0;

		private final Set<Feature> features;

		private Feature baseFeature;

		private int nMatched, testMatched;

		private boolean pharmPoint;

		public FeaturePoint() {
			features = new HashSet<Feature>();
		}

		/**
		 * @return the number of features in this point.
		 */
		@SuppressWarnings("unused")
		public int getSize() {
			return features.size();
		}

		/**
		 * Adds a feature to this point.
		 * 
		 * @param f
		 * @param updateCenter
		 *            set true to update centriod.
		 */
		public void addPoint(Feature f, boolean updateCenter) {
			assert !features.contains(f) : "feature already in set";

			double coordinate[] = f.getSavedCoordinate();
			features.add(f);
			if (differentMoleculeSet)
				assert checkSingleMolecules();
			double no = features.size();

			for (int i = 0; i < 3; i++) {
				sum[i] += coordinate[i];
				if (updateCenter)
					center[i] = sum[i] / no;
			}

		}

		/**
		 * Updates centriod coordinates.
		 */
		public void recalculateCenter() {
			double no = features.size();
			for (int i = 0; i < 3; i++) {
				center[i] = sum[i] / no;
			}
		}

		/**
		 * Removes a feature from this point.
		 * 
		 * @param f
		 * @param updateCenter
		 *            set true to recalculate centroid coordinates.
		 */
		public void removePoint(Feature f, boolean updateCenter) {
			assert features.contains(f) : "feature not in set";
			features.remove(f);
			double coordinate[] = f.getSavedCoordinate();
			double no = features.size();
			for (int i = 0; i < 3; i++) {
				sum[i] += coordinate[i];
				if (updateCenter)
					center[i] = sum[i] / no;
			}
		}

		/**
		 * @param f
		 * @return The distance from point centroid to this feature.
		 */
		public double sqrDistance(Feature f) {
			double coordinate[] = f.getSavedCoordinate();
			return Coord.sqrDistance(coordinate, center);
		}

		/**
		 * @param molecule
		 * @return true if the point already contains a feature from this
		 *         molecule
		 */
		public boolean containsMoleculeFeature(GaMolecule molecule) {
			for (Feature feature : features) {
				if (feature.getMolecule() == molecule)
					return true;
			}
			return false;
		}

		/**
		 * @return true if all features in this set are froim different
		 *         molecules.
		 */
		private boolean checkSingleMolecules() {
			for (Feature feature : features) {
				for (Feature other : features) {
					if (feature == other)
						continue;
					if (feature.getMolecule() == other.getMolecule())
						return false;
				}
			}
			return true;
		}

		/**
		 * Determines the score for this point. We try each feature in the point
		 * as the "base" feature and keep the one with the best score.
		 * 
		 * @return
		 */
		public double score() {

			double maxScore = -Double.MAX_VALUE;
			baseFeature = null;
			pharmPoint = false;
			nMatched = 0;
			score = .0;

			if (features.size() < 2)
				return .0;

			if (DEBUG)
				System.out.println("entering score function");

			// test each feature in turn.
			for (Feature testFeature : features) {
				// testFeature.setPharmPoint(false);

				double test = scoreFeature(testFeature);
				if (DEBUG)
					System.out.println("test " + test);
				// record best
				if (test > maxScore) {
					if (DEBUG)
						System.out.println("setting base feature to  " + testFeature);
					maxScore = test;
					baseFeature = testFeature;
					nMatched = testMatched;
				}

				// if we have two features the scores will be identical. Should
				// probably pick the feature with the highest group
				// probability.
				if (features.size() == 2)
					break;
			}

			if (baseFeature == null) {
				System.out.println("oops");
			}

			assert baseFeature != null;

			// check for pharmacophore point.
			if (nMatched > 1 || (nMolecules == 2 && nMatched == 1)) {
				assert !baseFeature.isPharmPoint();
				baseFeature.setPharmPoint(true);
				baseFeature.setnMatched(nMatched);
				pharmPoint = true;
			}

			// Normalize to maximum number of pairs
			maxScore = maxScore / ((double) nMolecules - 1);
			score = maxScore;

			return score;
		}

		/**
		 * Test a feature as the base feature.
		 * 
		 * @param testFeature
		 * @return score that we get from using this feature a the base feature.
		 */
		private double scoreFeature(Feature testFeature) {
			double testScore = .0;
			testMatched = 0;

			GaMolecule testMolecule = testFeature.getMolecule();

			// sum pair-wise score with all other features
			for (Feature otherFeature : features) {
				if (otherFeature == testFeature)
					continue;
				GaMolecule otherMolecule = otherFeature.getMolecule();

				// only one feature per molecule can contribute to a
				// pharmacophore point.
				if (otherMolecule == testMolecule)
					continue;

				// get molecule weighting
				double weight = (testMolecule.getWeight() + otherMolecule.getWeight()) / 2.0;
				// determine feature pair score.
				double featureScore = weight * testFeature.score(otherFeature);
				testScore += featureScore;

				// count any "pharmacophore" matches
				if (problem.isScalePharmacophore()
						&& testFeature.getGeometricScore() > problem
								.getGeometricWeightCriteria())
					testMatched++;
			}

			// detmine any pharmacophore scale up.
			if (testMatched > 2) {
				if (DEBUG)
					System.out.println("testMatched (2) " + testMatched);
				if (testMatched > 2)
					testScore *= FastMath.pow(testMatched, problem.getPharmFactor());
			}

			return testScore;
		}

		/**
		 * Removes previous feature and scoring information from this point.
		 * Seeds point with single feature. This should be the entry point for
		 * creating new points.
		 * 
		 * @param feature
		 */
		public void seed(Feature feature) {
			features.clear();
			for (int i = 0; i < 3; i++) {
				sum[i] = .0;
				center[i] = .0;
			}
			sum[3] = center[3] = 1.0;
			nMatched = 0;
			score = 0;
			pharmPoint = false;

			addPoint(feature, true);
		}

		/**
		 * @param feature
		 * @return ture if this point contains the feature.
		 */
		public boolean hasFeature(Feature feature) {
			return features.contains(feature);
		}

		/**
		 * @return the number of features in this point.
		 */
		public int getNFeatures() {
			return features.size();
		}

	}
}
