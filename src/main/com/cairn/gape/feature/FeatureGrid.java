package com.cairn.gape.feature;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Class to represent Grids of features for base molecule independent fitting.
 * Currently not implemented.
 * 
 * @author Gareth Jones
 * 
 */
public class FeatureGrid {

	// Want a separate grid for each feature.

	int featureSet;

	// Inter-point distance in grid
	static double gridPointSpacing;

	// use feature radius to set this- we only want to compare features that are
	// this distance apart.
	static double featureDistance;

	// free list of spare grid points
	ArrayList<FeatureGridPoint> freeGridPoints;

	// the grid is represented by a hashmap indexed by grid point
	HashMap<FeatureGridPointIndex, FeatureGridPoint> gridHash;

	/**
	 * Creates a grid for a particular type of feature (e,g, donor hydrogen or
	 * acceptor)
	 * 
	 * @param _featureSet
	 */
	FeatureGrid(int _featureSet) {
		featureSet = _featureSet;

		freeGridPoints = new ArrayList<FeatureGridPoint>();
		gridHash = new HashMap<FeatureGridPointIndex, FeatureGridPoint>();

	}

	/**
	 * Fress up a grid. Points are added to the free list so thay can be used
	 * again.
	 */
	void freeGridHash() {
		for (FeatureGridPoint point : gridHash.values()) {
			freeGridPoints.add(point);
			point.free();
		}
		gridHash.clear();
	}

	private FeatureGridPointIndex _currentIndex = new FeatureGridPointIndex();

	/**
	 * Gets the grid point at that index or creates a new grid point at that
	 * index.
	 * 
	 * @param x
	 * @param y
	 * @param z
	 * @return
	 */
	FeatureGridPoint getOrCreateGridPoint(int x, int y, int z) {
		_currentIndex.setIndices(x, y, z);

		// get from grid
		if (gridHash.containsKey(_currentIndex))
			return gridHash.get(_currentIndex);

		// or add to grid
		FeatureGridPoint point = null;
		int no = freeGridPoints.size();
		if (no > 0)
			point = freeGridPoints.remove(no - 1);
		else
			point = new FeatureGridPoint();

		point.setIndices(x, y, z);
		gridHash.put(point.index, point);

		return point;
	}

	/**
	 * Adds a feature to a grid. The number of grid points that a feature is
	 * added to depends on the point spacing and the maximum distance between
	 * features.
	 * 
	 * @param f
	 */
	void addFeatureToGrid(Feature f) {
		int gridDistance = (int) Math.round(featureDistance - gridPointSpacing);

		int x = getIndex(f.coordinate[0]);
		int y = getIndex(f.coordinate[0]);
		int z = getIndex(f.coordinate[0]);

		for (int i = -gridDistance; i <= gridDistance; i++)
			for (int j = -gridDistance; j <= gridDistance; j++)
				for (int k = -gridDistance; k <= gridDistance; k++) {

					int _x = x + i;
					int _y = y + j;
					int _z = z + k;

					FeatureGridPoint point = getOrCreateGridPoint(_x, _y, _z);
					point.addfeature(f);
				}
	}

	/**
	 * @param _v
	 * @return Grud point index for a coordinate.
	 */
	int getIndex(double _v) {
		return (int) Math.floor(_v / gridPointSpacing);
	}

	/**
	 * main method used for testing and debugging
	 * 
	 * @param args
	 */
	public static void main(String[] args) {

	}

	/**
	 * Represents a point in the grid
	 */
	class FeatureGridPoint {
		ArrayList<Feature> features;

		FeatureGridPointIndex index;

		/**
		 * Constructor
		 */
		FeatureGridPoint() {
			index = new FeatureGridPointIndex();
			features = new ArrayList<Feature>();
		}

		/**
		 * Clesr the feature list.
		 */
		void free() {
			features.clear();
		}

		/**
		 * Adds a feature to the point.
		 * 
		 * @param _feature
		 */
		void addfeature(Feature _feature) {
			features.add(_feature);
		}

		/**
		 * Sets the grid indices for this point.
		 * 
		 * @param x
		 * @param y
		 * @param z
		 */
		void setIndices(int x, int y, int z) {
			index.setIndices(x, y, z);
		}

		/**
		 * Cenerates the co-orinates for this grid point.
		 * 
		 * @param coord
		 */
		void getCoord(double coord[]) {
			coord[0] = index.x * gridPointSpacing;
			coord[1] = index.y * gridPointSpacing;
			coord[2] = index.z * gridPointSpacing;
			coord[3] = 1;
		}
	}

}
