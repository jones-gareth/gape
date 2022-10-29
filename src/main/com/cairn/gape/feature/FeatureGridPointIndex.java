package com.cairn.gape.feature;

/**
 * Class used to index points on the Feature Grid. Use for hashmap keys.
 * Currently not implemented.
 * 
 * @author Gareth Jones
 * 
 */
class FeatureGridPointIndex {

	// Grid point indices
	int x, y, z;

	FeatureGridPointIndex() {
	}

	/**
	 * Sets the grid point indices
	 * 
	 * @param _x
	 * @param _y
	 * @param _z
	 */
	void setIndices(int _x, int _y, int _z) {
		x = _x;
		y = _y;
		z = _z;
	}

	/*
	 * Returns true if indices match.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	public boolean equals(Object obj) {
		if (!(obj instanceof FeatureGridPointIndex))
			return false;

		FeatureGridPointIndex other = (FeatureGridPointIndex) obj;

		if (x == other.x && y == other.y && z == other.z)
			return true;
		return false;
	}

}
