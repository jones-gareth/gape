package com.cairn.gape.feature;

import java.util.ArrayList;
import java.util.List;

import com.cairn.gape.feature.Feature.FeatureType;

/**
 * This class encodes an array of a features that occur in a molecule. The
 * features all all of the same type (donor hydrogen, acceptor etc).
 * 
 * @author Gareth Jones
 * 
 */
public class FeatureMapping {
	// feature type id e.g. Feature.ACCEPTOR_ATOM
	private final FeatureType featureType;

	private final List<Feature> features;

	// set true if these features are to be used in chromosome mapping
	private final boolean mapping;

	/**
	 * Create the feature list.
	 * 
	 * @param set
	 * @param v
	 */
	public FeatureMapping(FeatureType set, List<Feature> v) {
		featureType = set;
		features = new ArrayList<Feature>(v);
		boolean mapping = true;
		if (features.size() == 0) {
			mapping = false;
		} else {
			mapping = features.get(0).isMappingFeature();
		}
		this.mapping = mapping;
	}

	public int getNFeatures() {
		return features.size();
	}

	public Feature getFeatures(int i) {
		return features.get(i);
	}

	public List<Feature> getFeatures() {
		return features;
	}

	public boolean isMapping() {
		return mapping;
	}

	/**
	 * @return the featureType
	 */
	public FeatureType getFeatureType() {
		return featureType;
	}

}
