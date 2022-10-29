package com.cairn.gape.ga;

import java.text.NumberFormat;

import com.cairn.gape.feature.Feature.FeatureType;
import com.cairn.gape.feature.UserFeatureSet;
import com.cairn.gape.utils.InfoMessageLogger;

/**
 * Interface to support different GAs. Allows other objects to be able to access
 * the main GA class. Some of these objects may not be defined by GAs- in which
 * case an exception should be thrown.
 * 
 * @author Gareth Jones
 * @see com.cairn.gape.ga.BaseSupervisor
 * 
 */
public interface GaSupervisor {
	/**
	 * Initializes the GA- runtime configuration is specified in the file
	 * 
	 * @param file
	 */
	void init(String file);

	/**
	 * Initializes the GA without any configuration
	 * 
	 */
	void init();

	/**
	 * @return a random number between 0 1and 1
	 */
	double normalRand();

	/**
	 * @param bottom
	 * @param top
	 * @return a random integer between top and bottom
	 */
	int randomInt(int bottom, int top);

	/**
	 * @return a random boolean
	 */
	boolean randomBoolean();

	/**
	 * GA settings are stored as key/value pairs. This returns true if the key
	 * is present.
	 * 
	 * @param key
	 * @return
	 */
	boolean hasKey(String key);

	/**
	 * @param key
	 * @return the integer value of the key
	 */
	int getIntValue(String key);

	/**
	 * @param key
	 * @return the boolean value of the key
	 */
	boolean getBooleanValue(String key);

	/**
	 * @param key
	 * @return the double value of the key
	 */
	double getDoubleValue(String key);

	/**
	 * @param key
	 * @return the string value of the key
	 */
	String getStringValue(String key);

	/**
	 * @param molName
	 * @return any activity associated wtih a molecule
	 */
	double getActivity(String molName);

	/**
	 * GAPE or smilar GA only.
	 * 
	 * @return The number of feature sets defined.
	 */
	int getnFeatureTypes();

	/**
	 * GAPE or similar only
	 * 
	 * @return the number of user feature sets defined.
	 */
	int getnUserFeatureTypes();

	/**
	 * GAPE or similar only,
	 * 
	 * @param no
	 * @return Returns a user feature set no (no is an index into user feature
	 *         sets - not alkl feature sets).
	 */
	UserFeatureSet getUserFeatureSet(FeatureType no);

	/**
	 * GAPE or similar only. Registers a new feature set.
	 * 
	 * @param userFeatureSet
	 * @return user feature set number
	 */
	FeatureType registerNextFeatureSetNo(UserFeatureSet userFeatureSet);

	/**
	 * @return a number format object for consistent formatting.
	 */
	NumberFormat getNumberFormat();

	/**
	 * @return information logger
	 */
	InfoMessageLogger getInfoMessageLogger();
}
