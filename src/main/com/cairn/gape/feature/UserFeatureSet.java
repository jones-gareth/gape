package com.cairn.gape.feature;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import org.apache.log4j.Logger;

import com.cairn.gape.feature.Feature.FeatureType;
import com.cairn.gape.ga.GaSupervisor;
import com.cairn.gape.molecule.GaMolecule;
import com.cairn.gape.utils.InfoMessageLogger;
import com.cairn.molecule.Atom;
import com.cairn.molecule.Molecule;

/**
 * Class for reading user-defined features and finding them in input structures.
 */
public class UserFeatureSet {
	private volatile FeatureType featureType;
	private volatile String featureSetName;
	private volatile double featureWeight, radius, alpha;
	private volatile boolean atomCentered;
	private volatile GaSupervisor problem;
	private volatile List<UserFeatureDefinition> userFeatureDefinitions;
	private static final Logger logger = Logger.getLogger(UserFeatureSet.class);
	static final boolean DEBUG = true;

	private UserFeatureSet() {
		;
	}

	public UserFeatureSet(GaSupervisor problem) {
		this.problem = problem;
	}

	/**
	 * Loads user defined Feature definitions.
	 * 
	 * @param filename
	 *            File containing definitions
	 */
	public void loadParameterFile(String filename) {
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(filename)));

			System.out.println("\nReading user feature definition from " + filename);
			loadParameterFile(in);
		} catch (IOException ex) {
			throw new RuntimeException("Failed to read " + filename + " " + ex);
		}
	}

	/**
	 * Loads user defined Feature definitions.
	 * 
	 * @param in
	 *            Reader to file containing definitions
	 */
	public void loadParameterFile(BufferedReader in) {
		userFeatureDefinitions = new ArrayList<UserFeatureDefinition>();

		try {
			boolean inFeatures = false;
			String line = null;
			if (problem != null)
				featureType = problem.registerNextFeatureSetNo(this);
			else
				featureType = FeatureType.USER_FEATURES1;

			while ((line = in.readLine()) != null) {
				if (line.startsWith("#") || line.startsWith(" ") || line.startsWith("\t"))
					continue;
				line = line.trim();
				if (line.equals(""))
					continue;
				if (inFeatures) {
					UserFeatureDefinition f = getFeature(line);
					userFeatureDefinitions.add(f);
				} else if (line.startsWith("name ")) {
					featureSetName = getValue(line);
				} else if (line.startsWith("weight ")) {
					String str = getValue(line);
					featureWeight = Double.valueOf(str).doubleValue();
				} else if (line.startsWith("radius ")) {
					String str = getValue(line);
					radius = Double.valueOf(str).doubleValue();
				} else if (line.startsWith("atom_centered ")) {
					String str = getValue(line);
					if (str.equalsIgnoreCase("yes") || str.equalsIgnoreCase("true"))
						atomCentered = true;
				} else if (line.startsWith("START_FEATURES"))
					inFeatures = true;
			}

			alpha = Feature.getAlpha(radius);
		} catch (IOException ex) {
			throw new RuntimeException("IOException " + ex);
		}

		logger.info("Loaded parameter file for feature set " + featureSetName);
		logger.info("feature weight " + featureWeight);
		logger.info("feature radius " + radius);
		logger.info("alpha " + alpha);
		logger.info("atom centered " + (atomCentered ? "yes" : "no"));
	}

	/**
	 * Returns the value for a line of the form 'key = value'.
	 * 
	 * @param line
	 *            Line to be processed.
	 */
	private String getValue(String line) {
		int equals = line.indexOf('=');
		// String key = line.substring(0, equals-1).trim();
		String value = line.substring(equals + 1).trim();
		return value;
	}

	/**
	 * Returns true is a feature is attractive (has a positive weight).
	 */
	public boolean isAttractive() {
		boolean attractive = featureWeight < 0 ? false : true;
		return attractive;
	}

	/**
	 * Reads a feature line
	 * 
	 * @param line
	 *            Line to be processed.
	 */
	private UserFeatureDefinition getFeature(String line) {
		StringTokenizer st = new StringTokenizer(line, "\t");
		String name = st.nextToken().trim();
		String str = st.nextToken().trim();
		double wt = Double.valueOf(str).doubleValue();
		String sln = st.nextToken().trim();
		UserFeatureDefinition f = new UserFeatureDefinition(this, name, wt, sln);
		return f;
	}

	/**
	 * Seaches molecule mol for features in this set
	 * 
	 * @param mol
	 *            Molecule to be searched
	 */
	public List<Feature> findUserFeatures(GaMolecule mol) {

		InfoMessageLogger infoMessageLogger = mol.getInfoMessageLogger();
		infoMessageLogger.infoMessageln(2, "Searching for features in user defined set "
				+ featureSetName);

		List<Feature> matches = new ArrayList<Feature>();
		// TODO can't an atom have multiple user features? And won't this clear
		// any other user feature type matches?
		if (atomCentered) {
			for (Atom atom : mol.getAtoms()) {
				atom.setUserFeature(null);
			}
		}

		for (UserFeatureDefinition userFeatureDefinition : userFeatureDefinitions) {
			List<UserFeature> userFeatures = userFeatureDefinition
					.findUserFeatureTypes(mol);
			if (!atomCentered) {
				matches.addAll(userFeatures);
			}
		}

		// Atom centered features are added to each atom.
		// TODO - can this be merged into the previous block?
		if (atomCentered) {
			for (Atom atom : mol.getAtoms()) {
				if (atom.getUserFeature() != null) {
					matches.add(atom.getUserFeature());
				}
			}
		}

		// Otherwise features are accumulated in the matches list
		int i = 0;
		for (Feature feature : matches) {
			UserFeature userFeature = (UserFeature) feature;
			userFeature.featureSetNo = i;
			infoMessageLogger.infoMessageln(2, feature.info() + " matches "
					+ userFeature.getUserFeatureType().getSln());
			i++;
		}

		return matches;
	}

	/**
	 * Loads in a user-defined feature set and matches it against one of more
	 * structures.
	 * 
	 * @param args
	 *            [] args[0] is the parameter file and args[1] is a MOL2 file of
	 *            strictures.
	 */
	public static void main(String args[]) {
		if (args.length != 2) {
			System.out.println("usage: UserFeatureSet <parameter_file> <structure_file>");
			System.exit(0);
		}

		UserFeatureSet set = new UserFeatureSet();
		set.loadParameterFile(args[0]);

		List<GaMolecule> molecules = GaMolecule.loadFiles(new String[] { args[1] },
				Molecule.FileType.MOL2, Molecule.Source.FILE);

		for (GaMolecule molecule : molecules) {
			System.out.println("Saerching " + molecule.getName());
			set.findUserFeatures(molecule);
			System.out.println();
		}

	}

	public boolean isAtomCentered() {
		return atomCentered;
	}

	public String getFeatureSetName() {
		return featureSetName;
	}

	public FeatureType getFeatureType() {
		return featureType;
	}

	public double getFeatureWeight() {
		return featureWeight;
	}

	public double getRadius() {
		return radius;
	}

	public double getAlpha() {
		return alpha;
	}

}
