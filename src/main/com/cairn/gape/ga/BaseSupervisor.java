package com.cairn.gape.ga;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.collections.MapUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Logger;

import com.cairn.common.utils.CommonUtils;
import com.cairn.gape.feature.Feature.FeatureType;
import com.cairn.gape.feature.UserFeatureSet;
import com.cairn.gape.utils.InfoMessageLogger;

/**
 * BaseSupervisor. Base class for implementing Genetic Algorithms
 * 
 * @author Gareth Jones
 * @version 0.1
 */
public class BaseSupervisor implements GaSupervisor {
	private static Logger logger;

	static {
		logger = Logger.getLogger(BaseSupervisor.class);
		// logger.setLevel(Level.DEBUG);
	}
	private java.util.Random randomGenerator;
	public static final String INFO_FILE = "ga.out";
	public static final String SEED_FILE = "ga.seed";
	private final HashMap<String, String> settings = new HashMap<String, String>();
	private HashMap<String, Double> activities;
	protected final InfoMessageLogger infoMessageLogger = new InfoMessageLogger();
	public static final boolean USE_SPRINTF = false;
	protected static final NumberFormat nf = NumberFormat.getInstance();
	protected String directory;
	protected boolean cancel = false;

	static {
		nf.setMaximumFractionDigits(2);
		nf.setGroupingUsed(false);
	}

	/**
	 * Constructor
	 */
	public BaseSupervisor() {
		String jobName = "GA_JOB_" + CommonUtils.getUniqLabel();
		Thread.currentThread().setName(jobName);
	}

	/**
	 * Sets up informational logging.
	 * 
	 * Redirects standard output to ga.stdoutLOG_FILE
	 */
	public void setupInfoOutput() {
		String logFile = fileName(INFO_FILE);
		System.out.println("See logfile " + logFile);
		if (!infoMessageLogger.setupInfoOutput(logFile))
			logger.error("Can't open " + logFile);

		if (hasKey("log_level")) {
			int logLevel = getIntValue("log_level");
			infoMessageLogger.setLogLevel(logLevel);
		}
		String date = java.text.DateFormat.getDateTimeInstance().format(new Date());
		infoMessageln("\nLog file created at " + date);
		infoMessageln("\nLog level " + infoMessageLogger.getLogLevel());
	}

	public void infoMessage(int minLevel, String message) {
		infoMessageLogger.infoMessage(minLevel, message);
	}

	public void infoMessage(String msg) {
		infoMessageLogger.infoMessage(msg);
	}

	public void infoMessageln(String msg) {
		infoMessageLogger.infoMessageln(msg);
	}

	public void infoMessageln(int minLevel, String msg) {
		infoMessageLogger.infoMessageln(minLevel, msg);
	}

	public void infoMessageln() {
		infoMessageLogger.infoMessageln();
	}

	/**
	 * @return the infoMessageLogger
	 */
	@Override
	public InfoMessageLogger getInfoMessageLogger() {
		return infoMessageLogger;
	}

	/**
	 * Loads the configuration file and sets up the random number generator.
	 * Calls {@link #getDefaults(String)}.
	 * 
	 * @param file
	 *            configuration file name
	 */
	@Override
	public void init(String file) {
		getDefaults(file);
		init();
	}

	/**
	 * Loads the configuration file from a resource and sets up the random
	 * number genersator.
	 * 
	 * @param resurce
	 */
	protected void initFromResource(String resource) {
		getDefaultsFromResource(resource);
		init();
	}

	/*
	 * Sets up random number generator and saves seed.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.ga.GaSupervisor#init()
	 */
	@Override
	public void init() {

		if (hasKey("directory")) {
			directory = getValue("directory");
			setupDirectory();
		}

		setupInfoOutput();
		logSettings();

		long seed;
		if (hasKey("seed_file")) {
			String seedFile = getStringValue("seed_file");
			seedFile = fileName(seedFile);
			try {
				BufferedReader in = new BufferedReader(new FileReader(new File(seedFile)));
				String line = in.readLine().trim();
				seed = Long.valueOf(line).longValue();
				in.close();
			} catch (FileNotFoundException ex) {
				throw new RuntimeException("can't read " + seedFile);
			} catch (IOException ex) {
				throw new RuntimeException("can't read " + seedFile);
			}
		} else {
			java.util.Random r = new java.util.Random();
			seed = r.nextLong();
		}

		logger.debug("Random number seed is " + seed);
		String seedFile = fileName(SEED_FILE);
		try {
			PrintWriter seedWriter = new PrintWriter(new FileOutputStream(new File(
					seedFile)));
			seedWriter.println(String.valueOf(seed));
			seedWriter.close();
		} catch (FileNotFoundException ex) {
			throw new RuntimeException("can't open " + seedFile);
		}
		randomGenerator = new java.util.Random(seed);
	}

	/**
	 * Sets up the random number generator. You don't need to call this if you
	 * use the init method. But if you just want to use the random number
	 * functions and not load a config file call this.
	 */
	public void setupRandomGenerator() {
		java.util.Random r = new java.util.Random();
		long seed = r.nextLong();
		randomGenerator = new java.util.Random(seed);
	}

	/**
	 * Returns a random number between 0 1and 1
	 */
	@Override
	public double normalRand() {
		return randomGenerator.nextDouble();
	}

	/**
	 * Returns a random integer within a specified (inclusive range).
	 * 
	 * @param bottom
	 *            bottom of the range
	 * @param top
	 *            top of the range
	 */
	@Override
	public int randomInt(int bottom, int top) {
		int d = randomGenerator.nextInt(top - bottom + 1);
		return bottom + d;
	}

	/**
	 * Return a random boolean.
	 */
	@Override
	public boolean randomBoolean() {
		return randomGenerator.nextBoolean();
	}

	/**
	 * Setups the program directory.
	 * 
	 */
	public void setupDirectory() {
		if (StringUtils.isEmpty(directory))
			return;
		File dir = new File(directory);
		if (!dir.exists()) {
			try {
				boolean ok = dir.mkdirs();
				if (!ok)
					throw new RuntimeException("Unable to create directory " + directory);
			} catch (SecurityException e) {
				throw new RuntimeException("Unable to create directory " + directory
						+ " security exception");
			}
		}

		if (!dir.isDirectory())
			throw new RuntimeException("Unable to use directory " + directory
					+ ": it exists and is not a directory");
		if (!dir.canRead())
			throw new RuntimeException("Unable to use directory " + directory
					+ ": cannot read it");
		if (!dir.canWrite())
			throw new RuntimeException("Unable to use directory " + directory
					+ ": cannot write to it");

		directory = dir.getAbsolutePath();
	}

	/**
	 * Adds the current directory to a filename
	 * 
	 * @param file
	 * @return
	 */
	public String fileName(String fileName) {
		if (StringUtils.isEmpty(fileName))
			return null;
		if (StringUtils.isEmpty(directory))
			return fileName;
		File file = new File(fileName);
		if (file.isAbsolute())
			return fileName;
		File dirFile = new File(directory, fileName);
		return dirFile.getPath();
	}

	/**
	 * Loads the defaults file.
	 * 
	 * @param file
	 *            defaults/configuration filename
	 */
	private void getDefaults(String file) {
		try {
			logger.debug("reading settings from file " + file);
			BufferedReader in = new BufferedReader(new FileReader(new File(file)));
			getDefaults(in);
		} catch (IOException ex) {
			throw new RuntimeException("Failed to read " + file + " " + ex);
		}
	}

	/**
	 * Loads defaults from a resource.
	 * 
	 * @param resource
	 */
	private void getDefaultsFromResource(String resource) {
		InputStream stream = getClass().getResourceAsStream(resource);
		if (stream == null)
			throw new RuntimeException("Unable to load resource " + resource);
		BufferedReader in = new BufferedReader(new InputStreamReader(stream));
		getDefaults(in);
	}

	/**
	 * Loads the defaults file.
	 * 
	 * @param in
	 *            BufferedReader for defaults/configuration filename
	 */
	public void getDefaults(BufferedReader in) {
		getDefaults(in, true);
	}

	/**
	 * Loads the defaults file.
	 * 
	 * @param in
	 *            BufferedReader for defaults/configuration filename
	 * @param log
	 *            set to log settings
	 */
	public void getDefaults(BufferedReader in, boolean log) {
		settings.clear();
		activities = new HashMap<String, Double>();
		getDefaults(in, log, settings, activities);
	}

	/**
	 * Loads the defaults file. Use this function if you don't want settings to
	 * go into class variables.
	 * 
	 * @param in
	 *            BufferedReader for defaults/configuration filename
	 * @param log
	 *            set to log settings
	 * @param hash
	 *            HashMap for settings
	 * @param act
	 *            HashMap for activity settings
	 */
	public void getDefaults(BufferedReader in, boolean log, Map<String, String> hash,
			Map<String, Double> act) {
		try {
			String line = null;
			while ((line = in.readLine()) != null) {
				if (line.startsWith("#") || line.startsWith(" ") || line.startsWith("\t"))
					continue;
				logger.debug("defaults line " + line);
				line = line.trim();
				if (line.equals(""))
					continue;
				else if (line.startsWith("activity "))
					addActivity(line, act, log);
				else
					addKey(line, hash, log);
			}
			in.close();
		} catch (IOException ex) {
			throw new RuntimeException("IOException " + ex);
		}

	}

	/**
	 * Processes a key pair line.
	 * 
	 * @param line
	 *            string containing key pair (format 'key = value')
	 * @param hash
	 *            HashMap to put key pair in
	 * @param log
	 *            set to log key pair
	 */
	private void addKey(String line, Map<String, String> hash, boolean log) {
		int equals = line.indexOf('=');
		assert equals != -1;
		String key = line.substring(0, equals - 1).trim();
		String value = line.substring(equals + 1).trim();
		hash.put(key, value);
	}

	/**
	 * Processes an activity line
	 * 
	 * @param line
	 *            string containing activity (format 'activity molname pki')
	 * @param hash
	 *            HashMap to put key pair in
	 * @param log
	 *            set to log key pair
	 */
	private void addActivity(String line, Map<String, Double> act, boolean log) {
		if (!line.startsWith("activity "))
			return;
		line = line.substring(9);
		int space = line.indexOf(' ');
		String actString = line.substring(0, space).trim();
		String molName = line.substring(space).trim();
		Double activity = Double.valueOf(actString);
		act.put(molName, activity);
	}

	private void logSettings() {
		if (infoMessageLogger.getLogLevel() <= 0)
			return;
		List<String> keys = new ArrayList<String>(settings.keySet());
		Collections.sort(keys);
		for (String key : keys)
			infoMessageln(key + " = " + settings.get(key));

		infoMessageln();
		if (MapUtils.isNotEmpty(activities)) {
			List<String> molNames = new ArrayList<String>(activities.keySet());
			Collections.sort(molNames);
			for (String molName : molNames)
				infoMessageln("activity " + activities.get(molName) + " " + molName);
		}
	}

	/**
	 * Retrieves molecular activity
	 * 
	 * @param molName
	 *            name of molecule
	 */
	@Override
	public double getActivity(String molName) {
		if (!activities.containsKey(molName))
			throw new RuntimeException("No activity for " + molName);
		double activity = activities.get(molName);
		return activity;
	}

	/**
	 * Returns a settings key value
	 * 
	 * @param key
	 *            the key of interest
	 */
	public String getValue(String key) {
		return settings.get(key);
	}

	/**
	 * Checks to see if a settings key is present.
	 * 
	 * @param key
	 *            the key of interest
	 */
	@Override
	public boolean hasKey(String key) {
		return settings.containsKey(key);
	}

	/**
	 * Returns a settings value as an integer
	 * 
	 * @param key
	 *            the key of interest
	 */
	@Override
	public int getIntValue(String key) {
		String v = getValue(key);
		if (v == null)
			throw new RuntimeException("key " + key + " not found");
		return Integer.valueOf(v).intValue();
	}

	/**
	 * Returns a settings value as a boolean
	 * 
	 * @param key
	 *            the key of interest
	 */
	@Override
	public boolean getBooleanValue(String key) {
		String v = getValue(key);
		if (v == null)
			throw new RuntimeException("key " + key + " not found");
		boolean val = false;
		if (v.toLowerCase().equals("yes") || v.toLowerCase().equals("true")
				|| v.equals("1"))
			val = true;
		else if (v.toLowerCase().equals("no") || v.toLowerCase().equals("false")
				|| v.equals("0"))
			val = false;
		else
			throw new RuntimeException("key " + key + " doesn't have a boolean value ["
					+ v + "]");
		return val;
	}

	/**
	 * Returns a settings value as a double
	 * 
	 * @param key
	 *            the key of interest
	 */
	@Override
	public double getDoubleValue(String key) {
		String v = getValue(key);
		if (v == null)
			throw new RuntimeException("key " + key + " not found");
		return Double.valueOf(v).doubleValue();
	}

	/**
	 * Returns a settings value as a String
	 * 
	 * @param key
	 *            the key of interest
	 */
	@Override
	public String getStringValue(String key) {
		String v = getValue(key);
		if (v == null)
			throw new RuntimeException("key " + key + " not found");
		return v;
	}

	/**
	 * Returns the class NumberFormat
	 */
	@Override
	public NumberFormat getNumberFormat() {
		return nf;
	}

	/**
	 * Throws error- method not defined
	 */
	@Override
	public int getnFeatureTypes() {
		throw new RuntimeException("Method getNFeatures not defined");
	}

	/**
	 * Throws error- method not defined
	 */
	@Override
	public int getnUserFeatureTypes() {
		throw new RuntimeException("Method getNUserFeatures not defined");
	}

	/**
	 * Throws error- method not defined
	 */
	@Override
	public UserFeatureSet getUserFeatureSet(FeatureType no) {
		throw new RuntimeException("Method getUserFeatureSet not defined");
	}

	/**
	 * Throws error- method not defined
	 */
	@Override
	public FeatureType registerNextFeatureSetNo(UserFeatureSet userFeatureSet) {
		throw new RuntimeException("Method getUserFeatureSet not defined");
	}

	/**
	 * Returns Logging level
	 */
	public int getLogLevel() {
		return infoMessageLogger.getLogLevel();
	}

	/**
	 * Interface for listening to progress of the algorithm
	 */
	public interface ProgressListener {
		public void reportProgress(ProgressReport progressReport);
	}

	protected List<ProgressListener> progressListeners = new ArrayList<ProgressListener>();

	public void addProgressListener(ProgressListener progressListener) {
		progressListeners.add(progressListener);
	}

	public void removeProgressListener(ProgressListener progressListener) {
		progressListeners.remove(progressListener);
	}

	/**
	 * Inform all listeners of progress
	 * 
	 * @param progressReport
	 */
	protected void tellListeners(ProgressReport progressReport) {

		for (ProgressListener progressListener : progressListeners)
			progressListener.reportProgress(progressReport);

	}

	/**
	 * Request that the application terminate.
	 */
	public void cancel() {
		cancel = true;
	}
}
