package com.cairn.gape.utils.config;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.StringWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.lang3.BooleanUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Logger;

import com.cairn.gape.ga.BaseSupervisor;

/**
 * A model for the configuration file
 * 
 * @author gjones
 * 
 */
public class ConfigModel {
	private ArrayList<ConfigSection> sections;
	private HashMap<String, ConfigEntry<?>> allEntries;

	private String defaultsResource = "/com/cairn/gape/utils/superposition.conf";
	private String header;

	private StructureFileConfigEntry structureFileEntry;

	private static Logger logger;
	static {
		logger = Logger.getLogger(ConfigModel.class);
		// logger.setLevel(Level.DEBUG);
	}

	/**
	 * Loads configuration from File f.
	 * 
	 * @param f
	 */
	public void loadFile(File f) {
		try {
			if (sections != null)
				clearOptionalValues();
			readConfigFile();
			BufferedReader in = new BufferedReader(new FileReader(f));
			logger.info("Setting values from file " + f);
			loadValues(in, false);
		} catch (FileNotFoundException ex) {
			logger.error("File " + f + " not found");
			return;
		}
		logger.info("File " + f + " loaded");
	}

	/**
	 * Loads settings from a string
	 * 
	 * @param str
	 */
	public void loadFromString(String str) {
		if (sections != null)
			clearOptionalValues();
		readConfigFile();
		BufferedReader in = new BufferedReader(new StringReader(str));
		logger.info("Setting values from string ");
		loadValues(in, false);
	}

	/**
	 * Loads in settings from a BufferedReader {@see BaseSupervisor}. Fill
	 * values on the config entry list and may also set defaults
	 * 
	 * @param in
	 */
	@SuppressWarnings("unchecked")
	private void loadValues(BufferedReader in, boolean setDefaults) {
		Map<String, String> stringValues = new HashMap<String, String>();
		// Activities stores molecule activities
		Map<String, Double> activities = new HashMap<String, Double>();

		BaseSupervisor bs = new BaseSupervisor();
		bs.getDefaults(in, false, stringValues, activities);

		for (ConfigEntry<?> entry : allEntries.values()) {
			String strValue = stringValues.get(entry.getName().toLowerCase());
			if (StringUtils.isEmpty(strValue)) {
				if (entry.isOptional()) {
					entry.setUse(false);
					entry.setValue(null);
					if (setDefaults)
						entry.setDefValue(null);
				} else
					logger.warn("no value for " + entry.getName());
				continue;
			}
			strValue = strValue.trim();
			Class<?> clazz = entry.getClazz();
			if (clazz == Boolean.class) {
				boolean value = strValue.equals("1") ? true : BooleanUtils
						.toBoolean(strValue);
				((ConfigEntry<Boolean>) entry).setValue(value);
				if (setDefaults)
					((ConfigEntry<Boolean>) entry).setDefValue(value);
			} else if (clazz == Integer.class) {
				int value = Integer.parseInt(strValue);
				((ConfigEntry<Integer>) entry).setValue(value);
				if (setDefaults)
					((ConfigEntry<Integer>) entry).setDefValue(value);
			} else if (clazz == Double.class) {
				double value = Double.parseDouble(strValue);
				((ConfigEntry<Double>) entry).setValue(value);
				if (setDefaults)
					((ConfigEntry<Double>) entry).setDefValue(value);
			} else if (clazz == String.class) {
				((ConfigEntry<String>) entry).setValue(strValue);
				if (setDefaults)
					((ConfigEntry<String>) entry).setDefValue(strValue);
			}
			if (entry.isOptional())
				entry.setUse(true);
		}

		for (ConfigSection section : sections)
			if (section instanceof ActivitySection)
				((ActivitySection) section).setActivities(activities);
	}

	/**
	 * Gets dialog configuration and default values from resource
	 * superposition.conf. The file is read twice- once to get all variable
	 * settings - the second time to create all the entries.
	 * 
	 */
	public boolean readConfigFile() {

		logger.info("loading config from " + defaultsResource);

		BufferedReader in = new BufferedReader(new InputStreamReader(getClass()
				.getResourceAsStream(defaultsResource)));
		sections = new ArrayList<ConfigSection>();
		allEntries = new HashMap<String, ConfigEntry<?>>();
		ConfigSection currentSection = null;
		boolean inSections = false;
		header = "";

		try {
			String line = in.readLine();
			structureFileEntry = null;

			while (line != null) {
				String origLine = line;
				line = line.trim();
				String next = in.readLine();

				// Look for a dialog definition string
				if (line.startsWith("[") && line.endsWith("]")) {
					inSections = true;
					line = line.substring(1, line.length() - 1);
					String info = null;
					if (next != null && StringUtils.isBlank(next)) {
						next = in.readLine();

						while (next != null && !StringUtils.isBlank(next)
								&& next.startsWith(" ")) {

							if (info != null)
								info += "\n" + next;
							else
								info = next;
							next = in.readLine();

						}
					}

					logger.debug("line " + line + " info " + info);

					int no = line.indexOf(' ');
					String name = line.substring(0, no);
					String rest = line.substring(no + 1);
					no = rest.indexOf(' ');
					String type = null, optionStr = null;
					if (no == -1) {
						type = rest;
					} else {
						type = rest.substring(0, no);
						optionStr = rest.substring(no + 1);
					}
					logger.debug("name '" + name + "' type '" + type + "' optionStr '"
							+ optionStr + "'");

					ConfigEntry<?> entry = null;
					if (type.equals("SECTION")) {
						// section header
						if (name.equals("ACTIVITIES"))
							currentSection = new ActivitySection(name, info);
						else
							currentSection = new ConfigSection(name, info);
						sections.add(currentSection);
					} else if (type.equals("BOOLEAN")) {
						// boolean entry
						entry = new ConfigEntry<Boolean>(name, info, optionStr,
								Boolean.class);
					} else if (type.equals("FILENAME")) {
						// filename entry
						entry = new FilenameConfigEntry(name, info, optionStr);
					} else if (type.equals("STRUCTUREFILE")) {
						// This entry corresponds to a structure file (SDF or
						// MOL2)
						entry = new StructureFileConfigEntry(name, info, optionStr);
						if (structureFileEntry != null)
							logger.error("More that one STRUCTUREFILE entry!");
						structureFileEntry = (StructureFileConfigEntry) entry;
					} else if (type.equals("STRING")) {
						// string entry
						entry = new ConfigEntry<String>(name, info, optionStr,
								String.class);
					} else if (type.equals("INTEGER")) {
						// integer entry
						entry = new ConfigEntry<Integer>(name, info, optionStr,
								Integer.class);
					} else if (type.equals("DOUBLE")) {
						// float number entry
						entry = new ConfigEntry<Double>(name, info, optionStr,
								Double.class);
					} else if (type.equals("MENU")) {
						entry = new MenuConfigEntry(name, info, optionStr);
					} else {
						System.out.println("unknown type " + type);
					}

					if (currentSection == null) {
						System.out.println("Current section not defined");
						return false;
					}
					if (entry != null) {
						String key = entry.name.toUpperCase();
						if (allEntries.containsKey(key))
							System.err.println(key + " already present");
						else {
							allEntries.put(key, entry);
							currentSection.addEntry(entry);
						}
					}
				}

				if (!inSections
						&& (StringUtils.isBlank(origLine) || origLine.startsWith(" ")))
					header += origLine + "\n";

				line = next;
			}
			in.close();

			// Get variable settings and activities
			logger.info("loading defaults and values from " + defaultsResource);
			in = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream(
					defaultsResource)));
			loadValues(in, true);

		} catch (IOException ex) {
			logger.error("IO error: " + ex);
			return false;
		}

		return true;
	}

	/**
	 * Unchecks all optional values
	 * 
	 */
	private void clearOptionalValues() {
		for (ConfigEntry<?> entry : allEntries.values())
			entry.clearOptional();
	}

	/**
	 * Prints variable values from the editor to a file
	 */
	public boolean save(String configFile) {
		Writer out = null;

		logger.debug(saveToString());

		try {
			out = new BufferedWriter(new FileWriter(configFile));
			if (!save(out)) {
				logger.error("Failed to save file " + configFile);
				return false;
			}
			out.close();
		} catch (IOException ex) {
			logger.error("IOException on save to configfile " + configFile, ex);
			return false;
		}
		logger.info("Saved file " + configFile);
		return true;
	}

	public boolean save(Writer out) {
		try {
			out.write(header);

			// for (int i=0; i<sections.size(); i++) {
			// ConfigSection section = (ConfigSection) sections.get(i);
			// section.printInfo(out);
			// }

			out.write("\n");
			out.write("\n");

			for (int i = 0; i < sections.size(); i++) {
				ConfigSection section = sections.get(i);
				// section.printEntries(out);
				section.printEntries(out, true);
			}

		} catch (IOException ex) {
			logger.error("IOException on save ", ex);
			return false;
		}
		return true;
	}

	/**
	 * Saves the config settings to a string.
	 * 
	 * @return
	 */
	public String saveToString() {
		StringWriter writer = new StringWriter();
		if (!save(new PrintWriter(writer))) {
			System.err.println("save error");
		}
		return writer.toString();
	}

	/**
	 * @return the sections
	 */
	public ArrayList<ConfigSection> getSections() {
		return sections;
	}

	/**
	 * @return the allEntries
	 */
	public HashMap<String, ConfigEntry<?>> getAllEntries() {
		return allEntries;
	}

	/**
	 * @return the header
	 */
	public String getHeader() {
		return header;
	}

	/**
	 * @return the structureFileEntry
	 */
	public StructureFileConfigEntry getStructureFileEntry() {
		return structureFileEntry;
	}

	/**
	 * @return the defaultsResource
	 */
	public String getDefaultsResource() {
		return defaultsResource;
	}

	/**
	 * @param defaultsResource
	 *            the defaultsResource to set
	 */
	public void setDefaultsResource(String defaultsResource) {
		if (!defaultsResource.startsWith("/"))
			defaultsResource = "/com/cairn/gape/utils/" + defaultsResource;
		this.defaultsResource = defaultsResource;
	}

	/**
	 * Test routine
	 * 
	 * @param args
	 */
	@SuppressWarnings("unchecked")
	public static void main(String args[]) {

		ConfigModel configModel = new ConfigModel();
		configModel.readConfigFile();

		for (String key : configModel.getAllEntries().keySet()) {
			ConfigEntry<?> entry = configModel.getAllEntries().get(key);
			System.out.println("Key " + key);
			if (key.equalsIgnoreCase("n_runs")) {
				((ConfigEntry<Integer>) entry).setValue(57);
			}
		}

		String info = configModel.saveToString();
		System.out.println("Info " + info);

		configModel.save("/tmp/superposition.conf");
	}
}
