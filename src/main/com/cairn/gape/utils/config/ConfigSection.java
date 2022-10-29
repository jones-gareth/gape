package com.cairn.gape.utils.config;

import java.io.Writer;
import java.util.ArrayList;

import org.apache.log4j.Logger;

public class ConfigSection {
	protected static Logger logger;
	static {
		logger = Logger.getLogger(ConfigSection.class);
		// logger.setLevel(Level.DEBUG);
	}
	protected String name;
	private final String description;
	protected ArrayList<ConfigEntry<?>> entries = new ArrayList<ConfigEntry<?>>();

	/**
	 * Constructor
	 * 
	 * @param n
	 *            section name
	 * @param d
	 *            description
	 */
	protected ConfigSection(String n, String d) {
		name = n;
		description = d;
	}

	/**
	 * Add a variable entry to this section
	 * 
	 * @param e
	 */
	public void addEntry(ConfigEntry<?> e) {
		entries.add(e);
	}

	public void printEntries(Writer out) {
		printEntries(out, false);
	}

	/**
	 * Prints description and values to a file.
	 * 
	 * @param out
	 * @param showDescription
	 */
	public void printEntries(Writer out, boolean showDescription) {
		try {
			out.write("\n");
			out.write(" [" + name + " SECTION]\n");
			out.write("\n");
			if (showDescription && description != null) {
				out.write(description);
				out.write("\n");
				out.write("\n");
			}

			for (int i = 0; i < entries.size(); i++) {
				final ConfigEntry<?> entry = entries.get(i);
				entry.printInfo(out);
				if (showDescription)
					entry.printHelp(out);
				out.write("\n");
				if (entry.isOptional() && !entry.isUse()) {
					if (showDescription)
						out.write(" Optional parameter " + entry.name.toLowerCase()
								+ " not set\n\n");
					continue;
				}
				if (entry.isUse())
					entry.printValue(out);
				out.write("\n");
			}
		} catch (Exception e) {
			logger.error("IO exception writing info ", e);
		}
	}

	/**
	 * uncheck all optional values in the section
	 * 
	 */
	public void clearOptionalValues() {
		for (int i = 0; i < entries.size(); i++) {
			final ConfigEntry<?> entry = entries.get(i);
			entry.clearOptional();
		}
	}

	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	/**
	 * @param name
	 *            the name to set
	 */
	public void setName(String name) {
		this.name = name;
	}

	/**
	 * @return the entries
	 */
	public ArrayList<ConfigEntry<?>> getEntries() {
		return entries;
	}

	/**
	 * @param entries
	 *            the entries to set
	 */
	public void setEntries(ArrayList<ConfigEntry<?>> entries) {
		this.entries = entries;
	}

	/**
	 * @return the description
	 */
	public String getDescription() {
		return description;
	}

}
