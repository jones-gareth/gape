package com.cairn.gape.utils.config;

import java.io.IOException;
import java.io.Writer;
import java.util.HashSet;
import java.util.Set;

import org.apache.log4j.Logger;

/**
 * Models an Entry in a configuration file
 * 
 * @author gjones
 * 
 * @param <V>
 */
public class ConfigEntry<V> {
	protected static Logger logger;
	static {
		logger = Logger.getLogger(ConfigEntry.class);
		// logger.setLevel(Level.DEBUG);
	}

	public enum ConfigEntryOptions {
		OPTIONAL
	};

	protected String help, name;
	protected V value, defValue;
	protected Class<V> clazz;
	protected Set<ConfigEntryOptions> options;
	protected boolean use = true;
	protected String optionsString;

	public ConfigEntry(String name, String help, String optionsStr, Class<V> clazz) {
		super();
		this.name = name;
		this.help = help;
		this.optionsString = optionsStr;
		this.options = stringToOptions(optionsStr);
		this.clazz = clazz;
	}

	/**
	 * Parses a string for config entry options.
	 * 
	 * @param optionsStr
	 * @return
	 */
	private static Set<ConfigEntryOptions> stringToOptions(String optionsStr) {
		Set<ConfigEntryOptions> options = new HashSet<ConfigEntryOptions>();
		if (optionsStr == null)
			return options;
		if (optionsStr.indexOf("OPTIONAL") != -1)
			options.add(ConfigEntryOptions.OPTIONAL);
		return options;
	}

	/**
	 * Sets an optional value to not being used
	 */
	public void clearOptional() {
		if (options.contains(ConfigEntryOptions.OPTIONAL))
			use = false;
	}

	/**
	 * Prints out information about this variable
	 * 
	 * @param out
	 */
	public void printInfo(Writer out) {
		try {
			String opt = optionsString != null ? " " + optionsString : "";
			out.write(" [" + name + " " + getType() + opt + "]\n");
		} catch (IOException e) {
			logger.error("IO exception writing info ", e);
		}
	}

	/**
	 * Prints out the variable value
	 * 
	 * @param out
	 */
	public void printValue(Writer out) {
		if (use) {
			try {
				String strValue = value == null ? "" : value.toString();
				if (value instanceof Boolean)
					strValue = (Boolean) value ? "yes" : "no";
				out.write(name.toLowerCase() + " = " + strValue + "\n");
			} catch (IOException e) {
				logger.error("IO exception writing value ", e);
			}
		}
	}

	/**
	 * Return the type name used in the config file
	 * 
	 * @return
	 */
	public String getType() {
		String className = clazz.getCanonicalName();
		String base = className.substring(className.lastIndexOf('.') + 1);
		return base.toUpperCase();
	}

	/**
	 * @return true if this entry is optional
	 */
	public boolean isOptional() {
		return options.contains(ConfigEntryOptions.OPTIONAL);
	}

	/**
	 * Print help string to file
	 * 
	 * @param out
	 */
	public void printHelp(Writer out) {
		if (help == null)
			return;
		try {
			out.write("\n");
			out.write(help);
			out.write("\n");
		} catch (IOException e) {
			logger.error("IO exception writing entry help", e);
		}
	}

	/**
	 * Reset value to default value
	 */
	public void resetValueToDefault() {
		value = defValue;
	}

	/**
	 * @return the value
	 */
	public V getValue() {
		return value;
	}

	/**
	 * @param value
	 *            the value to set
	 */
	public void setValue(V value) {
		this.value = value;
	}

	/**
	 * @return the help
	 */
	public String getHelp() {
		return help;
	}

	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	/**
	 * @return the clazz
	 */
	public Class<V> getClazz() {
		return clazz;
	}

	/**
	 * @return the options
	 */
	public Set<ConfigEntryOptions> getOptions() {
		return options;
	}

	/**
	 * @return the defValue
	 */
	public V getDefValue() {
		return defValue;
	}

	/**
	 * @param defValue
	 *            the defValue to set
	 */
	public void setDefValue(V defValue) {
		this.defValue = defValue;
	}

	/**
	 * @return the use
	 */
	public boolean isUse() {
		return use;
	}

	/**
	 * @param use
	 *            the use to set
	 */
	public void setUse(boolean use) {
		this.use = use;
	}

}
