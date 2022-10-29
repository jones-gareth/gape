package com.cairn.gape.utils.config;

/**
 * Models a filename config entry
 * 
 * @author gjones
 */
public class FilenameConfigEntry extends ConfigEntry<String> {

	public FilenameConfigEntry(String name, String help, String optionsString) {
		super(name, help, optionsString, String.class);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.utils.ConfigEntry#getType()
	 */
	@Override
	public String getType() {
		return "FILENAME";
	}

}
