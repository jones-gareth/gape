package com.cairn.gape.utils.config;

/**
 * Models a structure file config entry
 * 
 * @author gjones
 * 
 */
public class StructureFileConfigEntry extends ConfigEntry<String> {

	public StructureFileConfigEntry(String name, String help,
			String optionsString) {
		super(name, help, optionsString, String.class);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.utils.ConfigEntry#getType()
	 */
	@Override
	public String getType() {
		return "STRUCTUREFILE";
	}

}
