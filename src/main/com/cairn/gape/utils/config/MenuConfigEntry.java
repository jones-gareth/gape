package com.cairn.gape.utils.config;

/**
 * This class is used to represent a menu of items.
 * 
 * @author gjones
 * 
 */
public class MenuConfigEntry extends ConfigEntry<String> {
	private boolean molName = false;
	String items[];

	public MenuConfigEntry(String name, String help, String optionsString) {
		super(name, help, optionsString, String.class);
		optionStrToMenuItems(optionsString);
	}

	private void optionStrToMenuItems(String optionStr) {
		if (optionStr.indexOf("MOLNAME") != -1) {
			molName = true;
		} else {
			int start = optionStr.indexOf('(') + 1;
			int end = optionStr.indexOf(')');
			String menu = optionStr.substring(start, end).toLowerCase();
			items = menu.split("[, ]+");
		}
	}

	/**
	 * @return the molName
	 */
	public boolean isMolName() {
		return molName;
	}

	/**
	 * @return the items
	 */
	public String[] getItems() {
		return items;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.utils.ConfigEntry#getType()
	 */
	@Override
	public String getType() {
		return "MENU";
	}
}
