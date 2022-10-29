package com.cairn.gape.utils.config;

import javax.swing.JComponent;
import javax.swing.JPanel;

/**
 * All variable dialogs are sub-classes of this
 * 
 */
public abstract class ConfigEntryComponent<V> extends JPanel {
	protected ConfigEditor configEditor;
	protected ConfigEntry<V> configEntry;
	protected ConfigEntryComponent<?> otherHelp;
	protected JComponent editor;

	/*
	 * private static Logger logger; static { Utils.configureLogger(); logger =
	 * Logger.getLogger(ConfigEntryComponent.class); //
	 * logger.setLevel(Level.DEBUG); }
	 */
	/**
	 * Create the entry
	 * 
	 * @param n
	 * @param o
	 * @param h
	 */
	ConfigEntryComponent(ConfigEditor _editor, ConfigEntry<V> entry) {
		configEditor = _editor;
		this.configEntry = entry;
		buildEditor();
	}

	/**
	 * set value from defaults.
	 * 
	 */
	public void setDefault() {
		configEntry.setValue(configEntry.getDefValue());
	}

	/**
	 * Update the UI with a changed value from the model
	 */
	protected abstract void updateDialog();

	/**
	 * Update the Model with a changed value from the UI
	 */
	protected abstract boolean updateModel();

	public String getHelp() {
		return configEntry.getHelp();
	}

	public JComponent getHelpComponent() {
		JComponent helpTa = configEditor.infoText(configEntry.getHelp()
				.replaceAll("\n", ""));
		return helpTa;
	}

	/**
	 * Removes value if variable is optional
	 * 
	 */
	public void clearOptional() {
		if (!configEntry.isOptional())
			return;

		configEntry.clearOptional();
		setChecked(false);
	}

	public boolean isOptional() {
		return configEntry.isOptional();
	}

	public V getValue() {
		return configEntry.getValue();
	}

	/**
	 * Build a component for editing/setting the variable
	 * 
	 * @return
	 */
	public abstract void buildEditor();

	/**
	 * Used to enable editor for an optional variable
	 * 
	 * @param c
	 */
	public void setChecked(boolean checked) {
		configEntry.setUse(checked);
		editor.setEnabled(checked);
	}

	/**
	 * @return the configEntry
	 */
	public ConfigEntry<V> getConfigEntry() {
		return configEntry;
	}

	/**
	 * @return the editor
	 */
	public JComponent getEditor() {
		return editor;
	}

}
