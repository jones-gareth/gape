package com.cairn.gape.utils.config;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JTextField;

import org.apache.commons.lang3.StringUtils;

/**
 * Class to hold a string value
 */
class StringConfigEntryComponent extends ConfigEntryComponent<String> {
	private JTextField stringField;

	public StringConfigEntryComponent(ConfigEditor configEditor, ConfigEntry<String> entry) {
		super(configEditor, entry);
	}

	@Override
	public void buildEditor() {
		stringField = new JTextField(25);
		stringField.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent arg0) {
				updateModel();
			}
		});
		updateDialog();
		this.editor = stringField;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.utils.config.ConfigEntryComponent#updateDialog()
	 */
	@Override
	protected void updateDialog() {
		String val = configEntry.getValue();
		stringField.setText(val != null ? val : "");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.utils.config.ConfigEntryComponent#updateModel()
	 */
	@Override
	protected boolean updateModel() {
		String val = stringField.getText().trim();
		if (StringUtils.isBlank(val)) {
			this.configEditor.errorMessage("No value set for " + configEntry.getName());
			return false;
		}
		configEntry.setValue(val);
		return true;
	}

}