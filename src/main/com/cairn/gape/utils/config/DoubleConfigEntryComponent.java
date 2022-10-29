package com.cairn.gape.utils.config;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JTextField;

import org.apache.commons.lang3.StringUtils;

/**
 * class to hold a floating point number
 * 
 */
class DoubleConfigEntryComponent extends ConfigEntryComponent<Double> {
	private JTextField doubleField;

	public DoubleConfigEntryComponent(ConfigEditor configEditor, ConfigEntry<Double> entry) {
		super(configEditor, entry);
	}

	@Override
	public void buildEditor() {
		doubleField = new JTextField(25);
		doubleField.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent arg0) {
				updateModel();
			}
		});
		updateDialog();
		this.editor = doubleField;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.utils.config.ConfigEntryComponent#updateDialog()
	 */
	@Override
	protected void updateDialog() {
		Double val = configEntry.getValue();
		doubleField.setText(val != null ? String.valueOf(val) : "");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.utils.config.ConfigEntryComponent#updateModel()
	 */
	@Override
	protected boolean updateModel() {
		String val = doubleField.getText().trim();
		if (StringUtils.isEmpty(val)) {
			this.configEditor.errorMessage("No value set for " + configEntry.getName());
			return false;
		}
		try {
			double d = Double.parseDouble(val);
			configEntry.setValue(d);
		} catch (NumberFormatException ex) {
			this.configEditor.errorMessage("Value \"" + val + "\" set for "
					+ configEntry.getName() + " is not a double.");
			return false;
		}
		return true;
	}

}