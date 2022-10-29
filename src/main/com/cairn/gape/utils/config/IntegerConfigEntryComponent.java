package com.cairn.gape.utils.config;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JTextField;

import org.apache.commons.lang3.StringUtils;

/**
 * class to hold an integer variable
 * 
 */
class IntegerConfigEntryComponent extends ConfigEntryComponent<Integer> {
	private JTextField integerField;

	public IntegerConfigEntryComponent(ConfigEditor configEditor,
			ConfigEntry<Integer> entry) {
		super(configEditor, entry);
	}

	@Override
	public void buildEditor() {
		integerField = new JTextField(25);
		integerField.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent arg0) {
				updateModel();
			}
		});
		updateDialog();
		this.editor = integerField;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.utils.config.ConfigEntryComponent#updateDialog()
	 */
	@Override
	protected void updateDialog() {
		Integer val = configEntry.getValue();
		integerField.setText(val != null ? String.valueOf(val) : "");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.utils.config.ConfigEntryComponent#updateModel()
	 */
	@Override
	protected boolean updateModel() {
		String val = integerField.getText().trim();
		if (StringUtils.isEmpty(val)) {
			this.configEditor.errorMessage("No value set for " + configEntry.getName());
			return false;
		}
		try {
			int d = Integer.parseInt(val);
			configEntry.setValue(d);
		} catch (NumberFormatException ex) {
			this.configEditor.errorMessage("Value \"" + val + "\" set for "
					+ configEntry.getName() + " is not an integer.");
			return false;
		}
		return true;
	}

}