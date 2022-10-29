package com.cairn.gape.utils.config;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JCheckBox;

import org.apache.commons.lang3.BooleanUtils;

/**
 * This class is used to edit a boolean value
 * 
 */
class BooleanConfigEntryComponent extends ConfigEntryComponent<Boolean> {

	public BooleanConfigEntryComponent(ConfigEditor configEditor,
			ConfigEntry<Boolean> configEntry) {
		super(configEditor, configEntry);
	}

	@Override
	public void buildEditor() {
		JCheckBox editor = new JCheckBox();
		editor.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent arg0) {
				updateModel();
			}
		});
		this.editor = editor;
		updateDialog();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.utils.config.ConfigEntryComponent#updateDialog()
	 */
	@Override
	protected void updateDialog() {
		// the activities editor may not be present
		if (editor == null)
			return;
		((JCheckBox) editor).setSelected(BooleanUtils.isTrue(configEntry.getValue()));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.utils.config.ConfigEntryComponent#updateModel()
	 */
	@Override
	protected boolean updateModel() {
		boolean value = ((JCheckBox) editor).isSelected();
		configEntry.setValue(value);
		return true;
	}

}