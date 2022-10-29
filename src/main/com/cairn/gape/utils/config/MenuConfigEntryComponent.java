package com.cairn.gape.utils.config;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JComboBox;

import org.apache.commons.lang3.StringUtils;

/**
 * This class is used to represent a menu of items.
 * 
 */
class MenuConfigEntryComponent extends ConfigEntryComponent<String> implements
		StructureFileEntryComponent.StructureFileListener {
	private JComboBox<String> comboBox;

	// private boolean molName;

	public MenuConfigEntryComponent(ConfigEditor configEditor, ConfigEntry<String> entry) {
		super(configEditor, entry);
	}

	@Override
	public void buildEditor() {
		MenuConfigEntry entry = (MenuConfigEntry) configEntry;

		String[] items = null;
		StructureFileConfigEntry structureFileEntry = this.configEditor.configModel
				.getStructureFileEntry();
		if (entry.isMolName() && structureFileEntry != null) {
			if (StringUtils.isNotBlank(structureFileEntry.getValue()))
				items = ConfigEditor.getMoleculeNames(structureFileEntry.getValue());
			this.configEditor.structureFileEntryComponent.addListener(this);
		} else {
			items = ((MenuConfigEntry) configEntry).getItems();
		}

		if (items != null)
			comboBox = new JComboBox<String>(items);
		else
			comboBox = new JComboBox<String>();

		comboBox.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent arg0) {
				updateModel();

			}
		});
		updateDialog();
		this.editor = comboBox;
	}

	/**
	 * We can tie this menu to the list of molecules in a structure file.
	 */
	@Override
	public void tellListenerNames(String molNames[]) {
		comboBox.removeAllItems();
		boolean found = false;
		String value = configEntry.getValue();
		if (molNames != null)
			for (int i = 0; i < molNames.length; i++) {
				comboBox.addItem(molNames[i]);
				if (value != null && value.equalsIgnoreCase(molNames[i]))
					found = true;
			}
		if (found)
			comboBox.setSelectedItem(value);
		else {
			configEntry.setValue(null);
			updateDialog();
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.utils.config.ConfigEntryComponent#updateDialog()
	 */
	@Override
	protected void updateDialog() {
		String value = configEntry.getValue();
		if (value != null) {
			comboBox.setSelectedItem(value.toLowerCase());
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.utils.config.ConfigEntryComponent#updateModel()
	 */
	@Override
	protected boolean updateModel() {
		Object item = comboBox.getSelectedItem();
		if (item == null)
			return false;
		String value = ((String) item).toLowerCase();
		configEntry.setValue(value);
		return true;
	}

}