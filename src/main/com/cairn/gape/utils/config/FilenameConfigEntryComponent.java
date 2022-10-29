package com.cairn.gape.utils.config;

import java.awt.GridBagConstraints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.filechooser.FileFilter;

import com.cairn.common.utils.UtilLayout;

/**
 * This class is used to edit filenames
 * 
 */
class FilenameConfigEntryComponent extends ConfigEntryComponent<String> {

	private JTextField fileField;

	private JButton button;

	public FilenameConfigEntryComponent(ConfigEditor configEditor,
			ConfigEntry<String> entry) {
		super(configEditor, entry);
	}

	@Override
	public void buildEditor() {
		fileField = new JTextField(18);
		fileField.setEnabled(false);
		updateDialog();
		button = new JButton("Pick");
		button.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				File f = pickFile("Pick File for "
						+ configEntry.getName().replace('_', ' '));
				if (f != null) {
					String filename = f.getPath();
					fileField.setText(filename);
					if (updateModel()) {
						fileField.setText(filename);
					}
				}
			}
		});
		JPanel p = new JPanel();
		UtilLayout layout = new UtilLayout(p);
		layout.c.fill = GridBagConstraints.BOTH;
		layout.add(fileField, 0, 0);
		layout.c.fill = GridBagConstraints.NONE;
		layout.add(button, 1, 0);
		this.editor = p;
	}

	@Override
	public void setChecked(boolean v) {
		configEntry.setUse(v);
		fileField.setEnabled(v);
		button.setEnabled(v);
	}

	protected File pickFile(String title) {
		return pickFile(title, null);
	}

	protected File pickFile(String title, FileFilter f) {
		File dir = new File(".");
		JFileChooser chooser = new JFileChooser();
		chooser.setCurrentDirectory(dir);
		chooser.setDialogTitle(title);
		if (f != null)
			chooser.setFileFilter(f);
		int returnVal = chooser.showOpenDialog(this);
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			return chooser.getSelectedFile();
		}
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.utils.config.ConfigEntryComponent#updateDialog()
	 */
	@Override
	protected void updateDialog() {
		fileField.setText(configEntry.getValue());
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.utils.config.ConfigEntryComponent#updateModel()
	 */
	@Override
	protected boolean updateModel() {

		String fName = fileField.getText();
		if (fName == null || fName.equals("")) {
			this.configEditor.errorMessage("Filename for "
					+ configEntry.getName() + " not set");
			return false;
		}
		File f = new File(fName);
		if (!f.exists()) {
			this.configEditor.errorMessage("Filename " + f + " does not exist");
			return false;
		}
		if (!f.canRead()) {
			this.configEditor.errorMessage("Filename " + f + " not readable");
			return false;
		}

		configEntry.setValue(fName);
		return true;
	}

}