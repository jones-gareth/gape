package com.cairn.gape.utils.config;

import java.io.File;
import java.util.ArrayList;

import javax.swing.filechooser.FileFilter;

import org.apache.commons.lang3.StringUtils;

/**
 * This is an extension of the filename class to handle files containing
 * structures.
 * 
 */
class StructureFileEntryComponent extends FilenameConfigEntryComponent {

	private final ArrayList<StructureFileListener> listeners = new ArrayList<StructureFileListener>();

	String[] molNames;

	public StructureFileEntryComponent(ConfigEditor configEditor,
			ConfigEntry<String> entry) {
		super(configEditor, entry);
	}

	@Override
	public void setChecked(boolean v) {
		super.setChecked(v);
		if (v == true && StringUtils.isNotEmpty(configEntry.getValue())) {
			ConfigEditor.logger.debug("getting names");
			molNames = ConfigEditor.getMoleculeNames(configEntry.getValue());
		}
		if (v == false)
			molNames = null;
		tellListeners();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.utils.config.ConfigEditor.FilenameConfigEntry
	 * #updateDialog()
	 */
	@Override
	protected void updateDialog() {
		super.updateDialog();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.utils.config.ConfigEditor.FilenameConfigEntry
	 * #updateModel()
	 */
	@Override
	protected boolean updateModel() {
		return super.updateModel();
	}

	@Override
	public File pickFile(String title) {
		File f = super.pickFile(title, new StructureFileFilter());
		if (f == null)
			return null;
		molNames = ConfigEditor.getMoleculeNames(f.getPath());
		tellListeners();
		return f;
	}

	/**
	 * We can add listeners to this dialog, so that other components can be
	 * informed of structures in a file
	 * 
	 * @param listener
	 */
	public void addListener(StructureFileListener listener) {
		if (!listeners.contains(listener))
			listeners.add(listener);
	}

	/**
	 * Tell the listeners when we load a structure file and pass them the list
	 * of molecule names loaded
	 * 
	 */
	private void tellListeners() {
		for (int i = 0; i < listeners.size(); i++) {
			StructureFileListener listener = listeners.get(i);
			listener.tellListenerNames(molNames);
		}
	}

	/**
	 * Flilename filter for structure files.
	 * 
	 */
	class StructureFileFilter extends FileFilter {

		@Override
		public boolean accept(File f) {
			if (f.isDirectory())
				return true;
			String name = f.getName();
			if (name.toLowerCase().endsWith(".mol2"))
				return true;
			if (name.toLowerCase().endsWith(".mol"))
				return true;
			if (name.toLowerCase().endsWith(".sdf"))
				return true;
			return false;
		}

		@Override
		public String getDescription() {
			return "Structure (.mol2 .mol .sdf) files";
		}

	}

	/**
	 * Implement this to be informed of structures in a molecule file.
	 * 
	 */
	interface StructureFileListener {
		void tellListenerNames(String names[]);
	}

}