package com.cairn.gape.utils.config;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JEditorPane;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.UIManager;
import javax.swing.border.EtchedBorder;
import javax.swing.filechooser.FileFilter;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import com.cairn.common.utils.BasicWindowMonitor;
import com.cairn.gape.molecule.GaMolecule;

/**
 * This class is used to build the GAPE configuration file editor. Also see
 * documentation and comments in the example config files.
 * 
 * @author Gareth
 * 
 */
public class ConfigEditor extends JPanel {
	private File configFile;

	ConfigModel configModel;

	private ArrayList<ConfigSectionComponent> configSectionComponents;

	private HashMap<String, ConfigEntryComponent<?>> allConfigEntryComponents;

	StructureFileEntryComponent structureFileEntryComponent;

	private JTabbedPane tabbedPane;

	private Controls controls;

	private JTextArea infoPanel;

	private String messages;

	private static final int HELP_WIDTH = 50;

	static Logger logger;
	static {
		logger = Logger.getLogger(ConfigEditor.class);
		logger.setLevel(Level.DEBUG);
	}

	/**
	 * Simple constructor
	 * 
	 */
	public ConfigEditor() {
		common();
		init();
	}

	/**
	 * This constructor loads in a file
	 * 
	 * @param f
	 *            File object
	 */
	public ConfigEditor(File f) {
		common();
		init(f);
	}

	/**
	 * Given the name of a field returns the ConfigEntry for that field
	 * 
	 * @param name
	 * @return
	 */
	public ConfigEntryComponent<?> getEntry(String name) {
		return allConfigEntryComponents.get(name.toUpperCase());
	}

	/**
	 * Sets L&F
	 * 
	 */
	protected void common() {
		// On Linux I've had stack overflows in GTK look and feel
		// and also lots of initially blank applets.
		if (System.getProperty("os.name").equals("Linux")) {
			return;
		}

		try {
			// Set System L&F
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		} catch (Exception e) {
			System.err.println("Unable to set look and feel");
		}
	}

	/**
	 * Set up editor: loads defaults and displays
	 * 
	 */
	private void init() {
		configModel = new ConfigModel();
		configModel.readConfigFile();
		display();
		removeFile();
	}

	/**
	 * As init(), but defaults are loaded from a file.
	 * 
	 * @param f
	 */
	private void init(File f) {
		configModel = new ConfigModel();
		configModel.readConfigFile();
		configModel.loadFile(f);
		display();
		setFile(f);
	}

	/**
	 * Build the GUI
	 * 
	 */
	private void display() {
		setLayout(new BorderLayout());

		controls = new Controls();
		add(BorderLayout.NORTH, new JScrollPane(controls));
		buildUI();
		infoPanel = new JTextArea(3, 0);
		infoPanel.setEnabled(false);
		infoPanel.setText(messages);
		add(BorderLayout.SOUTH, new JScrollPane(infoPanel));
	}

	/**
	 * @return the configModel
	 */
	public ConfigModel getConfigModel() {
		return configModel;
	}

	/**
	 * Main method
	 * 
	 * @param args
	 */
	public static void main(String args[]) {
		// Properties p = new Properties(System.getProperties());
		// p.put("swing.metalTheme", "steel");
		// System.setProperties(p);

		ConfigEditor cf = null;
		if (args.length == 1)
			cf = new ConfigEditor(new File(args[0]));
		else if (args.length == 0) {
			cf = new ConfigEditor();
		} else {
			System.out.println("usage: ConfigEditor [config file]");
		}

		JFrame f = new JFrame("Configuration file editor");
		f.addWindowListener(new BasicWindowMonitor());
		f.setSize(950, 600);
		f.getContentPane().add(cf);
		f.setVisible(true);
	}

	/**
	 * Reset all values to defaults
	 */
	private void resetToDefaults() {

		for (ConfigEntryComponent<?> entryComponent : allConfigEntryComponents.values()) {
			ConfigEntry<?> entry = entryComponent.getConfigEntry();
			entry.resetValueToDefault();
			if (entry.isOptional())
				entry.setUse(false);
			entryComponent.updateDialog();
		}

	}

	/**
	 * Reset all values to model values
	 */
	private void resetToModelValues() {

		for (ConfigEntryComponent<?> entryComponent : allConfigEntryComponents.values()) {
			entryComponent.updateDialog();
		}

	}

	/**
	 * Makes sure the model is up to date - this may be necessary for textfields
	 * 
	 * @return
	 */
	private boolean dialogToModel() {
		for (ConfigEntryComponent<?> entryComponent : allConfigEntryComponents.values()) {
			ConfigEntry<?> entry = entryComponent.getConfigEntry();

			if (entry.isUse() && !entryComponent.updateModel()) {
				String name = entry.getName();
				logger.error("failed to update model for component " + name);
				errorMessage("Unable to save value for " + name);
				return false;
			}
		}
		return true;
	}

	/**
	 * Resets editor to defaults
	 * 
	 */
	private void newFile() {
		setDefaults();
		removeFile();
	}

	/**
	 * Resets editor to defaults- but keeps filename
	 * 
	 */
	private void setDefaults() {
		resetToDefaults();
		revalidate();
		// buildUI();
	}

	/**
	 * Sets the filename in the editor.
	 * 
	 * @param f
	 */
	private void setFile(File f) {
		configFile = f;
		// controls.fileField.setText(configFile.getAbsolutePath());
		controls.fileField.setText(configFile.getPath());
		controls.save.setEnabled(true);
	}

	private void save() {
		if (!dialogToModel())
			return;
		if (configFile == null) {
			errorMessage("No Filename chosen!");
			return;
		}
		configModel.save(configFile.getPath());
		infoMessage("Saved to " + configFile.getPath());
	}

	/**
	 * Clears the filename field in the editor
	 * 
	 */
	private void removeFile() {
		configFile = null;
		controls.fileField.setText("");
		controls.save.setEnabled(false);
	}

	/**
	 * Loads up settings from an existing file.
	 * 
	 */
	private void loadFile() {
		JFileChooser chooser = new JFileChooser();
		chooser.setCurrentDirectory(new File("."));
		chooser.setDialogTitle("Open Config File");
		chooser.setFileFilter(new ConfigFileFilter());
		int returnVal = chooser.showOpenDialog(this);

		File f = null;
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			f = chooser.getSelectedFile();
			// sections = null;
			if (tabbedPane != null)
				remove(tabbedPane);
			tabbedPane = null;
			configModel.loadFile(f);
			resetToModelValues();
			revalidate();
			setFile(f);
		} else {
			return;
		}
	}

	/**
	 * Picks a file then saves settings to it.
	 * 
	 */
	private void saveAs() {
		if (!dialogToModel())
			return;

		JFileChooser chooser = new JFileChooser();
		chooser.setCurrentDirectory(new File("."));
		chooser.setDialogTitle("Save Config File");
		chooser.setFileFilter(new ConfigFileFilter());
		int returnVal = chooser.showSaveDialog(this);

		File f = null;
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			f = chooser.getSelectedFile();
			setFile(f);
			configModel.save(f.getPath());
			infoMessage("Saved to " + configFile.getPath());
		} else {
			return;
		}
	}

	/**
	 * Filter for JFileChooser filename filter
	 */
	class ConfigFileFilter extends FileFilter {

		@Override
		public boolean accept(File f) {
			if (f.isDirectory())
				return true;
			String name = f.getName();
			if (name.toLowerCase().endsWith(".conf"))
				return true;
			return false;
		}

		@Override
		public String getDescription() {
			return "Config (.conf) files";
		}

	}

	/**
	 * Adds an informational message to the lower panel.
	 * 
	 * @param msg
	 */
	private void updateInfoPanel(String msg) {
		if (infoPanel != null)
			infoPanel.append("\n" + msg);
		if (messages != null)
			messages += "\nInfo: " + msg;
		else
			messages = msg;
	}

	/**
	 * Display an error message
	 * 
	 * @param msg
	 */
	void errorMessage(String msg) {
		JOptionPane.showMessageDialog(null, msg, "Error Message",
				JOptionPane.ERROR_MESSAGE);
		updateInfoPanel("Error: " + msg);
		System.out.println("Error: " + msg);

	}

	/**
	 * Display an informational message
	 * 
	 * @param msg
	 */
	private void infoMessage(String msg) {
		// JOptionPane.showMessageDialog
		// (null, msg,
		// "Information Message", JOptionPane.INFORMATION_MESSAGE);
		updateInfoPanel("Info: " + msg);
		System.out.println("Info: " + msg);
	}

	/**
	 * Creates a component for displaying help about a variable.
	 * 
	 * @param info
	 * @return
	 */
	protected JComponent infoText(String info) {
		JTextArea infoTa = new JTextArea(info.trim());
		infoTa.setOpaque(false);
		infoTa.setEnabled(false);
		infoTa.setLineWrap(true);
		infoTa.setWrapStyleWord(true);
		infoTa.setForeground(Color.BLACK);
		infoTa.setDisabledTextColor(Color.DARK_GRAY);
		infoTa.setColumns(HELP_WIDTH);
		return infoTa;
	}

	/**
	 * This class contains the main controls for the editor
	 * 
	 */
	private class Controls extends JPanel {
		private final JButton newB, load, saveAs, save, exit, help, defaults;

		private final JTextField fileField;

		private final JComboBox<String> defaultsBox;

		private Controls() {
			setLayout(new FlowLayout(FlowLayout.LEFT));

			newB = new JButton("New");
			newB.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					newFile();
				}
			});
			add(newB);

			defaults = new JButton("Default");
			defaults.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					setDefaults();
				}
			});
			add(defaults);

			defaultsBox = new JComboBox<>(new String[] { "GAPE", "PharmSearch", "GRIPS" });
			defaultsBox.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					String setting = (String) defaultsBox.getSelectedItem();
					System.out.println("Got " + setting);
					if (setting.equals("GAPE"))
						configModel.setDefaultsResource("superposition.conf");
					else if (setting.equals("PharmSearch"))
						configModel.setDefaultsResource("pharm_search.conf");
					else if (setting.equals("GRIPS"))
						configModel.setDefaultsResource("rigid_search.conf");
					configModel.readConfigFile();
					buildUI();
				}
			});
			add(defaultsBox);

			load = new JButton("Load");
			load.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					loadFile();
				}
			});
			add(load);

			save = new JButton("Save");
			save.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					save();
				}
			});
			add(save);

			saveAs = new JButton("SaveAs");
			saveAs.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					saveAs();
				}
			});
			add(saveAs);

			exit = new JButton("Exit");
			exit.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					System.exit(1);
				}
			});
			add(exit);

			fileField = new JTextField(20);
			add(fileField);
			fileField.setEnabled(false);
			if (configFile != null)
				// fileField.setText(configFile.getAbsolutePath());
				fileField.setText(configFile.getPath());
			JLabel l = new JLabel("File");
			add(l);

			help = new JButton("Help");
			help.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					new HelpWin("gape.html");
				}
			});
			add(help);

			setBorder(new EtchedBorder());
		}
	}

	/**
	 * Gets structure names from a molecule file
	 * 
	 * @param molFile
	 * @return
	 */
	static String[] getMoleculeNames(String molFile) {
		List<GaMolecule> molecules = null;

		try {
			molecules = GaMolecule.loadFiles(new String[] { molFile },
					GaMolecule.FileType.UNKNOWN, GaMolecule.Source.FILE);
		} catch (Exception ex) {
			logger.error("GaException: " + ex);
			return null;
		}

		List<String> names = molecules.stream().map(GaMolecule::getName)
				.collect(Collectors.toList());
		return names.toArray(new String[names.size()]);
	}

	/**
	 * @return the configFile
	 */
	public File getConfigFile() {
		return configFile;
	}

	/**
	 * Construct the user interface from the model.
	 */
	@SuppressWarnings("unchecked")
	private void buildUI() {
		configSectionComponents = new ArrayList<ConfigSectionComponent>();
		allConfigEntryComponents = new HashMap<String, ConfigEntryComponent<?>>();

		// first find any structure entry field
		if (configModel.getStructureFileEntry() != null)
			structureFileEntryComponent = new StructureFileEntryComponent(this,
					configModel.getStructureFileEntry());

		for (ConfigSection configSection : configModel.getSections()) {
			ConfigSectionComponent configSectionComponent = null;
			if (configSection.getName().equals("ACTIVITIES"))
				configSectionComponent = new ActivitySectionComponent(this, configSection);
			else
				configSectionComponent = new ConfigSectionComponent(this, configSection);

			for (ConfigEntry<?> entry : configSection.getEntries()) {

				Class<?> clazz = entry.getClazz();

				ConfigEntryComponent<?> component = null;
				if (entry instanceof StructureFileConfigEntry) {
					component = structureFileEntryComponent;
				} else if (entry instanceof FilenameConfigEntry) {
					component = new FilenameConfigEntryComponent(this,
							(ConfigEntry<String>) entry);
				} else if (entry instanceof MenuConfigEntry) {
					component = new MenuConfigEntryComponent(this,
							(ConfigEntry<String>) entry);
				} else if (clazz == Boolean.class) {
					component = new BooleanConfigEntryComponent(this,
							(ConfigEntry<Boolean>) entry);
				} else if (clazz == Integer.class) {
					component = new IntegerConfigEntryComponent(this,
							(ConfigEntry<Integer>) entry);
				} else if (clazz == Double.class) {
					component = new DoubleConfigEntryComponent(this,
							(ConfigEntry<Double>) entry);
				} else if (clazz == String.class) {
					component = new StringConfigEntryComponent(this,
							(ConfigEntry<String>) entry);
				}

				configSectionComponent.addConfigEntryComponent(component);
				allConfigEntryComponents.put(entry.getName(), component);
			}

			configSectionComponents.add(configSectionComponent);
		}

		if (tabbedPane != null)
			remove(tabbedPane);

		tabbedPane = new JTabbedPane();
		add(BorderLayout.CENTER, tabbedPane);

		for (int i = 0; i < configSectionComponents.size(); i++) {
			ConfigSectionComponent section = configSectionComponents.get(i);
			section.display();
			tabbedPane.addTab(section.getConfigSection().getName().replace('_', ' '),
					new JScrollPane(section));
		}

	}
}

class HelpWin extends JFrame {

	protected HelpWin(String resource) {
		setBounds(10, 10, 500, 600);
		// JEditorPane pane = new JEditorPane("text/html", html);
		JEditorPane pane = null;
		try {
			java.net.URL file = getClass().getResource(resource);
			pane = new JEditorPane(file);
		} catch (java.io.IOException ex) {
			System.err.println("Failed to open " + resource);
		}
		pane.setEditable(false);
		JScrollPane sp = new JScrollPane(pane);
		JPanel p = new JPanel(new FlowLayout(FlowLayout.LEFT));
		JButton b = new JButton("OK");
		b.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				dispose();
			}
		});
		p.add(b);
		getContentPane().setLayout(new BorderLayout());
		getContentPane().add(sp, BorderLayout.CENTER);
		getContentPane().add(p, BorderLayout.SOUTH);
		setVisible(true);
	}
}
