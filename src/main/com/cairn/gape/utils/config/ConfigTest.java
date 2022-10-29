package com.cairn.gape.utils.config;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.border.EtchedBorder;

import com.cairn.common.utils.BasicWindowMonitor;
import com.cairn.common.utils.UtilLayout;

/**
 * Test class for bulding limited Configuration GUI Given and editor and array
 * of fields displays an editor for only those fields.
 * 
 * Not sure that this works any more after refactoring to MVC
 * 
 * @author Gareth Jones
 * 
 */
@Deprecated
class ConfigTest {
	private final ConfigEditor editor;
	private final String fields[];

	private ConfigTest(ConfigEditor _editor, String _fields[]) {
		editor = _editor;
		fields = _fields;

		int nFields = fields.length;

		JPanel panel = new JPanel();
		panel.setBorder(new EtchedBorder());
		UtilLayout layout = new UtilLayout(panel);
		final JFrame f = new JFrame();
		f.setSize(400, 400);
		f.addWindowListener(new BasicWindowMonitor());

		int y = 0;
		for (int i = 0; i < nFields; i++) {
			ConfigEntryComponent<?> entry = editor.getEntry(fields[i]);
			if (entry == null) {
				System.err.println("No entry for " + fields[i]);
				System.exit(0);
			}

			// set optional entries checked. Use setChecked(true) not
			// setChecked() which sets checked only is the field has a value
			// If you don't do this the field value won't be set.
			if (entry.isOptional())
				entry.setChecked(true);

			JComponent component = entry.getEditor();
			layout.add(new JLabel(entry.getName()), 0, y);
			layout.add(component, 1, y++);

			// Add comments- help information
			layout.add(entry.getHelpComponent(), 0, y++, 2, 1);
		}

		JPanel p = new JPanel(new BorderLayout());
		p.add(BorderLayout.CENTER, new JScrollPane(panel));

		JPanel controlPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		controlPanel.setBorder(new EtchedBorder());

		// This button prints out the configuration to the standard output
		JButton toString = new JButton("To String");
		controlPanel.add(toString);
		toString.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				String str = editor.getConfigModel().saveToString();
				System.out.println(str);
			}
		});

		JButton ok = new JButton("OK");
		controlPanel.add(ok);
		ok.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				f.dispose();
			}
		});

		p.add(BorderLayout.SOUTH, controlPanel);

		f.getContentPane().add(p);

		f.setVisible(true);

	}

	/**
	 * Main test method. You can pass a list of fields or use the dafault list.
	 * 
	 * @param args
	 *            List of fields to edit
	 */
	public static void main(String args[]) {

		if (args == null || args.length == 0)
			args = new String[] { "POPSIZE", "N_ITERATIONS", "SEED_FILE",
					"BASE_MOLECULE_SELECTION", "FLATTEN_BONDS" };

		ConfigEditor editor = new ConfigEditor();
		new ConfigTest(editor, args);
	}

}
