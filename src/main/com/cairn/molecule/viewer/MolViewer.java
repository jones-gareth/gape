package com.cairn.molecule.viewer;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import com.cairn.common.utils.BasicWindowMonitor;
import com.cairn.molecule.Molecule;

/**
 * Application to create a simple molecular viewer in a frame.
 * 
 * @author Gareth Jones
 * 
 */
class MolViewer extends JFrame {
	private MolPanel molPanel;
	private static final int MAX_MOLS = 500;

	MolViewer(String files[], Molecule.FileType type, Molecule.Source source) {
		setTitle("Simple Molecule Viewer");
		addWindowListener(new BasicWindowMonitor());
		Utils.setLookAndFeel(this);

		List<DisplayMolecule> molecules = DisplayMolecule.loadFiles(files, type, source);
		if (molecules.size() > MAX_MOLS) {
			molecules = new ArrayList<>(molecules.subList(0, MAX_MOLS));
		}
		molPanel = new MolPanel(molecules);
		setSize(800, 500);
		int w = 500;
		int h = 500;
		molPanel.getMolOverlay().setPreferredSize(new Dimension(w, h));
		getContentPane().add(BorderLayout.CENTER, molPanel);
		int wc = molPanel.getMolControls().getPreferredSize().width;
		molPanel.setPreferredSize(new Dimension(w + wc + 20, h + 20));
		pack();
		setVisible(true);

		SwingUtilities.invokeLater(() -> molPanel.getMolOverlay().redraw());
	}

	public static void main(String argv[]) {

		Molecule.setAddHydrogensFlag(false);
		Molecule.setChargeFlag(false);
		Molecule.setSolvateFlag(false);

		if (argv.length > 0) {
			int no = 0;
			while (argv[no].startsWith("-")) {
				if (argv[no].equals("-i")) {
					Molecule.setAddHydrogensFlag(true);
					Molecule.setChargeFlag(true);
					Molecule.setSolvateFlag(true);
					continue;
				}
				if (argv[no].equals("-l")) {
					DisplayMolecule.setDisplayLonePairs(true);
					no++;
					continue;
				}
				if (argv[no].equals("-p")) {
					DisplayMolecule.setPrepare(true);
					no++;
					continue;
				}
			}

			new MolViewer(argv, Molecule.FileType.UNKNOWN, Molecule.Source.FILE);

		} else {
			System.err.println("Usage: MolViewer [-i] [-l] [-p] <structure_files>");
			System.exit(0);
		}
	}
}
