package com.cairn.molecule.viewer;

import java.awt.BorderLayout;
import java.util.List;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.border.BevelBorder;

import com.cairn.common.utils.UtilColor;

/**
 * A panel for displaying a 3D molecular viewer with controls.
 * 
 * @author Gareth Jones
 * 
 */
/**
 * @author Gareth Jones
 *
 */
class MolPanel extends JPanel {
	private final List<DisplayMolecule> molecules;
	private final MolOverlay molOverlay;
	private final MolControls molControls;

	public MolPanel(List<DisplayMolecule> m) {
		molecules = m;
		setBackground(UtilColor.BG);
		setLayout(new BorderLayout());
		JPanel p2 = new JPanel(new BorderLayout());
		p2.setBorder(new BevelBorder(BevelBorder.LOWERED));
		molOverlay = new MolOverlay(m);
		p2.add(BorderLayout.CENTER, molOverlay);
		add(p2, BorderLayout.CENTER);
		JPanel p = new JPanel();
		molControls = new MolControls(molOverlay);
		p.add(molControls);
		add(new JScrollPane(p), BorderLayout.WEST);
	}

	/**
	 * @return the molecules
	 */
	public List<DisplayMolecule> getMolecules() {
		return molecules;
	}

	/**
	 * @return the molOverlay
	 */
	public MolOverlay getMolOverlay() {
		return molOverlay;
	}

	/**
	 * @return the molControls
	 */
	public MolControls getMolControls() {
		return molControls;
	}

}
