package com.cairn.molecule.viewer;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.AbstractAction;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.border.BevelBorder;
import javax.swing.border.EtchedBorder;

import com.cairn.common.utils.UtilLayout;
import com.cairn.molecule.viewer.MolOverlay.Operation;
import com.cairn.molecule.viewer.MolView.BackgroundStyle;
import com.cairn.molecule.viewer.MolView.LabelType;

/**
 * Controls for a simple 3D molecular viewer.
 * 
 * @author Gareth Jones
 * 
 */
class MolControls extends JPanel {
	private MolOverlay molOverlay;
	private JComboBox<String> atomLabels, bondLabels;

	MolControls(final MolOverlay molOverlay) {
		this.molOverlay = molOverlay;
		molOverlay.getMolView().setAtomLabelType(LabelType.NONE);

		setBorder(new EtchedBorder());
		JPanel p = new JPanel();
		// p.setBorder(new EtchedBorder());
		p.setBorder(new BevelBorder(BevelBorder.RAISED));
		UtilLayout layout = new UtilLayout(p);
		layout.c.insets = new Insets(2, 2, 2, 2);
		int y = 0;

		molOverlay.getMolView().setHeavy(true);
		final JCheckBox heavyCb = new JCheckBox("Heavy", molOverlay.getMolView()
				.isHeavy());
		heavyCb.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent e) {
				boolean state = heavyCb.isSelected();
				molOverlay.getMolView().setHeavy(state);
				molOverlay.redraw();
			}
		});
		layout.add(heavyCb, 0, y);

		final JCheckBox ballAndStickCb = new JCheckBox("Stick", molOverlay.getMolView()
				.isBallAndStick());
		ballAndStickCb.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent e) {
				boolean state = ballAndStickCb.isSelected();
				molOverlay.getMolView().setBallAndStick(state);
				molOverlay.redraw();
			}
		});
		layout.add(ballAndStickCb, 1, y, 2, 1);
		y++;

		boolean bs = (molOverlay.getBackgroundStyle() == BackgroundStyle.BLACK) ? true
				: false;
		JCheckBox cb = new JCheckBox("Black BG", bs);
		layout.add(cb, 0, y);
		cb.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent e) {
				boolean state = ((JCheckBox) e.getSource()).isSelected();
				BackgroundStyle bg = state ? BackgroundStyle.BLACK
						: BackgroundStyle.WHITE;
				molOverlay.getMolView().setBackgroundStyle(bg);
				molOverlay.redraw();
			}
		});

		cb = new JCheckBox("Color by Type", molOverlay.isColorByType());
		layout.add(cb, 1, y++, 2, 1);
		cb.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent e) {
				boolean state = ((JCheckBox) e.getSource()).isSelected();
				molOverlay.setColorByType(state);
				molOverlay.redraw();
			}
		});

		JLabel l = new JLabel("Atoms");
		layout.add(l, 0, y);
		l = new JLabel("Bonds");
		layout.add(l, 1, y++);

		atomLabels = new JComboBox<>(new String[] { "none", "hetero", "type", "label",
				"id", "partial charge", "formal charge", "substructure" });
		atomLabels.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent e) {
				String selected = (String) atomLabels.getSelectedItem();
				if (selected.equals("hetero"))
					molOverlay.getMolView().setAtomLabelType(LabelType.HETERO);
				else if (selected.equals("type"))
					molOverlay.getMolView().setAtomLabelType(LabelType.TYPE);
				else if (selected.equals("label"))
					molOverlay.getMolView().setAtomLabelType(LabelType.LABEL);
				else if (selected.equals("id"))
					molOverlay.getMolView().setAtomLabelType(LabelType.ID);
				else if (selected.equals("none"))
					molOverlay.getMolView().setAtomLabelType(LabelType.NONE);
				else if (selected.equals("partial charge"))
					molOverlay.getMolView().setAtomLabelType(LabelType.PARTIAL_CHARGE);
				else if (selected.equals("formal charge"))
					molOverlay.getMolView().setAtomLabelType(LabelType.FORMAL_CHARGE);
				else if (selected.equals("substructure"))
					molOverlay.getMolView().setAtomLabelType(LabelType.SUBSTR);
				molOverlay.redraw();
			}
		});
		layout.add(atomLabels, 0, y);

		bondLabels = new JComboBox<>(new String[] { "none", "type", "id" });
		bondLabels.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent e) {
				String selected = (String) bondLabels.getSelectedItem();
				if (selected.equals("type"))
					molOverlay.getMolView().setBondLabelType(LabelType.TYPE);
				if (selected.equals("id"))
					molOverlay.getMolView().setBondLabelType(LabelType.ID);
				else if (selected.equals("none"))
					molOverlay.getMolView().setBondLabelType(LabelType.NONE);
				molOverlay.redraw();
			}
		});
		layout.add(bondLabels, 1, y);

		JButton b = new JButton("Reset");
		b.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				molOverlay.getMolOverlayBounds();
				molOverlay.redraw();
			}
		});
		layout.add(b, 2, y++);

		JPanel opPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		ButtonGroup opGrp = new ButtonGroup();
		final JRadioButton rotateRb = new JRadioButton("Rotate", true);
		rotateRb.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent e) {
				if (rotateRb.isSelected())
					molOverlay.setOperation(Operation.ROTATE);
			}
		});
		opPanel.add(rotateRb);
		opGrp.add(rotateRb);
		final JRadioButton moveRb = new JRadioButton("Move", false);
		moveRb.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent e) {
				if (moveRb.isSelected())
					molOverlay.setOperation(Operation.TRANSLATE);
			}
		});
		opPanel.add(moveRb);
		opGrp.add(moveRb);
		final JRadioButton zoomRb = new JRadioButton("Zoom", false);
		zoomRb.addItemListener(new ItemListener() {
			@Override
			public void itemStateChanged(ItemEvent e) {
				if (zoomRb.isSelected())
					molOverlay.setOperation(Operation.ZOOM);
			}
		});
		opPanel.add(zoomRb);
		opGrp.add(zoomRb);
		layout.add(opPanel, 0, y++, 3, 1);

		layout = new UtilLayout(this);
		layout.add(p, 0, 0);
		layout.add(getForcePanel(), 0, 1);
	}

	JPanel getForcePanel() {
		final JCheckBox checkBoxes[] = new JCheckBox[molOverlay.getMolecules().size()];
		final ItemListener[] listeners = new ItemListener[molOverlay.getMolecules()
				.size()];

		JPanel p2 = new JPanel(new FlowLayout(FlowLayout.LEFT));
		// p2.setBorder(new EtchedBorder());
		p2.setBorder(new BevelBorder(BevelBorder.RAISED));

		final JButton showAllButton = new JButton("Show All");
		showAllButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				for (int i = 0; i < checkBoxes.length; i++) {
					checkBoxes[i].removeItemListener(listeners[i]);
					checkBoxes[i].setSelected(true);
					checkBoxes[i].addItemListener(listeners[i]);
					molOverlay.setShow(i, true);
				}
				molOverlay.redraw();
			}
		});
		p2.add(showAllButton);
		final JButton hideAllButton = new JButton("Hide All");
		hideAllButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				for (int i = 0; i < checkBoxes.length; i++) {
					checkBoxes[i].removeItemListener(listeners[i]);
					checkBoxes[i].setSelected(false);
					checkBoxes[i].addItemListener(listeners[i]);
					molOverlay.setShow(i, false);
				}
				molOverlay.redraw();
			}
		});
		p2.add(hideAllButton);

		JPanel p = new JPanel();
		// p.setBorder(new EtchedBorder());
		p.setBorder(new BevelBorder(BevelBorder.RAISED));

		UtilLayout layout = new UtilLayout(p);
		layout.c.anchor = GridBagConstraints.CENTER;
		int y = 0;

		layout.add(new JLabel("No"), 0, y);
		layout.add(new JLabel("Name"), 1, y);
		layout.add(new JLabel("Show"), 2, y++);

		for (int j = 0; j < molOverlay.getMolecules().size(); j++) {
			final int i = j;
			int no = i + 1;
			JLabel l = new JLabel(String.valueOf(no));
			layout.add(l, 0, y);
			l = new JLabel(molOverlay.getMolecules().get(i).getName());
			layout.add(l, 1, y);

			final JCheckBox cb = new JCheckBox("", molOverlay.getShow().get(i));
			checkBoxes[i] = cb;
			listeners[i] = new ItemListener() {
				@Override
				public void itemStateChanged(ItemEvent e) {
					boolean state = cb.isSelected();
					molOverlay.setShow(i, state);
					molOverlay.redraw();
				}
			};
			cb.addItemListener(listeners[i]);
			layout.add(cb, 2, y);

			final JMenu menu = new JMenu("Operations");
			JMenuBar mb = new JMenuBar();
			mb.add(menu);
			mb.setBorder(new EtchedBorder());
			layout.add(mb, 3, y++);

			menu.add(new AbstractAction("Minimize") {
				@Override
				public void actionPerformed(ActionEvent e) {
					InteractiveTaff taff = new InteractiveTaff(molOverlay.getMolecules()
							.get(i), molOverlay);
					new Thread(taff).start();
				}
			});
			menu.add(new AbstractAction("Energy") {
				@Override
				public void actionPerformed(ActionEvent e) {
					InteractiveTaff taff = new InteractiveTaff(molOverlay.getMolecules()
							.get(i), molOverlay);
					taff.energyWindow();
				}
			});
			menu.add(new AbstractAction("Shake") {
				@Override
				public void actionPerformed(ActionEvent e) {
					molOverlay.getMolecules().get(i).shake();
					molOverlay.redraw();
				}
			});
			menu.add(new AbstractAction("Mol2") {
				@Override
				public void actionPerformed(ActionEvent e) {
					mol2Window(i);
				}
			});
		}

		JPanel p3 = new JPanel(new BorderLayout());
		p3.add(p, BorderLayout.CENTER);
		p3.add(p2, BorderLayout.NORTH);
		return p3;
	}

	void mol2Window(int molNo) {
		final JFrame f = new JFrame("Mol2 Format");
		f.setSize(400, 400);
		Mol2Window mol2 = new Mol2Window(molNo);
		f.getContentPane().add(BorderLayout.CENTER, mol2);
		JButton b = new JButton("OK");
		b.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				f.dispose();
			}
		});
		f.getContentPane().add(BorderLayout.SOUTH, b);
		f.setVisible(true);
	}

	class Mol2Window extends JPanel {

		Mol2Window(int molNo) {
			setLayout(new BorderLayout());
			setBorder(new EtchedBorder());
			String mol2 = molOverlay.getMolecules().get(molNo).sybylMol2String();
			JTextArea text = new JTextArea(mol2);
			text.setEditable(false);
			JScrollPane sp = new JScrollPane(text);
			add(BorderLayout.CENTER, sp);
		}
	}

}
