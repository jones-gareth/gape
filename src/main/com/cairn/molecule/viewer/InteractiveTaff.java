package com.cairn.molecule.viewer;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;

import org.apache.log4j.Logger;

import com.cairn.common.utils.UtilColor;
import com.cairn.molecule.Atom;
import com.cairn.molecule.Molecule;
import com.cairn.molecule.Taff;

/**
 * Performs interactive simplex minimization. The minimization is done atom by
 * atom.
 * 
 * @author Gareth Jones
 * 
 */
public class InteractiveTaff extends Taff implements Runnable {
	private final MolOverlay mo;
	private static final Logger logger = Logger.getLogger(InteractiveTaff.class);

	public InteractiveTaff(DisplayMolecule mol) {
		super(mol);
		mo = null;
	}

	public InteractiveTaff(DisplayMolecule mol, MolOverlay mo) {
		super(mol);
		this.mo = mo;
	}

	@Override
	public void run() {
		minimizeAtoms();
	}

	public void energyWindow() {
		molEnergy();
		new StatusWindow();
	}

	@Override
	protected void redraw() {
		((DisplayMolecule) getMol()).recalculateRingCenters();
		mo.redraw();
	}

	@Override
	protected synchronized void minimizeAtoms() {
		StatusWindow sw = null;
		if (mo != null)
			sw = new StatusWindow();
		molEnergy();
		int no = getnCycles();
		TaffSimplexAtom simplex = new TaffSimplexAtom();
		double last = Double.MIN_VALUE;
		for (int i = 1; i <= no; i++) {
			Molecule mol = getMol();
			for (Atom test : mol.getAtoms()) {
				synchronized (mol) {
					setCurrentAtom(test);
					simplex.initSimplex();
					simplex.optimize();
					simplex.result();
				}
			}
			molEnergy();
			if (mo != null)
				redraw();
			if (sw != null)
				sw.reportCycle(i);
			logger.info("Cycle " + i + " Energy is " + getEnergy());
			if (Math.abs(last - getEnergy()) < getTolerance() && i > no / 2) {
				System.out.println("Converged");
				break;
			}
			last = getEnergy();
		}
		molEnergy();
		if (mo != null) {
			redraw();
			sw.rc.repaint();
		}

		if (logger.isDebugEnabled()) {
			logger.debug("VDW energy " + geteVdw());
			logger.debug("Bond Stretch " + geteBond());
			logger.debug("Angle Bend " + geteAng());
			logger.debug("Out of Plane " + geteOop());
			logger.debug("Energy is " + getEnergy());
		}
	}

	class StatusWindow extends JFrame {
		ReportCanvas rc;

		StatusWindow() {
			super("Energy");
			setResizable(false);
			addWindowListener(new BasicFrameMonitor());
			getContentPane().setLayout(new BorderLayout());
			rc = new ReportCanvas();
			getContentPane().add(rc, BorderLayout.CENTER);
			// setPreferredSize(new Dimension(250, 150));
			JButton b = new JButton("OK");
			b.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					dispose();
				}
			});
			getContentPane().add(b, BorderLayout.SOUTH);
			pack();
			setVisible(true);
		}

		void reportCycle(int no) {
			rc.reportCycle(no);
		}
	}

	class ReportCanvas extends JPanel {
		int width = 250;
		int height = 150;

		ReportCanvas() {
			setPreferredSize(new Dimension(width, height));
		}

		void reportCycle(int no) {
			Graphics g = this.getGraphics();
			if (g == null)
				return;
			FontMetrics fm = g.getFontMetrics();
			// ((Frame)getParent()).show();
			g.clearRect(0, 0, width, height);
			g.setColor(UtilColor.BG);
			g.fillRect(0, 0, width, height);
			g.setColor(Color.black);
			String msg = "Cycle " + no + " Energy " + format(getEnergy());
			int h = fm.getHeight();
			g.drawString(msg, 10, h * 2);
			g.dispose();
			// repaint();
		}

		@Override
		public void paintComponent(Graphics g) {
			super.paintComponent(g);
			FontMetrics fm = g.getFontMetrics();
			// ((Frame)getParent()).show();
			g.clearRect(0, 0, width, height);
			g.setColor(UtilColor.BG);
			g.fillRect(0, 0, width, height);
			g.setColor(Color.black);
			int h = fm.getHeight();
			int w = fm.stringWidth("Out of Plane") + 20;
			g.drawString("Total Energy", 10, h * 2);
			g.drawString(format(getEnergy()), w, h * 2);
			g.drawString("VDW", 10, h * 3);
			g.drawString(format(geteVdw()), w, h * 3);
			g.drawString("Bond Stretch", 10, h * 4);
			g.drawString(format(geteBond()), w, h * 4);
			g.drawString("Angle Bend", 10, h * 5);
			g.drawString(format(geteAng()), w, h * 5);
			g.drawString("Out of Plane", 10, h * 6);
			g.drawString(format(geteOop()), w, h * 6);
			g.drawString("Torsional", 10, h * 7);
			g.drawString(format(geteTor()), w, h * 7);
		}

		String format(double number) {
			return String.format("%.3f", number);
		}
	}

}
