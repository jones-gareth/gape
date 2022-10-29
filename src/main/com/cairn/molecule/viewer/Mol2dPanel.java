package com.cairn.molecule.viewer;

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Image;

import javax.swing.JPanel;

import com.cairn.molecule.viewer.MolView.BackgroundStyle;

public class Mol2dPanel extends JPanel {
	Mol2dGraphics mol2dGraphics;

	Image buf;
	double angleX = .0, angleY = .0;

	Mol2dPanel() {
	}

	public Mol2dPanel(DisplayMolecule m) {
		mol2dGraphics = new Mol2dGraphics(m);
	}

	public boolean getHeavy() {
		return mol2dGraphics.getHeavy();
	}

	public void setHeavy(boolean h) {
		mol2dGraphics.setHeavy(h);
	}

	public boolean getHetero() {
		return mol2dGraphics.getHetero();
	}

	public void setHetero(boolean h) {
		mol2dGraphics.setHetero(h);
	}

	public void setBackgroundStyle(BackgroundStyle c) {
		mol2dGraphics.setBackgroundStyle(c);
	}

	public BackgroundStyle getBackgroundStyle() {
		return mol2dGraphics.getBackgroundStyle();
	}

	void redraw() {
		Graphics g = getGraphics();
		offImage();
		g.drawImage(buf, 0, 0, this);
		g.dispose();
		this.validate();
	}

	@Override
	public void paintComponent(Graphics g) {
		if (buf == null)
			offImage();
		g.drawImage(buf, 0, 0, this);
	}

	void offImage() {
		Dimension d = getSize();
		int width = mol2dGraphics.width;
		int height = mol2dGraphics.height;
		if (buf == null || height != d.height || width != d.width) {
			height = d.height;
			width = d.width;
			buf = createImage(width, height);
			if (buf == null) {
				System.err.println("No Image");
				return;
			}
			mol2dGraphics.setSize(new Dimension(width, height));
		}
		// Matrix trans = new Matrix();
		// trans.identity();
		Graphics g = buf.getGraphics();
		mol2dGraphics.plot(g);
		g.dispose();
	}
}
