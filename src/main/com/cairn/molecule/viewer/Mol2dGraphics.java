package com.cairn.molecule.viewer;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Image;

import com.cairn.common.utils.Coord;
import com.cairn.molecule.viewer.DisplayMolecule.Bounds;
import com.cairn.molecule.viewer.MolView.BackgroundStyle;
import com.cairn.molecule.viewer.MolView.LabelType;

public class Mol2dGraphics implements Viewport.Screen {
	int x_coord, y_coord;
	int width, height;
	Image buf;
	DisplayMolecule molecule;
	double xMin, xMax, yMin, yMax, zMin, zMax;
	MolView vp;
	boolean hetero = false, heavy = true;
	BackgroundStyle backgroundStyle = BackgroundStyle.WHITE;
	Color background = Color.white;

	public Mol2dGraphics(DisplayMolecule m) {
		molecule = m;
		getMol2dPanelBounds();
		x_coord = 0;
		y_coord = 0;
	}

	public Mol2dGraphics(DisplayMolecule m, int x, int y, int w, int h) {
		this(m);
		width = w;
		height = h;
		x_coord = x;
		y_coord = y;
	}

	public boolean getHeavy() {
		return heavy;
	}

	public void setHeavy(boolean h) {
		heavy = h;
	}

	public boolean getHetero() {
		return hetero;
	}

	public void setHetero(boolean h) {
		hetero = h;
	}

	public void setBackgroundStyle(BackgroundStyle c) {
		backgroundStyle = c;
	}

	public BackgroundStyle getBackgroundStyle() {
		return backgroundStyle;
	}

	void getMol2dPanelBounds() {
		molecule.createMoleculeBounds();
		molecule.createDisplayCoordinates();
		Bounds bounds = molecule.getBounds();
		xMin = bounds.getxMin();
		yMin = bounds.getyMin();
		zMin = bounds.getzMin();
		xMax = bounds.getxMax();
		yMax = bounds.getyMax();
		zMax = bounds.getzMax();
		double diff = xMax - xMin;
		if ((yMax - yMin) > diff)
			diff = yMax - yMin;
		if (diff < .1)
			diff = .1;
		double x = (xMax + xMin) / 2;
		double y = (yMax + yMin) / 2;
		double z = (zMax + zMin) / 2;
		double scale = .55;
		xMin = -diff * scale;
		xMax = diff * scale;
		yMin = -diff * scale;
		yMax = diff * scale;
		zMin = -diff * scale;
		zMax = diff * scale;
		double trans[][] = new double[4][4];
		Coord.identity(trans);
		trans[3][0] = -x;
		trans[3][1] = -y;
		trans[3][2] = -z;
		molecule.resetDisplay();
		molecule.translateDisplay(trans);
		vp = new MolView(this, xMax, xMin, yMax, yMin, zMax, zMin);

		molecule.setCoords(molecule.getDisplay());
		molecule.recalculateRingCenters();
	}

	@Override
	public Dimension getSize() {
		return new Dimension(width, height);
	}

	public void setSize(Dimension size) {
		width = size.width;
		height = size.height;
	}

	@Override
	public Color getBackground() {
		return background;
	}

	@Override
	public void setBackground(Color c) {
		background = c;
	}

	public void plot(Graphics g) {
		g.translate(x_coord, y_coord);
		if (hetero)
			vp.setAtomLabelType(LabelType.HETERO);
		if (heavy)
			vp.setHeavy(true);
		vp.setBackgroundStyle(backgroundStyle);
		vp.setBounds();
		g.setColor(vp.getBackgroundColor());
		g.fillRect(0, 0, width, height);
		vp.g = g;
		vp.g.setFont(new Font("SansSerif", Font.PLAIN, 8));
		molecule.display(vp, null);
		g.translate(-x_coord, -y_coord);
	}

}
