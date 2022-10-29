package com.cairn.molecule.viewer;

import java.awt.BasicStroke;
import java.awt.Canvas;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;

import javax.swing.JPanel;

import com.cairn.common.utils.Coord;

public class Viewport {
	public Canvas screen;
	public JPanel screen2;
	public Screen screen3;
	// public Matrix transform;
	public double maxX, minX, maxY, minY, maxZ, minZ;
	public double width, height;
	public Graphics g;
	public boolean preserveRatio = true;

	public interface Screen {
		Dimension getSize();

		Color getBackground();

		void setBackground(Color c);
	}

	private Viewport(double maxX, double minX, double maxY, double minY, double maxZ,
			double minZ) {
		this.maxX = maxX;
		this.minX = minX;
		this.maxY = maxY;
		this.minY = minY;
		this.maxZ = maxZ;
		this.minZ = minZ;
	}

	public Viewport(Canvas c, double maxX, double minX, double maxY, double minY,
			double maxZ, double minZ) {
		this(maxX, minX, maxY, minY, maxZ, minZ);
		screen = c;
	}

	public Viewport(JPanel c, double maxX, double minX, double maxY, double minY,
			double maxZ, double minZ) {
		this(maxX, minX, maxY, minY, maxZ, minZ);
		screen2 = c;
	}

	public Viewport(Screen c, double maxX, double minX, double maxY, double minY,
			double maxZ, double minZ) {
		this(maxX, minX, maxY, minY, maxZ, minZ);
		screen3 = c;
	}

	/**
	 * Get the PreserveRatio value.
	 * 
	 * @return the PreserveRatio value.
	 */
	public boolean isPreserveRatio() {
		return preserveRatio;
	}

	/**
	 * Set the PreserveRatio value.
	 * 
	 * @param newPreserveRatio
	 *            The new PreserveRatio value.
	 */
	public void setPreserveRatio(boolean newPreserveRatio) {
		this.preserveRatio = newPreserveRatio;
	}

	synchronized public void zoom(double scale) {
		double sum = maxX + minX;
		double diff = maxX - minX;
		maxX = (sum + diff * scale) / 2;
		minX = (sum - diff * scale) / 2;
		sum = maxY + minY;
		diff = maxY - minY;
		maxY = (sum + diff * scale) / 2;
		minY = (sum - diff * scale) / 2;
		sum = maxZ + minZ;
		diff = maxZ - minZ;
		maxZ = (sum + diff * scale) / 2;
		minZ = (sum - diff * scale) / 2;
	}

	synchronized public void translate(double x, double y) {
		maxX += x;
		minX += x;
		maxY += y;
		minY += y;
	}

	public void setBounds() {
		Dimension d = null;
		if (screen != null)
			d = screen.getSize();
		if (screen2 != null)
			d = screen2.getSize();
		if (screen3 != null)
			d = screen3.getSize();
		width = d.width;
		height = d.height;

		if (preserveRatio) {
			double xDiff = maxX - minX;
			double yDiff = maxY - minY;
			double r = xDiff / yDiff;
			double newWidth = width;
			double newHeight = height;

			// Preserve ratios
			if (r > 1)
				newHeight = width * r;
			else
				newWidth = height / r;

			width = newWidth;
			height = newHeight;
		}
		// System.out.println("with "+width+" height "+height);
	}

	public int[] toScreen(double point[]) {
		// if (transform != null) point = transform.transPoint(point);
		return toScreen(point[0], point[1]);
	}

	public int[] toScreen(double x, double y) {
		int coords[] = new int[2];
		double sx = (x - minX) * width / (maxX - minX) + .5;
		double sy = height - (y - minY) * height / (maxY - minY) + .5;
		coords[0] = (int) sx;
		coords[1] = (int) sy;
		return coords;
	}

	public double[] toViewport(int _x, int _y) {
		double x = _x;
		double y = _y;

		x = x * (maxX - minX) / width + minX;
		y = (height - y) * (maxY - minY) / height + minY;
		return new double[] { x, y };
	}

	public double[] toScreenD(double point[]) {
		// if (transform != null) point = transform.transPoint(point);
		return toScreenD(point[0], point[1]);
	}

	public double[] toScreenD(double x, double y) {
		double coords[] = new double[2];
		double sx = (x - minX) * width / (maxX - minX) + .5;
		double sy = height - (y - minY) * height / (maxY - minY) + .5;
		coords[0] = sx;
		coords[1] = sy;
		return coords;
	}

	public void setColor(Color c) {
		try {
			g.setColor(c);
		} catch (java.lang.IllegalStateException ex) {
			;
		} catch (java.lang.NullPointerException ex) {
			;
		}
	}

	public void setColor(Color c, double v[]) {
		double z = v[2];
		double mid = (maxZ + minZ) / 2;
		// System.out.println("max "+max_z+" min "+min_z);
		if (z > mid) {
			double diff = (z - mid) / (maxZ - mid);
			// System.out.println ("+diff "+diff);
			for (double i = 0.25; i < diff; i += .5)
				c = c.brighter();
		} else {
			double diff = (mid - z) / (mid - minZ);
			// System.out.println ("-diff "+diff);
			for (double i = 0.25; i < diff; i += .5)
				c = c.darker();
		}
		setColor(c);
	}

	public void drawLine(double p1[], double p2[]) {
		int[] sp1 = toScreen(p1);
		int[] sp2 = toScreen(p2);
		try {
			g.drawLine(sp1[0], sp1[1], sp2[0], sp2[1]);
		} catch (java.lang.IllegalStateException ex) {
			;
		} catch (java.lang.NullPointerException ex) {
			;
		}
	}

	public void drawPolygon(double vs[][], boolean fill) {
		int xPoints[] = new int[vs.length];
		int yPoints[] = new int[vs.length];
		for (int i = 0; i < vs.length; i++) {
			int[] sp = toScreen(vs[i]);
			xPoints[i] = sp[0];
			yPoints[i] = sp[1];
		}
		try {
			if (fill)
				g.fillPolygon(xPoints, yPoints, vs.length);
			else
				g.drawPolygon(xPoints, yPoints, vs.length);
		} catch (java.lang.IllegalStateException ex) {
			;
		} catch (java.lang.NullPointerException ex) {
			;
		}
	}

	public void drawThickLine(double p1[], double p2[], double w) {
		Graphics2D g = (Graphics2D) this.g;
		double sw = w * width / (maxX - minX);

		double[] sp1 = toScreenD(p1);
		double[] sp2 = toScreenD(p2);
		java.awt.geom.Line2D.Double line = new java.awt.geom.Line2D.Double(sp1[0],
				sp1[1], sp2[0], sp2[1]);
		g.setStroke(new BasicStroke((float) sw));
		try {
			g.draw(line);
		} catch (java.lang.IllegalStateException ex) {
			;
		} catch (java.lang.NullPointerException ex) {
			;
		}
	}

	public void drawDashedLine(double p1[], double p2[], double size) {
		double[] start = new double[4];
		double[] end = new double[4];
		double[] dash = new double[4];
		double[] line = new double[4];
		double[] s = new double[4];
		double[] e = new double[4];

		Coord.subtract(p2, p1, line);
		double step = 0;
		Coord.copy(line, dash);
		Coord.setLength(dash, size);
		while (Coord.mag(start) < Coord.mag(line)) {
			Coord.copy(dash, end);
			Coord.setLength(end, step + size);
			Coord.add(p1, start, s);
			if (Coord.mag(end) > Coord.mag(line)) {
				Coord.add(p1, line, e);
				drawLine(s, e);
			} else {
				Coord.add(p1, end, e);
				drawLine(s, e);
			}
			step += 3 * size;
			Coord.copy(dash, start);
			Coord.setLength(start, step);
		}
	}

	public void drawLine(double x1, double y1, double x2, double y2) {
		int[] sp1 = toScreen(x1, y1);
		int[] sp2 = toScreen(x2, y2);
		try {
			g.drawLine(sp1[0], sp1[1], sp2[0], sp2[1]);
		} catch (java.lang.IllegalStateException ex) {
			;
		} catch (java.lang.NullPointerException ex) {
			;
		}
	}

	public void drawString(String l, double p[]) {
		int[] sp = toScreen(p);
		try {
			g.drawString(l, sp[0], sp[1]);
		} catch (java.lang.IllegalStateException ex) {
			;
		} catch (java.lang.NullPointerException ex) {
			;
		}
	}

	public void drawCenteredString(String l, double p[], Color c) {
		int[] sp = toScreen(p);
		FontMetrics fm = g.getFontMetrics();
		int w = fm.stringWidth(l);
		int h = fm.getHeight();
		int a = fm.getAscent();
		int x = sp[0] - w / 2;
		int y = sp[1] + a / 2;
		int y2 = sp[1] - a / 2;

		try {
			if (screen != null)
				g.setColor(screen.getBackground());
			if (screen2 != null)
				g.setColor(screen2.getBackground());
			if (screen3 != null)
				g.setColor(screen3.getBackground());
			g.fillRect(x, y2, w, h);
			g.setColor(c);
			g.drawString(l, x, y);
		} catch (java.lang.IllegalStateException ex) {
			;
		} catch (java.lang.NullPointerException ex) {
			;
		}
	}

	public void drawString(String l, double x, double y) {
		int[] sp = toScreen(x, y);
		try {
			g.drawString(l, sp[0], sp[1]);
		} catch (java.lang.IllegalStateException ex) {
			;
		} catch (java.lang.NullPointerException ex) {
			;
		}
	}

	private double d[] = new double[3];
	private double spV[] = new double[3];
	private double opV[] = new double[3];

	public void fillOval(double p[], double width, double height) {
		d[0] = width / 2;
		d[1] = height / 2;
		d[2] = 0;
		Coord.subtract(p, d, spV);
		Coord.add(p, d, opV);
		int[] sp = toScreen(spV);
		int[] op = toScreen(opV);
		int w = op[0] - sp[0];
		int h = sp[1] - op[1];
		// System.out.println("w "+w+" h "+h);
		try {
			g.fillOval(sp[0], op[1], w, h);
		} catch (java.lang.IllegalStateException ex) {
			;
		} catch (java.lang.NullPointerException ex) {
			;
		}
	}
}
