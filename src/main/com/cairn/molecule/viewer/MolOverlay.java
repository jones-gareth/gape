package com.cairn.molecule.viewer;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JPanel;

import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;
import com.cairn.molecule.viewer.DisplayMolecule.Bounds;
import com.cairn.molecule.viewer.MolView.BackgroundStyle;

/**
 * A 3D molecular viewer
 * 
 * @author Gareth Jones
 * 
 */
class MolOverlay extends JPanel implements MouseListener, MouseMotionListener {
	private static final Logger logger = Logger.getLogger(MolOverlay.class);
	static {
		// logger.setLevel(Level.DEBUG);
	}

	private int width, height;
	private Image buf;
	private final List<DisplayMolecule> molecules;
	private final List<Boolean> show;
	private boolean moving = false, drawing = false;
	private double xMin, xMax, yMin, yMax, zMin, zMax;
	private double angleX = .0, angleY = .0;

	private final double xAxis[] = new double[] { 5, 0, 0, 1 };
	private final double yAxis[] = new double[] { 0, 5, 0, 1 };
	private final double zAxis[] = new double[] { 0, 0, 5, 1 };
	private final double origin[] = new double[] { 0, 0, 0, 1 };
	private int lastX, lastY;
	private MolView molView;

	public enum Operation {
		ROTATE, ZOOM, TRANSLATE
	};

	private Operation operation = Operation.ROTATE;
	private BackgroundStyle backgroundStyle = BackgroundStyle.BLACK;
	private boolean colorByType = true;

	private final Color colors[] = new Color[] { Color.white, Color.green, Color.orange,
			Color.pink, Color.magenta, Color.orange, Color.yellow };

	public MolOverlay(List<DisplayMolecule> m) {
		molecules = m;
		show = new ArrayList<>();
		for (int i = 0; i < m.size(); i++) {
			show.add(true);
		}
		getMolOverlayBounds();
		this.addMouseListener(this);
		this.addMouseMotionListener(this);
		setColorByType(colorByType);
	}

	void setColorByType(boolean b) {
		colorByType = b;
		int no = 0;
		for (DisplayMolecule molecule : molecules) {
			molecule.setColorByType(b);
			if (!colorByType)
				molecule.setColor(colors[no % colors.length]);
			no++;
		}
	}

	void getMolOverlayBounds() {
		boolean boundsSet = false;

		for (DisplayMolecule molecule : molecules) {
			if (molecule.getnAtoms() == 0)
				continue;
			molecule.createMoleculeBounds();
			Bounds bounds = molecule.getBounds();
			if (!boundsSet) {
				xMin = bounds.getxMin();
				yMin = bounds.getyMin();
				zMin = bounds.getzMin();
				xMax = bounds.getxMax();
				yMax = bounds.getyMax();
				zMax = bounds.getzMax();
				boundsSet = true;
			} else {
				if (xMin > bounds.getxMin())
					xMin = bounds.getxMin();
				if (xMax < bounds.getxMax())
					xMax = bounds.getxMax();
				if (yMin > bounds.getyMin())
					yMin = bounds.getyMin();
				if (yMax < bounds.getyMax())
					yMax = bounds.getyMax();
				if (zMin > bounds.getzMin())
					zMin = bounds.getzMin();
				if (zMax < bounds.getzMax())
					zMax = bounds.getzMax();
			}
		}
		logger.debug("Bounds: " + xMin + " " + yMin + " " + zMin + " " + xMax + " "
				+ yMax + " " + zMax);
		double diff = xMax - xMin;
		if ((yMax - yMin) > diff)
			diff = yMax - yMin;
		if ((zMax - zMin) > diff)
			diff = zMax - zMin;
		if (diff < 2.5)
			diff = 2.5;
		double x = (xMax + xMin) / 2;
		double y = (yMax + yMin) / 2;
		double z = (zMax + zMin) / 2;
		xMin = -diff * .65;
		xMax = +diff * .65;
		yMin = -diff * .65;
		yMax = diff * .65;
		zMin = -diff * .65;
		zMax = diff * .65;

		double trans[][] = new double[4][4];
		Coord.identity(trans);
		trans[3][0] = -x;
		trans[3][1] = -y;
		trans[3][2] = -z;
		for (DisplayMolecule molecule : molecules) {
			molecule.setDisplay(trans);
		}

		logger.debug("Bounds: " + xMin + " " + yMin + " " + zMin + " " + xMax + " "
				+ yMax + " " + zMax);
		molView = new MolView(this, xMax, xMin, yMax, yMin, zMax, zMin);
		molView.setBackgroundStyle(backgroundStyle);
	}

	void redraw() {
		Graphics g = getGraphics();
		offImage();
		if (g != null) {
			g.drawImage(buf, 0, 0, this);
			g.dispose();
		}
		this.validate();
	}

	@Override
	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		if (buf == null)
			offImage();
		g.drawImage(buf, 0, 0, this);
	}

	void offImage() {
		if (drawing)
			return;
		drawing = true;
		Graphics g = null;
		try {
			Dimension d = getSize();
			logger.debug("Screen width = " + d.width + " height " + d.height);
			if (d.width == 0 || d.height == 0) {
				return;
			}
			if (buf == null || height != d.height || width != d.width) {
				height = d.height;
				width = d.width;
				buf = createImage(width, height);
				if (buf == null) {
					logger.error("No Image");
					return;
				}
			}
			g = buf.getGraphics();
			molView.g = g;

			double yTrans[][] = new double[4][4];
			double xTrans[][] = new double[4][4];
			double trans[][] = new double[4][4];

			double a1 = (angleY * Math.PI) / 180;
			Coord.determineRotation(origin, yAxis, a1, yTrans);
			Coord.transPointInPlace(yTrans, xAxis);
			double a2 = (angleX * Math.PI) / 180;
			Coord.determineRotation(origin, xAxis, a2, xTrans);
			Coord.transPointInPlace(xTrans, yAxis);
			Coord.product(yTrans, xTrans, trans);
			molView.setBounds();
			g.setColor(getBackground());
			g.fillRect(0, 0, width, height);
			for (int i = 0; i < molecules.size(); i++) {
				if (show.get(i))
					molecules.get(i).display(molView, trans);
				else
					molecules.get(i).translateDisplay(trans);
			}

			xAxis[0] = 5;
			xAxis[1] = 0;
			xAxis[2] = 0;
			yAxis[0] = 0;
			yAxis[1] = 5;
			yAxis[2] = 0;
			zAxis[0] = 0;
			zAxis[1] = 0;
			zAxis[2] = 5;

			angleX = angleY = 0;
		} finally {
			drawing = false;
			if (g != null) {
				g.dispose();
			}
		}
	}

	@Override
	public void mousePressed(MouseEvent e) {
		lastX = e.getX();
		lastY = e.getY();
	}

	@Override
	public void mouseDragged(MouseEvent e) {
		if (moving)
			return;
		moving = true;
		int x = e.getX();
		int y = e.getY();

		if (operation == Operation.ROTATE) {
			angleY += (double) (x - lastX) * 2;
			angleX += (double) (y - lastY) * 2;
			if (angleX < 0.0)
				angleX += 360.0;
			if (angleX > 360.0)
				angleX -= 360.0;
			if (angleY < 0.0)
				angleY += 360.0;
			if (angleY > 360.0)
				angleY -= 360.0;
			// System.out.println("x "+angleX+" y "+angleY);
		} else if (operation == Operation.TRANSLATE) {
			double[] viewportXY = molView.toViewport(x, y);
			double[] lastViewportXY = molView.toViewport(lastX, lastY);
			molView.translate(lastViewportXY[0] - viewportXY[0], lastViewportXY[1]
					- viewportXY[1]);
		} else if (operation == Operation.ZOOM) {
			if (x < lastX) {
				molView.zoom(1.1);
			} else if (lastX < x) {
				molView.zoom(.95);
			}
		}
		redraw();
		lastX = x;
		lastY = y;
		moving = false;
	}

	@Override
	public void mouseReleased(MouseEvent e) {
		;
	}

	@Override
	public void mouseClicked(MouseEvent e) {
		;
	}

	@Override
	public void mouseEntered(MouseEvent e) {
		;
	}

	@Override
	public void mouseExited(MouseEvent e) {
		;
	}

	@Override
	public void mouseMoved(MouseEvent e) {
		;
	}

	/**
	 * @return the width
	 */
	@Override
	public int getWidth() {
		return width;
	}

	/**
	 * @return the molecules
	 */
	public List<DisplayMolecule> getMolecules() {
		return molecules;
	}

	/**
	 * @return the show
	 */
	public List<Boolean> getShow() {
		return show;
	}

	/**
	 * @return the moving
	 */
	public boolean isMoving() {
		return moving;
	}

	/**
	 * @return the xMin
	 */
	public double getxMin() {
		return xMin;
	}

	/**
	 * @return the xMax
	 */
	public double getxMax() {
		return xMax;
	}

	/**
	 * @return the yMin
	 */
	public double getyMin() {
		return yMin;
	}

	/**
	 * @return the yMax
	 */
	public double getyMax() {
		return yMax;
	}

	/**
	 * @return the zMin
	 */
	public double getzMin() {
		return zMin;
	}

	/**
	 * @return the zMax
	 */
	public double getzMax() {
		return zMax;
	}

	/**
	 * @return the xAxis
	 */
	public double[] getxAxis() {
		return xAxis;
	}

	/**
	 * @return the yAxis
	 */
	public double[] getyAxis() {
		return yAxis;
	}

	/**
	 * @return the zAxis
	 */
	public double[] getzAxis() {
		return zAxis;
	}

	/**
	 * @return the origin
	 */
	public double[] getOrigin() {
		return origin;
	}

	/**
	 * @param width
	 *            the width to set
	 */
	public void setWidth(int width) {
		this.width = width;
	}

	/**
	 * @param show
	 *            the show to set
	 */
	public void setShow(int no, boolean show) {
		this.show.set(no, show);
	}

	/**
	 * @param moving
	 *            the moving to set
	 */
	public void setMoving(boolean moving) {
		this.moving = moving;
	}

	/**
	 * @param xMin
	 *            the xMin to set
	 */
	public void setxMin(double xMin) {
		this.xMin = xMin;
	}

	/**
	 * @param xMax
	 *            the xMax to set
	 */
	public void setxMax(double xMax) {
		this.xMax = xMax;
	}

	/**
	 * @param yMin
	 *            the yMin to set
	 */
	public void setyMin(double yMin) {
		this.yMin = yMin;
	}

	/**
	 * @param yMax
	 *            the yMax to set
	 */
	public void setyMax(double yMax) {
		this.yMax = yMax;
	}

	/**
	 * @param zMin
	 *            the zMin to set
	 */
	public void setzMin(double zMin) {
		this.zMin = zMin;
	}

	/**
	 * @param zMax
	 *            the zMax to set
	 */
	public void setzMax(double zMax) {
		this.zMax = zMax;
	}

	/**
	 * @return the backgroundStyle
	 */
	public BackgroundStyle getBackgroundStyle() {
		return backgroundStyle;
	}

	/**
	 * @param backgroundStyle
	 *            the backgroundStyle to set
	 */
	public void setBackgroundStyle(BackgroundStyle backgroundStyle) {
		this.backgroundStyle = backgroundStyle;
	}

	/**
	 * @return the colorByType
	 */
	public boolean isColorByType() {
		return colorByType;
	}

	/**
	 * @return the operation
	 */
	public Operation getOperation() {
		return operation;
	}

	/**
	 * @param operation
	 *            the operation to set
	 */
	public void setOperation(Operation operation) {
		this.operation = operation;
	}

	/**
	 * @return the molView
	 */
	protected MolView getMolView() {
		return molView;
	}

}
