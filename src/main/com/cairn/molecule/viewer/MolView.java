package com.cairn.molecule.viewer;

import java.awt.Canvas;
import java.awt.Color;

import javax.swing.JPanel;

import com.cairn.common.utils.UtilColor;

/**
 * A viewport and some display settings for the molecular viewer.
 * 
 * @author Gareth Jones
 * 
 */
public class MolView extends Viewport {
	public enum LabelType {
		NONE, TYPE, LABEL, ID, HETERO, PARTIAL_CHARGE, FORMAL_CHARGE, SUBSTR,
	};

	private LabelType atomLabelType = LabelType.NONE;
	private LabelType bondLabelType = LabelType.NONE;
	private boolean heavy = false;
	private boolean ballAndStick = false;

	public enum BackgroundStyle {
		BLACK, WHITE
	};

	private BackgroundStyle backgroundStyle = BackgroundStyle.WHITE;

	public void setBackgroundStyle(BackgroundStyle b) {
		backgroundStyle = b;
		if (b == BackgroundStyle.WHITE) {
			if (screen != null)
				screen.setBackground(Color.white);
			if (screen2 != null)
				screen2.setBackground(Color.white);
			if (screen3 != null)
				screen3.setBackground(Color.white);
		} else if (b == BackgroundStyle.BLACK) {
			if (screen != null)
				screen.setBackground(Color.black);
			if (screen2 != null)
				screen2.setBackground(Color.black);
			if (screen3 != null)
				screen3.setBackground(Color.black);
		}
	}

	public MolView(Canvas c, double maxX, double minX, double maxY, double minY,
			double maxZ, double minZ) {
		super(c, maxX, minX, maxY, minY, maxZ, minZ);
	}

	public MolView(JPanel c, double maxX, double minX, double maxY, double minY,
			double maxZ, double minZ) {
		super(c, maxX, minX, maxY, minY, maxZ, minZ);
	}

	public MolView(Screen c, double maxX, double minX, double maxY, double minY,
			double maxZ, double minZ) {
		super(c, maxX, minX, maxY, minY, maxZ, minZ);
	}

	private Color checkBackgroundStyle(Color c) {
		if (backgroundStyle == BackgroundStyle.WHITE && c == Color.white)
			c = Color.black;
		else if (backgroundStyle == BackgroundStyle.BLACK && c == Color.black)
			c = Color.white;
		else if (backgroundStyle == BackgroundStyle.WHITE && c == UtilColor.PART_WHITE)
			c = UtilColor.PART_BLACK;
		else if (backgroundStyle == BackgroundStyle.BLACK && c == UtilColor.PART_BLACK)
			c = UtilColor.PART_WHITE;
		return c;
	}

	@Override
	public void setColor(Color c) {
		super.setColor(checkBackgroundStyle(c));
	}

	@Override
	public void setColor(Color c, double v[]) {
		super.setColor(checkBackgroundStyle(c), v);
	}

	@Override
	public void drawCenteredString(String l, double p[], Color c) {
		super.drawCenteredString(l, p, checkBackgroundStyle(c));
	}

	Color getBackgroundColor() {
		if (backgroundStyle == BackgroundStyle.WHITE)
			return Color.white;
		if (backgroundStyle == BackgroundStyle.BLACK)
			return Color.black;
		return null;
	}

	/**
	 * @return the atomLabelType
	 */
	protected LabelType getAtomLabelType() {
		return atomLabelType;
	}

	/**
	 * @return the bondLabelType
	 */
	protected LabelType getBondLabelType() {
		return bondLabelType;
	}

	/**
	 * @return the heavy
	 */
	protected boolean isHeavy() {
		return heavy;
	}

	/**
	 * @return the ballAndStick
	 */
	protected boolean isBallAndStick() {
		return ballAndStick;
	}

	/**
	 * @return the backgroundStyle
	 */
	protected BackgroundStyle getBackgroundStyle() {
		return backgroundStyle;
	}

	/**
	 * @param atomLabelType
	 *            the atomLabelType to set
	 */
	protected void setAtomLabelType(LabelType atomLabelType) {
		this.atomLabelType = atomLabelType;
	}

	/**
	 * @param bondLabelType
	 *            the bondLabelType to set
	 */
	protected void setBondLabelType(LabelType bondLabelType) {
		this.bondLabelType = bondLabelType;
	}

	/**
	 * @param heavy
	 *            the heavy to set
	 */
	protected void setHeavy(boolean heavy) {
		this.heavy = heavy;
	}

	/**
	 * @param ballAndStick
	 *            the ballAndStick to set
	 */
	protected void setBallAndStick(boolean ballAndStick) {
		this.ballAndStick = ballAndStick;
	}

}
