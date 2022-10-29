package com.cairn.gape.utils;

import org.apache.commons.math3.util.FastMath;

import com.cairn.common.utils.Coord;

/**
 * Class to represent a 3D Gaussian. Largly superceeded by GaussianList to avoid
 * overhead of object creation.
 * 
 * G = n * exp(-alpha(r-rN)
 * 
 * @author Gareth Jones
 * 
 */
public class Gaussian {

	private final double n; // Normalization factor

	private final double alpha; // Width exponent

	private final double[] rN; // Center

	private Gaussian parent1, parent2;

	private double volume;

	private static final boolean DEBUG = false;

	private static double maxAlpha = 0, minN = 0;

	public Gaussian(double _n, double _alpha, double[] _rN) {
		n = _n;
		alpha = _alpha;
		rN = _rN;
		volume();
		if (DEBUG) {
			System.out.print("Gausian N " + n + " alpha " + alpha + " vol " + volume
					+ " center " + Coord.info(rN));
		}
	}

	public double volume() {
		if (volume > 0)
			return volume;
		double f = Math.PI / alpha;
		double vol = n * f * Math.sqrt(f);
		if (DEBUG)
			System.out.println("Gaussian vol " + vol + " n " + n);
		volume = vol;
		return vol;
	}

	private final double rNScale[] = new double[4];

	private final double rNBScale[] = new double[4];

	public Gaussian overlap(Gaussian b) {
		double nB = b.n;
		double alphaB = b.alpha;
		double rNB[] = b.rN;

		double rABSqr = Coord.sqrDistance(rN, rNB);
		double newN = n * nB
				* FastMath.exp(-((alpha * alphaB) / (alpha + alphaB)) * rABSqr);
		double newAlpha = alphaB + alpha;

		Coord.copy(rN, rNScale);
		Coord.copy(rNB, rNBScale);
		Coord.multiply(rNScale, alpha);
		Coord.multiply(rNBScale, alphaB);
		double newRN[] = new double[4];
		Coord.add(rNScale, rNBScale, newRN);
		Coord.multiply(newRN, 1 / (alpha + alphaB));

		Gaussian g = new Gaussian(newN, newAlpha, newRN);
		g.parent1 = this;
		g.parent2 = b;

		if (DEBUG) {
			if (newN < minN)
				minN = newN;
			if (newAlpha > maxAlpha)
				maxAlpha = newAlpha;
		}

		return g;
	}

	protected double getAlpha() {
		return alpha;
	}

	protected double getN() {
		return n;
	}

	protected double[] getRN() {
		return rN;
	}

	protected Gaussian getParent1() {
		return parent1;
	}

	protected Gaussian getParent2() {
		return parent2;
	}

}
