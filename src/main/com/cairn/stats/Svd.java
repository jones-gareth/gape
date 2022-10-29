package com.cairn.common.stats;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

/**
 * A simple wrapper around commons SVD.
 * 
 * @author Gareth Jones
 * 
 */
public class Svd {
	private final double[][] matrix;
	private SingularValueDecomposition svd;

	/**
	 * @param matrix
	 *            input matrix m by n
	 */
	public Svd(double[][] matrix) {
		this.matrix = matrix;
	}

	/**
	 * Perform SVD. Decomposes input matrix of m by n into U . Epsilon . VT
	 */
	public void solve() {
		Array2DRowRealMatrix inputMatrix = new Array2DRowRealMatrix(matrix);
		svd = new SingularValueDecomposition(inputMatrix);
	}

	/**
	 * @return the singular values
	 */
	public double[] getSingularValues() {
		return svd.getSingularValues();
	}

	/**
	 * Returns U the m by m orthogonal matrix
	 * 
	 * @return
	 */
	public double[][] getU() {
		RealMatrix u = svd.getU();
		return u.getData();
	}

	/**
	 * Returns the m by n matrix Epsilon. Epsilon is zero except for the
	 * diagonal which contains singular values.
	 * 
	 * @return
	 */
	public double[][] getEpsilon() {
		return svd.getS().getData();
	}

	/**
	 * Returns V the n by n orthogonal matrix
	 * 
	 * @return
	 */
	public double[][] getV() {
		RealMatrix v = svd.getV();
		return v.getData();
	}

	public double[][] getUT() {
		RealMatrix ut = svd.getUT();
		return ut.getData();
	}

	public double[][] getVT() {
		RealMatrix vt = svd.getVT();
		return vt.getData();
	}

}
