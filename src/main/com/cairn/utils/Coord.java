package com.cairn.common.utils;

import java.util.List;

import org.apache.commons.math3.util.FastMath;

import com.cairn.common.stats.Svd;

/**
 * Class of static methods for dealing with coordinate arrays. A coordinate is
 * of type double[4], the first three values being x,y,z and the last double is
 * normally 1.0. The user is responsible for keeping value 4 = 1. The matrix
 * class contains much of this in a OO class. These static methods were created
 * for efficiency.
 * 
 * Now should be thread safe
 * 
 * Some Java lang.Math methods use lang.StrictMath which is very slow. There is
 * a native JNI call for acos(). Commons Math FastMath now seems to offer the
 * best performance.
 * 
 * @author Gareth Jones
 * 
 */
public class Coord {

	static final boolean DEBUG = false;

	private Coord() {
		;
	}

	/**
	 * Sets all values in a vector to 0.
	 * 
	 * @param p
	 */
	public static final void zero(double p[]) {
		p[0] = p[1] = p[2] = .0;
	}

	private static ThreadLocal<double[]> temp = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[4];
		}
	};

	/**
	 * Return a String representation of a point
	 * 
	 * @param p
	 * @return
	 */
	public static final String toString(double p[]) {
		return "[" + p[0] + "," + p[1] + "," + p[2] + "]";
	}

	/**
	 * Converts a String representation (created by toString) back to double[4].
	 * 
	 * @param s
	 * @param p
	 * 
	 * @see #toString(double[])
	 */
	public static final void toPoint(String str, double p[]) {
		String vals[] = str.split(",");
		for (int i = 0; i < 3; i++) {
			String s = vals[i];
			s = s.replace('[', ' ');
			s = s.replace(']', ' ');
			s = s.trim();
			p[i] = Double.valueOf(s).doubleValue();
		}
		p[3] = 1.0;
	}

	/**
	 * Vector substraction of b from a. The result is placed in temporary
	 * storage, which may be overwritten in later calls to Coord functions.
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static final double[] subtract(double a[], double b[]) {
		subtract(a, b, temp.get());
		return temp.get();
	}

	/**
	 * Vector substraction of b from a. The result is placed in c.
	 * 
	 * @param a
	 * @param b
	 * @param c
	 */
	public static final void subtract(double a[], double b[], double c[]) {
		c[0] = a[0] - b[0];
		c[1] = a[1] - b[1];
		c[2] = a[2] - b[2];
	}

	/**
	 * Vector substraction of b from a. The result is placed in a.
	 * 
	 * @param a
	 * @param b
	 */
	public static final void subtractInPlace(double a[], double b[]) {
		a[0] = a[0] - b[0];
		a[1] = a[1] - b[1];
		a[2] = a[2] - b[2];
	}

	/**
	 * Multiplies x,y,z values of p by scalar m
	 * 
	 * @param p
	 * @param m
	 */
	public static final void multiply(double p[], double m) {
		p[0] = p[0] * m;
		p[1] = p[1] * m;
		p[2] = p[2] * m;
	}

	/**
	 * Determines the mid-point of vectors a and b. The result is placed in
	 * temporary storage, which may be overwritten in later calls to Coord
	 * functions.
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static final double[] midPoint(double a[], double b[]) {
		midPoint(a, b, temp.get());
		return temp.get();
	}

	/**
	 * Determines the mid-point of vectors a and b. The result is placed in m.
	 * 
	 * @param a
	 * @param b
	 * @param m
	 */
	public static final void midPoint(double a[], double b[], double m[]) {
		m[0] = (a[0] + b[0]) / 2.0;
		m[1] = (a[1] + b[1]) / 2.0;
		m[2] = (a[2] + b[2]) / 2.0;
	}

	/**
	 * Determines the mid-point of vectors a and b. The result is placed back in
	 * a.
	 * 
	 * @param a
	 * @param b
	 */
	public static final void midPointInPlace(double a[], double b[]) {
		a[0] = (a[0] + b[0]) / 2.0;
		a[1] = (a[1] + b[1]) / 2.0;
		a[2] = (a[2] + b[2]) / 2.0;
	}

	/**
	 * Returns the magnitude of p.
	 * 
	 * @param p
	 * @return
	 */
	public static final double mag(double p[]) {
		return Math.sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
	}

	/**
	 * Returns the squared-distance between a and b.
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static final double sqrDistance(double[] a, double b[]) {
		return (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1])
				+ (a[2] - b[2]) * (a[2] - b[2]);
	}

	/**
	 * Returns the distance between a and b
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static final double distance(double a[], double b[]) {
		return Math.sqrt(sqrDistance(a, b));
	}

	/**
	 * Adds vectors a and b. The result is placed in temporary storage, which
	 * may be overwritten in later calls to Coord functions.
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static final double[] add(double[] a, double b[]) {
		add(a, b, temp.get());
		return temp.get();
	}

	/**
	 * Adds vectors a and b. The result is placed in c.
	 * 
	 * @param a
	 * @param b
	 * @param c
	 */
	public static final void add(double[] a, double b[], double c[]) {
		c[0] = a[0] + b[0];
		c[1] = a[1] + b[1];
		c[2] = a[2] + b[2];
	}

	/**
	 * Adds vectors a and b. The result is placed back in a
	 * 
	 * @param a
	 * @param b
	 */
	public static final void addInPlace(double[] a, double b[]) {
		a[0] = a[0] + b[0];
		a[1] = a[1] + b[1];
		a[2] = a[2] + b[2];
	}

	/**
	 * Adds vectors a and b. The result is placed in temporary storage, which
	 * may be overwritten in later calls to Coord functions.
	 * 
	 * @param x
	 * @param y
	 * @return
	 */
	public static final double[] vectorProduct(double x[], double y[]) {
		vectorProduct(x, y, temp.get());
		return temp.get();
	}

	/**
	 * Determine the vector product of vectors x and y. The product is placed in
	 * vp.
	 * 
	 * @param x
	 * @param y
	 * @param vp
	 */
	public static final void vectorProduct(double x[], double y[], double vp[]) {
		vp[0] = x[1] * y[2] - x[2] * y[1];
		vp[1] = x[2] * y[0] - x[0] * y[2];
		vp[2] = x[0] * y[1] - x[1] * y[0];
	}

	/**
	 * Scales vector p (in place) to a normal vector.
	 * 
	 * @param p
	 */
	public static final void unit(double p[]) {
		double m = mag(p);
		p[0] = p[0] / m;
		p[1] = p[1] / m;
		p[2] = p[2] / m;
	}

	/**
	 * Returns the scalar product of a and b
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static final double scalarProduct(double a[], double b[]) {
		return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	}

	/**
	 * Determines cos(c) where c is the angle between vectors a and b.
	 * 
	 * @param a
	 * @param b
	 * @return
	 * 
	 * @see #scalarProduct(double[], double[])
	 */
	public static final double cosine(double a[], double b[]) {
		double sp = scalarProduct(a, b);
		double mag = mag(a);
		double bMag = mag(b);
		double cosine = sp / (mag * bMag);
		return cosine;
	}

	/**
	 * Determines sin(c) where c is the angle between vectors a and b.
	 * 
	 * @param a
	 * @param b
	 * @return
	 * 
	 * @see #vectorProduct(double[], double[])
	 */
	public static final double sin(double a[], double b[]) {
		double vp = mag(vectorProduct(a, b));
		double mag = mag(a);
		double bMag = mag(b);
		double sin = vp / (mag * bMag);
		return sin;
	}

	/**
	 * Returns th angle beteen vectors a and b
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static final double angle(double a[], double b[]) {
		double cosine = cosine(a, b);
		if (cosine > 1.0)
			return 0;
		if (cosine < -1.0)
			return Math.PI;
		if (DEBUG)
			System.out.println("cosine " + cosine);
		return FastMath.acos(cosine);
	}

	/**
	 * The distance that you have to go out on b to draw a line parallel to c in
	 * the plane of the two vectors at a distance d from this.
	 * 
	 * @param c
	 * @param b
	 * @param d
	 * @return
	 */
	public static final double parallelDistance(double c[], double b[], double d) {
		double cosine = cosine(c, b);
		double a = FastMath.acos(cosine);
		a = Math.PI / 2 - a;
		double d2 = d / Math.cos(a);
		if (d2 < .0)
			d2 = -d2;
		if (DEBUG)
			System.out.println("parallel dist " + d2);
		return d2;
	}

	/**
	 * Returns the distance (d2) that you go out on a when you go out a distance
	 * d on b such that d2 is the hypotenuse of a right angle triangle.
	 * 
	 * @param a
	 * @param b
	 * @param d
	 * @return
	 */
	public static final double perpendicularDistance(double[] a, double b[], double d) {
		double cosine = cosine(a, b);
		double d2 = d / cosine;
		if (d2 < .0)
			d2 = -d2;
		if (DEBUG)
			System.out.println("perpendicular dist " + d2);
		return d2;
	}

	/**
	 * Sets the lenght (in place) of vector a to mag.
	 * 
	 * @param a
	 * @param mag
	 */
	public static final void setLength(double a[], double mag) {
		unit(a);
		a[0] = a[0] * mag;
		a[1] = a[1] * mag;
		a[2] = a[2] * mag;
	}

	/**
	 * Almost the same as toString, but returns 4th coordinate.
	 * 
	 * @param p
	 * @return
	 * 
	 * @see #toString(double[])
	 */
	public static final String info(double p[]) {
		return "[" + p[0] + "," + p[1] + "," + p[2] + "," + p[3] + "]";
	}

	/**
	 * Print our information about the string.
	 * 
	 * @param p
	 */
	public static final void print(double p[]) {
		System.out.println(info(p));
	}

	/**
	 * Copies a to b
	 * 
	 * @param a
	 * @param b
	 */
	public static final void copy(double a[], double b[]) {
		b[0] = a[0];
		b[1] = a[1];
		b[2] = a[2];
		b[3] = 1;
	}

	/**
	 * Determines the center of points and puts in c.
	 * 
	 * @param points
	 * @param c
	 */
	public static final void center(double points[][], double c[]) {
		c[0] = c[1] = c[2] = .0;
		for (int i = 0; i < points.length; i++) {
			c[0] += points[i][0];
			c[1] += points[i][1];
			c[2] += points[i][2];
		}
		double n = points.length;
		c[0] /= n;
		c[1] /= n;
		c[2] /= n;
	}

	/**
	 * Determines the center of points and puts in c.
	 * 
	 * @param points
	 * @param c
	 */
	public static final void center(List<double[]> points, double[] c) {
		c[0] = c[1] = c[2] = .0;
		for (double[] point : points) {
			c[0] += point[0];
			c[1] += point[1];
			c[2] += point[2];
		}
		double n = points.size();
		c[0] /= n;
		c[1] /= n;
		c[2] /= n;
	}

	static private ThreadLocal<double[]> i = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[3];
		}
	};
	static private ThreadLocal<double[]> j = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[3];
		}
	};
	static private ThreadLocal<double[]> k = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[3];
		}
	};
	static private ThreadLocal<double[]> n1 = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[3];
		}
	};
	static private ThreadLocal<double[]> n2 = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[3];
		}
	};
	static private ThreadLocal<double[]> tvp = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[3];
		}
	};

	/**
	 * Returns the torsion angle betwwen 4 points.
	 * 
	 * @param a
	 * @param b
	 * @param c
	 * @param d
	 * @return
	 */
	public static final double torsion(double a[], double b[], double c[], double d[]) {
		double i[] = Coord.i.get();
		double j[] = Coord.j.get();
		double k[] = Coord.k.get();
		double n1[] = Coord.n1.get();
		double n2[] = Coord.n2.get();
		double tvp[] = Coord.tvp.get();

		subtract(b, a, i);
		subtract(c, b, j);
		subtract(d, c, k);
		// find triple vector product
		vectorProduct(i, j, n1);
		vectorProduct(j, k, n2);
		vectorProduct(n1, n2, tvp);
		double sp = scalarProduct(n2, n1);
		double sp2 = scalarProduct(tvp, j);
		double mag = mag(n1) * mag(n2);
		double cosVal = sp / mag;
		if (cosVal > 1.0)
			cosVal = 1.0;
		if (cosVal < -1.0)
			cosVal = -1.0;
		double angle = FastMath.acos(cosVal);
		if (sp2 < 0.0)
			return -angle;
		return angle;
	}

	/**
	 * Returns true if the distance between the two points is less than tol.
	 * 
	 * @param a
	 * @param b
	 * @param tol
	 * @return
	 */
	public static final boolean withinDistance(double a[], double b[], double tol) {
		double sqrDist = sqrDistance(a, b);
		if (sqrDist < tol * tol)
			return true;
		return false;
	}

	// Matrices are treated as type double[][]

	/**
	 * Transposes matirx m and places the result in t.
	 * 
	 * @param m
	 * @param t
	 */
	public static final void transpose(double m[][], double t[][]) {
		for (int i = 0; i < m.length; i++)
			for (int j = 0; j < m[i].length; j++)
				t[j][i] = m[i][j];
	}

	static private ThreadLocal<double[][]> trans = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};

	/**
	 * Transponses matrix m and returns the result. The result is placed in
	 * temporary storage which may be overwritten in later calls. Works for 4x4
	 * matrix.
	 * 
	 * @param m
	 */
	public static final void transpose(double m[][]) {
		double trans[][] = Coord.trans.get();
		for (int i = 0; i < m.length; i++)
			for (int j = 0; j < m[i].length; j++)
				trans[j][i] = m[i][j];
		copy(trans, m);
	}

	/**
	 * Copies matrix m to c.
	 * 
	 * @param m
	 * @param c
	 */
	public static final void copy(double m[][], double c[][]) {
		for (int i = 0; i < m.length; i++)
			for (int j = 0; j < m[i].length; j++)
				c[i][j] = m[i][j];
	}

	/**
	 * Performs scalar multiplication of matrix. The matrix is modified in
	 * place.
	 * 
	 * @param m
	 * @param scale
	 */
	public static final void scale(double m[][], double scale) {
		for (int i = 0; i < m.length; i++)
			for (int j = 0; j < m[i].length; j++)
				m[i][j] *= scale;

	}

	/**
	 * Performs matrix addition. The sum is placed in m.
	 * 
	 * @param m
	 * @param m2
	 */
	public static final void add(double m[][], double m2[][]) {
		for (int i = 0; i < m.length; i++)
			for (int j = 0; j < m[i].length; j++)
				m[i][j] += m2[i][j];
	}

	/**
	 * Performs matrix subtraction. The difference m-m2 is placed in m.
	 * 
	 * @param m
	 * @param m2
	 */
	public static final void subtract(double m[][], double m2[][]) {
		for (int i = 0; i < m.length; i++)
			for (int j = 0; j < m[i].length; j++)
				m[i][j] -= m2[i][j];
	}

	/**
	 * Sets all elements of a matrix to zero.
	 * 
	 * @param m
	 */
	public static final void zero(double m[][]) {
		for (int i = 0; i < m.length; i++)
			for (int j = 0; j < m[i].length; j++)
				m[i][j] = 0;
	}

	/**
	 * Sets m to an identity matrix.
	 * 
	 * @param m
	 */
	public static final void identity(double m[][]) {
		zero(m);
		for (int i = 0; i < m.length; i++)
			m[i][i] = 1;
	}

	/**
	 * Calculates the matrix product of m1 and m2 and places the result in p.
	 * 
	 * @param m1
	 * @param m2
	 * @param p
	 */
	public static final void product(double m1[][], double m2[][], double p[][]) {
		zero(p);
		for (int k = 0; k < m2[0].length; k++)
			for (int i = 0; i < m1.length; i++)
				for (int s = 0; s < m1[0].length; s++)
					p[i][k] += m1[i][s] * m2[s][k];
	}

	/**
	 * Calculates the weighted matrix product of ma and mb and places the result
	 * in p.
	 * 
	 * @param ma
	 * @param mb
	 * @param weights
	 * @param p
	 */
	public static final void weightedProduct(double ma[][], double mb[][],
			double weights[], double p[][]) {
		zero(p);
		for (int k = 0; k < mb[0].length; k++)
			for (int i = 0; i < ma.length; i++)
				for (int s = 0; s < ma[0].length; s++)
					p[i][k] += weights[s] * ma[i][s] * mb[s][k];
	}

	/**
	 * Determines the determinant of matrix m (order 2 or 3 matrices only)
	 * 
	 * @param m
	 * @return
	 */
	public static final double determinant(double m[][]) {
		if (m.length == 2) {
			return m[0][0] * m[1][1] - m[0][1] * m[1][0];
		} else if (m.length == 3) {
			double det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1]
					* (m[1][0] * m[2][2] - m[2][0] * m[1][2]) + m[0][2]
					* (m[1][0] * m[2][1] - m[2][0] * m[1][1]);
			return det;
		}
		return .0;
	}

	/**
	 * Multiplies vector p by matrix m and returns the result in t.
	 * 
	 * @param m
	 * @param p
	 * @param t
	 */
	public static final void transPoint(double m[][], double p[], double t[]) {
		int l = p.length;
		for (int k = 0; k < l; k++) {
			t[k] = 0;
			for (int s = 0; s < l; s++)
				t[k] += p[s] * m[s][k];
		}
	}

	// 4x4 Matrix only!!
	/**
	 * Mutliplies vector p by matrix m and stores the result in p. For 4x4
	 * matrices and vectors of length 4 only.
	 * 
	 * @param m
	 * @param p
	 */
	public static final void transPointInPlace(double m[][], double p[]) {
		double x = 0, y = 0, z = 0, w = 0;
		for (int s = 0; s < 4; s++) {
			x += p[s] * m[s][0];
			y += p[s] * m[s][1];
			z += p[s] * m[s][2];
			w += p[s] * m[s][3];
		}
		p[0] = x;
		p[1] = y;
		p[2] = z;
		p[3] = w;
	}

	static private ThreadLocal<double[][]> rotPoint = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};

	/**
	 * Performs rotation of point p3 around the axis defined by points p1 and
	 * p2. The new point is stored in p. 4x4 system only (obviously).
	 * 
	 * @param p1
	 * @param p2
	 * @param p3
	 * @param angle
	 * @param p
	 * 
	 * @see #determineRotation(double[], double[], double, double[][])
	 */
	public static final void rotatePoint(double p1[], double p2[], double p3[],
			double angle, double p[]) {
		double[][] rotPoint = Coord.rotPoint.get();
		determineRotation(p1, p2, angle, rotPoint);
		p3[3] = 1.0;
		transPoint(rotPoint, p3, p);
	}

	/**
	 * Information string for matrices.
	 * 
	 * @param m
	 * @return
	 */
	public static final String info(double m[][]) {
		StringBuilder rtn = new StringBuilder();
		rtn.append("[");
		for (int i = 0; i < m.length; i++) {
			rtn.append(" [");
			for (int j = 0; j < m[i].length; j++) {
				rtn.append(String.valueOf(m[i][j]));
				if (j != m[i].length - 1)
					rtn.append(",");
			}
			rtn.append("]");
		}
		rtn.append("]");
		return rtn.toString();
	}

	/**
	 * Prints out matrix info
	 * 
	 * @param m
	 */
	public static final void print(double m[][]) {
		System.out.println(info(m));
	}

	// 4x4 Matrix only!!

	static private ThreadLocal<double[][]> rot = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};
	static private ThreadLocal<double[][]> xRot = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};
	static private ThreadLocal<double[][]> yRot = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};
	static private ThreadLocal<double[][]> invTrans = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};
	static private ThreadLocal<double[][]> invXrot = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};
	static private ThreadLocal<double[][]> invYrot = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};
	static private ThreadLocal<double[][]> rotation = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};
	static private ThreadLocal<double[][]> product1 = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};
	static private ThreadLocal<double[][]> product2 = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};
	static private ThreadLocal<double[]> diff = new ThreadLocal<double[]>() {
		@Override
		protected double[] initialValue() {
			return new double[4];
		}
	};

	/**
	 * Determines the matrix transformation to perform a rotation (of angle
	 * angle) about an arbitrary axis (defined by p1 and p2). The result is
	 * returned in matrix.
	 * 
	 * Ref: Newman & Sproull, Principles of Interactive computer Graphics,
	 * McGraw Hill, 1981 pp 346-348.
	 * 
	 * @param p1
	 * @param p2
	 * @param angle
	 * @param matrix
	 */
	public static final void determineRotation(double p1[], double p2[], double angle,
			double matrix[][]) {

		double[][] rot = Coord.rot.get();
		double[][] xRot = Coord.xRot.get();
		double[][] yRot = Coord.yRot.get();
		double[][] invTrans = Coord.invTrans.get();
		double[][] invXrot = Coord.invXrot.get();
		double[][] invYrot = Coord.invYrot.get();
		double[][] rotation = Coord.rotation.get();
		double[][] product1 = Coord.product1.get();
		double[][] product2 = Coord.product2.get();
		double[] diff = Coord.diff.get();
		double[][] trans = Coord.trans.get();

		double a, b, c, v, x, y, z;

		x = p1[0];
		y = p1[1];
		z = p1[2];
		subtract(p1, p2, diff);
		unit(diff);
		a = diff[0];
		b = diff[1];
		c = diff[2];

		// trans: translates point one to the origin. xRot and yRot:
		// rotation about x and y axis such that point2 is aligned on
		// the z-axis. rot: rotation of angle radians about the
		// z-axis. invTrans, invXrot and invYrot: inverse
		// translations and rotations.

		identity(trans);
		trans[3][0] = -x;
		trans[3][1] = -y;
		trans[3][2] = -z;
		identity(invTrans);
		invTrans[3][0] = x;
		invTrans[3][1] = y;
		invTrans[3][2] = z;
		v = Math.sqrt(b * b + c * c);
		identity(xRot);
		identity(invXrot);
		if (v > .0) {
			xRot[1][1] = xRot[2][2] = invXrot[1][1] = invXrot[2][2] = c / v;
			xRot[1][2] = invXrot[2][1] = b / v;
			xRot[2][1] = invXrot[1][2] = -b / v;
		}
		identity(yRot);
		identity(invYrot);
		yRot[0][0] = yRot[2][2] = invYrot[0][0] = invYrot[2][2] = v;
		yRot[0][2] = invYrot[2][0] = a;
		yRot[2][0] = invYrot[0][2] = -a;
		rot[2][2] = rot[3][3] = 1;
		rot[0][0] = rot[1][1] = Math.cos(angle);
		rot[0][1] = -1.0 * Math.sin(angle);
		rot[1][0] = Math.sin(angle);
		product(xRot, yRot, product1);
		product(product1, rot, product2);
		product(product2, invYrot, product1);
		product(product1, invXrot, rotation);
		product(trans, rotation, product1);
		product(product1, invTrans, matrix);
	}

	/**
	 * Uses SVD to least squares fit the xCoords to yCoords. A matrix transform
	 * is returned. The transform can then be applied to the xCoords to best fit
	 * them to yCoords. The transform includes a correction to prevent
	 * inversion.
	 * 
	 * @param xCoords
	 * @param yCoords
	 * @param matrix
	 */
	public static final void leastSquaresFit(double xCoords[][], double yCoords[][],
			double matrix[][]) {
		leastSquaresFit(xCoords, yCoords, null, matrix);
	}

	static private ThreadLocal<double[][]> yPrimeX = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[3][3];
		}
	};
	static private ThreadLocal<double[][]> uPrime = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[3][3];
		}
	};
	static private ThreadLocal<double[][]> rotMat = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[3][3];
		}
	};
	static private ThreadLocal<double[][]> in = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};
	static private ThreadLocal<double[][]> out = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};
	static private ThreadLocal<double[][]> product = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};
	static private ThreadLocal<double[][]> rotate = new ThreadLocal<double[][]>() {
		@Override
		protected double[][] initialValue() {
			return new double[4][4];
		}
	};

	/**
	 * Uses SVD to least squares fit the xCoords to yCoords with weighted
	 * points. A matrix transform is returned. The transform can then be applied
	 * to the xCoords to best fit them to yCoords. The transform includes a
	 * correction to prevent inversion.
	 * 
	 * @param xCoords
	 * @param yCoords
	 * @param weights
	 * @param matrix
	 */
	public static final void leastSquaresFit(double xCoords[][], double yCoords[][],
			double weights[], double matrix[][]) {

		double yPrimeX[][] = Coord.yPrimeX.get();
		double uPrime[][] = Coord.uPrime.get();
		double rotMat[][] = Coord.rotMat.get();
		double v[][] = null;// Coord.v.get();
		double w[] = null; // Coord.w.get();
		double in[][] = Coord.in.get();
		double out[][] = Coord.out.get();
		double product[][] = Coord.product.get();
		double rotate[][] = Coord.rotate.get();

		int nPoints = 0;
		for (int i = 0; i < xCoords.length; i++)
			if (xCoords[i] != null)
				nPoints++;
		if (DEBUG)
			System.out.println("No Points " + nPoints);

		// First create X and Y'
		double yPrime[][] = new double[3][nPoints];
		double x[][] = new double[nPoints][3];

		double x0 = .0, x1 = .0, x2 = .0, y0 = .0, y1 = .0, y2 = .0, weightSum = .0;
		int no = 0;

		for (int i = 0; i < xCoords.length; i++) {
			if (xCoords[i] == null)
				continue;
			double wt = 1.0;
			if (weights != null)
				wt = weights[i];
			double xP[] = xCoords[i];
			double yP[] = yCoords[i];

			x0 += xP[0] * wt;
			x1 += xP[1] * wt;
			x2 += xP[2] * wt;

			y0 += yP[0] * wt;
			y1 += yP[1] * wt;
			y2 += yP[2] * wt;

			weightSum += wt;

			x[no][0] = xP[0];
			x[no][1] = xP[1];
			x[no][2] = xP[2];
			yPrime[0][no] = yP[0];
			yPrime[1][no] = yP[1];
			yPrime[2][no] = yP[2];

			no++;
		}

		x0 /= weightSum;
		x1 /= weightSum;
		x2 /= weightSum;
		y0 /= weightSum;
		y1 /= weightSum;
		y2 /= weightSum;

		if (DEBUG) {
			System.out.println("X");
			print(x);
			System.out.println("Y Prime");
			print(yPrime);
		}

		for (int i = 0; i < nPoints; i++) {
			x[i][0] -= x0;
			x[i][1] -= x1;
			x[i][2] -= x2;
			yPrime[0][i] -= y0;
			yPrime[1][i] -= y1;
			yPrime[2][i] -= y2;
		}

		if (DEBUG) {
			System.out.println("X Centered");
			print(x);
			System.out.println("Y Prime Centered");
			print(yPrime);
		}

		// Find Y'X
		if (weights != null)
			weightedProduct(yPrime, x, weights, yPrimeX);
		else
			product(yPrime, x, yPrimeX);

		// Do SVD to get the rotation matrix

		if (DEBUG) {
			System.out.println("Product matrix prior to svd");
			print(yPrimeX);
		}

		Svd svd = new Svd(yPrimeX);
		svd.solve();
		w = svd.getSingularValues();
		v = svd.getV();
		double[][] u = svd.getU();

		if (DEBUG) {
			System.out.println("U post svd");
			print(u);
			System.out.println("V post svd");
			print(v);
			System.out.println("Singular values");
			System.out.println(w[0] + " " + w[1] + " " + w[2]);
		}

		// Rotation Matrix is VU'
		transpose(u, uPrime);
		product(v, uPrime, rotMat);
		double det = determinant(rotMat);
		if (DEBUG)
			System.out.println("Determinant of U'V is " + det);

		if (det < .0) {
			// If the determinant is 1 then the matrix is only
			// rotation. Otherwise we look through the singular values
			// (they are not ordered) and change the sign of the row
			// in U' that corresponds to the smallest singular value.

			int minRow = 0;
			double minSig = w[0];
			for (int i = 1; i < 3; i++)
				if (w[i] < minSig) {
					minSig = w[i];
					minRow = i;
				}

			// Modify U'
			uPrime[minRow][0] *= -1;
			uPrime[minRow][1] *= -1;
			uPrime[minRow][2] *= -1;

			product(v, uPrime, rotMat);

			if (DEBUG) {
				det = determinant(rotMat);
				System.out.println("Determinant of modified U'V is " + det);
			}
		}

		if (DEBUG) {
			System.out.println("Rotation matrix (U'V)");
			print(rotMat);
		}

		// Transformation in a 4x4 tranformation matrix
		zero(rotate);
		rotate[3][3] = 1.0;

		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				rotate[i][j] = rotMat[i][j];

		// In centers X
		// Out moves X into Y
		in[0][0] = in[1][1] = in[2][2] = in[3][3] = 1;
		in[3][0] = -x0;
		in[3][1] = -x1;
		in[3][2] = -x2;
		out[0][0] = out[1][1] = out[2][2] = out[3][3] = 1;
		out[3][0] = y0;
		out[3][1] = y1;
		out[3][2] = y2;

		product(in, rotate, product);
		// matrix is full transformation
		product(product, out, matrix);

		if (DEBUG) {
			System.out.println("Final transformation ");
			print(matrix);
		}
	}

}
