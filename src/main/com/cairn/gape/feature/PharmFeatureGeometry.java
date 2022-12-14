package com.cairn.gape.feature;

import com.cairn.common.utils.Coord;

/**
 * This class is used to store the Geometry definition of a pharmacophore
 * feature. This exists in extension to the Feature class as different feature
 * types share the same geometry (e.g. Sp1 acceptors and h-donor are both
 * defined by a vector.
 * 
 * @author Gareth Jones
 * 
 */
public abstract class PharmFeatureGeometry {
	enum Geometry {
		SPHERE, VECTOR, ARC, CONE, MULTI_VECTOR
	}

	protected Geometry geometry;

	final static boolean DEBUG = false;

	/**
	 * Convert the geometry to a string.
	 * 
	 * @return
	 */
	abstract public String summary();

	/**
	 * Returns a feature geometry by parsing a summary string. The string should
	 * have been created by the summary method. So you can save the geometry to
	 * a string with summary, then convert it back using this method.
	 * 
	 * @param summary
	 * @return
	 * 
	 * @see #summary()
	 */
	static PharmFeatureGeometry parseString(String summary) {
		if (summary.startsWith("{VECTOR"))
			return new VectorPharmFeatureGeometry(summary);
		if (summary.startsWith("{MULTI_VECTOR"))
			return new MultiVectorPharmFeatureGeometry(summary);
		else if (summary.startsWith("{SPHERE"))
			return new SpherePharmFeatureGeometry(summary);
		else if (summary.startsWith("{ARC"))
			return new ArcPharmFeatureGeometry(summary);
		else if (summary.startsWith("{CONE"))
			return new ConePharmFeatureGeometry(summary);
		throw new RuntimeException("PharmFeatureGeometry: parseString: can't parse "
				+ summary);
	}

	/**
	 * Gets the number of points used by this feature.
	 * 
	 * @return
	 */
	abstract public int getNPoints();

	/**
	 * Returns the coordinates of a point
	 * 
	 * @param no
	 * @return
	 */
	abstract public double[] getPoint(int no);

	/**
	 * Performs an approximate validation using vector lengths.
	 * 
	 * @return
	 */
	abstract public boolean check();

	public Geometry getGeometry() {
		return geometry;
	}

	/**
	 * Checks that the distance between any two points in the parmacophore is
	 * not crazy.
	 * 
	 * @param point1
	 * @param point2
	 * @return
	 */
	public boolean checkVectorLength(double[] point1, double point2[]) {
		double len = Coord.distance(point1, point2);
		// the maximum distance should be 6, which is the normal distance for
		// aromatic rings.
		return len < 7;
	}
}

/**
 * Represents a pharmacophore feature that is a vector. For example, a donor
 * hydrogen or an sp1 acceptor.
 * 
 */
class VectorPharmFeatureGeometry extends PharmFeatureGeometry {
	double points[][];

	/**
	 * Constructor - recreates class from a string generated by summary method
	 * 
	 * @param summary
	 * 
	 * @see #summary()
	 */
	VectorPharmFeatureGeometry(String summary) {
		summary = summary.trim();
		String orig = summary;

		if (!summary.startsWith("{VECTOR "))
			throw new RuntimeException("VectorPharmFeatureGeometry: can't parse summary");
		if (!summary.endsWith("}"))
			throw new RuntimeException("VectorPharmFeatureGeometry: can't parse summary");
		summary = summary.substring(8, summary.length() - 1);

		String coordStr[] = summary.split(" ");
		double point1[] = new double[4];
		double point2[] = new double[4];
		Coord.toPoint(coordStr[0], point1);
		Coord.toPoint(coordStr[1], point2);

		points = new double[][] { point1, point2 };
		geometry = Geometry.VECTOR;

		String check = summary();
		if (DEBUG) {
			System.out.println("Orig:  " + orig);
			System.out.println("check: " + check);
		}
		if (!check.equals(orig))
			throw new RuntimeException(
					"VectorPharmFeatureGeometry: can't reporduce original");

		assert check() : "Invalid pharmacophore vector";
	}

	/**
	 * Constructor. The feature directionality is from point 1 to point 2.
	 *
	 * @param point1
	 * @param point2
	 */
	VectorPharmFeatureGeometry(double point1[], double point2[]) {
		points = new double[][] { point1, point2 };
		geometry = Geometry.VECTOR;

		assert check() : "Invalid pharmacophore vector";
	}

	@Override
	public String summary() {
		return "{VECTOR " + Coord.toString(points[0]) + " " + Coord.toString(points[1])
				+ "}";
	}

	@Override
	public int getNPoints() {
		return 2;
	}

	@Override
	public double[] getPoint(int no) {
		return points[no];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.feature.PharmFeatureGeometry#check()
	 */
	@Override
	public boolean check() {
		return checkVectorLength(points[0], points[1]);
	}

}

/**
 * Represents a pharmacophore feature that comprises multiple vectors, each
 * vector sharing one point. For example, an acceptor with two or three lone
 * pairs. For n vectors class ecnodes n+1 points. The first point is the center
 * and the next n are the end points of the vectors
 */
class MultiVectorPharmFeatureGeometry extends PharmFeatureGeometry {
	double points[][];
	int nPoints;

	/**
	 * Constructor - recreates class from a string generated by summary method
	 * 
	 * @param summary
	 * 
	 * @see #summary()
	 */
	MultiVectorPharmFeatureGeometry(String summary) {
		summary = summary.trim();
		String orig = summary;

		if (!summary.startsWith("{MULTI_VECTOR "))
			throw new RuntimeException("VectorPharmFeatureGeometry: can't parse summary");
		if (!summary.endsWith("}"))
			throw new RuntimeException("VectorPharmFeatureGeometry: can't parse summary");
		summary = summary.substring(14, summary.length() - 1);

		String coordStr[] = summary.split(" ");
		int nVectors = Integer.parseInt(coordStr[0]);
		nPoints = nVectors + 1;
		points = new double[nPoints][];
		double center[] = new double[4];
		Coord.toPoint(coordStr[1], center);
		points[0] = center;

		for (int i = 0; i < nVectors; i++) {
			double point[] = new double[4];
			Coord.toPoint(coordStr[i + 2], point);
			points[i + 1] = point;
		}

		geometry = Geometry.MULTI_VECTOR;

		String check = summary();
		if (DEBUG) {
			System.out.println("Orig:  " + orig);
			System.out.println("check: " + check);
		}
		if (!check.equals(orig))
			throw new RuntimeException(
					"VectorPharmFeatureGeometry: can't reporduce original");

		assert check() : "Invalid pharmacophore vector";
	}

	/**
	 * Constructor. Each vector goes from the center out.
	 * 
	 * @param point1
	 * @param point2
	 */
	MultiVectorPharmFeatureGeometry(int nVectors, double center[], double _points[][]) {
		nPoints = nVectors + 1;
		points = new double[nPoints][];
		points[0] = center;
		for (int i = 0; i < nVectors; i++) {
			points[i + 1] = _points[i];
		}

		geometry = Geometry.MULTI_VECTOR;
		assert check() : "Invalid pharmacophore vector";
	}

	@Override
	public String summary() {
		int nVectors = nPoints - 1;
		String rtn = "{MULTI_VECTOR " + nVectors;
		for (int i = 0; i < nPoints; i++)
			rtn += " " + Coord.toString(points[i]);
		rtn += "}";
		return rtn;
	}

	@Override
	public int getNPoints() {
		return nPoints;
	}

	@Override
	public double[] getPoint(int no) {
		return points[no];
	}

	/**
	 * @return the number of vectors
	 */
	public int getNVectors() {
		return nPoints - 1;
	}

	/**
	 * @return the center of the set of vectors
	 */
	public double[] getCenter() {
		return points[0];
	}

	/**
	 * @param no
	 * @return the end point of vector no.
	 */
	public double[] getEnd(int no) {
		return points[no + 1];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.feature.PharmFeatureGeometry#check()
	 */
	@Override
	public boolean check() {
		for (int i = 1; i < points.length; i++) {
			if (!checkVectorLength(points[0], points[i])) {
				return false;
			}
		}
		return true;
	}

}

/**
 * Represents a pharmacophore feature that has no directionality. For example,
 * certain acceptors and hydrophobic centers.
 * 
 */
class SpherePharmFeatureGeometry extends PharmFeatureGeometry {
	double radius;

	double point[];

	/**
	 * Constructor - recreates class from a string generated by summary method
	 * 
	 * @param summary
	 * 
	 * @see #summary()
	 */
	SpherePharmFeatureGeometry(String summary) {
		summary = summary.trim();
		String orig = summary;
		if (!summary.startsWith("{SPHERE "))
			throw new RuntimeException("SpherePharmFeatureGeometry: can't parse "
					+ summary);
		if (!summary.endsWith("}"))
			throw new RuntimeException("SpherePharmFeatureGeometry: can't parse "
					+ summary);
		summary = summary.substring(8, summary.length() - 1);

		String coordStr[] = summary.split(" ");
		point = new double[4];
		Coord.toPoint(coordStr[0], point);

		if (!coordStr[1].equals("RADIUS"))
			throw new RuntimeException("SpherePharmFeatureGeometry: can't parse "
					+ summary);
		radius = Double.parseDouble(coordStr[2]);
		geometry = Geometry.SPHERE;

		String check = summary();
		if (DEBUG) {
			System.out.println("Orig:  " + orig);
			System.out.println("check: " + check);
		}
		if (!check.equals(orig))
			throw new RuntimeException(
					"SpherePharmFeatureGeometry: can't reporduce original");
	}

	/**
	 * Construcutor for a feature center and radius,
	 * 
	 * @param _point
	 * @param _radius
	 */
	SpherePharmFeatureGeometry(double _point[], double _radius) {
		point = _point;
		radius = _radius;
		geometry = Geometry.SPHERE;
	}

	@Override
	public String summary() {
		return "{SPHERE " + Coord.toString(point) + " RADIUS " + radius + "}";
	}

	@Override
	public int getNPoints() {
		return 1;
	}

	@Override
	public double[] getPoint(int no) {
		if (no == 0)
			return point;
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.feature.PharmFeatureGeometry#check()
	 */
	@Override
	public boolean check() {
		return true;
	}

}

/**
 * Represents a feature that is a arc- typically an sp2 acceptor with two lone
 * pairs.
 * 
 */
class ArcPharmFeatureGeometry extends PharmFeatureGeometry {
	double points[][];

	/**
	 * Constructor - recreates class from a string generated by summary method
	 * 
	 * @param summary
	 * 
	 * @see #summary()
	 */
	ArcPharmFeatureGeometry(String summary) {
		summary = summary.trim();
		String orig = summary;
		if (!summary.startsWith("{ARC "))
			throw new RuntimeException("ArcPharmFeatureGeometry: can't parse summary");
		if (!summary.endsWith("}"))
			throw new RuntimeException("ArcPharmFeatureGeometry: can't parse summary");
		summary = summary.substring(5, summary.length() - 1);

		String coordStr[] = summary.split(" ");
		double point1[] = new double[4];
		double point2[] = new double[4];
		double point3[] = new double[4];
		Coord.toPoint(coordStr[0], point1);
		Coord.toPoint(coordStr[1], point2);
		Coord.toPoint(coordStr[2], point3);

		points = new double[][] { point1, point2, point3 };
		geometry = Geometry.ARC;

		String check = summary();
		if (DEBUG) {
			System.out.println("Orig:  " + orig);
			System.out.println("check: " + check);
		}
		if (!check.equals(orig))
			throw new RuntimeException(
					"ArcPharmFeatureGeometry: can't reporduce original");

		assert check() : "Invalid pharmacophore vector";
	}

	/**
	 * The feature center is at point1 and point2 and point3 define the extent
	 * of the arc.
	 * 
	 * @param point1
	 * @param point2
	 * @param point3
	 */
	ArcPharmFeatureGeometry(double point1[], double point2[], double point3[]) {
		points = new double[][] { point1, point2, point3 };
		geometry = Geometry.ARC;

		assert check() : "Invalid pharmacophore vector";
	}

	@Override
	public String summary() {
		return "{ARC " + Coord.toString(points[0]) + " " + Coord.toString(points[1])
				+ " " + Coord.toString(points[2]) + "}";
	}

	@Override
	public int getNPoints() {
		return 3;
	}

	@Override
	public double[] getPoint(int no) {
		return points[no];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.feature.PharmFeatureGeometry#check()
	 */
	@Override
	public boolean check() {
		return checkVectorLength(points[0], points[1])
				&& checkVectorLength(points[1], points[2]);
	}

}

/**
 * Represents a feature that is a solid cone- typically an sp3 acceptor with
 * three lone pairs.
 * 
 */
class ConePharmFeatureGeometry extends PharmFeatureGeometry {
	double points[][];

	/**
	 * Constructor - recreates class from a string generated by summary method
	 * 
	 * @param summary
	 * 
	 * @see #summary()
	 */
	ConePharmFeatureGeometry(String summary) {
		summary = summary.trim();
		String orig = summary;

		if (!summary.startsWith("{CONE "))
			throw new RuntimeException("ConePharmFeatureGeometry: can't parse summary");
		if (!summary.endsWith("}"))
			throw new RuntimeException("ConePharmFeatureGeometry: can't parse summary");
		summary = summary.substring(6, summary.length() - 1);

		String coordStr[] = summary.split(" ");
		double point1[] = new double[4];
		double point2[] = new double[4];
		double point3[] = new double[4];
		Coord.toPoint(coordStr[0], point1);
		Coord.toPoint(coordStr[1], point2);
		Coord.toPoint(coordStr[2], point3);

		points = new double[][] { point1, point2, point3 };
		geometry = Geometry.CONE;

		String check = summary();
		if (DEBUG) {
			System.out.println("Orig:  " + orig);
			System.out.println("check: " + check);
		}
		if (!check.equals(orig))
			throw new RuntimeException(
					"ConePharmFeatureGeometry: can't reporduce original");

		assert check() : "Invalid pharmacophore vector";
	}

	/**
	 * point1 is the center of the cone and point2 and point3 are two points on
	 * the outside edge of the cone base such that the center of the circular
	 * base of the cone lies between points 2 and 3.
	 * 
	 * @param point1
	 * @param point2
	 * @param point3
	 */
	ConePharmFeatureGeometry(double point1[], double point2[], double point3[]) {
		points = new double[][] { point1, point2, point3 };
		geometry = Geometry.CONE;

		assert check() : "Invalid pharmacophore vector";

	}

	@Override
	public String summary() {
		return "{CONE " + Coord.toString(points[0]) + " " + Coord.toString(points[1])
				+ " " + Coord.toString(points[2]) + "}";
	}

	@Override
	public int getNPoints() {
		return 3;
	}

	@Override
	public double[] getPoint(int no) {
		return points[no];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.feature.PharmFeatureGeometry#check()
	 */
	@Override
	public boolean check() {
		return checkVectorLength(points[0], points[1])
				&& checkVectorLength(points[1], points[2]);
	}
}
