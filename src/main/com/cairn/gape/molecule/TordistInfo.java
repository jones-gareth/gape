package com.cairn.gape.molecule;

import java.util.List;

import org.apache.commons.math3.util.FastMath;
import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;
import com.cairn.molecule.Atom;
import com.cairn.molecule.Molecule;
import com.cairn.molecule.RotatableBond;

/**
 * Class for mapping torsional distributions to rotatable Bonds
 * 
 * @author Gareth Jones
 * 
 */
public class TordistInfo {
	private volatile Atom atom1, atom2, atom3, atom4;

	private volatile int count[], nCuts;

	private volatile double weight, prob[], startCuts[], endCuts[], total, refAngle,
			torAngle;

	private volatile boolean reverse;

	private volatile Molecule molecule;

	// private RotatableBond rotatableBond;

	private volatile TorsionalDistributions distributions;

	static final boolean DEBUG = false;

	private static final Logger logger = Logger.getLogger(TordistInfo.class);

	public TordistInfo() {
		;
	}

	public TordistInfo(RotatableBond rb, TorsionalDistribution td, Molecule m) {
		molecule = m;
		// rotatableBond = rb;
		distributions = td.getDistributions();
		copyTordistInfo(td);
	}

	/**
	 * Copy relevant information from the matching torsional distribution'
	 * 
	 * @param tordist
	 */
	private void copyTordistInfo(TorsionalDistribution tordist) {
		atom1 = tordist.getAtom1().getMatchedAtom();
		atom2 = tordist.getAtom2().getMatchedAtom();
		atom3 = tordist.getAtom3().getMatchedAtom();
		atom4 = tordist.getAtom4().getMatchedAtom();

		if (logger.isDebugEnabled()) {
			int a1 = atom1.getNo() + 1;
			int a2 = atom2.getNo() + 1;
			int a3 = atom3.getNo() + 1;
			int a4 = atom4.getNo() + 1;
			logger.debug("Full torsion [" + a1 + " " + a2 + " " + a3 + " " + a4
					+ "] weight " + tordist.getWeight());
		}

		int nBins = tordist.getDistributions().getNBins();
		count = new int[nBins];
		prob = new double[nBins];
		for (int i = 0; i < count.length; i++)
			count[i] = tordist.getCount(i);
		weight = tordist.getWeight();

		distributions = tordist.getDistributions();
	}

	/**
	 * Initialise the structure.
	 */
	public void setupTordistInfo() {

		logger.debug(TorsionalDistribution.histogramInfo(count));

		// find the largest peak
		double max = 0;
		for (int i = 0; i < count.length; i++) {
			double val = count[i];
			if (val > max)
				max = val;
		}

		// apply any high energy cutoff, see HTML documentation for
		// the maths
		if (distributions.isRemoveHighEnergy()) {
			double ratio = FastMath.exp(distributions.getDeltaE() / .5898);
			if (DEBUG)
				System.out.println("Ratio to screen out torsions of energy "
						+ distributions.getDeltaE() + " is " + ratio);

			for (int i = 0; i < count.length; i++) {
				double val = count[i];
				double r = max / val;
				if (r > ratio)
					count[i] = 0;
			}
		}

		double sum = .0;
		for (int i = 0; i < count.length; i++) {
			double val = count[i];
			sum += val;
		}

		// Generate probabilities, this isn't complete now, I don't
		// know how to convert between the distribution and a
		// probability density for torsional angles.
		if (distributions.isEnergyProbs()) {
			double logsum = .0;
			for (int i = 0; i < count.length; i++) {
				if (count[i] > 0) {
					double val = count[i];
					logsum += Math.log(val);
				}
			}
			for (int i = 0; i < count.length; i++) {
				if (count[i] > 0) {
					double val = count[i];
					prob[i] = Math.log(val) / logsum;
				}
			}
		} else {
			for (int i = 0; i < count.length; i++) {
				double val = count[i];
				prob[i] = val / sum;
			}
		}

		// Convert the torsional histogram to a number of "cuts", with
		// start and end positions.
		startCuts = new double[count.length];
		endCuts = new double[count.length];
		boolean inside = false;
		double step = distributions.getGranularity();

		int i = 0;
		for (; i < count.length; i++) {
			if (count[i] > 0) {
				if (!inside) {
					inside = true;
					startCuts[nCuts] = i * step - 180.0;
				}
				total += step;
			} else {
				if (inside) {
					inside = false;
					endCuts[nCuts] = i * step - 180.0;
					nCuts++;
				}
			}
		}
		if (inside) {
			inside = false;
			endCuts[nCuts] = i * step - 180.0;
			nCuts++;
		}

		if (DEBUG) {
			System.out.print(TorsionalDistribution.histogramInfo(count));
			for (int j = 0; j < nCuts; j++)
				System.out.println("[" + startCuts[j] + " " + endCuts[j] + "]");
			System.out.println("Total " + total);
			System.out.println();
		}
		refAngle = getTorAngle();
		// Find a reference offset angle for this match distribution.
		// When the angle of rotation about the bond is 0, then the
		// torsional dirtribution angle is the reference angle.
		if (DEBUG) {
			double deg = refAngle * 180.0 / Math.PI;
			System.out.println("Reference torsion angle is " + deg);
		}
	}

	public double getTorAngle() {
		List<double[]> coords = molecule.getCoords();
		double point1[] = coords.get(atom1.getNo());
		double point2[] = coords.get(atom2.getNo());
		double point3[] = coords.get(atom3.getNo());
		double point4[] = coords.get(atom4.getNo());
		return Coord.torsion(point1, point2, point3, point4);
	}

	/**
	 * Decode frac is a number between 0 and 1. This function finds the
	 * distribution angle associated with this fraction. The reference angle for
	 * the rotatable bond is then used to return the angle of rotation that
	 * should be applied to the bond.
	 * 
	 * @param decodeFrac
	 * @return
	 */
	public double getTordistRotationAngle(double decodeFrac) {
		torAngle = total * decodeFrac;
		// Work through cuts to find torsional distribution angle
		for (int i = 0; i < nCuts; i++) {
			torAngle += startCuts[i];
			if (torAngle <= endCuts[i])
				break;
			else
				torAngle -= endCuts[i];
		}

		torAngle *= Math.PI / 180.0;
		// Use reference angle to determine angle of rotation.
		double angle = torAngle - refAngle;
		return angle;
	}

	public Atom getAtom1() {
		return atom1;
	}

	public Atom getAtom2() {
		return atom2;
	}

	public Atom getAtom3() {
		return atom3;
	}

	public Atom getAtom4() {
		return atom4;
	}

	public boolean isReverse() {
		return reverse;
	}

	public void setReverse(boolean reverse) {
		this.reverse = reverse;
	}

	/**
	 * @return the weight
	 */
	public double getWeight() {
		return weight;
	}

	/**
	 * @return the refAngle
	 */
	public double getRefAngle() {
		return refAngle;
	}

}
