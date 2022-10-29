package com.cairn.common.stats;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.log4j.Logger;

/**
 * A base class for performing simplex minimization using commons math.
 * 
 * @author Gareth Jones
 * 
 */
public abstract class Simplex {
	private NelderMeadSimplex nelderMeadSimplex;
	private SimplexOptimizer simplexOptimizer;
	private double steps[], startingPoint[];
	private int maxEvals = 1000;
	private double[] minimum;
	private double[] currentValues;
	private static final Logger logger = Logger.getLogger(Simplex.class);

	static {
		// logger.setLevel(Level.DEBUG);
	}

	/**
	 * Performs optimization. The step sizes and starting point must have been
	 * previously set.
	 * 
	 */
	public void optimize() {
		// see NR in C 10.4 for information on Nelder Mead simplex
		nelderMeadSimplex = new NelderMeadSimplex(steps);
		simplexOptimizer = new SimplexOptimizer(1e-10, 1e-30);

		InitialGuess initialGuess = new InitialGuess(startingPoint);
		MaxEval maxEval = new MaxEval(maxEvals);
		MultivariateFunction mvFunction = new MultivariateFunction() {

			@Override
			public double value(double[] point) throws IllegalArgumentException {
				return calculateValue(point);
			}
		};
		ObjectiveFunction objectiveFunction = new ObjectiveFunction(mvFunction);
		boolean failed = false;
		currentValues = new double[steps.length];
		try {
			PointValuePair solution = simplexOptimizer.optimize(nelderMeadSimplex,
					maxEval, objectiveFunction, GoalType.MINIMIZE, initialGuess);
			minimum = solution.getPoint();
		} catch (TooManyEvaluationsException e) {
			logger.debug("failed to converge in " + simplexOptimizer.getEvaluations()
					+ " evaluations");
			failed = true;
		} catch (ConvergenceException e) {
			logger.debug("Convergence exception in " + simplexOptimizer.getEvaluations()
					+ " evaluations: " + e);
			failed = true;
		}
		if (failed) {
			minimum = currentValues;
		}
		int nEvals = simplexOptimizer.getEvaluations();
		int nIterations = simplexOptimizer.getIterations();

		logger.debug(
				"Performed " + nEvals + " evaluations, " + nIterations + " iterations");
	}

	/**
	 * @return The minimum point.
	 */
	public double[] getMinimum() {
		return minimum;
	}

	/**
	 * Stores the current position, so we can at least return where we've go to
	 * if there is an exception.
	 * 
	 * @param position
	 * @return
	 */
	private double calculateValue(double[] position) {
		for (int i = 0; i < position.length; i++) {
			currentValues[i] = position[i];
		}
		return value(position);
	}

	/**
	 * Subclasses need to implement this for their objective function.
	 * 
	 * @param position
	 * @return
	 */
	protected abstract double value(double[] position);

	/**
	 * @param steps
	 *            the steps to set
	 */
	public void setSteps(double[] steps) {
		this.steps = steps;
	}

	/**
	 * @param startingPoint
	 *            the starting point to set
	 */
	public void setStartingPoint(double[] startingPoint) {
		this.startingPoint = startingPoint;
	}

	/**
	 * @param maxEvals
	 *            the maxEvals to set
	 */
	public void setMaxEvals(int maxEvals) {
		this.maxEvals = maxEvals;
	}

}
