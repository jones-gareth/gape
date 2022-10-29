package com.cairn.gape.chromosome;

import java.text.NumberFormat;

import org.apache.commons.math3.util.FastMath;

/**
 * Test chromosome for binary string GA. Function binary f6 from Lawrence Davis,
 * editor. Handbook of Genetic Algorithms. Van Nostrand Reinhold, New York,
 * 1991.
 * 
 * This function reaches its global maximum of 1.0 at the point (0.0,0.0). For
 * this problem, the GA's chromosomes are strings of 44 bits; 22 bits represent
 * each of x and y in the range -100.0 to +100.0.
 * 
 * @author Gareth Jones
 * 
 */
public class BinaryGaChromosome extends BinaryStringChromosome {
	private boolean hasFitness = false;

	// Use 22 bits to encode each double.
	static final int N_BITS_PER_DOUBLE = 22;

	// Maximum value
	static final double max = FastMath.pow(2, N_BITS_PER_DOUBLE);

	static private final NumberFormat nf = NumberFormat.getInstance();
	static {
		nf.setMinimumFractionDigits(5);
	}

	private double x, y, fitness;

	/**
	 * Constructor. Set size to store 2 variables.
	 */
	public BinaryGaChromosome() {
		nBits = N_BITS_PER_DOUBLE * 2;
	}

	/**
	 * @param pos
	 * @return double value encoded on binary string starting at pos
	 */
	private double getDoubleValue(int pos) {
		// posToInt gets an integer value at a position- to check a bit use
		// "bitSet"
		int intVal = posToInt(pos, N_BITS_PER_DOUBLE);
		double frac = (intVal) / (max);
		return frac * 200 - 100;
	}

	/*
	 * Decodes the chromosome and generates x and y. Fitness is binary f6
	 * 
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.ga.BinaryStringChromosome#getFitness()
	 */
	@Override
	public double getFitness() {
		if (hasFitness)
			return fitness;

		x = getDoubleValue(0);
		y = getDoubleValue(N_BITS_PER_DOUBLE);
		double sqr = x * x + y * y;
		double a = Math.sin(Math.sqrt(sqr));
		double b = 1.0 + 0.001 * sqr;

		fitness = 0.5 - (a * a - 0.5) / (b * b);

		// System.out.println(fitnessInfo());
		hasFitness = true;
		return fitness;
	}

	/*
	 * Returns information about the chromosome
	 * 
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.ga.BinaryStringChromosome#fitnessInfo()
	 */
	@Override
	public String fitnessInfo() {
		return geneInfo() + "x " + nf.format(x) + " y " + nf.format(y) + " Fit "
				+ nf.format(fitness) + "\n";
	}

	/*
	 * Returns an empty chromosome, used for the migration operator.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.IntegerStringChromosome#getChromosome()
	 */
	@Override
	public Chromosome getChromosome() {
		BinaryGaChromosome c = new BinaryGaChromosome();
		c.createEmpty();
		c.setProblem(problem);
		return c;
	}

	/*
	 * Implement this if you need to support "failed" chromsomes
	 * 
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.ga.BinaryStringChromosome#ok()
	 */
	@Override
	public boolean ok() {
		return true;
	}

	/*
	 * Implement this is you want to define niching.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.ga.BinaryStringChromosome#sameNiche(com.cairn
	 * .gape.ga.Chromosome)
	 */
	@Override
	public boolean sameNiche(Chromosome c) {
		return false;
	}

	/*
	 * This doesn't do anything here. We can use it to recycle Chromsome objects
	 * if object churn becomes an issue.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.ga.BinaryStringChromosome#freeChromosome()
	 */
	@Override
	public void freeChromosome() {
		;
	}

	/*
	 * We need this for the migration and scaling operator. Just recalculates
	 * the fitness.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.ga.BinaryStringChromosome#rebuild()
	 */
	@Override
	public double rebuild() {
		hasFitness = false;
		return getFitness();
	}

}
