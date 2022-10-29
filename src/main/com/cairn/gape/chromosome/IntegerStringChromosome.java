package com.cairn.gape.chromosome;

import com.cairn.gape.ga.GaSupervisor;

/**
 * Class to represent integer string chromosomes.
 * 
 * @author Gareth Jones
 * 
 */
public class IntegerStringChromosome implements Chromosome {
	protected int length = 0, ranges[], values[] = null;

	// set true to allow "null" values in string. -1 is used for the null value.
	protected boolean allowNulls[], allowNullDefault = true;

	protected GaSupervisor problem;

	protected boolean allowSwitch = false;

	protected IntegerStringChromosome() {
		;
	}

	/**
	 * Constructor
	 * 
	 * @param p
	 *            GA problem
	 * @param l
	 *            length of chromosome
	 * @param r
	 *            ranges. r[p] is the maximum value for position p.
	 */
	public IntegerStringChromosome(GaSupervisor p, int l, int r[]) {
		length = l;
		ranges = r;
		problem = p;
	}

	// This has never been used
	@Override
	public Chromosome create(GaSupervisor problemHandle) {
		this.problem = problemHandle;
		createEmpty();
		return this;
	}

	// This has never been used
	@Override
	public Chromosome create(GaSupervisor problemHandle, Object initialData) {

		this.problem = problemHandle;
		createEmpty();

		int[] data = null;
		try {
			data = (int[]) initialData;
		} catch (ClassCastException ex) {
			throw new IllegalStateException(
					"IntegerStringChromosome: create initial data must be of type int[]");
		}

		if (data.length != length)
			throw new IllegalStateException(
					"IntegerStringChromosome: create initial data must be of length "
							+ length);

		for (int i = 0; i < length; i++) {
			values[i] = data[i];
		}

		return this;
	}

	@Override
	public double getFitness() {
		throw new RuntimeException("method not defined");
	}

	@Override
	public String fitnessInfo() {
		return "\n";
	}

	@Override
	public double distance(Chromosome c) {
		throw new RuntimeException("method not defined");
	}

	@Override
	public boolean sameNiche(Chromosome c) {
		throw new RuntimeException("method not defined");
	}

	@Override
	public boolean ok() {
		throw new RuntimeException("method not defined");
	}

	@Override
	public double rebuild() {
		throw new RuntimeException("method not defined");
	}

	public Chromosome getChromosome(Chromosome c) {
		throw new RuntimeException("method not defined");
	}

	@Override
	public void freeChromosome() {
		throw new RuntimeException("method not defined");
	}

	@Override
	public int getChromosomeId() {
		throw new RuntimeException("method not defined");
	}

	@Override
	public Chromosome getChromosome() {
		throw new RuntimeException("method not defined");
	}

	/**
	 * Initializes integer string with random values.
	 */
	public void initialize() {

		if (values == null)
			createEmpty();
		for (int i = 0; i < length; i++) {
			int bottom = allowNulls == null ? (allowNullDefault ? -1 : 0)
					: (allowNulls[i] ? -1 : 0);
			if (ranges[i] - 1 == bottom)
				values[i] = bottom;
			else
				values[i] = problem.randomInt(bottom, ranges[i] - 1);
		}
	}

	/**
	 * Allocates storage
	 */
	public void createEmpty() {
		values = new int[length];
	}

	/**
	 * Performs mutation. The per-integer mutation probability is 1/length. The
	 * operation is repeated until a mutation occurs. An integer is just mutated
	 * to another allowed value with equal probability.
	 */
	public void mutate() {
		double pMutate = 1.0 / length;
		boolean mutated = false;

		for (int i = 0; i < length; i++) {
			boolean allowNull = allowNulls == null ? allowNullDefault : allowNulls[i];
			if (allowNull && ranges[i] == 0)
				continue;
			int bottom = allowNull ? -1 : 0;
			if (problem.normalRand() < pMutate) {
				int newVal = ranges[i] - 2 == 0 ? bottom : problem.randomInt(bottom,
						ranges[i] - 2);
				if (newVal >= values[i])
					newVal++;
				values[i] = newVal;
				mutated = true;
			}
		}
		if (!mutated)
			mutate();
	}

	/**
	 * Performs two point crossover.
	 * 
	 * @param parent2
	 * @param child1
	 * @param child2
	 */
	public void crossover(IntegerStringChromosome parent2,
			IntegerStringChromosome child1, IntegerStringChromosome child2) {
		twoPointCrossover(parent2, child1, child2);
	}

	/**
	 * Performs two point crossover.
	 * 
	 * @param parent2
	 * @param child1
	 * @param child2
	 */
	public void twoPointCrossover(IntegerStringChromosome parent2,
			IntegerStringChromosome child1, IntegerStringChromosome child2) {

		// choose cross point sites
		int site1 = problem.randomInt(0, length);
		int site2 = problem.randomInt(0, length - 1);
		if (site2 >= site1)
			site2++;
		else {
			int n = site1;
			site1 = site2;
			site2 = n;
		}

		int p1[] = values;
		int p2[] = parent2.values;
		int c1[] = null, c2[] = null;

		boolean switchFlg = false;
		if (allowSwitch)
			switchFlg = problem.randomBoolean();

		if (switchFlg) {
			c1 = child2.values;
			c2 = child1.values;
		} else {
			c1 = child1.values;
			c2 = child2.values;
		}

		// create children before 1st cross point
		int i = 0;
		for (; i < site1; i++) {
			c1[i] = p1[i];
			c2[i] = p2[i];
		}

		// children between cross points
		for (; i < site2; i++) {
			c1[i] = p2[i];
			c2[i] = p1[i];
		}

		// children after 2nd cross point
		for (; i < length; i++) {
			c1[i] = p1[i];
			c2[i] = p2[i];
		}
	}

	/**
	 * Performs one point crossover.
	 * 
	 * @param parent2
	 * @param child1
	 * @param child2
	 */
	public void onePointCrossover(IntegerStringChromosome parent2,
			IntegerStringChromosome child1, IntegerStringChromosome child2) {
		// choose site
		int site = problem.randomInt(0, length - 1);

		// do child1 then child2 */
		int p1[] = values;
		int p2[] = parent2.values;
		int c1[] = null, c2[] = null;

		boolean switchFlg = false;
		if (allowSwitch)
			switchFlg = problem.randomBoolean();

		if (switchFlg) {
			c1 = child2.values;
			c2 = child1.values;
		} else {
			c1 = child1.values;
			c2 = child2.values;
		}

		// create children before cross point
		int i = 0;
		for (; i < site; i++) {
			c1[i] = p1[i];
			c2[i] = p2[i];
		}

		// children after cross point
		for (; i < length; i++) {
			c1[i] = p1[i];
			c2[i] = p2[i];
		}
	}

	/**
	 * Assumes that the strings primarily comprise dummy values. The children
	 * then comprise all non dummy values from the parents. Where both parents
	 * have non-dummy values set at the same position child1 takes the value
	 * from parent1 and child2 takes the value from parent2.
	 * 
	 * @param parent2
	 * @param child1
	 * @param child2
	 */
	public void fullMixing(IntegerStringChromosome parent2,
			IntegerStringChromosome child1, IntegerStringChromosome child2) {
		int p1[] = values;
		int p2[] = parent2.values;
		int c1[] = child1.values;
		int c2[] = child2.values;

		for (int i = 0; i < length; i++) {
			if (p1[i] == -1 && p2[i] == -1)
				c1[i] = c2[i] = -1;
			else if (p1[i] != -1 && p2[i] == -1)
				c1[i] = c2[i] = p1[i];
			else if (p1[i] == -1 && p2[i] != -1)
				c1[i] = c2[i] = p2[i];
			else {
				c1[i] = p1[i];
				c2[i] = p2[i];
			}
		}
	}

	/**
	 * Assumes that the strings primarily comprise dummy values. Proceed as for
	 * 1 point crossover on integer strings. However if the child is given a
	 * dummy value when the other parent has a non-dummy value at a position,
	 * then that value is copied to the child.
	 * 
	 * @param parent2
	 * @param child1
	 * @param child2
	 */
	public void fullMixingAndCross(IntegerStringChromosome parent2,
			IntegerStringChromosome child1, IntegerStringChromosome child2) {
		int p1[] = values;
		int p2[] = parent2.values;
		int c1[] = null;
		int c2[] = null;

		// select crossover sites
		int site = problem.randomInt(0, length);

		boolean switchFlg = false;
		if (allowSwitch)
			switchFlg = problem.randomBoolean();

		if (switchFlg) {
			c1 = child2.values;
			c2 = child1.values;
		} else {
			c1 = child1.values;
			c2 = child2.values;
		}

		// create child before cross point
		int i = 0;
		for (; i < site; i++) {
			if (p1[i] == -1 && p2[i] != -1)
				c1[i] = c2[i] = p2[i];
			else if (p1[i] != -1 && p2[i] == -1)
				c1[i] = c2[i] = p1[i];
			else {
				c2[i] = p2[i];
				c1[i] = p1[i];
			}
		}

		// child after cross point
		for (; i < length; i++) {
			if (p1[i] == -1 && p2[i] != -1)
				c1[i] = c2[i] = p2[i];
			else if (p1[i] != -1 && p2[i] == -1)
				c1[i] = c2[i] = p1[i];
			else {
				c2[i] = p1[i];
				c1[i] = p2[i];
			}
		}

	}

	/*
	 * Copies this integer string to another chromosome.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.ga.Chromosome#copyGene(com.cairn.gape.ga.Chromosome
	 * )
	 */
	@Override
	public void copyGene(Chromosome c2) {
		IntegerStringChromosome c = (IntegerStringChromosome) c2;
		if (c.values == null || c.length != length) {
			c.values = new int[length];
			c.length = length;
		}
		for (int i = 0; i < length; i++)
			c.values[i] = values[i];
	}

	/*
	 * Return true if both integer strings are equal.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.ga.Chromosome#equals(com.cairn.gape.ga.Chromosome
	 * )
	 */
	@Override
	public boolean equals(Chromosome c2) {
		IntegerStringChromosome c = (IntegerStringChromosome) c2;
		for (int i = 0; i < length; i++)
			if (c.values[i] != values[i])
				return false;
		return true;
	}

	/*
	 * Formats the integer string.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.ga.Chromosome#geneInfo()
	 */
	@Override
	public String geneInfo() {
		String rtn = "";
		for (int i = 0; i < length; i++)
			rtn += values[i] + " ";
		return rtn + "\n";
	}

	public int[] getValues() {
		return values;
	}

	public int getValues(int i) {
		return values[i];
	}

	public void setValues(int values[]) {
		this.values = values;
	}

	public void setValues(int i, int value) {
		values[i] = value;
	}

	public int[] getRanges() {
		return ranges;
	}

	public void setRanges(int[] ranges) {
		this.ranges = ranges;
	}

	public void setAllowNulls(boolean[] allowNulls) {
		this.allowNulls = allowNulls;
	}

	/**
	 * @return the allowNullDefault
	 */
	public boolean isAllowNullDefault() {
		return allowNullDefault;
	}

	/**
	 * @param allowNullDefault
	 *            the allowNullDefault to set
	 */
	public void setAllowNullDefault(boolean allowNullDefault) {
		this.allowNullDefault = allowNullDefault;
	}

	/**
	 * @return the allowNulls
	 */
	public boolean[] getAllowNulls() {
		return allowNulls;
	}

}
