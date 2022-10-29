package com.cairn.gape.chromosome;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Logger;

import com.cairn.gape.ga.GaSupervisor;

/**
 * @author Gareth Jones
 * 
 */
public abstract class BinaryAndIntegerChromosome implements Chromosome {

	protected volatile BinaryStringChromosome binaryString;

	protected volatile IntegerStringChromosome integerString;

	private static final ThreadLocal<Boolean> allocated = ThreadLocal
			.withInitial(() -> false);

	private static ThreadLocal<String> threadName = new ThreadLocal<>();

	private final int chromosomeNo;

	private static AtomicInteger counter = new AtomicInteger(0);

	protected static final List<BinaryAndIntegerChromosome> chromosomeList = new ArrayList<>();

	protected boolean free;

	protected GaSupervisor problem;

	protected double binaryOperationProbability = 0.5;

	private static final Logger logger = Logger
			.getLogger(BinaryAndIntegerChromosome.class);

	/**
	 * Empty constructor for a new chromosome. Increases chromosome count.
	 * 
	 */
	public BinaryAndIntegerChromosome() {
		chromosomeNo = counter.incrementAndGet();
	}

	@Override
	public void copyGene(Chromosome c2) {

		BinaryAndIntegerChromosome c = (BinaryAndIntegerChromosome) c2;
		if (binaryString != null)
			binaryString.copyGene(c.binaryString);
		integerString.copyGene(c.integerString);

	}

	/*
	 * Returns a spare chromosome from the pool. This class maintains a pool of
	 * chromosomes to prevent the overheard of object creation and deletion.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.chromosome.Chromosome#create(com.cairn.gape.
	 * ga.GaSupervisor)
	 */
	@Override
	public Chromosome create(GaSupervisor problemHandle) {

		this.problem = problemHandle;
		if (!isAllocated()) {
			allocate(problemHandle);
			allocated.set(true);
		}
		BinaryAndIntegerChromosome c = getFreedChromosome();
		c.initialize();
		return c;
	}

	@Override
	public Chromosome create(GaSupervisor problemHandle, Object initialData) {
		throw new RuntimeException(
				"Method is create(GaSupervisor, initialData is not available for SuperpositionChromosome");
	}

	@Override
	public abstract double distance(Chromosome c);

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.chromosome.Chromosome#equals(com.cairn.gape.
	 * chromosome.Chromosome) Two chromosomes are equal if their binary and
	 * integer stings match.
	 */
	@Override
	public boolean equals(Chromosome c2) {

		BinaryAndIntegerChromosome c = (BinaryAndIntegerChromosome) c2;
		if (!integerString.equals(c.integerString))
			return false;
		if (binaryString == null)
			return true;
		if (!binaryString.equals(c.binaryString))
			return false;
		return true;
	}

	@Override
	public abstract String fitnessInfo();

	/**
	 * Discards a chomosome. The object is kept on a freelist for reuse.
	 */
	@Override
	public void freeChromosome() {
		synchronized (chromosomeList) {
			int no = chromosomeList.size();
			logger.debug("freeing chromosome no " + no + " id " + chromosomeNo);
			chromosomeList.add(this);
			free = true;
		}
	}

	/**
	 * Gets a chromosome off the free list.
	 * 
	 * @return
	 */
	static public BinaryAndIntegerChromosome getFreedChromosome() {
		synchronized (chromosomeList) {
			int no = chromosomeList.size() - 1;
			BinaryAndIntegerChromosome c = chromosomeList.remove(no);
			c.free = false;
			logger.debug("getting chromosome " + no + " ID " + c.chromosomeNo);
			return c;
		}
	}

	@Override
	public String geneInfo() {
		String rtn = "";
		if (binaryString != null)
			rtn += binaryString.geneInfo();
		rtn += integerString.geneInfo();
		return rtn;
	}

	static boolean isAllocated() {
		String name = Thread.currentThread().getName();
		// check that thread id has not been reused
		if (!StringUtils.equals(name, threadName.get())) {
			threadName.set(name);
			allocated.set(false);
			return false;
		}

		return allocated.get();
	}

	public static void setAllocated(boolean allocated) {
		BinaryAndIntegerChromosome.allocated.set(allocated);
	}

	/**
	 * Returns the chromosome number
	 */
	@Override
	public int getChromosomeId() {
		return chromosomeNo;
	}

	/**
	 * Mutates this chromosome. performs mutation on the binary or integer
	 * string with equal probability.
	 */
	public void mutate() {
		if (binaryString != null && binaryString.getNBits() > 0
				&& problem.normalRand() < binaryOperationProbability) {
			binaryString.mutate();
		} else {
			integerString.mutate();
		}
	}

	/**
	 * Crossover operator. Either does crossover on the binary string or full
	 * mixing on the integer string with equal probability
	 * 
	 * @param parent2
	 * @param child1
	 * @param child2
	 */
	public void crossover(BinaryAndIntegerChromosome parent2,
			BinaryAndIntegerChromosome child1, BinaryAndIntegerChromosome child2) {

		this.copyGene(child1);
		parent2.copyGene(child2);

		if (binaryString != null && binaryString.getNBits() > 0
				&& problem.normalRand() < binaryOperationProbability) {
			binaryString.crossover(parent2.binaryString, child1.binaryString,
					child2.binaryString);
		} else {
			integerString.fullMixing(parent2.integerString, child1.integerString,
					child2.integerString);
		}
	}

	@Override
	public Chromosome getChromosome() {
		return getFreedChromosome();
	}

	@Override
	public abstract double getFitness();

	@Override
	public abstract boolean ok();

	@Override
	public abstract double rebuild();

	@Override
	public abstract boolean sameNiche(Chromosome c);

	/**
	 * Initializes binary and integer strings to random values.
	 * 
	 */
	protected void initialize() {
		if (binaryString != null)
			binaryString.initialize();
		integerString.initialize();
	}

	/**
	 * Creates a pool of chromosomes in ChromosomeList
	 * 
	 * @param problemHandle
	 */
	public abstract void allocate(GaSupervisor problemHandle);

	/**
	 * @return the binaryString
	 */
	public BinaryStringChromosome getBinaryString() {
		return binaryString;
	}

	/**
	 * @return the integerString
	 */
	public IntegerStringChromosome getIntegerString() {
		return integerString;
	}

	/**
	 * @return the problem
	 */
	public GaSupervisor getProblem() {
		return problem;
	}

}
