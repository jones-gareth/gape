package com.cairn.gape;

import com.cairn.gape.chromosome.BinaryGaChromosome;
import com.cairn.gape.chromosome.BinaryStringCrossover;
import com.cairn.gape.chromosome.BinaryStringMutation;
import com.cairn.gape.ga.BaseSupervisor;
import com.cairn.gape.ga.IslandModel;
import com.cairn.gape.ga.Operator;
import com.cairn.gape.ga.Population;

/**
 * @author Gareth Jones
 * 
 *         Test function for binary string GA. Function binary f6 from Lawrence
 *         Davis, editor. Handbook of Genetic Algorithms. Van Nostrand Reinhold,
 *         New York, 1991.
 * 
 *         This function reaches its global maximum of 1.0 at the point
 *         (0.0,0.0). For this problem, the GA's chromosomes are strings of 44
 *         bits; 22 bits represent each of x and y in the range -100.0 to
 *         +100.0.
 * 
 * @see BinaryGaChromosome
 * 
 */
public class BinaryGa extends BaseSupervisor {
	// private static final Logger logger= Logger.getLogger(BinaryGa.class);

	static {
		// logger.setLevel(Level.DEBUG);
	}
	// total number of GA operations (crossover mutation and migration) applied.
	private final int nIterations = 50000;

	// number of islands in model
	private final int nIslands = 5;

	// population size- no of individuals in an island
	private final int popSize = 100;

	// selection pressure for rank-based normalized selection
	private final double selectPressure = 1.01;

	// this is the chromosome class for the binary f6 function. You can put your
	// own class here
	String chromosomeClass = "com.cairn.gape.chromosome.BinaryGaChromosome";

	/**
	 * @param args
	 */

	public static void main(String[] args) {
		BinaryGa ga = new BinaryGa();
		ga.run();
	}

	/**
	 * Runs the GA
	 */
	public void run() {
		init();

		// If you just want to use a single population use
		// LinkedPopLinearSel rather than the IslandModel.

		// Population pop = new LinkedPopLinearSel(this);

		// IslandModel handles GA population. 5 is % chance of migration
		Population pop = new IslandModel(this, nIslands, 5);

		// set selection pressure
		pop.setSelectPressure(selectPressure);

		// set population size
		pop.setSize(popSize);

		// create population - pass class name for chromosome
		pop.create(chromosomeClass, this);
		// out.print(pop.popInfo());

		// create genetic operators of crossover and mutation
		Operator mutationOp = new BinaryStringMutation(this, chromosomeClass);
		Operator crossoverOp = new BinaryStringCrossover(this, chromosomeClass);
		pop.addOperator(mutationOp, 45);
		pop.addOperator(crossoverOp, 45);

		// run the GA
		for (int i = 0; i < nIterations; i++)
			pop.iterate();

		// out.print(pop.popInfo());
		// best solution so far
		BinaryGaChromosome best = (BinaryGaChromosome) pop.getBest();
		infoMessageln("\nBest: " + best.fitnessInfo());

	}
}
