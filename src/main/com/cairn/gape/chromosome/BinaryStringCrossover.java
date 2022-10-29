package com.cairn.gape.chromosome;

import com.cairn.gape.ga.GaSupervisor;
import com.cairn.gape.ga.Operator;

/**
 * Implements crossover on a Binary string.
 * 
 * @author Gareth Jones
 *
 */
public class BinaryStringCrossover implements Operator {
	private GaSupervisor gaSupervisor;

	private Class<?> chromosomeClass;

	private static final boolean DEBUG = false;

	private Chromosome children[] = new Chromosome[2];

	public BinaryStringCrossover(GaSupervisor gaSupervisor, String chromosomeClassString) {
		this.gaSupervisor = gaSupervisor;

		chromosomeClass = null;
		try {
			// get chromsome class
			chromosomeClass = Class.forName(chromosomeClassString);
		} catch (ClassNotFoundException ex) {
			throw new RuntimeException("No chromosome class " + chromosomeClassString);
		}

	}

	@Override
	public Chromosome[] apply(Chromosome[] parents) {
		BinaryStringChromosome child1, child2, parent1, parent2;

		try {
			child1 = (BinaryStringChromosome) chromosomeClass.newInstance();
			child2 = (BinaryStringChromosome) chromosomeClass.newInstance();
		} catch (IllegalAccessException ex1) {
			throw new RuntimeException("Illegal access for creating chromosome "
					+ chromosomeClass);
		} catch (InstantiationException ex2) {
			throw new RuntimeException("Exception on instantiating " + chromosomeClass);
		}

		child1.createEmpty();
		child1.setProblem(gaSupervisor);
		child2.createEmpty();
		child2.setProblem(gaSupervisor);

		parent1 = (BinaryStringChromosome) parents[0];
		parent2 = (BinaryStringChromosome) parents[1];

		if (DEBUG) {
			System.out.println("Doing Crossover");
			System.out.print(parent1.geneInfo());
			System.out.print(parent2.geneInfo());
			System.out.println("-->");
		}
		parent1.crossover(parent2, child1, child2);
		if (DEBUG) {
			System.out.print(child1.geneInfo());
			System.out.print(child2.geneInfo());
		}

		children[0] = child1;
		children[1] = child2;
		return children;
	}

	@Override
	public int nChildren() {
		return 2;
	}

	@Override
	public int nParents() {
		return 2;
	}

}
