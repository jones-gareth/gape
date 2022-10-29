package com.cairn.gape.chromosome;

import com.cairn.gape.ga.GaSupervisor;
import com.cairn.gape.ga.Operator;

/**
 * Implements mutation on a Binary string.
 * 
 * @author Gareth Jones
 *
 */
public class BinaryStringMutation implements Operator {
	private GaSupervisor gaSupervisor;

	private Class<?> chromosomeClass;

	private static final boolean DEBUG = false;

	private Chromosome children[] = new Chromosome[1];

	public BinaryStringMutation(GaSupervisor gaSupervisor, String chromosomeClassString) {
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
		BinaryStringChromosome parent, child;

		try {
			child = (BinaryStringChromosome) chromosomeClass.newInstance();
		} catch (IllegalAccessException ex1) {
			throw new RuntimeException("Illegal access for creating chromosome "
					+ chromosomeClass);
		} catch (InstantiationException ex2) {
			throw new RuntimeException("Exception on instantiating " + chromosomeClass);
		}

		child.createEmpty();
		child.setProblem(gaSupervisor);
		parent = (BinaryStringChromosome) parents[0];
		parent.copyGene(child);
		child.mutate();

		if (DEBUG) {
			System.out.println("Mutating");
			System.out.print(parent.geneInfo());
			System.out.println("-->");
			System.out.print(child.geneInfo());
		}

		children[0] = child;
		return children;
	}

	@Override
	public int nChildren() {
		return 1;
	}

	@Override
	public int nParents() {
		return 1;
	}

}
