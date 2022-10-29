package com.cairn.gape.chromosome;

import org.apache.log4j.Logger;

import com.cairn.gape.Superposition;
import com.cairn.gape.ga.Operator;

/**
 * SuperpositionCrossover. Crossover operator class for GAPE- uses full mixing
 * 
 * @author Gareth Jones
 * @version 0.1
 */
public class SuperpositionCrossover implements Operator {
	private final Chromosome children[] = new Chromosome[2];
	private static final Logger logger = Logger.getLogger(SuperpositionChromosome.class);

	public SuperpositionCrossover(Superposition p) {
	}

	@Override
	public int nParents() {
		return 2;
	}

	@Override
	public int nChildren() {
		return 2;
	}

	@Override
	public Chromosome[] apply(Chromosome parents[]) {
		BinaryAndIntegerChromosome child1, child2, parent1, parent2;

		child1 = BinaryAndIntegerChromosome.getFreedChromosome();
		child2 = BinaryAndIntegerChromosome.getFreedChromosome();

		parent1 = (BinaryAndIntegerChromosome) parents[0];
		parent2 = (BinaryAndIntegerChromosome) parents[1];

		if (logger.isDebugEnabled()) {
			logger.debug("Crossover from " + parent1.getChromosomeId() + ","
					+ parent2.getChromosomeId() + " to " + child1.getChromosomeId() + ","
					+ child2.getChromosomeId());
			logger.debug(parent1.geneInfo());
			logger.debug(parent2.geneInfo() + "-->");
		}
		parent1.crossover(parent2, child1, child2);
		if (logger.isDebugEnabled()) {
			logger.debug(child1.geneInfo());
			logger.debug(child2.geneInfo());
		}

		children[0] = child1;
		children[1] = child2;
		return children;
	}
}
