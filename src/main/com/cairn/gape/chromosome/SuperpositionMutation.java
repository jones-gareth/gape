package com.cairn.gape.chromosome;

import org.apache.log4j.Logger;

import com.cairn.gape.Superposition;
import com.cairn.gape.ga.Operator;

/**
 * SuperpositionMuation. Mutation operator class for GAPE.
 * 
 * @author Gareth Jones
 * @version 0.1
 */
public class SuperpositionMutation implements Operator {
	private final Chromosome children[] = new Chromosome[1];
	private static final Logger logger = Logger.getLogger(SuperpositionMutation.class);

	public SuperpositionMutation(Superposition p) {
	}

	@Override
	public int nParents() {
		return 1;
	}

	@Override
	public int nChildren() {
		return 1;
	}

	@Override
	public Chromosome[] apply(Chromosome parents[]) {
		BinaryAndIntegerChromosome parent, child;

		child = BinaryAndIntegerChromosome.getFreedChromosome();

		parent = (BinaryAndIntegerChromosome) parents[0];
		parent.copyGene(child);
		child.mutate();

		if (logger.isDebugEnabled()) {
			logger.debug("Mutating from " + parent.getChromosomeId() + " to "
					+ child.getChromosomeId());
			logger.debug("Mutating");
			logger.debug(parent.geneInfo() + "-->");
			logger.debug(child.geneInfo());
		}

		children[0] = child;
		return children;
	}
}
