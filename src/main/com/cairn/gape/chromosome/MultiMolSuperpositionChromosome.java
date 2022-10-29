package com.cairn.gape.chromosome;

import com.cairn.gape.MultiMolSuperposition;
import com.cairn.gape.Superposition;
import com.cairn.gape.molecule.MultiGaMolecule;

/**
 * The chromosome for use in the multi-conformer version of gape.
 * 
 * @author Gareth Jones
 * 
 */
public class MultiMolSuperpositionChromosome extends SuperpositionChromosome {
	private IntegerStringChromosome conformationChromosome;
	private double binaryConformerOperatorProbability;

	public MultiMolSuperpositionChromosome() {
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @seecom.cairn.gape.chromosome.SuperpositionChromosome#
	 * copyReferenceCoordinates(int)
	 * 
	 * Create coordinates from conformer coordinates.
	 */
	@Override
	protected void copyReferenceCoordinates(int moleculeNo) {
		MultiMolSuperposition multiMolSuperposition = (MultiMolSuperposition) problem;
		MultiGaMolecule molecule = (MultiGaMolecule) molecules.get(moleculeNo);
		int nConformers = molecule.nConformers();
		int conformerNo = 0;
		if (nConformers > 1) {
			int position = multiMolSuperposition.getConformerEntryPoint(moleculeNo);
			conformerNo = conformationChromosome.getValues(position);
			assert conformerNo >= 0 && conformerNo < nConformers : "Invalid conformer number";
		}

		molecule.setConformer(conformerNo);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.chromosome.SuperpositionChromosome#createChromosome()
	 */
	@Override
	public SuperpositionChromosome createChromosome() {
		return new MultiMolSuperpositionChromosome();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.chromosome.SuperpositionChromosome#createEmpty(com
	 * .cairn.gape.Superposition)
	 */
	@Override
	public void createEmpty(Superposition s) {
		super.createEmpty(s);

		MultiMolSuperposition superposition = (MultiMolSuperposition) s;
		if (superposition.getConformerStringLength() > 0) {
			conformationChromosome = new IntegerStringChromosome(superposition,
					superposition.getConformerStringLength(),
					superposition.getConformerStringRanges());
			conformationChromosome.setAllowNullDefault(false);
		}

		if (binaryString != null) {
			binaryConformerOperatorProbability = ((double) binaryString.nBytes)
					/ ((double) (binaryString.nBytes + conformationChromosome.length));
		}

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.chromosome.BinaryAndIntegerChromosome#mutate()
	 */
	@Override
	public void mutate() {
		if (hasConformationalData() && problem.normalRand() < 0.5) {
			if (problem.normalRand() < binaryConformerOperatorProbability)
				binaryString.mutate();
			else
				conformationChromosome.mutate();
		} else {
			integerString.mutate();
		}
	}

	/**
	 * @return true if we have conformational information
	 */
	private boolean hasConformationalData() {
		if ((binaryString != null && binaryString.nBytes > 0)
				|| (conformationChromosome != null && conformationChromosome.length > 0))
			return true;
		return false;
	}

	/**
	 * Crossover operator. Either does crossover on the binary string or full
	 * mixing and crossover on the integer string with equal probability
	 * 
	 * @param parent2
	 * @param child1
	 * @param child2
	 */
	@Override
	public void crossover(BinaryAndIntegerChromosome p2, BinaryAndIntegerChromosome c1,
			BinaryAndIntegerChromosome c2) {

		MultiMolSuperpositionChromosome parent2 = (MultiMolSuperpositionChromosome) p2;
		MultiMolSuperpositionChromosome child1 = (MultiMolSuperpositionChromosome) c1;
		MultiMolSuperpositionChromosome child2 = (MultiMolSuperpositionChromosome) c2;

		this.copyGene(child1);
		parent2.copyGene(child2);

		if (hasConformationalData() && problem.normalRand() < 0.5) {
			if (problem.normalRand() < binaryConformerOperatorProbability)
				binaryString.crossover(parent2.binaryString, child1.binaryString,
						child2.binaryString);
			else
				conformationChromosome.crossover(parent2.conformationChromosome,
						child1.conformationChromosome, child2.conformationChromosome);
		} else {
			integerString.fullMixing(parent2.integerString, child1.integerString,
					child2.integerString);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.chromosome.BinaryAndIntegerChromosome#copyGene(com
	 * .cairn.gape.chromosome.Chromosome)
	 */
	@Override
	public void copyGene(Chromosome c2) {
		super.copyGene(c2);

		MultiMolSuperpositionChromosome c = (MultiMolSuperpositionChromosome) c2;
		conformationChromosome.copyGene(c.conformationChromosome);

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.chromosome.BinaryAndIntegerChromosome#initialize()
	 */
	@Override
	protected void initialize() {
		super.initialize();
		if (conformationChromosome != null)
			conformationChromosome.initialize();
	}

}
