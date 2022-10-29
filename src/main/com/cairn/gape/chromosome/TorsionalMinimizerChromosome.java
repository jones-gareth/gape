package com.cairn.gape.chromosome;

import com.cairn.gape.TorsionalMinimizer;
import com.cairn.gape.molecule.GaMolecule;

public class TorsionalMinimizerChromosome extends BinaryStringChromosome implements
		Chromosome {
	private TorsionalMinimizer tm;
	private boolean hasFitness = false;
	private double fitness;

	public TorsionalMinimizerChromosome() {
		;
	}

	public TorsionalMinimizerChromosome(TorsionalMinimizer tm) {
		this.tm = tm;
		this.problem = tm;
	}

	@Override
	public Chromosome getChromosome(Chromosome c) {
		throw new RuntimeException("method not defined");
	}

	@Override
	public void freeChromosome() {

	}

	@Override
	public void createEmpty() {
		this.tm = (TorsionalMinimizer) problem;
		nBits = this.tm.getMolecule().getnRotatableBonds() * 8;
		super.createEmpty();
	}

	@Override
	public double getFitness() {
		if (hasFitness)
			return fitness;
		GaMolecule m = tm.getMolecule();
		m.copyReferenceCorordinates();
		m.generateConformation(this);
		// fitness = -m.taff.molEnergy();
		fitness = -m.conformationalEnergy();
		hasFitness = true;
		return fitness;
	}

	public void printFitness() {
		System.out.println("Fit " + fitness);
	}

	@Override
	public double rebuild() {
		hasFitness = false;
		return getFitness();
	}

	@Override
	public double distance(Chromosome c) {
		return -1;
	}

	@Override
	public boolean sameNiche(Chromosome c) {
		return false;
	}

	@Override
	public boolean ok() {
		return true;
	}

}
