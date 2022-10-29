package com.cairn.gape.chromosome;

import com.cairn.gape.ga.BaseSupervisor;
import com.cairn.gape.ga.GaSupervisor;

/**
 * Interface for defining GA Chromosome object
 * 
 * @author Gareth Jones
 * @version 0.1
 */
public interface Chromosome {

	/**
	 * Creates a new chromosome and initializes it. To be used in Population
	 * creation.
	 * 
	 * @param problemHandle
	 *            GA control object typically a subclass of
	 *            {@link BaseSupervisor}
	 */
	public Chromosome create(GaSupervisor problemHandle);

	/**
	 * Creates a new chromosome and initializes it from an Object. To be used in
	 * Population creation. See individual chromosome classes to see what Object
	 * needs to be or if this method is supported.
	 * 
	 * @param problemHandle
	 * @param initialData
	 * @return
	 * @throws GaException
	 */
	public Chromosome create(GaSupervisor problemHandle, Object initialData);

	/**
	 * Free this chromosome. The interface is designed so that chromosome
	 * objects cat be recycled (as opposed to having to constantly create and
	 * garbage collect objects).
	 */
	public void freeChromosome();

	/**
	 * Return or calculate the chromosome's fitness.
	 */
	public double getFitness();

	/**
	 * Rebuilds the chromosome. This ensures that all data structures are
	 * conststent with this chromosome. For example there may be molecular
	 * co-ordinates in a shared molecule object that need to be updated.
	 */
	public double rebuild();

	/**
	 * Return true if two chromosmes are equal at the Genetic level.
	 */
	public boolean equals(Chromosome c);

	/**
	 * Return a distance measure for two chromosomes
	 */
	public double distance(Chromosome c);

	/**
	 * Return true if the two chromosomes share a niche.
	 */
	public boolean sameNiche(Chromosome c);

	/**
	 * Return false if the chromosome is invalid. For example the chromosome
	 * cannot be decoded.
	 */
	public boolean ok();

	/**
	 * The genetic information in Chromosome c is copied to this chromosome.
	 * 
	 * @param c
	 *            the Chromosome to copy from.
	 */
	public void copyGene(Chromosome c);

	/**
	 * Print a summary of the fitness.
	 */
	public String fitnessInfo();

	/**
	 * Print a summary of the gene.
	 */
	public String geneInfo();

	/**
	 * Print out a numeric identifier for the Chromosome. This could be a
	 * sequence number. The identifier should be unique to an object- but should
	 * be retained if the Chromosome object is recycled.
	 */
	public int getChromosomeId();

	/**
	 * Returns an empty Chromosome object for use in a genetic operator.
	 */
	public Chromosome getChromosome();
}
