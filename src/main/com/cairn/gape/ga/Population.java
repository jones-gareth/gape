package com.cairn.gape.ga;

import com.cairn.gape.chromosome.Chromosome;

/**
 * Interface to define GA population.
 * 
 * @author Gareth Jones
 * 
 */
public interface Population {

	/**
	 * @return Population size
	 */
	int getSize();

	/**
	 * Set population size
	 * 
	 * @param size
	 */
	void setSize(int size);

	/**
	 * Sets the niche size. Only this number of chromosomes will be able to
	 * share the niche (assuming the population supports niching).
	 * 
	 * @param nicheSize
	 */
	/**
	 * @param nicheSize
	 */
	void setNicheSize(int nicheSize);

	/**
	 * Used to set selection pressure. Selection pressure is the ratio of
	 * probabilities that the best chromosome will be selected as a parent,
	 * relative to the worst chromosome. Not all populations will support this.
	 * 
	 * @param selectPressure
	 */
	void setSelectPressure(double selectPressure);

	/**
	 * @return an informational string about the population. This should return
	 *         detailed information including a description of each chromosme in
	 *         the population.
	 */
	String popInfo();

	/**
	 * Adds a chromosome to the population.
	 * 
	 * @param c
	 * @return
	 * */
	boolean addChromosome(Chromosome c);

	/**
	 * @return a parent chromosome for use in GA operator.
	 * 
	 */
	Chromosome selectParent();

	/**
	 * Adds a chromosome to the population. The population should also remove a
	 * chromosome so that the size of the population is unchanged.
	 * 
	 * @param c
	 * @return
	 */
	Chromosome replaceChromosome(Chromosome c);

	/**
	 * Applies a genetic operator to the population. Any children are returned.
	 * 
	 * @param o
	 * @return
	 * 
	 * @see Operator
	 */
	Chromosome[] applyOperator(Operator o);

	/**
	 * Selects an operator. Typically roulette wheel selecion is used to pick an
	 * operator.
	 * 
	 * @return
	 */
	Operator selectOperator();

	/**
	 * Registers an operator with the population.
	 * 
	 * @param o
	 *            the operator
	 * @param weight
	 *            a weight to use in operator selection.
	 */
	void addOperator(Operator o, int weight);

	/**
	 * Creates a new population. Arguments define the chromosome class to
	 * generate and a handle to the problem.
	 * 
	 * @param chromosomeClass
	 * @param problemHandle
	 */
	void create(String chromosomeClass, GaSupervisor problemHandle);

	/**
	 * Creates a new population. Arguments define the chromosome class to
	 * generate, a handle to the problem and initial data to seed the
	 * population.
	 * 
	 * @param chromosomeClass
	 * @param problemHandle
	 * @param initialData
	 */
	void create(String chromosomeClass, GaSupervisor problemHandle, Object initialData);

	/**
	 * Performs an iteration- in a steady-state GA typically involves picking an
	 * operator, applying it and updating population with any children.
	 * 
	 */
	void iterate();

	/**
	 * Frees up memory associated with the population.
	 * 
	 */
	void freePop();

	/**
	 * @return A simple (one line) string describing the current state of the
	 *         population- this is typically information about the best
	 *         chromosome in the population.
	 */
	String info();

	/**
	 * @return The most fit chromosome in the population.
	 */
	Chromosome getBest();

	/**
	 * Rebuild all chromosomes in the population. This may be required it we do
	 * things like annealing parameters through a GA run. For example GAPE
	 * adjusts feature overlap parameters- when the parameters change the
	 * population needs to be rebuilt.
	 * 
	 */
	void rebuild();
}
