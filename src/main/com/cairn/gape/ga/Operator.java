package com.cairn.gape.ga;

import com.cairn.gape.chromosome.Chromosome;

/**
 * Interface for genetic operators (crossover, mutation..)
 * 
 * @author Gareth Jones
 * 
 */
public interface Operator {
	/**
	 * @return the number of parents the operator requires
	 * @throws GaException
	 */
	int nParents();

	/**
	 * @return the number of child chromsomes the operator creates
	 * @throws GaException
	 */
	int nChildren();

	/**
	 * Apply a genetic operator. Parents are passed in an argument and an array
	 * of child chromosomes are returned.
	 * 
	 * @param c
	 *            array of parents
	 * @return array of children
	 * @throws GaException
	 */
	Chromosome[] apply(Chromosome c[]);

}
