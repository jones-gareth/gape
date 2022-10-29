package com.cairn.gape.ga;

import org.apache.commons.lang3.StringUtils;

import com.cairn.gape.chromosome.Chromosome;

/**
 * Wrapper class used to create linked lists of chromosome. Useful for creating
 * populations soted by fitness - makes parent selection for rank-based GAs
 * easy. JDK 1.4 means that this could be replace by newer collections.
 * 
 * @author Gareth Jones
 * 
 */
class ChromosomeList {
	Chromosome chromosome;

	// links for list
	ChromosomeList next, prev;

	static final boolean DEBUG = false;

	// free list
	private static ThreadLocal<ChromosomeList[]> freeChromosomeList = new ThreadLocal<ChromosomeList[]>();

	private static ThreadLocal<Integer> nAllocated = new ThreadLocal<Integer>();

	private static ThreadLocal<Boolean> allocated = new ThreadLocal<Boolean>() {
		@Override
		protected Boolean initialValue() {
			return false;
		}
	};

	private static ThreadLocal<String> threadName = new ThreadLocal<String>();

	static synchronized boolean isAllocated() {
		String name = Thread.currentThread().getName();
		// check that thread id has not been reused
		if (!StringUtils.equals(name, threadName.get())) {
			threadName.set(name);
			allocated.set(false);
			return false;
		}

		return allocated.get();
	}

	/**
	 * Create a Chromosome List Cache- should be of size total popsize + 1. We
	 * don't want to be constantly creating and garbage collecting objects.
	 * 
	 * @param size
	 */
	static synchronized void createListCache(int size) {
		if (isAllocated())
			return;
		ChromosomeList[] freeChromosomes = new ChromosomeList[size];
		for (int i = 0; i < freeChromosomes.length; i++)
			freeChromosomes[i] = new ChromosomeList();
		freeChromosomeList.set(freeChromosomes);
		nAllocated.set(size);
		allocated.set(true);
	}

	/**
	 * Constructor
	 */
	private ChromosomeList() {
	}

	/**
	 * @return the next unused chromosome list
	 */
	static ChromosomeList getChromosomeList() {
		int no = nAllocated.get();
		no--;
		nAllocated.set(no);
		return freeChromosomeList.get()[no];
	}

	/**
	 * Returns a chromosome list, with the chromosome added. The chromosome's
	 * fitness is determined.
	 * 
	 * @param c
	 * @return the next unused chomosome list.
	 */
	static ChromosomeList getChromosomeList(Chromosome c) {
		ChromosomeList cL = getChromosomeList();
		cL.chromosome = c;
		c.getFitness();
		return cL;
	}

	/**
	 * Put this chromsome list on the free list- called when we no longer use
	 * this in a population.
	 */
	void freeChromosomeList() {
		int no = nAllocated.get();
		freeChromosomeList.get()[no] = this;
		no++;
		nAllocated.set(no);
	}

	/**
	 * Adds the chromsome before this list. Creates a new ChromosomeList to
	 * store the chromosome and ensures all links are built correctly.
	 * 
	 * @param c
	 */
	void addBefore(Chromosome c) {
		ChromosomeList member = getChromosomeList(c);
		member.prev = prev;
		member.next = this;
		if (prev != null)
			prev.next = member;
		prev = member;
	}

	/**
	 * Adds the chromsome after this list. Creates a new ChromosomeList to store
	 * the chromosome and ensures all links are built correctly.
	 * 
	 * @param c
	 */
	void addAfter(Chromosome c) {
		ChromosomeList member = getChromosomeList(c);
		member.prev = this;
		member.next = next;
		if (next != null)
			next.prev = member;
		next = member;
	}

	/**
	 * Inserts a chromsome recursively, going down the links until the next link
	 * has a fitness greater than this cromosome's or we come to the end of the
	 * linked list. Call this on the head of the list.
	 * 
	 * @param c
	 */
	void insert(Chromosome c) {
		double f = c.getFitness();
		if (f > chromosome.getFitness())
			addBefore(c);
		else if (next == null)
			addAfter(c);
		else
			next.insert(c);
	}

	/**
	 * Removes this from the linked list. Make sure links are properly built.
	 * 
	 */
	void unlink() {
		if (DEBUG)
			System.out.println("UnLinking " + chromosome.getChromosomeId());
		if (prev != null)
			prev.next = next;
		if (next != null)
			next.prev = prev;
		next = prev = null;
		chromosome = null;
		freeChromosomeList();
	}

	/**
	 * Recursive function that goes doen the list counting the number of links.
	 * Call this on the list head.
	 * 
	 * @return size of the linked list
	 */
	int checkSize() {
		if (next == null)
			return 1;
		return 1 + next.checkSize();
	}
}
