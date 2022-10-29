package com.cairn.gape.ga;

import org.apache.log4j.Logger;

import com.cairn.gape.chromosome.Chromosome;
import com.cairn.gape.utils.InfoMessageLogger;

/**
 * Implementation of steady state no-duplicates GA with rank-based linear
 * selection pressure. The chromosomes are ordered in a linked list.
 * 
 * @author Gareth Jones
 * 
 */
public class LinkedPopLinearSel implements Population {

	private static Logger logger;
	private static boolean logDebug;

	static {
		logger = Logger.getLogger(LinkedPopLinearSel.class);
		// logger.setLevel(Level.DEBUG);
		logDebug = logger.isDebugEnabled();
	}
	private int size = 100;

	private double selectPressure = 1.001;

	private final double selectStart = 5000;

	private double selectDecrement, totalFitness, bestFitness;

	private GaSupervisor problemHandle;

	private ChromosomeList startPop, endPop;

	private final NicheMatch nicheMatch;

	private final Operator operators[];

	int nOperators, operatorWeights[], totalOperatorWeight, nicheSize = 5, popNo,
			nOperations, nAdded, nNiche, nFail, nDuplicates = 0;

	private java.text.NumberFormat nf;

	protected InfoMessageLogger infoMessageLogger = new InfoMessageLogger();

	/**
	 * Constructor- allocates some storage
	 */
	private LinkedPopLinearSel() {
		nicheMatch = new NicheMatch();
		operatorWeights = new int[100];
		operators = new Operator[100];
	}

	/**
	 * Constructor - sets GA problem
	 * 
	 * @param ph
	 */
	public LinkedPopLinearSel(GaSupervisor ph) {
		this();
		problemHandle = ph;
		nf = ph.getNumberFormat();
		infoMessageLogger = ph.getInfoMessageLogger();
	}

	/**
	 * @return population size
	 */
	@Override
	public int getSize() {
		return size;
	}

	/**
	 * Set population size
	 * 
	 * @param v
	 */
	@Override
	public void setSize(int v) {
		ChromosomeList.createListCache(v + 1);
		this.size = v;
	}

	/**
	 * @return niche size
	 */
	public int getNicheSize() {
		return nicheSize;
	}

	/**
	 * Sets niche size. One this number of chromosomes will be allowed in a
	 * niche.
	 */
	@Override
	public void setNicheSize(int v) {
		this.nicheSize = v;
	}

	/**
	 * @return selection pressure
	 */
	public double getSelectPressure() {
		return selectPressure;
	}

	/**
	 * Sets selection pressure- this is the ratio of the probability that the
	 * best chromosome is selected as a parent to the probability that the least
	 * fit chromosome is chosen
	 * 
	 * @param v
	 */
	@Override
	public void setSelectPressure(double v) {
		this.selectPressure = v;
	}

	/**
	 * Sets up variables for performing roulette wheel parent selection on
	 * linear ranked fitness values.
	 * 
	 */
	private void setupSelectPressure() {
		double p = size - 1;

		if (selectPressure > 2.0)
			throw new RuntimeException(
					"set_sel_press: The selection pressure must be less than 2.0\n");

		double x = selectStart * (1 - selectPressure);
		double d = (selectPressure / 2.0 - 1.0) * p;
		selectDecrement = x / d;

		double currentFitness;
		totalFitness = currentFitness = selectStart;
		for (int i = 1; i < size; i++) {
			currentFitness = selectStart + i * selectDecrement;
			totalFitness += currentFitness;
		}

		if (logger.isDebugEnabled())
			logger.debug("Rank fitness base is " + selectStart + " increment is "
					+ selectDecrement + " PopSize is " + size + " select pressure is "
					+ selectPressure + " total fitness is " + totalFitness);

	}

	/**
	 * Searches the population to find any chromosome that is an exact match to
	 * c.
	 * 
	 * @param c
	 * @return
	 */
	private ChromosomeList findExactMatch(Chromosome c) {
		for (ChromosomeList cL = endPop; cL != null; cL = cL.prev) {
			if (c.equals(cL.chromosome)) {
				nDuplicates++;
				return cL;
			}
		}
		return null;
	}

	/**
	 * Class for handing niching
	 */
	private class NicheMatch {
		private boolean hasNicheMatch, betterThanMatch;

		private int nMatched;

		private ChromosomeList match;

		/**
		 * Searches the population for a niche match to the chromsome. The
		 * actual match is the least fit individual in the niche. We aslo record
		 * if the chromsome is more fit than the niche. If so we can replace the
		 * worst member of the niche by this chromosome.
		 * 
		 * @param c
		 */
		public void findNicheMatch(Chromosome c) {
			hasNicheMatch = false;
			betterThanMatch = false;
			match = null;
			nMatched = 0;
			// Work from the bottom to make sure niche checks against the
			// least fit.
			for (ChromosomeList cL = endPop; cL != null; cL = cL.prev) {
				Chromosome test = cL.chromosome;
				if (c.sameNiche(test)) {
					if (match == null) {
						match = cL;
						if (test.getFitness() < c.getFitness())
							betterThanMatch = true;
					}
					nMatched++;
					if (nMatched == nicheSize) {
						hasNicheMatch = true;
						nNiche++;
						return;
					}
				}
			}
			match = null;
		}

	}

	/**
	 * Removes a chromsome from the linked list. Correctly maintains the start
	 * and end of the list
	 * 
	 * @param cL
	 * @see ChromosomeList#unlink()
	 */
	private void removeLink(ChromosomeList cL) {
		if (logger.isDebugEnabled())
			logger.debug("Removing Link " + cL.chromosome.getChromosomeId() + " pop "
					+ popNo);
		if (cL.chromosome != null)
			cL.chromosome.freeChromosome();
		if (startPop == cL) {
			startPop = cL.next;
			startPop.prev = null;
		}
		if (endPop == cL) {
			endPop = cL.prev;
			endPop.next = null;
		}
		cL.unlink();
	}

	/**
	 * Inserts a chromsome into the linked list. Correctly maintains the start
	 * and end of the list
	 * 
	 * @param c
	 * @see ChromosomeList#insert(Chromosome);
	 */
	private void insert(Chromosome c) {
		if (logger.isDebugEnabled()) {
			logger.debug("Inserting " + c.getChromosomeId() + " pop " + popNo);
			logger.debug("Pop Size " + checkSize());
		}
		startPop.insert(c);
		if (startPop.prev != null)
			startPop = startPop.prev;
		if (endPop.next != null)
			endPop = endPop.next;
		if (logger.isDebugEnabled())
			logger.debug("Pop Size " + checkSize());
	}

	/**
	 * @return the size of the population by counting the linked list. Used for
	 *         debugging
	 * @see ChromosomeList#checkSize()
	 */
	private int checkSize() {
		return startPop.checkSize();
	}

	/*
	 * Adds a chromosome to a building population. Check for duplicates and
	 * niche matches. Returns true only if the chromsome is added to the
	 * population. (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.Population#addChromosome(com.cairn.gape.Chromosome
	 * )
	 */
	@Override
	public boolean addChromosome(Chromosome c) {

		// check to see chromsome is decoded ok
		if (!c.ok()) {
			nFail++;
			if (logDebug)
				logger.debug("addChromosome: Failed: freeing");
			c.freeChromosome();
			return false;
		}

		// check for duplicates
		if (findExactMatch(c) != null) {
			if (logDebug)
				logger.debug("addChromosome: Exact match: freeing");
			c.freeChromosome();
			return false;
		}

		// check for niche matches
		nicheMatch.findNicheMatch(c);
		if (nicheMatch.hasNicheMatch) {
			if (nicheMatch.betterThanMatch) {
				// niche match, but better than worst member of niche
				if (logDebug)
					logger.debug("addChromosome: Niche match: replacing");
				// in this case remove the niche match
				removeLink(nicheMatch.match);
			} else {
				if (logDebug)
					logger.debug("addChromosome: Niche match: freeing");
				c.freeChromosome();
				return false;
			}
		}

		// add chromsome
		insert(c);
		return true;
	}

	/*
	 * Adds a chromosome to a population in staeady state. The chromosome should
	 * replace another member. Check for duplicates and niche matches. Return
	 * replaced chromsome or null if we don't add the new
	 * chromsome.(non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.Population#replaceChromosome(com.cairn.gape.
	 * Chromosome)
	 */
	@Override
	public Chromosome replaceChromosome(Chromosome c) {
		ChromosomeList match = null;

		if (logDebug)
			logger.debug("Replacing " + c.getChromosomeId() + " pop " + popNo);

		// unable to decode chromosome
		if (!c.ok()) {
			nFail++;
			c.freeChromosome();
			return null;
		}
		if (findExactMatch(c) != null) {
			// duplicate- discard
			c.freeChromosome();
			return null;
		}

		nicheMatch.findNicheMatch(c);
		if (nicheMatch.hasNicheMatch) {
			// check for niche match
			if (nicheMatch.betterThanMatch) {
				// replace worst member in niche if new chromsome is fitter
				match = nicheMatch.match;
			} else {
				// else discard new chromsome
				c.freeChromosome();
				return null;
			}
		} else {
			// if we doon't hit a niche replace the worst chromsome in the
			// currect population
			match = endPop;
		}

		removeLink(match);
		nAdded++;

		insert(c);
		return match.chromosome;
	}

	/*
	 * Implement roulette-wheel parent selection on the population to select a
	 * parent. Rank-based fitnesses are use with a user-defined selection
	 * pressure. (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Population#selectParent()
	 * 
	 * @see #setSelectPressure(double)
	 */
	@Override
	public Chromosome selectParent() {
		double n = problemHandle.normalRand() * totalFitness;
		double total = 0, currentFitness = 0;
		int i = 0;
		for (ChromosomeList cL = startPop; cL != null; cL = cL.next) {
			currentFitness = selectStart + i++ * selectDecrement;
			total += currentFitness;
			if (total >= n)
				return cL.chromosome;
		}
		throw new RuntimeException("selectParent: parent selection error");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.Population#addOperator(com.cairn.gape.Operator,
	 * int)
	 */
	@Override
	public void addOperator(Operator op, int weight) {
		operators[nOperators] = op;
		operatorWeights[nOperators] = weight;
		totalOperatorWeight += weight;
		nOperators++;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Population#selectOperator()
	 */
	@Override
	public Operator selectOperator() {
		int wt = problemHandle.randomInt(0, totalOperatorWeight - 1);
		Operator op = null;
		int currentWt = 0;
		for (int i = 0; i < nOperators; i++) {
			currentWt += operatorWeights[i];
			if (wt < currentWt) {
				op = operators[i];
				break;
			}
		}
		return op;
	}

	static private final ThreadLocal<Chromosome[]> parents = new ThreadLocal<Chromosome[]>() {
		@Override
		protected Chromosome[] initialValue() {
			return new Chromosome[10];
		}
	};

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.Population#applyOperator(com.cairn.gape.Operator
	 * )
	 */
	@Override
	public Chromosome[] applyOperator(Operator op) {
		if (logDebug)
			logger.debug("Applying Operator");
		int nParents = op.nParents();

		// Pick parents
		for (int i = 0; i < nParents;) {
			Chromosome c = selectParent();
			boolean ok = true;
			for (int j = 0; j < i; j++)
				if (c == parents.get()[j])
					ok = false;
			if (ok) {
				parents.get()[i] = c;
				i++;
			}
		}
		// apply operator and get children
		Chromosome children[] = op.apply(parents.get());
		for (int i = 0; i < op.nChildren(); i++) {
			if (children[i] != null) {
				if (logDebug)
					logger.debug("Replacing");
				// add a child into the population
				replaceChromosome(children[i]);
			}
		}

		nOperations++;
		return children;
	}

	@Override
	public void create(String chromosomeClass, GaSupervisor h) {
		create(chromosomeClass, h, null);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Population#create(java.lang.String,
	 * com.cairn.gape.GaSupervisor)
	 */
	@Override
	public void create(String chromosomeClass, GaSupervisor h, Object initialData) {

		problemHandle = h;
		Class<?> cls = null;
		try {
			// get chromsome class
			cls = Class.forName(chromosomeClass);
		} catch (ClassNotFoundException ex) {
			throw new RuntimeException("No chromosome class " + chromosomeClass);
		}

		// extract initial data array if present.
		Object initialChromosomeData[] = null;
		if (initialData != null) {
			try {
				initialChromosomeData = (Object[]) initialData;
			} catch (ClassCastException ex) {
				throw new RuntimeException(
						"LinkedPopLinearSel: initial data must be of type Object[]");
			}
		}

		int n = 0;
		int nInitial = 0;

		while (n < size) {
			Chromosome c;
			try {
				// create a new chromosome instance
				c = (Chromosome) cls.newInstance();
			} catch (IllegalAccessException ex1) {
				throw new RuntimeException("Illegal access for creating chromosome "
						+ chromosomeClass);
			} catch (InstantiationException ex2) {
				throw new RuntimeException("Exception on instantiating "
						+ chromosomeClass);
			}

			// Throw exception if for some reason we can't generate the
			// population
			if (nFail > 1000 * size)
				throw new RuntimeException("Too many failures in creating population");
			if (nNiche > 1000 * size)
				throw new RuntimeException(
						"Too many niche matches in creating population");
			if (nDuplicates > 1000 * size)
				throw new RuntimeException("Too many duplicates in creating population");

			// initialize the chromsome
			if (initialChromosomeData != null && initialChromosomeData.length > nInitial) {
				c = c.create(problemHandle, initialChromosomeData[nInitial]);
				nInitial++;
			} else {
				c = c.create(problemHandle);
			}

			// check to see if we can decode chromsome
			if (!c.ok()) {
				nFail++;
				c.freeChromosome();
				continue;
			}

			// add to new population
			if (startPop == null) {
				startPop = endPop = ChromosomeList.getChromosomeList(c);
				n = 1;
			} else {
				if (addChromosome(c) && !nicheMatch.hasNicheMatch) {
					n++;
					if (logger.isDebugEnabled())
						logger.debug("Pop Size " + checkSize() + " n " + n);
				}
			}
		}

		// finish creating population
		bestFitness = startPop.chromosome.getFitness();
		infoMessageLogger.infoMessage(bestInfo());
		logger.debug("Pop Size " + checkSize());

		// set up parent selection parameters
		setupSelectPressure();
	}

	/**
	 * @return an Informational string about the best individual in the
	 *         population
	 */
	public String bestInfo() {
		Chromosome best = startPop.chromosome;
		infoMessageLogger.infoMessage("Op " + nOperations + " ");
		return best.fitnessInfo();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Population#getBest()
	 */
	@Override
	public Chromosome getBest() {
		return startPop.chromosome;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Population#info()
	 */
	@Override
	public String info() {
		String rtn = "";
		if (BaseSupervisor.USE_SPRINTF) {
			rtn += String.format(
					"Op %5d Fit %7.1f Added %5d Dups %5d Niche %5d Fail %4d\n",
					nOperations, getBest().getFitness(), nAdded, nDuplicates, nNiche,
					nFail);
		} else {
			rtn += "Op " + nOperations + " Fit " + nf.format(getBest().getFitness())
					+ " Added " + nAdded + " Dups " + nDuplicates + " nNiche " + nNiche
					+ " Fail " + nFail + "\n";
		}
		return rtn;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Population#iterate()
	 */
	@Override
	public void iterate() {
		Operator op = selectOperator();
		applyOperator(op);

		double fit = startPop.chromosome.getFitness();
		// check to see if we have a new best solution
		if (fit > bestFitness) {
			bestFitness = fit;
			infoMessageLogger.infoMessage(bestInfo());
		}

		if (nOperations > 0 && nOperations % 50 == 0)
			infoMessageLogger.infoMessage(info());
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Population#popInfo()
	 */
	@Override
	public String popInfo() {
		String rtn = "";
		for (ChromosomeList cL = startPop; cL != null; cL = cL.next) {
			Chromosome c = cL.chromosome;
			double fitness = c.getFitness();
			rtn += fitness + "\t";
			rtn += c.geneInfo();
		}
		return rtn;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Population#rebuild()
	 */
	@Override
	public void rebuild() {
		// cL hold the start of the current population
		ChromosomeList cL = startPop;
		if (logger.isDebugEnabled())
			logger.debug("Rebuilding pop " + popNo + " size " + checkSize());

		// remove link references
		startPop = endPop = null;

		// rebuild population redermining fitnesses of all chromomsomes.
		while (cL != null) {
			Chromosome c = cL.chromosome;
			// get new fitness score for chromosome
			c.rebuild();
			if (startPop == null) {
				// new link references
				startPop = endPop = ChromosomeList.getChromosomeList(c);
			} else {
				insert(c);
			}

			// move down the old list
			ChromosomeList next = cL.next;
			cL.freeChromosomeList();
			cL.next = cL.prev = null;
			cL = next;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Population#freePop()
	 */
	@Override
	public void freePop() {
		ChromosomeList cL = startPop;
		startPop = endPop = null;

		while (cL != null) {
			Chromosome c = cL.chromosome;
			c.freeChromosome();
			ChromosomeList next = cL.next;
			cL.freeChromosomeList();
			cL.next = cL.prev = null;
			cL = next;
		}
	}

}
