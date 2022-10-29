package com.cairn.gape.ga;

import org.apache.log4j.Logger;

import com.cairn.gape.chromosome.Chromosome;
import com.cairn.gape.utils.InfoMessageLogger;

// Island Model

// Try to implement Population as much as possible

/**
 * 
 * Implementation of the Isalnd Model. Not all population mehods can be
 * implemented. The Island model consists for a set of GA populations (each
 * population being an island). Genetic operators can allow crossover of
 * mutation with populations or migration may take place between islands
 * 
 * @see com.cairn.gape.ga.LinkedPopLinearSel
 * 
 * @author Gareth Jones
 * 
 */
public class IslandModel implements Population {
	private static Logger logger;
	private static boolean logDebug;

	static {
		logger = Logger.getLogger(IslandModel.class);
		// logger.setLevel(Level.DEBUG);
		logDebug = logger.isDebugEnabled();
	}
	double selectPressure = 1.001, bestFitness;

	private final GaSupervisor problemHandle;

	private int nIslands = 5, size = 100;

	private final int migrateWt;

	private int nMigrate, totalOperatorWt, currentPopNo, nOperations, bestPopNo;

	// islands
	private final LinkedPopLinearSel populations[];

	private Chromosome best;

	private final InfoMessageLogger infoMessageLogger;

	/**
	 * Create an island model
	 * 
	 * @param ph
	 *            GA
	 * @param n
	 *            number of islands
	 * @param wt
	 *            weight of migration operator
	 */
	public IslandModel(GaSupervisor ph, int n, int wt) {
		nIslands = n;
		problemHandle = ph;
		infoMessageLogger = ph.getInfoMessageLogger();
		migrateWt = wt;

		populations = new LinkedPopLinearSel[n];
		for (int i = 0; i < populations.length; i++) {
			populations[i] = new LinkedPopLinearSel(ph);
			populations[i].popNo = i;
		}
	}

	// Size is the size of each population:

	/*
	 * This returns the size of each population (as opposed to the total number
	 * of individuals in the model, which is this multiplied by the number of
	 * islands) (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Population#getSize()
	 */
	@Override
	public int getSize() {
		return size;
	}

	/*
	 * Sets the size of each population.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Population#setSize(int)
	 */
	@Override
	public void setSize(int s) {
		int v = s * populations.length + 1;
		ChromosomeList.createListCache(v);
		size = s;
		for (int i = 0; i < populations.length; i++)
			populations[i].setSize(s);
	}

	/*
	 * Sets the nichesize for each island population.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Population#setNicheSize(int)
	 */
	@Override
	public void setNicheSize(int s) {
		for (int i = 0; i < populations.length; i++)
			populations[i].setNicheSize(s);
	}

	/*
	 * Sets the selection pressult in each island population.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Population#setSelectPressure(double)
	 */
	@Override
	public void setSelectPressure(double p) {
		selectPressure = p;
		for (int i = 0; i < populations.length; i++)
			populations[i].setSelectPressure(p);
	}

	/*
	 * Calls popInfo for each Island (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Population#popInfo()
	 */
	@Override
	public String popInfo() {
		String rtn = "";
		for (int i = 0; i < populations.length; i++) {
			int no = i + 1;
			rtn += "Population " + no;
			rtn += populations[i].popInfo();
		}
		return rtn;
	}

	// These population methods are not implemented in the Island model

	@Override
	public boolean addChromosome(Chromosome c) {
		throw new RuntimeException("method not implemented by Island Model");
	}

	@Override
	public Chromosome selectParent() {
		throw new RuntimeException("method not implemented by Island Model");
	}

	@Override
	public Chromosome replaceChromosome(Chromosome c) {
		throw new RuntimeException("method not implemented by Island Model");
	}

	@Override
	public Chromosome[] applyOperator(Operator o) {
		throw new RuntimeException("method not implemented by Island Model");
	}

	@Override
	public Operator selectOperator() {
		throw new RuntimeException("method not implemented by Island Model");
	}

	/*
	 * Adds the operator to each of the islands. (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.Population#addOperator(com.cairn.gape.Operator,
	 * int)
	 */
	@Override
	public void addOperator(Operator o, int weight) {
		totalOperatorWt += weight;
		for (int i = 0; i < populations.length; i++)
			populations[i].addOperator(o, weight);

	}

	/*
	 * Call create for each of the island populations. (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Population#create(java.lang.String,
	 * com.cairn.gape.GaSupervisor)
	 */
	@Override
	public void create(String chromosomeClass, GaSupervisor problemHandle) {
		for (int i = 0; i < populations.length; i++) {
			int no = i + 1;
			infoMessageLogger.infoMessageln("Creating Population " + no);
			populations[i].create(chromosomeClass, problemHandle);
		}

		getBest();
	}

	/*
	 * Calls create of each of the islands. The initialData should be of type
	 * Object[nIslands][]. initialDate[i] will be passed to island i.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.ga.Population#create(java.lang.String,
	 * com.cairn.gape.ga.GaSupervisor, java.lang.Object)
	 */
	@Override
	public void create(String chromosomeClass, GaSupervisor problemHandle,
			Object initialData) {

		Object initialPopData[] = null;
		try {
			initialPopData = (Object[]) initialData;
		} catch (ClassCastException ex) {
			throw new RuntimeException(
					"IsalandModel: initial data must be of type Object[]");
		}

		if (initialPopData.length != nIslands)
			throw new RuntimeException(
					"IsalandModel: initial data of type Object[] must be of length nIslands");

		for (int i = 0; i < populations.length; i++) {
			int no = i + 1;
			infoMessageLogger.infoMessageln("Creating Population " + no);
			populations[i].create(chromosomeClass, problemHandle, initialPopData[i]);
		}
		getBest();
	}

	/*
	 * Piopulations are selected in turn and the iterate method applied.
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Population#iterate()
	 */
	@Override
	public void iterate() {
		LinkedPopLinearSel pop = populations[currentPopNo];
		int wt = problemHandle.randomInt(0, totalOperatorWt + migrateWt - 1);
		if (wt < migrateWt && nIslands > 1) {
			migrate();
		} else {
			Operator op = pop.selectOperator();
			pop.applyOperator(op);
		}
		if (pop.getBest().getFitness() > bestFitness) {
			getBest();
			infoMessageLogger.infoMessage(bestInfo());
		}
		nOperations++;
		currentPopNo++;
		if (currentPopNo == nIslands)
			currentPopNo = 0;
	}

	/*
	 * Calls rebuild for each of the island populations (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Population#rebuild()
	 */
	@Override
	public void rebuild() {
		infoMessageLogger.infoMessageln("Rebuilding ..");
		for (int i = 0; i < populations.length; i++) {
			populations[i].rebuild();
		}
		getBest();
		infoMessageLogger.infoMessage(bestInfo());
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
			rtn += String.format("Op %5d nMigrate %5d :\n", nOperations, nMigrate);
			for (int i = 0; i < populations.length; i++) {
				int no = i + 1;
				rtn += String.format("Pop %2d ", no);
				rtn += populations[i].info();
			}
		} else {
			rtn += "Op " + nOperations + " nMigrate " + nMigrate + " :\n";
			for (int i = 0; i < populations.length; i++) {
				int no = i + 1;
				rtn += "Pop " + no + " " + populations[i].info();
			}
		}
		return rtn;
	}

	/**
	 * @return Information about the best individual across all islands
	 */
	public String bestInfo() {
		int no = bestPopNo + 1;
		String rtn = "";
		if (BaseSupervisor.USE_SPRINTF) {
			rtn += String.format("Pop %2d Op %4d ", no, nOperations);
		} else {
			rtn += "Pop " + no + " Op " + nOperations + " ";
		}
		rtn += best.fitnessInfo();
		return rtn;
	}

	/*
	 * Finds the best chromsome by searching all island populations.
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Population#getBest()
	 */
	@Override
	public Chromosome getBest() {
		best = populations[0].getBest();
		bestPopNo = 0;
		for (int i = 1; i < populations.length; i++) {
			Chromosome test = populations[i].getBest();
			if (test.getFitness() > best.getFitness()) {
				best = test;
				bestPopNo = i;
			}
		}
		bestFitness = best.getFitness();
		return best;
	}

	/*
	 * Calls freePop in each of the island populations. (non-Javadoc)
	 * 
	 * @see com.cairn.gape.Population#freePop()
	 */
	@Override
	public void freePop() {
		for (int i = 0; i < populations.length; i++)
			populations[i].freePop();
	}

	/**
	 * Migrates (copies) an individual from one island population to another
	 * (randomly chosen) island.
	 * 
	 */
	void migrate() {
		int otherPopNo = problemHandle.randomInt(0, nIslands - 2);
		if (otherPopNo >= currentPopNo)
			otherPopNo += 1;
		if (logDebug)
			logger.debug("Migrating from " + otherPopNo + " to " + currentPopNo);
		Chromosome parent = populations[otherPopNo].selectParent();
		Chromosome child = parent.getChromosome();
		if (logDebug)
			logger.debug("Migrating from " + parent.getChromosomeId() + " to "
					+ child.getChromosomeId());
		parent.copyGene(child);
		child.rebuild();
		populations[currentPopNo].replaceChromosome(child);
		nMigrate++;
	}

}
