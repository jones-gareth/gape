package com.cairn.gape;

import org.apache.log4j.Logger;

import com.cairn.gape.chromosome.Chromosome;
import com.cairn.gape.chromosome.TorsionalMinimizerChromosome;
import com.cairn.gape.ga.BaseSupervisor;
import com.cairn.gape.ga.LinkedPopLinearSel;
import com.cairn.gape.ga.Operator;
import com.cairn.gape.ga.Population;
import com.cairn.gape.molecule.GaMolecule;
import com.cairn.gape.molecule.TorsionalDistributions;

/**
 * A simple GA to perform torsional minimization of a molecule
 * 
 * @author Gareth Jones
 * 
 */
public class TorsionalMinimizer extends BaseSupervisor {
	private volatile GaMolecule molecule;

	private volatile Population pop;

	private volatile TorsionalMinimizerChromosome best;

	private volatile TorsionalDistributions tordist;

	/**
	 * Runs the GA as an application.
	 * 
	 * @param args
	 */
	public static void main(String args[]) {
		if (args.length != 2) {
			System.out.println("Usage: TorsionalMinimizer <conf file> <molecule file>");
			System.exit(0);
		}
		GaMolecule molecule = new GaMolecule();
		molecule.setFindFeatures(false);
		molecule.loadFile(args[1]);
		TorsionalMinimizer tm = new TorsionalMinimizer();
		tm.minimize(args[0], molecule);
	}

	/**
	 * Initializes the GA
	 * 
	 * @param configFile
	 * @param m
	 */
	private void setup(String configFile, GaMolecule m) {
		init(configFile);
		if (hasKey("tordist_file")) {
			tordist = new TorsionalDistributions();
			tordist.readTordistFile(getStringValue("tordist_file"));
		}
		molecule = m;
		molecule.setProblem(this);
		molecule.setup();
		if (tordist != null)
			tordist.moleculeMatchTordist(molecule);
	}

	/**
	 * Uses a GA to search the torsional space of a Molecule to minimize it.
	 * 
	 * @param configFile
	 * @param m
	 */
	public void minimize(String configFile, GaMolecule m) {
		minimize(configFile, m, null);

	}

	/**
	 * Uses a GA to search the torsional space of a Molecule to minimize it. You
	 * can supply chromosome initialization data for this routine.
	 * 
	 * @param configFile
	 * @param m
	 * @param initData
	 */
	public void minimize(String configFile, GaMolecule m, Object initData) {
		setup(configFile, m);

		int nIterations = getIntValue("n_iterations");
		int popSize = getIntValue("popsize");
		double selectPressure = getDoubleValue("select_pressure");

		pop = new LinkedPopLinearSel(this);
		pop.setSelectPressure(selectPressure);
		pop.setSize(popSize);
		pop.create("com.cairn.gape.chromosome.TorsionalMinimizerChromosome", this,
				initData);
		infoMessage(pop.popInfo());

		pop.addOperator(new TorsionalMinimizerCrossover(this), 10);
		pop.addOperator(new TorsionalMinimizerMutation(this), 10);

		for (int i = 0; i < nIterations; i++)
			pop.iterate();

		infoMessage(pop.popInfo());
		best = (TorsionalMinimizerChromosome) pop.getBest();
		best.rebuild();
		String file = "GA Minimized " + molecule.getBaseName() + ".mol2";
		infoMessageln("Writing solution to " + file);
		molecule.writeSybylMol2File(file, "Minimized GA Molecule");
	}

	public GaMolecule getMolecule() {
		return molecule;
	}

}

/**
 * Crossover operator for Torsional Minimizer.
 * 
 */
class TorsionalMinimizerCrossover implements Operator {
	private final TorsionalMinimizer problem;

	private final Chromosome children[] = new Chromosome[2];

	private static final Logger logger = Logger
			.getLogger(TorsionalMinimizerCrossover.class);

	TorsionalMinimizerCrossover(TorsionalMinimizer p) {
		problem = p;
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
		TorsionalMinimizerChromosome child1, child2, parent1, parent2;

		child1 = new TorsionalMinimizerChromosome(problem);
		child2 = new TorsionalMinimizerChromosome(problem);

		child1.createEmpty();
		child2.createEmpty();
		parent1 = (TorsionalMinimizerChromosome) parents[0];
		parent2 = (TorsionalMinimizerChromosome) parents[1];

		if (logger.isDebugEnabled()) {
			logger.debug("Doing Crossover");
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

/**
 * Mutation operator for torsional minimizer.
 * 
 */
class TorsionalMinimizerMutation implements Operator {
	private final TorsionalMinimizer problem;

	private final Chromosome children[] = new Chromosome[1];

	private static final boolean DEBUG = false;

	TorsionalMinimizerMutation(TorsionalMinimizer p) {
		problem = p;
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
		TorsionalMinimizerChromosome parent, child;

		child = new TorsionalMinimizerChromosome(problem);

		child.createEmpty();
		parent = (TorsionalMinimizerChromosome) parents[0];
		parent.copyGene(child);
		child.mutate();

		if (DEBUG) {
			System.out.println("Mutating");
			System.out.print(parent.geneInfo());
			System.out.println("-->");
			System.out.print(child.geneInfo());
		}

		children[0] = child;
		return children;
	}
}
