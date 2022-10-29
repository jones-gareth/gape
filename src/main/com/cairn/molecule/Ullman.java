package com.cairn.molecule;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.log4j.Logger;

/**
 * Java version of Ullman Algorithm for subgraph isomorphism.
 * 
 * @author Gareth Jones
 * 
 */
public class Ullman {
	private static final Logger logger = Logger.getLogger(Ullman.class);

	static {
		// logger.setLevel(Level.DEBUG);
	}

	private volatile boolean ullMat[][], aMat[][], bMat[][];

	private volatile int pAlpha, pBeta, isomorphisms;

	// molA is the query and molB the target
	private Molecule molA, molB;

	private volatile List<Atom> atomsA, atomsB;

	public enum MatchType {
		MATCH, FRAG, FRAGDU, PATTERN, FRAG_HEAVY
	};

	private MatchType matchType = MatchType.FRAG;

	private volatile boolean heavy = false, subGraph = true;

	// make aromatic and sp2 atoms equivalent
	private volatile boolean ignoreAromatic = false;

	// ignore hydrogens when setting up adjacency matrix- for matches only
	private volatile boolean ignoreHydrogensInMatch = true;

	// connect on all neighbours (or just all not-dummy neighbours)
	private volatile boolean connectAllNeighbours;

	// interface class for handling matches
	UllmanCallback ullCallback;

	// STATIC_STORAGE is an attempt to reduce overhead by reusing allocated
	// arrays between searches- actually seems to make things a little worse.
	// Also not thread safe!
	private static final boolean STATIC_STORAGE = false;

	// match on elemental types rather than specific atom types
	private boolean matchElementalTypes = false;

	private static final ThreadLocal<boolean[][]> adjA = new ThreadLocal<boolean[][]>();
	private static final ThreadLocal<boolean[][]> adjB = new ThreadLocal<boolean[][]>();
	private static final ThreadLocal<Integer> adjASize = ThreadLocal.withInitial(() -> 0);
	private static final ThreadLocal<Integer> adjBSize = ThreadLocal.withInitial(() -> 0);

	/**
	 * Creates Ullman depth first for two molecules
	 * 
	 * @param a
	 * @param b
	 */
	public Ullman(Molecule a, Molecule b) {
		molA = a;
		molB = b;
	}

	/**
	 * Creates Ullman depth first for two molecules, setting the match type.
	 * 
	 * @param a
	 * @param b
	 * @param t
	 */
	public Ullman(Molecule a, Molecule b, MatchType t) {
		this(a, b);
		matchType = t;
	}

	/**
	 * @param ignoreAromatic
	 *            the ignoreAromatic to set
	 */
	public void setIgnoreAromatic(boolean ignoreAromatic) {
		this.ignoreAromatic = ignoreAromatic;
	}

	/**
	 * @return the ignoreAromatic
	 */
	public boolean isIgnoreAromatic() {
		return ignoreAromatic;
	}

	/**
	 * @return the ignoreHydrogensInMatch
	 */
	public boolean isIgnoreHydrogensInMatch() {
		return ignoreHydrogensInMatch;
	}

	/**
	 * @param ignoreHydrogensInMatch
	 *            the ignoreHydrogensInMatch to set
	 */
	public void setIgnoreHydrogensInMatch(boolean ignoreHydrogensInMatch) {
		this.ignoreHydrogensInMatch = ignoreHydrogensInMatch;
	}

	/**
	 * Setup the adjacancy matrix for a molecule.
	 * 
	 * @param molecule
	 * @return
	 */
	private boolean[][] createAdjacencyMatrix(Molecule molecule, boolean second) {
		List<Atom> atoms = null;
		int no;
		if (heavy) {
			atoms = molecule.getHeavyAtoms();
			no = atoms.size();
		} else {
			atoms = molecule.getAtoms();
			no = molecule.getnAtoms();
		}

		boolean matrix[][] = null;

		if (STATIC_STORAGE) {
			if (!second) {
				if (no > adjASize.get()) {
					adjASize.set(no);
					adjA.set(new boolean[no][no]);
				}
			} else {
				if (no > adjBSize.get()) {
					adjBSize.set(no);
					adjB.set(new boolean[no][no]);
				}
			}

			matrix = second ? adjB.get() : adjA.get();
		} else
			matrix = new boolean[no][no];

		for (int i = 0; i < no; i++)
			for (int j = i + 1; j < no; j++) {
				if (atoms.get(i).isNotDummyNeighbour(atoms.get(j)))
					matrix[i][j] = matrix[j][i] = true;
				else if (STATIC_STORAGE)
					matrix[i][j] = matrix[j][i] = false;
			}
		return matrix;
	}

	private static final ThreadLocal<boolean[][]> ullMat2 = new ThreadLocal<boolean[][]>();
	private static final ThreadLocal<Integer> ullMat2Size = new ThreadLocal<Integer>() {
		@Override
		protected Integer initialValue() {
			return 0;
		}
	};

	/**
	 * Setup the ullman matrix.
	 * 
	 */
	private void createUllmanMatrix() {

		if (heavy) {
			atomsA = molA.getHeavyAtoms();
			atomsB = molB.getHeavyAtoms();
			pAlpha = atomsA.size();
			pBeta = atomsB.size();
		} else {
			atomsA = molA.getAtoms();
			atomsB = molB.getAtoms();
			pAlpha = molA.getnAtoms();
			pBeta = molB.getnAtoms();
		}

		if (atomsA.size() > atomsB.size())
			throw new RuntimeException("Query is larger than target");

		boolean matrix[][] = null;

		if (STATIC_STORAGE) {
			int size = pAlpha > pBeta ? pAlpha : pBeta;
			if (size > ullMat2Size.get()) {
				ullMat2.set(new boolean[size][size]);
				ullMat2Size.set(size);
			}
			matrix = ullMat2.get();
		} else
			matrix = new boolean[pAlpha][pBeta];

		// Set matrix[i][j] if atom i from the query and atom j from
		// the target have the same atom type and atom j has the same or
		// more neighbours as atom i.

		for (int i = 0; i < pAlpha; i++) {
			for (int j = 0; j < pBeta; j++) {

				// Connect on all neighbors or just not-dummy neighbors

				Atom atomA = atomsA.get(i);
				Atom atomB = atomsB.get(j);

				int na = connectAllNeighbours ? atomA.getnNeighbours()
						: atomA.getnNotDummyNeighbours();
				int nb = connectAllNeighbours ? atomB.getnNeighbours()
						: atomB.getnNotDummyNeighbours();

				if (heavy && ignoreHydrogensInMatch) {
					na -= atomA.getHydrogenCount();
					nb -= atomB.getHydrogenCount();
				}

				if (matchType == MatchType.PATTERN
						&& ((QueryAtom) atomA).getnSlnHydrogens() != null) {
					na += ((QueryAtom) atomA).getnSlnHydrogens();
				}

				AtomType type1 = atomA.getType();
				AtomType type2 = atomB.getType();

				if (STATIC_STORAGE)
					matrix[i][j] = false;

				// na == nb for graph isomorphism. na <= nb for subgraph
				// isomorphism

				// can use either matchType for more exact matches or
				// matchElementalType for atom only matches

				boolean match = matchElementalTypes
						? type1.matchElementalType(type2, ignoreAromatic)
						: type1.matchType(type2, ignoreAromatic);

				if (!subGraph && match && na == nb) {
					matrix[i][j] = true;
					if (atomA.isInRing() != atomB.isInRing())
						matrix[i][j] = false;
				} else if (subGraph && match && na <= nb) {
					matrix[i][j] = true;
					if (atomA.isInRing() && !atomB.isInRing())
						matrix[i][j] = false;
				}

			}
		}

		ullMat = matrix;
	}

	/**
	 * @return the matchElementalTypes
	 */
	public boolean isMatchElementalTypes() {
		return matchElementalTypes;
	}

	/**
	 * @param matchElementalTypes
	 *            the matchElementalTypes to set
	 */
	public void setMatchElementalTypes(boolean matchElementalTypes) {
		this.matchElementalTypes = matchElementalTypes;
	}

	/**
	 * Subclass this to do something useful
	 * 
	 */
	public interface UllmanCallback {
		/**
		 * This is called when we find a match.
		 * 
		 * @param match
		 *            match[i][j] is set when atom i in molA maps to atom j in
		 *            molB
		 */
		// void callback(boolean match[][]);

		/**
		 * This is called when we find a match.
		 * 
		 * @param atomMapping
		 *            A mapping of matching atoms in the two molecules.
		 */
		void callback(Map<Atom, Atom> atomMapping);
	}

	/**
	 * Prints out match - used in test callback
	 * 
	 * @param m
	 */
	private void printMatch(Map<Atom, Atom> matchMap) {
		System.out.println("Isomorphism " + isomorphisms);
		for (Entry<Atom, Atom> entry : matchMap.entrySet()) {
			System.out.println(
					"Atom " + entry.getKey().getNo() + "-->" + entry.getValue().getNo());
		}
	}

	/**
	 * Simple test application that searches a query molecule(s) against a
	 * target molecule(s).
	 * 
	 * @param args
	 */
	public static void main(String args[]) {
		if (args.length < 2) {
			System.out.println("usage: ullman <query.mol2> <target.mol2> [MATCH|FRAG]");
			System.exit(0);
		}

		List<Molecule> queryMols = Molecule.loadFiles(new String[] { args[0] });
		List<Molecule> targetMols = Molecule.loadFiles(new String[] { args[1] });
		assert queryMols.size() == targetMols
				.size() : "Mismatch between query and target numbers";

		MatchType matchType = MatchType.FRAG;
		if (args.length == 3) {
			if (args[2].toUpperCase().equals("MATCH")) {
				matchType = MatchType.MATCH;
				System.out.println("EXACT MATCH");
			} else if (args[2].toUpperCase().equals("FRAG")) {
				matchType = MatchType.FRAG;
				System.out.println("SUBSTRUCTURE MATCH");
			}
		}

		for (int no = 0; no < queryMols.size(); no++) {
			Molecule query = queryMols.get(no);
			Molecule target = targetMols.get(no);
			System.out.println("Matching " + no + " query " + query.getName()
					+ " against target " + target.getName());
			Ullman ull = new Ullman(query, target, matchType);
			ull.doUllman();
			System.out.println(
					"Match " + no + ": " + ull.getIsomorphisms() + " isomorphisms found");
		}
	}

	/**
	 * doUllman with printMatch callback
	 * 
	 * @throws UllmanException
	 */
	public void doUllman() {
		doUllman(new UllmanCallback() {
			@Override
			public void callback(Map<Atom, Atom> matchMap) {
				printMatch(matchMap);
			}
		});
	}

	private int _maxPBeta;

	private boolean _k[];

	/**
	 * This is the entry point to the Ullman algorithm. Mola is the query
	 * molecule and molb the target. Callback.callback is called when a match is
	 * found.
	 * 
	 * @param callback
	 * @throws UllmanException
	 */
	public synchronized void doUllman(UllmanCallback callback) {
		ullCallback = callback;

		// Initilialization
		if (matchType == MatchType.MATCH) {
			// ULL_MATCH just matches on heavy atoms only
			heavy = true;
			subGraph = false;
			connectAllNeighbours = false;
		}

		else if (matchType == MatchType.FRAG || matchType == MatchType.FRAGDU) {
			// match on all atoms
			heavy = false;
			subGraph = true;
			// if ULL_FRAGDU is required count dummy atoms in neighbour
			if (matchType == MatchType.FRAGDU)
				connectAllNeighbours = true;
			else
				connectAllNeighbours = false;
		}

		else if (matchType == MatchType.FRAG_HEAVY) {
			// subgraph matching on heavy atoms only
			heavy = true;
			subGraph = true;
			connectAllNeighbours = false;
		}

		else if (matchType == MatchType.PATTERN) {
			// rules for SLN matching using molecule from ga.MolPattern
			heavy = false;
			subGraph = true;
			connectAllNeighbours = false;
		}

		aMat = createAdjacencyMatrix(molA, false);
		bMat = createAdjacencyMatrix(molB, true);
		createUllmanMatrix();

		boolean k[] = null;

		if (STATIC_STORAGE) {
			if (pBeta > _maxPBeta) {
				_k = new boolean[pBeta];
				_maxPBeta = pBeta;
			}
			k = _k;
			for (int i = 0; i < pBeta; i++)
				k[i] = false;
		} else
			k = new boolean[pBeta];

		isomorphisms = 0;

		/* do the business */
		depthFirst(0, ullMat, k);

		logger.debug(isomorphisms + " Isomorphisms found");

	}

	/**
	 * Copy matrix old to n.
	 * 
	 * @param n
	 * @param old
	 */
	private void copyMatrix(boolean n[][], boolean old[][], int nRows, int nCols) {
		for (int i = 0; i < nRows; i++)
			for (int j = 0; j < nCols; j++)
				n[i][j] = old[i][j];

	}

	/**
	 * Creates an atom to atom mapping for a match from the correspondence
	 * matrix. In the returned map a key is a query atom and a value is the
	 * target atom the query is mapped to. Map is ordered on query atoms
	 * 
	 * @param m
	 * @return
	 */
	public Map<Atom, Atom> getMatchMapping(boolean m[][]) {

		Map<Atom, Atom> match = new LinkedHashMap<>();
		for (int i = 0; i < pAlpha; i++) {
			Atom atomA = atomsA.get(i);
			for (int j = 0; j < pBeta; j++) {
				if (m[i][j]) {
					Atom atomB = atomsB.get(j);
					match.put(atomA, atomB);
				}
			}
		}
		return match;
	}

	// depthFirst(), refine(), rowAndColumn(),
	// corresponds(), allFalse()

	// are java versions of published pseudocode. If you want to know
	// more about these then have a look at some of Pete Willetts
	// publications

	private void depthFirst(int d, boolean m[][], boolean f[]) {

		logger.debug("DepthFirst depth " + d);

		boolean finish = false;
		boolean temp1[][] = new boolean[pAlpha][pBeta];

		copyMatrix(temp1, m, pAlpha, pBeta);
		int node = 0;

		for (node = 0; node < pBeta; node++) {

			if (m[d][node] && !f[node]) {
				// m at this depth and node has not been visited */

				logger.debug("checking node " + node + " at depth " + d);

				for (int j = 0; j < pBeta; j++)
					m[d][j] = false;

				m[d][node] = true;
				f[node] = true;

				finish = refine(m);

				if (!finish) {
					if (d == (pAlpha - 1)) {
						isomorphisms++;
						ullCallback.callback(getMatchMapping(m));
					} else {
						depthFirst(d + 1, m, f);
					}
				}

				f[node] = false;

				copyMatrix(m, temp1, pAlpha, pBeta);
			}
		}
	}

	private boolean refine(boolean m[][]) {
		logger.debug("refining ..");
		boolean noChange = false;

		while (!noChange) {
			noChange = true;

			for (int i = 0; i < pAlpha; i++) {
				for (int j = 0; j < pBeta; j++) {

					if (m[i][j] && !corresponds(m, i, j)) {
						logger.debug("Setting m at " + i + ", " + j + " false");
						m[i][j] = false;
						noChange = false;
						if (allFalse(m, i)) {
							return true;
						}
					}
				}
			}
		}
		return false;
	}

	private boolean rowAndColumn(boolean m[][], int x, int j) {
		for (int y = 0; y < pBeta; y++)
			if (m[x][y] && bMat[y][j])
				return true;

		return false;
	}

	boolean corresponds(boolean m[][], int i, int j) {
		for (int x = 0; x < pAlpha; x++)
			if (aMat[i][x] && !rowAndColumn(m, x, j)) {
				logger.debug("Returning false for i " + i + ", j " + j + ", x " + x);
				return false;
			}

		return true;
	}

	boolean allFalse(boolean m[][], int row) {
		for (int col = 0; col < pBeta; col++)
			if (m[row][col])
				return false;

		return true;
	}

	/**
	 * @return the molA
	 */
	public Molecule getMolA() {
		return molA;
	}

	/**
	 * @return the molB
	 */
	public Molecule getMolB() {
		return molB;
	}

	/**
	 * @return the isomorphisms
	 */
	public int getIsomorphisms() {
		return isomorphisms;
	}

}
