package com.cairn.gape.utils;

import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.StringTokenizer;

/**
 * Uses public domain source code from
 * ftp://dimacs.rutgers.edu/pub/challenge/graph/solvers/dfmax.c
 * 
 * dfmax.c Semi-Exhaustive Greedy Independent Set from Matula and Johri
 * 
 * dfmax.c comments start:
 * -------------------------------------------------------
 * 
 * written October 22, 1983 by DSJ modified February 20, 1987 for checkpointing
 * modified August 1988 to find maximum independent sets modified April 1993 for
 * new data structures, bounding modified September 1993 for dimacs .b input
 * 
 * ------------------------------------------------------- dfmax.c comments end
 * 
 * This code is an adaptation of the Carraghan and Pardalos algorithm. Adapted
 * for Java by GJ. Tidied up globals and objectfied- changed fortran style
 * indices for input and output arrays. Java is about 50% slower than C on test
 * cases (JDK 1.5.0 vs Cygwin). Code reproduces results for test files.
 * 
 * @author Gareth Jones
 * @version 0.1
 */
public class Clique {

	static final boolean REORDER = true, DEBUG = false;

	// graph input parameters

	// Number of vertices
	private int nrVert;

	// Edges
	private boolean edge[][];

	// size at which exhaustive search begins
	private int setlim = 1;

	// array in which vertices reside and are moved
	private int vertex[];

	// vertex degrees with resp. to uncolored vertices
	private int degree[];

	// non-zero entries are members of current set
	private int set[];

	// non-zero entries constitute currently best ind. set
	private int bestset[];

	// size of current best set
	private int bestsize;

	private boolean verbose = false;

	public Clique() {
	}

	/**
	 * Sets up clique detection
	 * 
	 * @param e
	 *            edge matrix or undirected and unweighted graph.
	 * @param nrVert
	 *            number of vertices
	 */
	public void setEdges(boolean e[][], int nrVert) {
		edge = e;
		init(nrVert);
	}

	private int maxSize;

	/**
	 * Sets up storage
	 * 
	 * @param size
	 *            Number of vertices
	 */
	public void init(int size) {
		nrVert = size;

		if (size <= maxSize)
			return;
		maxSize = size;
		vertex = new int[size + 1];
		degree = new int[size + 1];
		set = new int[size];
		bestset = new int[size];
	}

	/**
	 * Main routine runs the same tests as dfmax.c (and produces similar
	 * output!). See README and data files at
	 * ftp://dimacs.rutgers.edu/pub/challenge/graph/solvers/
	 * 
	 * @param args
	 *            Data file name and optional minimum clique size.
	 */
	public static void main(String args[]) {
		Clique clique = new Clique();

		// read input
		if (args.length < 1 || args.length > 2) {
			System.out.println("Usage: dfmax <filename> [setlim]");
			System.exit(1);
		}
		DataInputStream inputfile = null;
		try {
			inputfile = new DataInputStream(new FileInputStream(args[0]));
		} catch (FileNotFoundException ex) {
			System.out.println("Can't open " + args[0] + " " + ex);
			System.exit(0);
		}

		try {
			clique.readGraph(inputfile);
		} catch (IOException ex) {
			ex.printStackTrace(System.err);
			// System.err.println(ex);
			System.exit(0);
		}

		/* initialize algorithmic parameters */
		clique.setlim = 1;
		if (args.length > 1) {
			clique.setlim = Integer.valueOf(args[1]).intValue();
			System.out.println("setlim=" + clique.setlim);
		}

		System.out.println("Loaded data");

		clique.verbose = false;
		clique.solve();
		clique.printBest();

		System.out.println("done");

	}

	/**
	 * Prints out largest clique found
	 */
	@Deprecated
	public void printBest() {
		System.out.println("Best:");
		for (int i = 0; i < bestsize; i++) {
			System.out.print(" " + String.valueOf(bestset[i] + 1));
		}
		System.out.println();
	}

	/**
	 * @return information about the largest clique
	 */
	public String bestInfo() {
		StringBuilder sb = new StringBuilder();
		sb.append("Best:");
		for (int i = 0; i < bestsize; i++) {
			sb.append(" ").append(String.valueOf(bestset[i] + 1));
		}
		return sb.toString();
	}

	/**
	 * Sets up vertex and degree arrays. Mainly original code from dfmax.c
	 */
	void setVertices() {
		if (REORDER) {

			for (int i = 1; i <= nrVert; i++) {
				degree[i] = 0;
				for (int j = 1; j <= nrVert; j++)
					if (!edge[i - 1][j - 1])
						degree[i]++;
			}

			int dmax = -1;
			int cand = 0;

			for (int i = 1; i <= nrVert; i++)
				if (degree[i] > dmax) {
					dmax = degree[i];
					cand = i;
				}
			vertex[nrVert] = cand;

			int newcand = 0;
			for (int j = nrVert - 1; j >= 1; j--) {
				degree[cand] = -9;
				dmax = -1;
				for (int i = 1; i <= nrVert; i++) {
					if (!edge[cand - 1][i - 1])
						degree[i]--;
					if (degree[i] > dmax) {
						dmax = degree[i];
						newcand = i;
					}
				}
				vertex[j] = cand = newcand;
			}
		} else {
			for (int j = 1; j <= nrVert; j++)
				vertex[j] = j;
		}

	}

	/**
	 * Finds the maximum clique and returns it's size. The clique is in the
	 * bestset array.
	 */
	public int solve() {
		setVertices();
		bestsize = 0;
		bestsize = maxind(nrVert, setlim, vertex, 1);
		return bestsize;
	}

	/**
	 * Maximum independent set code from dfmax.c.
	 */
	public int maxind(int top, int goal, int array[], int depth) {
		int newarray[] = new int[nrVert];
		int i, w, z;
		int best, restbest, newgoal;
		int canthrow;

		if (top <= 1) {
			if (top == 0)
				depth--;
			if (depth > bestsize) {
				bestsize = depth;
				if (top == 1)
					set[bestsize - 1] = array[top] - 1;
				for (i = 0; i < bestsize; i++)
					bestset[i] = set[i];
				if (DEBUG)
					check();
				if (verbose)
					System.out.println("Size = " + bestsize + " found");
			}
			return (top);
		}
		best = 1;
		newgoal = goal - 1;
		if (newgoal <= 1)
			newgoal = 1;
		outer: for (i = top; i >= goal; i--) {
			int pnew = 0;
			w = array[i];
			set[depth - 1] = w - 1;
			canthrow = i - goal;
			int pold = 1;
			while (pold < i) {
				z = array[pold++];
				if (edge[z - 1][w - 1]) {
					newarray[++pnew] = z;
				} else {
					if (canthrow == 0)
						continue outer;
					canthrow--;
				}
			}
			restbest = maxind(pnew, newgoal, newarray, depth + 1);
			if (restbest >= newgoal) {
				best = newgoal = restbest + 1;
				goal = best + 1;
			}
			if (top == nrVert && verbose) {
				System.out.println("N = " + i + " best = " + best);
			}
		}
		return (best);
	}

	byte bitmap[][];

	char masks[] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };

	/**
	 * Routine for extracting edges from binary data files.
	 */
	private boolean getEdge(int i, int j) {
		if (i < j) {
			int k = i;
			i = j;
			j = k;
		}

		int bit = 7 - (j & 0x00000007);
		int byte_ = j >>> 3;

		char mask = masks[bit];
		return ((bitmap[i][byte_] & mask) == mask);
	}

	/**
	 * Reads a binary test file. Dimacs format (whatever that is!).
	 * 
	 * @param input
	 *            Stream to data file.
	 */
	public void readGraph(DataInputStream input) throws IOException {
		String line = input.readLine();

		int length = Integer.valueOf(line).intValue();
		char preambleArray[] = new char[length];

		// getting nrVert and Nr_edge from the preamble string,
		// containing Dimacs format "p ??? num num"

		for (int i = 0; i < length; i++)
			preambleArray[i] = (char) input.readByte();

		String preamble = new String(preambleArray);
		System.out.println("Preamble " + preamble);

		StringTokenizer tk = new StringTokenizer(preamble);
		tk.nextToken();
		tk.nextToken();
		nrVert = Integer.valueOf(tk.nextToken()).intValue();
		int nrEdges = Integer.valueOf(tk.nextToken()).intValue();
		System.out.println("N vertices " + nrVert + " N edges" + nrEdges);

		edge = new boolean[nrVert + 1][nrVert + 1];
		init(nrVert);
		boolean extra = nrVert % 8 > 0 ? true : false;
		int s = nrVert / 8;
		if (extra)
			s++;
		bitmap = new byte[nrVert][s];

		try {
			for (int i = 0; i < nrVert; i++) {
				input.readFully(bitmap[i], 0, (i + 8) / 8);
			}
		} catch (IOException ex) {
			System.err.println(ex);
			System.exit(0);
		}

		input.close();

		for (int i = 0; i < nrVert; i++)
			for (int j = 0; j < nrVert; j++)
				if (getEdge(i, j))
					edge[i][j] = true;

	}

	/**
	 * Checks clique in bestset is actually a clique.
	 */
	void check() {
		// checks that the vertices with color K are an independent set
		for (int i = 0; i < bestsize; i++)
			for (int j = i + 1; j < bestsize; j++)
				if (!edge[bestset[i]][bestset[j]])
					System.out.println("Edge in set\n");
	}

	public int getSetlim() {
		return setlim;
	}

	public void setSetlim(int setlim) {
		this.setlim = setlim;
	}

	public boolean isVerbose() {
		return verbose;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	// these getters return the size and the nodes of the maximum clique

	public int[] getBestset() {
		return bestset;
	}

	public int getBestset(int i) {
		return bestset[i];
	}

	public int getBestsize() {
		return bestsize;
	}

}