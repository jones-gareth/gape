package com.cairn.molecule;

import java.util.ArrayList;
import java.util.List;

/**
 * Solves the all-pairs shortest path problem using Floyd-Warshall algorithm.
 * 
 * Based on code from Paul Watson, but added predecessor array for enumerating
 * shortest path.
 * 
 * @author Gareth Jones
 * 
 */
public class ShortestPath {
	// path lengths
	private int pathLengths[][];

	// predecessors
	private int predecessors[][];

	// current molecule
	private Molecule molecule;

	public ShortestPath() {
		;
	}

	/**
	 * Set molecule and determine shorstest path lengths.
	 * 
	 * @param molecule
	 */
	public void setMolecule(Molecule molecule) {
		this.molecule = molecule;
		shortestPath();
	}

	/**
	 * Solves the all-pairs shortest path problem using Floyd-Warshall
	 * algorithm.
	 */
	private void shortestPath() {

		int nAtoms = molecule.getnAtoms();
		pathLengths = new int[nAtoms][nAtoms];
		predecessors = new int[nAtoms][nAtoms];

		// Initialize path lengths
		for (int i = 0; i < nAtoms; i++) {
			for (int j = 0; j < nAtoms; j++) {
				if (i == j)
					pathLengths[i][j] = 0;
				else
					pathLengths[i][j] = nAtoms;
			}
		}

		// Starting lengths of 1 from bond table.
		for (Bond bond : molecule.getBonds()) {
			int atom1 = bond.getAtom1().getNo();
			int atom2 = bond.getAtom2().getNo();
			pathLengths[atom1][atom2] = 1;
			pathLengths[atom2][atom1] = 1;
			predecessors[atom1][atom2] = atom1;
			predecessors[atom2][atom1] = atom2;
		}

		// Algorithm loop
		for (int k = 0; k < nAtoms; k++) {
			for (int i = 0; i < nAtoms; i++) {
				for (int j = 0; j < nAtoms; j++) {
					if ((pathLengths[i][k] + pathLengths[k][j]) < pathLengths[i][j]) {
						pathLengths[i][j] = pathLengths[i][k] + pathLengths[k][j];
						predecessors[i][j] = predecessors[k][j];
					}
				}
			}
		}
	}

	/**
	 * @param i
	 * @param j
	 * @return The shortest path between two atoms.
	 */
	public int getPathLength(int i, int j) {
		return pathLengths[i][j];
	}

	/**
	 * @return The array of path lengths
	 */
	public int[][] getPathLengths() {
		return pathLengths;
	}

	/**
	 * @param i
	 * @param j
	 * @return The fully enumerated shortest path between two atoms.
	 */
	public List<Integer> getPath(int i, int j) {
		List<Integer> path = new ArrayList<Integer>();
		if (i == j)
			return path;
		path.add(j);
		getPath(i, j, path);
		assert path.size() == (getPathLength(i, j) + 1);
		return path;
	}

	/**
	 * Uses the predecessor array to work back one step on the shortest path.
	 * 
	 * @param start
	 * @param current
	 * @param path
	 */
	private void getPath(int start, int current, List<Integer> path) {
		int predecessor = predecessors[start][current];
		path.add(predecessor);
		if (predecessor == start)
			return;
		getPath(start, predecessor, path);
	}

	/**
	 * @return a summary of all shortest paths in the molecule.
	 */
	public String pathInfo() {
		StringBuffer string = new StringBuffer();
		int nAtoms = molecule.getnAtoms();
		for (int i = 0; i < nAtoms; i++) {
			for (int j = i + 1; j < nAtoms; j++) {
				int pathLength = getPathLength(i, j);
				List<Integer> path = getPath(i, j);
				string.append(i + 1).append(",").append(j + 1).append(": ")
						.append(pathLength).append(" [");
				boolean start = true;
				for (int no : path) {
					if (!start)
						string.append(",");
					start = false;
					string.append(no + 1);
				}
				string.append("]\n");
			}
		}

		return string.toString();
	}

	/**
	 * Prints a summary of all shortest paths in the molecule.
	 * 
	 * @see #pathInfo()
	 */
	public void printInfo() {
		System.out.println(pathInfo());
	}

	/**
	 * Test method. Reads structure files and displays shortest path
	 * information.
	 * 
	 * @param args
	 */
	public static void main(String args[]) {
		if (args.length == 0) {
			System.out.println("usage: ShortestPath <molecule files..>");
			System.exit(0);
		}

		List<Molecule> molecules = Molecule.loadFiles(args);

		ShortestPath shortestPath = new ShortestPath();

		for (Molecule molecule : molecules) {
			System.out.print("\n\n" + molecule.getName() + "\n\n");
			shortestPath.setMolecule(molecule);
			shortestPath.printInfo();
		}
	}

}
