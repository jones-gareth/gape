package com.cairn.molecule;

import java.util.ArrayList;

import org.apache.log4j.Logger;

//
//  File: sssr.c
//
//  Gareth Jones 1/02
//
//  This file contains sssr (smalllest set of smallest rings).  It's
//  based on an algorithm for Zamorra (REF), which is very fortran like
//  and pretty messy.  Also I didn't encode phase III, which is a ring
//  totally enclosed in other rings.
//
//  Having said all this everything always works, even on proteins.
//
public class Sssr {
	private final Molecule m;
	private AtomNode graph[];
	private RingSystem ringSystems[];
	private int nRingSystems;
	private boolean used[], usedEdges[];
	private int treePosition;
	private Ring rings[];
	private int nRings;

	private static final Logger logger = Logger.getLogger(Sssr.class);
	private ArrayList<Ring> ringArray;

	public Sssr(Molecule mol) {
		m = mol;
	}

	public void doSssr() {
		for (Atom atom : m.getAtoms()) {
			atom.getRings().clear();
		}
		for (Bond bond : m.getBonds()) {
			bond.getRings().clear();
		}

		int nAtoms = m.getnAtoms();
		graph = new AtomNode[nAtoms];
		for (int i = 0; i < nAtoms; i++)
			graph[i] = new AtomNode(m.getAtom(i), m);
		findRingSystems();
		ringArray = new ArrayList<Ring>();

		for (int i = 0; i < nRingSystems; i++)
			ringSystems[i].sssrPhases(this);

		rings = new Ring[ringArray.size()];
		rings = ringArray.toArray(rings);
		nRings = ringArray.size();

		if (logger.isDebugEnabled()) {
			logger.debug("SSSR----------");
			for (Ring ring : ringArray) {
				logger.debug("Ring " + ring.info());
			}
		}

		m.clearRings();
		for (int i = 0; i < nRings; i++) {
			m.addRing(rings[i]);
		}
	}

	void findRingSystems() {
		int nAtoms = m.getnAtoms();
		used = new boolean[nAtoms];
		usedEdges = new boolean[m.getnBonds()];
		nRingSystems = 0;
		ArrayList<RingSystem> ringSystemsArray = new ArrayList<RingSystem>();
		for (int i = 0; i < nAtoms; i++) {
			if (used[i])
				continue;
			used[i] = true;
			AtomNode node = graph[i];
			if (node.nConnections == 0)
				continue;
			treePosition = 1;
			RingSystem currentRs = new RingSystem(m, node.node);
			ringSystemsArray.add(currentRs);
			buildRingSystem(node, currentRs.getAtoms());
			currentRs.setnAtoms(treePosition);
			if (logger.isDebugEnabled()) {
				logger.debug("Found Ring System ");
				for (int j = 0; j < currentRs.getnAtoms(); j++)
					logger.debug(String.valueOf(currentRs.getAtoms()[j] + 1) + " ");
			}
			nRingSystems++;
		}

		ringSystems = new RingSystem[nRingSystems];
		ringSystems = ringSystemsArray.toArray(ringSystems);
		for (int i = 0; i < nRingSystems; i++)
			addEdgesToRingSystem(ringSystems[i]);
	}

	void buildRingSystem(AtomNode node, int tree[]) {
		for (int i = 0; i < node.nConnections; i++) {
			int test = node.connections[i];
			if (!used[test]) {
				used[test] = true;
				tree[treePosition] = test;
				treePosition++;
				buildRingSystem(graph[test], tree);
			}
		}
	}

	void addEdgesToRingSystem(RingSystem rs) {
		int nAtoms = rs.getnAtoms();
		for (int i = 0; i < nAtoms; i++) {
			AtomNode node = graph[rs.getAtoms()[i]];
			for (int j = 0; j < node.nConnections; j++) {
				int edge = node.edges[j];
				if (usedEdges[edge])
					continue;
				usedEdges[edge] = true;
				rs.addEdge(edge);
			}
		}
	}

	/**
	 * @return the a ring system
	 */
	public RingSystem getRingSystem(int i) {
		return ringSystems[i];
	}

	/**
	 * @return the nRingSystems
	 */
	public int getnRingSystems() {
		return nRingSystems;
	}

	/**
	 * @return the nRings
	 */
	int getnRings() {
		return nRings;
	}

	/**
	 * @param nRings
	 *            the nRings to set
	 */
	void setnRings(int nRings) {
		this.nRings = nRings;
	}

	void addRing(Ring ring) {
		ringArray.add(ring);
	}

	Ring getRing(int no) {
		return ringArray.get(no);
	}

	/**
	 * @return the graph
	 */
	AtomNode[] getGraph() {
		return graph;
	}

	static class AtomNode {
		private int node;
		private int nConnections;
		private int connections[];
		private int edges[];
		private Atom atom;

		AtomNode(Atom a, Molecule m) {
			atom = a;
			node = a.getNo();
			Atom neighbours[] = new Atom[12];
			nConnections = m.getRingNeighbours(atom, neighbours);
			connections = new int[nConnections];
			edges = new int[nConnections];
			for (int i = 0; i < nConnections; i++) {
				connections[i] = neighbours[i].getNo();
				Bond b = m.getBond(atom, neighbours[i]);
				edges[i] = b.getNo();
			}
		}

		/**
		 * @return the node
		 */
		int getNode() {
			return node;
		}

		/**
		 * @return the nConnections
		 */
		int getnConnections() {
			return nConnections;
		}

		/**
		 * @return the connections
		 */
		int[] getConnections() {
			return connections;
		}

		/**
		 * @return the edges
		 */
		int[] getEdges() {
			return edges;
		}

	}
}
