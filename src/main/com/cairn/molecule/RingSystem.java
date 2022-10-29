package com.cairn.molecule;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.log4j.Logger;

import com.cairn.molecule.Sssr.AtomNode;

/**
 * Ring system class for SSSR
 * 
 * @author Gareth Jones
 * 
 */
class RingSystem {
	private int nAtoms;
	private int nEdges;
	private int edges[];
	private int atoms[];
	private boolean usage[];
	private int edgeUsage[];
	private int best[];
	private int ringSize;
	private final Molecule molecule;
	private AtomNode graph[];
	private Sssr sssr;

	private static final Logger logger = Logger.getLogger(RingSystem.class);

	RingSystem(Molecule m, int start) {
		nAtoms = nEdges = 0;
		atoms = new int[m.getnAtoms()];
		edges = new int[m.getnBonds()];
		atoms[0] = start;
		this.molecule = m;
	}

	private final List<Ring> phase3Rings = new ArrayList<>();

	// int nPhase3Rings = 0;

	void sssrPhases(Sssr sssr) {
		this.sssr = sssr;
		graph = sssr.getGraph();
		int first = -1, second = -1;
		usage = new boolean[molecule.getnAtoms()];
		edgeUsage = new int[nEdges];
		best = new int[nAtoms];
		int nFound = 0;
		int nSssr = nEdges - nAtoms + 1;
		// Phase 1
		logger.debug("Phase 1");
		while (true) {
			int max = -1;
			first = -1;
			for (int i = 0; i < nAtoms; i++)
				if (!usage[atoms[i]]) {
					AtomNode node = graph[atoms[i]];
					if (node.getnConnections() > max) {
						max = node.getnConnections();
						first = atoms[i];
					}
				}
			if (first == -1)
				break;
			boolean foundRing = findSmallRing(first, -1, false, false);
			if (foundRing) {
				nFound++;
				addRing();
			} else {
				throw new RuntimeException("doSssr: failed to find phase 1 ring");
			}
		}
		// Phase 2
		logger.debug("Phase 2");
		int nAtoms = molecule.getnAtoms();
		for (int i = 0; i < nAtoms; i++)
			usage[i] = false;
		while (true) {
			int startEdgePos = -1;
			for (int i = 0; i < nEdges; i++)
				if (edgeUsage[i] == 0) {
					startEdgePos = i;
					break;
				}
			if (startEdgePos == -1)
				break;
			int edge = edges[startEdgePos];
			Bond bond = molecule.getBond(edge);
			first = bond.getAtom1().getNo();
			second = bond.getAtom2().getNo();
			boolean foundRing = findSmallRing(first, second, true, false);
			if (foundRing) {
				nFound++;
				addRing();
			} else {
				throw new RuntimeException("doSssr: failed to find phase 2 ring.");
			}
		}
		if (nFound == nSssr)
			return;

		// Phase 3 rings
		logger.debug("Phase 3");
		int nPhase3Rings = 0;
		while (true) {
			int startEdgePos = -1;
			for (int i = 0; i < nEdges; i++)
				if (edgeUsage[i] == 1) {
					startEdgePos = i;
					break;
				}
			if (startEdgePos == -1)
				break;
			logger.debug("Start Edge Pos " + startEdgePos);
			edgeUsage[startEdgePos] = 2;
			int edge = edges[startEdgePos];
			Bond bond = molecule.getBond(edge);
			first = bond.getAtom1().getNo();
			second = bond.getAtom2().getNo();
			boolean foundRing = findSmallRing(first, second, true, true);
			if (foundRing) {
				Ring r = getFoundRing();
				phase3Rings.add(r);
				nPhase3Rings++;
			}
		}

		if (nFound + nPhase3Rings < nSssr)
			throw new RuntimeException("doSssr: phase 3 failed to find enough rings");

		// Pick smallest phase3 rings
		logger.debug("Picking smallest Phase 3 ring");
		Collections.sort(phase3Rings,
				(o1, o2) -> Integer.compare(o1.getAtoms().size(), o2.getAtoms().size()));

		// Because of the relaxation of phase III staring points it's
		// possible the rings we find here are already present.

		int nAlreadyFound = nFound;
		for (int i = nAlreadyFound; i < nFound + nPhase3Rings; i++) {
			Ring r = phase3Rings.get(i - nAlreadyFound);

			// Check to see if ring is already
			boolean present = false;
			for (int j = 0; j < nAlreadyFound; j++) {
				if (sssr.getRing(j).sameRing(r)) {
					logger.debug("Phase 3 ring " + r.info() + " already in SSSR");
					present = true;
					break;
				}
			}
			if (present)
				continue;

			sssr.addRing(r);
			nFound++;
			logger.debug("Adding Phase 3 ring " + r.info());
			if (nFound == nSssr)
				break;
		}
		if (nFound < nSssr)
			throw new RuntimeException(
					"doSssr: phase 3 addition failed to find enough rings");

		sssr.setnRings(nSssr);
		logger.debug("SSSR size " + nSssr);
	}

	Ring getFoundRing() {
		logger.debug("Found Ring ");
		for (int i = 0; i < ringSize; i++) {
			usage[best[i]] = true;
			int pos = -1;
			if (i > 0) {
				pos = edgePosition(best[i], best[i - 1]);
			} else {
				pos = edgePosition(best[0], best[ringSize - 1]);
			}
			edgeUsage[pos]++;
		}

		Ring r = new Ring(ringSize, best, molecule);
		logger.debug(r.info());
		return r;
	}

	void addRing() {
		Ring r = getFoundRing();
		sssr.addRing(r);
	}

	int edgePosition(int a1, int a2) {
		AtomNode node = graph[a2];
		int edge = -1;
		int nConnections = node.getnConnections();
		for (int i = 0; i < nConnections; i++)
			if (node.getConnections()[i] == a1) {
				edge = node.getEdges()[i];
				break;
			}
		for (int i = 0; i < nEdges; i++)
			if (edges[i] == edge)
				return i;
		return -1;
	}

	boolean findSmallRing(int first, int second, boolean useSecond, boolean phase3)
			throws RuntimeException {
		int nAtoms = molecule.getnAtoms();
		int xAttachments[] = new int[nAtoms];
		int path[] = new int[nAtoms];
		boolean use[] = new boolean[nAtoms];
		for (int i = 0; i < nAtoms; i++) {
			use[i] = false;
			xAttachments[i] = -1;
		}
		int size = 0;
		ringSize = nAtoms + 1;
		path[0] = first;
		use[first] = true;
		int current = first;
		boolean ringFound = false;
		if (useSecond) {
			path[1] = second;
			xAttachments[first] = graph[first].getnConnections();
			size = 1;
			use[second] = true;
			current = second;
		}
		// Alpha !!
		while (true) {
			boolean flag = false;
			xAttachments[current]++;
			AtomNode node = graph[current];
			if (xAttachments[current] >= node.getnConnections()) {
				if (size == 0)
					break;
				flag = true; // scan
			} else {
				int attached = node.getConnections()[xAttachments[current]];
				if (size > 0 && attached == path[size - 1])
					continue;

				// Zamorra defines inside faces/rings with all atoms
				// having at least three cyclic connections. However,
				// analysis of the iResearch library shows that there
				// are many inside faces where some atoms have only
				// two connections - everything seems to work if we
				// take out this restriction.

				// if (phase3 && graph[attached].nConnections < 3)
				// continue;

				if (phase3 && edgeUsage[edgePosition(current, attached)] >= 2)
					continue;
				else if (use[attached]) {
					if (attached != first) // beta
						continue;
					else { // found ring
						if (size + 1 < ringSize) {
							// criteria for equal rings here
							ringSize = size + 1; // better
							for (int i = 0; i < ringSize; i++)
								best[i] = path[i];
							ringFound = true;
						}
					}
					flag = true; // scan
				} else {
					size++;
					path[size] = attached;
					use[attached] = true;
					current = attached;
				}
			}
			if (flag) { // scan or backtrak
				xAttachments[current] = -1;
				use[current] = false;
				size--;
				current = path[size];
				node = graph[current];
				if (size == 0)
					if ((xAttachments[current] + 1) >= node.getnConnections())
						break;
			}
		}
		return ringFound;
	}

	/**
	 * @return the nAtoms
	 */
	int getnAtoms() {
		return nAtoms;
	}

	/**
	 * @return the atoms
	 */
	int[] getAtoms() {
		return atoms;
	}

	/**
	 * @param nAtoms
	 *            the nAtoms to set
	 */
	void setnAtoms(int nAtoms) {
		this.nAtoms = nAtoms;
	}

	void addEdge(int edge) {
		edges[nEdges++] = edge;
	}

}
