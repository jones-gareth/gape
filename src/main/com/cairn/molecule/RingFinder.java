package com.cairn.molecule;

/**
 * Uses a spanning tree algorithm to find all cyclic atoms and bonds in a
 * structure.
 * 
 * @author Gareth Jones
 * 
 */
public class RingFinder {
	private final Molecule molecule;
	private boolean spanningTreeSet[];
	private int atomTree[];
	private int bondTab[][], bondPos;
	private static final boolean IN_SET = true, OUT_SET = false;

	public RingFinder(Molecule mol) {
		molecule = mol;
	}

	public void findRings() {
		for (Atom atom : molecule.getAtoms()) {
			atom.setNotInRing();
		}
		for (Bond bond : molecule.getBonds()) {
			bond.setNotInRing();
		}

		int nAtoms = molecule.getnAtoms();
		spanningTreeSet = new boolean[nAtoms];
		atomTree = new int[nAtoms];
		bondTab = new int[2][nAtoms];
		int connectedAtoms[] = new int[12];
		for (int i = 0; i < nAtoms; i++)
			spanningTreeSet[i] = IN_SET;
		int atomNo = 0;
		while (atomNo < nAtoms) {
			if (spanningTreeSet[atomNo] == IN_SET) {
				spanningTreeSet[atomNo] = OUT_SET;
				atomTree[0] = atomNo;
				int nConnected = getConnectedAtoms(atomNo, connectedAtoms);
				for (int i = 0; i < nConnected; i++) {
					if (spanningTreeSet[connectedAtoms[i]] == IN_SET)
						growTree(connectedAtoms[i], atomNo, 1);
				}
			}
			atomNo++;
		}
	}

	void growTree(int thisAtom, int previousAtom, int depth) {
		atomTree[depth] = thisAtom;
		spanningTreeSet[thisAtom] = OUT_SET;
		int connectedAtoms[] = new int[12];
		int nConnected = getConnectedAtoms(thisAtom, connectedAtoms);
		for (int i = 0; i < nConnected; i++) {
			int nextAtom = connectedAtoms[i];
			if (spanningTreeSet[nextAtom] == IN_SET)
				growTree(nextAtom, thisAtom, depth + 1);
			else if (nextAtom != previousAtom) {
				boolean duplicateBond = false;
				if (bondPos > 0)
					for (int posn = 0; posn < bondPos; posn++)
						if ((nextAtom == bondTab[0][posn] && thisAtom == bondTab[1][posn])
								|| (nextAtom == bondTab[1][posn] && thisAtom == bondTab[0][posn])) {
							duplicateBond = true;
							break;
						}
				if (!duplicateBond) {
					bondTab[0][bondPos] = nextAtom;
					bondTab[1][bondPos++] = thisAtom;
					setRingFlag(nextAtom, thisAtom);
					int node1 = nextAtom;
					int node2 = thisAtom;
					int backPos = depth;
					while (node2 != nextAtom && backPos > 0) {
						node1 = node2;
						node2 = atomTree[--backPos];
						setRingFlag(node1, node2);
					}
				}
			}
		}
	}

	void setRingFlag(int atom1No, int atom2No) {
		Atom atom1 = molecule.getAtom(atom1No);
		Atom atom2 = molecule.getAtom(atom2No);
		atom1.setInRing(true);
		atom2.setInRing(true);
		Bond bond = molecule.getBond(atom1, atom2);
		bond.setInRing(true);
	}

	int getConnectedAtoms(int atomNo, int connectedAtoms[]) {
		Atom atom = molecule.getAtom(atomNo);
		int no = 0;
		for (Atom neighbourAtom : atom.getNotDummyNeighbours()) {
			// Note that we need this test for
			// com.cairn.gape.PatternMatch queries
			// TODO - not sure about this
			if (neighbourAtom == null)
				break;
			connectedAtoms[no] = neighbourAtom.getNo();
			no++;
		}
		return no;
	}
}
