package com.cairn.gape.molecule;

import com.cairn.molecule.Atom;
import com.cairn.molecule.Atom.AtomFragmentType;
import com.cairn.molecule.AtomType;
import com.cairn.molecule.Fragments;

//TODO see if this can be replaced by SLN pattern matching

/**
 * Node in a Torsional distribution.
 * 
 * Torsional distribution routines
 * 
 * @author Gareth Jones
 * 
 */
public class TordistNode {
	private boolean isSybType;

	private AtomType sybType;

	public enum ElementalType {
		ET_NONE, ET_C, ET_P, ET_S, ET_N, ET_O
	};

	private int nHydrogens = -1, nNeighbours = 0, nStart, nEnd;

	private AtomFragmentType specialFrag = AtomFragmentType.ATM_NONE;

	private ElementalType elType = ElementalType.ET_NONE;

	public enum ParentLinkage {
		PL_ANY, PL_SINGLE, PL_DOUBLE, PL_AROM
	};

	private ParentLinkage parentLinkage = ParentLinkage.PL_ANY;

	private TordistNode neighbours[], parent = null;

	private Atom matchedAtom;

	private double weight;

	private TorsionalDistributions distributions;

	TordistNode(TorsionalDistributions d) {
		distributions = d;
		neighbours = new TordistNode[8];
	}

	/**
	 * ProcessNode reads in a node. The function is called recursively: each
	 * node is of the form "node [neighbours]", where neighbours is of the form
	 * "node1 node2 .."
	 * 
	 * @param entry
	 * @return
	 * 
	 */
	boolean processNode(String entry) {

		/* assign the atomic definition of the node */
		if (!assignNodeAtomId(entry))
			return false;

		nNeighbours = 0;
		nEnd = -1;

		/* find and process neighbours */
		String nEntry = null;
		while ((nEntry = findNextNeighbour(entry)) != null) {
			TordistNode nNode = new TordistNode(distributions);
			nNode.parent = this;
			neighbours[nNeighbours] = nNode;
			/* each neighbour is a node */
			if (nNode.processNode(nEntry))
				nNeighbours++;
		}

		return true;
	}

	// findNextNeighbour

	/**
	 * Find the next neighbour in the list of neighbours. Remember that
	 * neighbours may contain nodes with neighbours and so on.
	 * 
	 * @param entry
	 * @return
	 * 
	 */
	String findNextNeighbour(String entry) {
		int lastEnd = nEnd;
		char chars[] = entry.toCharArray();
		int cnt = 0;

		if (lastEnd == -1) {
			int max = chars.length;
			for (; cnt < max; cnt++)
				if (chars[cnt] == '(')
					break;
			if (cnt == max)
				return null;
		} else
			cnt = lastEnd;

		cnt++;
		while (chars[cnt] == ' ' || chars[cnt] == '\t')
			cnt++;
		if (chars[cnt] == ')')
			return null;
		nStart = cnt;

		// in atom id
		while (chars[cnt] != ' ' && chars[cnt] != '\t' && chars[cnt] != ')') {
			cnt++;
			if (cnt == entry.length())
				throw new RuntimeException(
						"findNextNeighbour: end of string while parsing: " + entry);
		}

		// outside atom id
		if (chars[cnt] == ')') {
			nEnd = cnt - 1;
			String s = entry.substring(nStart, cnt);
			return s;
		}
		while (chars[cnt] == ' ' || chars[cnt] == '\t') {
			cnt++;
			if (cnt == entry.length())
				throw new RuntimeException(
						"findNextNeighbour: end of string while parsing: " + entry);
		}

		if (chars[cnt] == '[') {
			// in special fragment string
			while (chars[cnt] != ']')
				cnt++;
			cnt++;
			while (chars[cnt] == ' ' || chars[cnt] == '\t') {
				cnt++;
				if (cnt == entry.length())
					throw new RuntimeException(
							"findNextNeighbour: end of string while parsing: " + entry);
			}
		}

		if (chars[cnt] != '(') {
			nEnd = cnt - 1;
			String s = entry.substring(nStart, cnt);
			return s;
		}

		int depth = 1;
		while (depth > 0) {
			cnt++;
			if (cnt == entry.length())
				throw new RuntimeException(
						"findNextNeighbour: end of string while parsing: " + entry);
			if (chars[cnt] == ')')
				depth--;
			if (chars[cnt] == '(')
				depth++;
		}

		nEnd = cnt;
		String s = entry.substring(nStart, cnt + 1);
		return s;
	}

	/**
	 * Determines the atomic definition of the node
	 * 
	 * @param entry
	 * @return
	 * 
	 */
	boolean assignNodeAtomId(String entry) {

		String atom = getNextWord(entry);
		char chars[] = atom.toCharArray();
		int cnt = 0;

		// Look for bond linkage with parent
		switch (chars[cnt]) {
		case '~':
			parentLinkage = ParentLinkage.PL_AROM;
			cnt++;
			break;
		case '=':
			parentLinkage = ParentLinkage.PL_DOUBLE;
			cnt++;
			break;
		case '-':
			parentLinkage = ParentLinkage.PL_SINGLE;
			cnt++;
			break;
		}

		/* cehck to see if node just set hydrogen count in parent */
		if (Character.isDigit(chars[cnt]) && chars[cnt + 1] == 'H') {
			int nh = Character.digit(chars[cnt], 10);
			if (parent == null)
				throw new RuntimeException(
						"assign_node_atom_id: no parent node to assign "
								+ "hydrogen count to:\nreading " + entry);
			parent.nHydrogens = nh;
			// No new node as we just change parent
			return false;
		}

		atom = atom.substring(cnt);
		// Find elemement type
		if ((elType = getElementType(atom)) != ElementalType.ET_NONE) {
			isSybType = false;
		}
		// or sybyl atom type
		else {
			isSybType = true;
			sybType = AtomType.sybType(atom);
		}

		// see if a fragment match is required
		specialFrag = getSpecialFrag(entry);

		return true;
	}

	/**
	 * Determines the weight of a node
	 * 
	 * @return
	 */
	double getWeight() {
		// find weight of atomic definition of node
		weight = .0;
		if (isSybType)
			weight += distributions.getSybTypeWt();
		else
			weight += distributions.getElTypeWt();
		if (parentLinkage != ParentLinkage.PL_ANY)
			weight += distributions.getLinkageWt();
		if (specialFrag != AtomFragmentType.ATM_NONE)
			weight += distributions.getFragWt();
		if (nHydrogens > -1)
			weight += distributions.getHCntWt();

		double wt = weight;
		// Add in weight of neighbours recursively
		for (int i = 0; i < nNeighbours; i++)
			wt += neighbours[i].getWeight();

		return wt;
	}

	/**
	 * Find any special fragments that apply to the node.
	 * 
	 * @param entry
	 * @return
	 * 
	 */
	AtomFragmentType getSpecialFrag(String entry) {
		int cnt = 0;
		char chars[] = entry.toCharArray();
		while (chars[cnt] != '[') {
			cnt++;
			if (cnt == chars.length)
				return AtomFragmentType.ATM_NONE;
		}
		cnt++;

		int start = cnt;
		while (chars[cnt] != ']') {
			cnt++;
			if (cnt == chars.length)
				throw new RuntimeException("getSpecialFrag: can't parse " + entry);
		}
		int end = cnt;
		String special = entry.substring(start, end);

		AtomFragmentType frag = Fragments.fragNameToId(special);
		if (frag == AtomFragmentType.ATM_NONE)
			throw new RuntimeException("getSpecialFrag: don't recognise fragment "
					+ special);
		return frag;
	}

	/**
	 * Assign the element type
	 * 
	 * @param atom
	 * @return
	 */
	ElementalType getElementType(String atom) {
		if (atom.length() != 1)
			return ElementalType.ET_NONE;
		switch (atom.charAt(0)) {
		case 'C':
			return ElementalType.ET_C;
		case 'P':
			return ElementalType.ET_P;
		case 'S':
			return ElementalType.ET_S;
		case 'N':
			return ElementalType.ET_N;
		case 'O':
			return ElementalType.ET_O;
		default:
			return ElementalType.ET_NONE;
		}
	}

	/**
	 * Copy the next word in line
	 * 
	 * @param line
	 * @return
	 * 
	 */
	String getNextWord(String line) {
		char chars[] = line.toCharArray();
		int cnt = 0;

		while (chars[cnt] == ' ' || chars[cnt] == '\t')
			cnt++;

		int start = cnt;
		if (chars[cnt] == ']' || chars[cnt] == ']' || chars[cnt] == '('
				|| chars[cnt] == ')' || chars[cnt] == '\n')
			throw new RuntimeException("getNextWord: can't get next word from line "
					+ line);

		while (chars[cnt] != ' ' && chars[cnt] != '\t' && chars[cnt] != '\n'
				&& chars[cnt] != '|' && chars[cnt] != ']' && chars[cnt] != ']'
				&& chars[cnt] != '(' && chars[cnt] != ')') {
			cnt++;
			if (cnt == chars.length)
				break;
		}
		String rtn = line.substring(start, cnt);
		return rtn;
	}

	/**
	 * Return description of a node
	 * 
	 * @return
	 * 
	 */
	String info() {
		String rtn = " ";
		// atomic definition
		switch (parentLinkage) {
		case PL_ANY:
			break;
		case PL_SINGLE:
			rtn += "-";
			break;
		case PL_DOUBLE:
			rtn += "=";
			break;
		case PL_AROM:
			rtn += "~";
			break;
		default:
			throw new RuntimeException("print: unknown parent linkage\n");
		}
		if (isSybType)
			rtn += sybType.getName();
		else
			rtn += elementTypeInfo(elType);
		rtn += " ";
		if (specialFrag != AtomFragmentType.ATM_NONE) {
			rtn += "[" + specialFragInfo() + "] ";
		}

		// print out neighbours
		if (nHydrogens >= 0 || nNeighbours > 0) {
			rtn += "(";
			if (nHydrogens >= 0)
				rtn += nHydrogens + "H";
			for (int i = 0; i < nNeighbours; i++)
				rtn += neighbours[i].info();
			rtn += ") ";
		}

		return rtn;
	}

	/**
	 * Describe an element
	 * 
	 * @param type
	 * @return
	 * 
	 */
	String elementTypeInfo(ElementalType type) {
		switch (type) {
		case ET_NONE:
			return "none";
		case ET_C:
			return "C";
		case ET_P:
			return "P";
		case ET_S:
			return "S";
		case ET_N:
			return "N";
		case ET_O:
			return "O";
		default:
			throw new RuntimeException("printElementType: unknown element type");
		}
	}

	/**
	 * describe a fragment
	 * 
	 * @return @
	 */
	String specialFragInfo() {
		String name = Fragments.idToFragName(specialFrag);
		if (name != null)
			return name;
		else
			throw new RuntimeException("printSpecialFrag: unknown fragment type");
	}

	/**
	 * Return true if node matches atom_no.
	 * 
	 * @param parent
	 * @param atom
	 * @param parentAtom
	 * @param td
	 * @return
	 * @throws GaException
	 */
	boolean matchNode(TordistNode parent, Atom atom, Atom parentAtom,
			TorsionalDistribution td) {

		// Atom is already matched to another node
		if (td.isAtomNodesUsed(atom.getNo()) && parent != null)
			return false;

		// Check to see if atom types match
		if (!matchAtomType(atom.getType()))
			return false;

		// get the environment for hydrogens and number of neighbours
		if (nHydrogens != -1) {
			if (nHydrogens != atom.getHydrogenCount())
				return false;
		}
		if (atom.getnNeighbours() < nNeighbours)
			return false;

		// Check linkage to parent
		if (parentLinkage != ParentLinkage.PL_ANY && parent != null) {
			switch (parentLinkage) {
			case PL_SINGLE:
				if (!atom.isSingleNeighbour(parentAtom))
					return false;
				break;
			case PL_DOUBLE:
				if (!atom.isDoubleNeighbour(parentAtom))
					return false;
				break;
			case PL_AROM:
				if (!atom.isAromaticNeighbour(parentAtom))
					return false;
				break;
			default:
				throw new IllegalArgumentException("Unknown linkage " + parentLinkage);
			}
		}

		// Check fragments
		if (specialFrag != AtomFragmentType.ATM_NONE) {
			if (!atom.inFragment(specialFrag))
				return false;
		}

		// node matches
		td.setAtomNodesUsed(atom.getNo(), true);

		// now check neighbours
		if (matchNeighbours(atom, 0, td)) {
			// all neighbours match, mark this node as matched to atom
			matchedAtom = atom;
			return true;
		} else {
			// failed to match node, if this isn't a toplevel node
			// backtrack.
			if (parent != null)
				td.setAtomNodesUsed(atom.getNo(), false);
			return false;
		}
	}

	/**
	 * Returns true if the atoms bound to parent_no match the neighbours in
	 * parent node. This is done recursively, with each neighbor being matched
	 * in turn, the current neighbout is neighbour_no. *
	 * 
	 * @param parent
	 * @param neighbourNo
	 * @param td
	 * @return
	 * @throws GaException
	 */
	boolean matchNeighbours(Atom parent, int neighbourNo, TorsionalDistribution td) {

		if (neighbourNo == nNeighbours)
			return true; // all neighbours have been matched

		// loop through all the atoms bonded to parent atom
		for (Atom parentNeighbour : parent.getNotDummyNeighbours()) {
			// match one neighbour node
			if (!neighbours[neighbourNo].matchNode(this, parentNeighbour, parent, td))
				continue; // failed: try the next atom

			// matched this neighbour: try to match the others
			if (matchNeighbours(parent, neighbourNo + 1, td))
				return true; // All neighbours matched
			// Failed, backtrack
			neighbours[neighbourNo].markAtomNodesUnused(td);
		}

		return false; // Failed
	}

	/**
	 * Matches node against atom type.
	 * 
	 * @param atype
	 * @return
	 * @throws GaException
	 */
	boolean matchAtomType(AtomType atype) {
		// match sybyl types
		if (isSybType)
			return sybType.matchType(atype);

		// match element types
		switch (elType) {
		case ET_C:
			return atype.isCarbonType();
		case ET_P:
			return atype.isPhosphorousType();
		case ET_S:
			return atype.isSulphurType();
		case ET_N:
			return atype.isNitrogenType();
		case ET_O:
			return atype.isOxygenType();
		default:
			throw new RuntimeException("match_atom_type: unknown elemental type");
		}
	}

	/**
	 * A match to a node has failed. Make sure all the atom involved are no
	 * longer marked as used.
	 * 
	 * @param td
	 */
	void markAtomNodesUnused(TorsionalDistribution td) {

		// The two atoms defining the rotatable bond can not be unmapped
		if (matchedAtom.getNo() != td.getNo2Atom().getNo()
				&& matchedAtom.getNo() != td.getNo3Atom().getNo())
			td.setAtomNodesUsed(matchedAtom.getNo(), false);
		for (int i = 0; i < nNeighbours; i++)
			neighbours[i].markAtomNodesUnused(td);
	}

	public Atom getMatchedAtom() {
		return matchedAtom;
	}

	public TordistNode[] getNeighbours() {
		return neighbours;
	}

	public int getNNeighbours() {
		return nNeighbours;
	}

	protected void setNeighbours(TordistNode[] neighbours) {
		this.neighbours = neighbours;
	}

	protected void setNeighbours(int i, TordistNode neighbour) {
		this.neighbours[i] = neighbour;
	}

	protected void setNNeighbours(int neighbours) {
		nNeighbours = neighbours;
	}

	protected void addNeighbour(TordistNode neighbour) {
		neighbours[nNeighbours++] = neighbour;
	}

	protected void removeLastNeighbour() {
		nNeighbours--;
	}

	protected void setParentLinkage(ParentLinkage parentLinkage) {
		this.parentLinkage = parentLinkage;
	}

	public ParentLinkage getParentLinkage() {
		return parentLinkage;
	}

}
