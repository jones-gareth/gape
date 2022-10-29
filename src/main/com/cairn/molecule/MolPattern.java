package com.cairn.molecule;

import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

/**
 * MolPattern- a SLN-compliant pattern for specifying substructure queries
 * 
 * @author Gareth Jones
 * 
 */
public class MolPattern {
	private static final Logger logger = Logger.getLogger(MolPattern.class);

	private String pattern;

	static final int MAX_ATOMS = 100;

	// The pattern is comprises of a number of tokens- each token is an atom or
	// extra bond (ring closure or branch) in
	// an sln string.
	private PatternToken tokens[];

	private int nTokens, position, branchLevel;

	private char charArray[];

	// Molecule class generated from SLN/pattern
	private Molecule patternMol;

	/**
	 * Simple constructor
	 */
	MolPattern() {
		;
	}

	/**
	 * Constructor from SLN string
	 * 
	 * @param p
	 *            SLN
	 */
	private MolPattern(String p) {
		pattern = p;
	}

	/**
	 * Constructs and parses a pattern from and SLN string.
	 * 
	 * @param pattern
	 *            SLN
	 * @return
	 */
	public static final MolPattern generateMolPattern(String pattern) {
		MolPattern p = new MolPattern(pattern);
		p.parsePattern();
		return p;
	}

	/**
	 * Test routine- converts an SLN to a MOL2 file.
	 * 
	 * @param args
	 */
	public static void main(String args[]) {
		if (args.length != 2) {
			System.err.println("Usage: MolPattern <pattern> <mol2file>");
			System.exit(0);
		}

		MolPattern p = generateMolPattern(args[0]);
		p.patternToMolecule(true);
		p.patternMol.build();
		p.patternMol.writeSybylMol2File(args[1], "Generated from SLN");

	}

	/**
	 * Parse an SLN into an array of Tokens.
	 * 
	 */
	private void parsePattern() {
		tokens = new PatternToken[MAX_ATOMS];
		charArray = pattern.toCharArray();
		position = 0;

		int no = 0;
		while (position < pattern.length()) {
			tokens[no] = new PatternToken();
			logger.debug("Getting token " + no + " position " + position + " at "
					+ pattern.substring(position));
			tokens[no].no = no;
			tokens[no].getToken();
			no++;
		}
		nTokens = no;

	}

	/**
	 * Converts a parsed (tokenized) pattern into a Molecule
	 * 
	 * @see #patternToMolecule(boolean)
	 */
	public void patternToMolecule() {
		patternToMolecule(false);
	}

	/**
	 * Converts a parsed (tokenized) pattern into a Molecule
	 * 
	 * @param explicitHydrogens
	 *            set this to include connected hydrogens in the constructed
	 *            Molecule. Otherwise it will be a heavy atom model.
	 */
	private void patternToMolecule(boolean explicitHydrogens) {

		patternMol = new Molecule();
		patternMol.setName("PatternMolecule");

		for (int i = 0; i < nTokens; i++) {
			tokens[i].toPatternMol();
		}

		if (explicitHydrogens) {
			int no = patternMol.getnAtoms();
			for (int i = 0; i < no; i++) {
				QueryAtom atom = (QueryAtom) patternMol.getAtom(no);
				if (atom.getnSlnHydrogens() != null) {
					for (int j = 0; j < atom.getnSlnHydrogens(); j++) {
						Atom hydAtom = new Atom(patternMol, patternMol.getnAtoms(),
								AtomType.Type.H);
						patternMol.addAtom(hydAtom, new double[] { 0, 0, 0, 1.0 });
						Bond bond = new Bond(patternMol.getnBonds(), hydAtom, atom,
								BondType.Type.SINGLE);
						patternMol.addBond(bond);
					}
				}
			}
		}

		logger.debug("nAtoms " + patternMol.getnAtoms());
		logger.debug("nBonds " + patternMol.getnBonds());
	}

	/**
	 * Counts the bonds in the pattern
	 * 
	 * @return
	 */
	@SuppressWarnings("unused")
	private int countBonds() {
		// Each token has a bond to the next token (except for the last token!)
		// - there are also tokens for ring closures and branches
		return nTokens - 1;
	}

	/**
	 * Counts the total number of hydrogens in the parsed pattern.
	 * 
	 * @return
	 */
	@SuppressWarnings("unused")
	private int countHydrogens() {
		int nH = 0;

		for (int i = 0; i < nTokens; i++)
			if (tokens[i].hasHydrogens)
				nH += tokens[i].nHydrogens;
		return nH;
	}

	/**
	 * Counts the total number of atoms in the parsed pattern.
	 * 
	 * @return
	 */
	@SuppressWarnings("unused")
	private int countAtoms() {
		int nAtoms = 0;

		for (int i = 0; i < nTokens; i++)
			if (tokens[i].type != null)
				nAtoms++;
		return nAtoms;
	}

	/**
	 * @return the patternMol
	 */
	public Molecule getPatternMol() {
		return patternMol;
	}

	/**
	 * Class to represent an element of an SLN string. Is an atom with bond,
	 * ring closure bond or branch bond.
	 * 
	 */
	private class PatternToken {

		private boolean singleBond, doubleBond, aromaticBond, tripleBond, openBranch,
				closeBranch, hasRingClosure, hasLabel, hasHydrogens, fullyConnected,
				noHydrogens, hasIsList, hasNotList;
		private Boolean inRing;

		private AtomType type;

		private QueryAtom atom;

		private PatternToken matchingToken, previousNeighbour, ringClosureToken;

		private int level, ringClosure, label, nHydrogens, nBonds, no, nConnected;
		private Integer formalCharge;

		private double partialCharge;

		private List<AtomType> isList, notList;

		static final char SINGLE = '-', DOUBLE = '=', TRIPLE = '#', AROMATIC = ':';

		/**
		 * Gets the next token from the current position on the SLN
		 * 
		 */
		private void getToken() {
			getOpenBranch();
			getBond();
			getAtomType();
			if (type == null) {
				getRingClosure();
				if (!hasRingClosure)
					throw new RuntimeException(
							"MolPattern: getToken: failed to get atom type or ring closure at "
									+ pattern.substring(position));
			}
			getCloseBranch();
			if (type != null && no > 0) {
				previousNeighbour = getPrevious();
				logger.debug("Previous " + previousNeighbour.no);
			}
		}

		/**
		 * Find a branch opening
		 */
		private void getOpenBranch() {
			if (position >= charArray.length)
				return;

			if (charArray[position] == '(') {
				branchLevel++;
				level = branchLevel;
				openBranch = true;
				logger.debug("Opening Branch level " + level);
				position++;
			}
		}

		/**
		 * Find a branch closing
		 */
		private void getCloseBranch() {
			if (position >= charArray.length)
				return;

			if (charArray[position] != ')')
				return;

			level = branchLevel;
			logger.debug("Closing Branch level " + level);

			// find matching token
			for (int i = no; i >= 0; i--) {
				if (tokens[i].openBranch && tokens[i].level == level) {
					matchingToken = tokens[i];
					tokens[i].matchingToken = this;
					break;
				}
			}

			logger.debug("Matching token " + matchingToken.no);
			branchLevel--;
			closeBranch = true;
			position++;
		}

		/**
		 * Gets the previous token. Skips across branches.
		 * 
		 * @return
		 */
		private PatternToken getPrevious() {
			PatternToken previous = tokens[no - 1];
			if (previous.closeBranch)
				return previous.matchingToken.getPrevious();
			return previous;
		}

		/**
		 * Gets the bond type
		 */
		private void getBond() {
			boolean hasBond = false;
			boolean inBond = true;
			nBonds = 0;

			while (inBond) {
				inBond = false;

				if (charArray[position] == SINGLE) {
					position++;
					nBonds++;
					hasBond = true;
					singleBond = true;
					logger.debug("Single Bond");
					inBond = true;
				}
				if (charArray[position] == DOUBLE) {
					position++;
					nBonds++;
					hasBond = true;
					doubleBond = true;
					logger.debug("Double Bond");
					inBond = true;
				}
				if (charArray[position] == TRIPLE) {
					position++;
					nBonds++;
					hasBond = true;
					tripleBond = true;
					logger.debug("Triple Bond");
					inBond = true;
				}
				if (charArray[position] == AROMATIC) {
					position++;
					nBonds++;
					hasBond = true;
					aromaticBond = true;
					logger.debug("Aromatic Bond");
					inBond = true;
				}
			}

			if (!hasBond) {
				singleBond = true;
				nBonds++;
				logger.debug("Implied Single Bond");
			}

		}

		/**
		 * Determines atom type.
		 * 
		 */
		private void getAtomType() {
			AtomType types[] = AtomType.types;
			type = null;
			for (int i = 0; i < types.length; i++) {
				if (pattern.startsWith(types[i].getName(), position)) {
					if (type == null
							|| type.getName().length() < types[i].getName().length())
						type = types[i];
				}
			}

			if (type != null) {
				logger.debug("Matching atom " + type.getName());
				position += type.getName().length();
				getAtomLabel();
				getHydrogenCount();
			}
		}

		/**
		 * Gets atom label and attributes (ring labels, hydrogen count, fully
		 * connected, is lists etc).
		 * 
		 */
		private void getAtomLabel() {
			if (position >= charArray.length)
				return;

			if (charArray[position] != '[') {
				return;
			}
			position++;
			int start = position;

			// here we get the atom label - so we're looking for
			// [<label>] or [<label>:<rest>] where label is an integer

			int end = pattern.indexOf(']', position);
			int endLabel = pattern.indexOf(':', position);
			if (endLabel > end)
				endLabel = -1;
			if (endLabel == -1
					&& Character.getType(charArray[position]) == Character.DECIMAL_DIGIT_NUMBER)
				endLabel = end;

			if (end == -1)
				throw new RuntimeException("MolPattern: getAtomLabel: parse error at "
						+ pattern.substring(position));

			if (endLabel != -1) {
				String labelStr = pattern.substring(position, endLabel);
				try {
					label = Integer.valueOf(labelStr).intValue();
				} catch (NumberFormatException ex) {
					throw new RuntimeException("MolPattern: getAtomLabel: label "
							+ labelStr + " is not an integer at "
							+ pattern.substring(position));
				}
				hasLabel = true;
				logger.debug("Got atom label " + label);
			}

			if (endLabel != -1)
				start = endLabel + 1;
			position = end + 1;

			logger.debug("Start " + start + " end " + end);
			if (start >= end)
				return;

			// Here we're looking for tokens like
			// [<label>:<token1>;<token2>..] or [<token1>;<token2>..]

			String rest = pattern.substring(start, end);
			logger.debug("Rest is " + rest);

			int p = 0;
			while (p != -1) {

				int n = rest.indexOf(';', p + 1);
				String token = n == -1 ? rest.substring(p) : rest.substring(p, n);

				p = n;
				if (n != -1)
					p++;

				logger.debug("Got token " + token);
				if (token.equals("f")) {
					fullyConnected = true;
					logger.debug("fully connected");
				} else if (token.equals("r")) {
					inRing = true;
					logger.debug("in ring");
				} else if (token.equals("!r")) {
					inRing = false;
					logger.debug("not in ring");
				} else if (token.equals("!h")) {
					noHydrogens = true;
					logger.debug("no hydrogens");
				} else if (token.startsWith("is=")) {
					logger.debug("is List");
					hasIsList = true;
					String list = token.substring(3);
					isList = getAtomList(list);
				} else if (token.startsWith("not=")) {
					logger.debug("not List");
					hasNotList = true;
					String list = token.substring(4);
					notList = getAtomList(list);
				} else if (token.startsWith("nconn=")) {
					String nStr = token.substring(6);
					int no = Integer.valueOf(nStr).intValue();
					logger.debug("nconn = " + no);
					nConnected = no;
				} else if (token.startsWith("charge=")) {
					String nStr = token.substring(7);
					int no = Integer.valueOf(nStr).intValue();
					logger.debug("formal charge = " + no);
					formalCharge = no;
				} else if (token.startsWith("fcharge=")) {
					String nStr = token.substring(8);
					partialCharge = Double.parseDouble(nStr);
					logger.debug("partial charge = " + partialCharge);
				}

				else {
					logger.warn("Unknown token " + token + " in mol pattern");
				}
			}

		}

		/**
		 * Determines hydrogen count.
		 * 
		 */
		private void getHydrogenCount() {
			if (position >= charArray.length)
				return;

			if (charArray[position] != 'H')
				return;
			hasHydrogens = true;
			position++;
			if (position >= charArray.length) {
				nHydrogens = 1;
			} else {
				if (Character.getType(charArray[position]) == Character.DECIMAL_DIGIT_NUMBER) {
					nHydrogens = charArray[position] - '0';
					position++;
				} else
					nHydrogens = 1;
			}
			logger.debug("Got " + nHydrogens + " Hydrogens");
		}

		/**
		 * Finds a ring closure and the matching ring opening.
		 * 
		 * @throws RuntimeException
		 */
		private void getRingClosure() {
			if (charArray[position] != '@')
				return;
			position++;
			int start = position;
			while (position < charArray.length
					&& Character.getType(charArray[position]) == Character.DECIMAL_DIGIT_NUMBER) {
				position++;
			}
			if (position == start)
				throw new RuntimeException(
						"MolPattern: ringClosure: failed to get ring close label at "
								+ pattern.substring(position));
			String labelStr = pattern.substring(start, position);
			ringClosure = Integer.valueOf(labelStr).intValue();
			logger.debug("Got ring closure label " + ringClosure);
			hasRingClosure = true;

			// find matching ring opening
			for (int i = no - 1; i >= 0; i--) {
				if (tokens[i].hasLabel && tokens[i].label == ringClosure) {
					ringClosureToken = tokens[i];
					tokens[i].ringClosureToken = this;
					break;
				}
			}

			logger.debug("Matching ring closure token " + ringClosureToken.no);

		}

		/**
		 * Converts the information in this token to the equivalent in a
		 * Molecule class.
		 */
		private void toPatternMol() {

			logger.debug("Converting token " + no + " to molecule");

			if (type != null) {
				atom = new QueryAtom(patternMol, patternMol.getnAtoms(), type.getType());
				if (hasHydrogens) {
					atom.setnSlnHydrogens(nHydrogens);
				}

				if (fullyConnected) {
					atom.setFullyConnected(true);
				}
				if (noHydrogens) {
					atom.setNoHydrogens(true);
				}
				if (hasIsList) {
					atom.setIs(isList);
				}
				if (hasNotList) {
					atom.setNot(notList);
				}
				if (nConnected > 0) {
					atom.setnConnected(nConnected);
				}
				if (formalCharge != null) {
					atom.setFormalCharge(formalCharge);
				}
				if (inRing != null) {
					atom.setQueryInRing(inRing);
				}
				if (partialCharge != .0) {
					atom.setPartialCharge(partialCharge);
				}

				patternMol.addAtom(atom, new double[] { 0, 0, 0, 1.0 });
			} else {
				assert hasRingClosure;
				atom = ringClosureToken.atom;
			}

			if (no == 0)
				return;

			QueryBond bond = null;
			Atom atom2 = null;

			if (hasRingClosure) {
				atom = tokens[no - 1].atom;
				atom2 = ringClosureToken.atom;
			} else {
				atom2 = previousNeighbour.atom;
			}

			int bondNo = patternMol.getnBonds();
			if (singleBond)
				bond = new QueryBond(bondNo, atom, atom2, BondType.Type.SINGLE);
			else if (doubleBond)
				bond = new QueryBond(bondNo, atom, atom2, BondType.Type.DOUBLE);
			else if (tripleBond)
				bond = new QueryBond(bondNo, atom, atom2, BondType.Type.TRIPLE);
			else if (aromaticBond)
				bond = new QueryBond(bondNo, atom, atom2, BondType.Type.AR);

			if (nBonds > 1) {
				bond.setHasMultipleBondTypes(true);
				if (singleBond)
					bond.setSingle(true);
				if (doubleBond)
					bond.setDouble(true);
				if (tripleBond)
					bond.setTriple(true);
				if (aromaticBond)
					bond.setAromatic(true);
			}

			patternMol.addBond(bond);
			bondNo++;

		}

		/**
		 * Determines a list of atom types for an is= or not= atom type
		 * specification
		 * 
		 * @param list
		 * @return
		 */
		private List<AtomType> getAtomList(String list) {
			String strTypes[] = list.split(",");
			List<AtomType> types = new ArrayList<>();
			for (int i = 0; i < strTypes.length; i++) {
				AtomType type = AtomType.sybType(strTypes[i]);
				logger.debug("List atom type " + type.getName());
				types.add(type);
			}
			return types;
		}
	}

}
