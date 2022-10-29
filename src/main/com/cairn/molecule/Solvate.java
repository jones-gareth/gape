package com.cairn.molecule;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import org.apache.log4j.Logger;

import com.cairn.gape.utils.InfoMessageLogger;

/**
 * Class for solvating a molecule. Sln patterns are stored in a text file. This
 * class matches pattern against atoms and protonates bases and removes
 * hydrogens from acids as appropriate.
 * 
 * @author Gareth Jones
 * 
 */
public class Solvate {
	private volatile Molecule molecule;

	private boolean acidMatches[], baseMatches[];

	private static final Logger logger = Logger.getLogger(Solvate.class);
	static {
		// logger.setLevel(Level.DEBUG);
	}

	/**
	 * Empty constructor.
	 */
	public Solvate() {
		;
	}

	/**
	 * Constructor- set's molecule.
	 * 
	 * @param _molecule
	 */
	public Solvate(Molecule _molecule) {
		setMolecule(_molecule);
	}

	/**
	 * Sets current molecule.
	 * 
	 * @param _molecule
	 */
	public void setMolecule(Molecule _molecule) {
		molecule = _molecule;
		acidMatches = new boolean[molecule.getnAtoms()];
		baseMatches = new boolean[molecule.getnBonds()];
	}

	/**
	 * Does the work- solvates current molecule.
	 * 
	 */
	public void solvate() {
		// We want to do solvation twice- the second pass can apply any rule
		// corrections

		for (int pass = 0; pass < 2; pass++) {
			logger.debug("Pass " + pass);
			SolvationGroup groups[] = SolvationGroup.solvationGroups.get();
			boolean baseMatched = false;
			boolean acidMatched = false;

			for (int i = 0; i < groups.length; i++) {
				boolean match = groups[i].matchMolecule(this);
				if (match) {
					if (groups[i].type == SolvationGroup.Type.BASE)
						baseMatched = true;
					else if (groups[i].type == SolvationGroup.Type.ACID)
						acidMatched = true;
				}
			}

			if (baseMatched) {
				// matched at least one base- add hydrogens
				AddHydrogens addHydrogens = new AddHydrogens(molecule);
				addHydrogens.addHydrogens();
			} else if (acidMatched)
				// matched at least one acid- update molecule after removal of
				// hydrogens.
				molecule.update();

			// no need to do pass 2
			if (!baseMatched && !acidMatched)
				break;
		}
	}

	/**
	 * Test routine- solvates a molecule file.
	 */
	public static void main(String[] args) {

		if (logger.isDebugEnabled()) {
			for (int i = 0; i < SolvationGroup.solvationGroups.get().length; i++) {
				SolvationGroup solvationGroup = SolvationGroup.solvationGroups.get()[i];
				logger.debug(solvationGroup.info());
			}
		}

		if (args.length == 0) {
			System.out.println("Usage: Solvate <mol file>");
			System.exit(0);
		}

		String file = args[0];

		try {

			List<Molecule> mols = Molecule.loadFiles(new String[] { file });

			for (Molecule mol : mols) {
				mol.assignAtomTypes();
				Solvate solvate = new Solvate(mol);
				solvate.solvate();
			}

			Molecule.write(mols, "solvate_" + file, "Solvated");

		} catch (Exception ex) {
			System.err.println("Exception " + ex);
		}
	}

	/**
	 * Class for representing a solvation group. Contains routines to find
	 * groups and add or remove hydrogens.
	 * 
	 * @author Gareth Jones
	 * 
	 */
	static class SolvationGroup {
		private String name, sln;

		private double pK;

		enum Type {
			ACID, BASE
		};

		private Type type;

		private MolPattern pattern;

		private AtomPatternMatch patternMatch;

		static final ThreadLocal<SolvationGroup[]> solvationGroups = ThreadLocal
				.withInitial(() -> loadSolvationGroups());

		private SolvationGroup() {
			;
		}

		/**
		 * Group constructor
		 * 
		 * @param _type
		 * @param _name
		 * @param _sln
		 * @param _pK
		 */
		private SolvationGroup(Type _type, String _name, String _sln, double _pK) {
			name = _name;
			sln = _sln;
			pK = _pK;
			type = _type;

			pattern = MolPattern.generateMolPattern(sln);
			patternMatch = new AtomPatternMatch(pattern);

		}

		/**
		 * Returns informational string
		 * 
		 * @return
		 */
		String info() {
			String rtn = name + " ";
			rtn += type == Type.ACID ? "Acid" : "Base";
			rtn += " " + sln;
			rtn += " [pK " + pK + "]";
			return rtn;
		}

		/**
		 * Searches the current molecule in solvate and finds any matching
		 * groups. Matching bases have formal charge set so that hydrogens can
		 * be added and matching acids have hydrogens removed.
		 * 
		 * @param solvate
		 * @return
		 */
		boolean matchMolecule(Solvate solvate) {
			Molecule molecule = solvate.molecule;

			InfoMessageLogger infoMessageLogger = molecule.getInfoMessageLogger();

			patternMatch.match(molecule);
			ArrayList<Atom> matches = patternMatch.getMatches();
			if (matches.size() == 0)
				return false;
			boolean match = false;

			if (type == Type.BASE) {
				for (Atom atom : matches) {
					// only match an atom once.
					if (solvate.baseMatches[atom.getNo()])
						continue;
					solvate.baseMatches[atom.getNo()] = true;
					atom.setFormalCharge(1);
					match = true;
					String message = "atom " + atom.info() + " matches " + info();
					infoMessageLogger.infoMessageln(3, message);
					logger.debug(message);
				}
			}

			else if (type == Type.ACID) {
				for (Atom atom : matches) {
					if (solvate.acidMatches[atom.getNo()])
						continue;
					solvate.acidMatches[atom.getNo()] = true;
					atom.setFormalCharge(-1);
					Optional<Atom> findH = atom.getNotDummyNeighbours().stream()
							.filter(a -> a.getAtomType() == AtomType.Type.H).findFirst();
					if (!findH.isPresent())
						throw new RuntimeException(
								"SolvationGroup:matchMolecule: can't find Acid hydrogen");
					Atom hydrogen = findH.get();
					molecule.deleteAtom(hydrogen, false);
					match = true;
					String message = "atom " + atom.info() + " matches " + info();
					infoMessageLogger.infoMessageln(3, message);
					logger.debug(message);
				}
			}

			return match;
		}

		/**
		 * Loads solvation group definitions from a text resource file
		 * (solvate.txt).
		 * 
		 * File format is group name on one line then type, sln and pK on the
		 * second.
		 * 
		 * @return
		 */
		private static SolvationGroup[] loadSolvationGroups() {
			ArrayList<SolvationGroup> solvationGroups = new ArrayList<SolvationGroup>();

			BufferedReader in = new BufferedReader(new InputStreamReader(
					(new SolvationGroup()).getClass().getResourceAsStream("solvate.txt")));

			try {
				String line = in.readLine();

				while (line != null) {
					line = line.trim();
					if (line.equals("") || line.startsWith("#")) {
						line = in.readLine();
						continue;
					}

					// first line is name
					String name = line;

					// second line is info
					line = in.readLine();
					if (line == null)
						throw new RuntimeException("empty definition for " + name);
					String vals[] = line.split("\\s+");

					Type type = null;
					if (vals[0].equalsIgnoreCase("base"))
						type = Type.BASE;
					else if (vals[0].equalsIgnoreCase("acid"))
						type = Type.ACID;
					else
						throw new RuntimeException("bad definition " + line);

					String sln = vals[1];
					double pK = Double.parseDouble(vals[2]);
					logger.debug("def type " + type + " name " + name + " sln " + sln
							+ " pK " + pK);
					solvationGroups.add(new SolvationGroup(type, name, sln, pK));

					line = in.readLine();
				}
			} catch (RuntimeException ex) {
				System.out.println(ex.toString());
				System.exit(0);
			} catch (IOException ex) {
				System.out.println(ex.toString());
				System.exit(0);

			}
			return solvationGroups.toArray(new SolvationGroup[solvationGroups.size()]);

		}

	}

}
