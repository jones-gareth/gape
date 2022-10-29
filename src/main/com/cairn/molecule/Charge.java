package com.cairn.molecule;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

/**
 * Class for charging a molecule. Sln patterns are stored in a text file. This
 * class matches patterns against atoms and sets formal and partial charges as
 * appropriate.
 * 
 * @author Gareth Jones
 * 
 */
public class Charge {

	private static Logger logger = Logger.getLogger(Charge.class);
	static {
		// logger.setLevel(Level.DEBUG);
	}
	private Molecule molecule;

	private boolean matches[];

	/**
	 * Empty constructor.
	 */
	public Charge() {
		;
	}

	/**
	 * Constructor- set's molecule.
	 * 
	 * @param _molecule
	 */
	public Charge(Molecule _molecule) {
		setMolecule(_molecule);
	}

	/**
	 * Sets current molecule.
	 * 
	 * @param _molecule
	 */
	public void setMolecule(Molecule _molecule) {
		molecule = _molecule;
		matches = new boolean[molecule.getnAtoms()];
	}

	/**
	 * Does the work- charges current molecule.
	 * 
	 */
	public void charge() {
		ChargeGroup groups[] = ChargeGroup.chargeGroups.get();

		for (int i = 0; i < groups.length; i++) {
			groups[i].matchMolecule(this);
		}
	}

	/**
	 * Test routine- charges a molecule file.
	 */
	public static void main(String[] args) {

		if (logger.isDebugEnabled()) {
			for (int i = 0; i < ChargeGroup.chargeGroups.get().length; i++) {
				ChargeGroup chargeGroup = ChargeGroup.chargeGroups.get()[i];
				logger.debug(chargeGroup.info());
			}
		}

		if (args.length == 0) {
			System.out.println("Usage: Charge <mol file>");
			System.exit(0);
		}

		String file = args[0];

		try {

			List<Molecule> mols = Molecule.loadFiles(new String[] { file });
			for (Molecule mol : mols) {
				mol.assignAtomTypes();
				Charge charge = new Charge(mol);
				charge.charge();
			}

			Molecule.write(mols, "charge_" + file, "Charged");

		} catch (Exception ex) {
			System.err.println("Exception " + ex);
		}
	}

	/**
	 * Class to represent a charged group
	 * 
	 * @author Gareth Jones
	 * 
	 */
	private static class ChargeGroup extends PatternMatch {
		String sln;

		int formalCharge;

		double partialCharge;

		boolean multipleFormalChargeMatch;

		Charge charge;

		static final ThreadLocal<ChargeGroup[]> chargeGroups = new ThreadLocal<ChargeGroup[]>() {
			@Override
			protected ChargeGroup[] initialValue() {
				return loadChargeGroups();
			}

		};

		/**
		 * Empty constructor
		 */
		private ChargeGroup() {
			;
		}

		/**
		 * Creates a charged group
		 * 
		 * @param _sln
		 * @param _formalCharge
		 * @param _partialCharge
		 */
		private ChargeGroup(String _sln, int _formalCharge, double _partialCharge,
				boolean _multipleFormalChargeMatch) {
			super(_sln);
			sln = _sln;
			formalCharge = _formalCharge;
			partialCharge = _partialCharge;
			multipleFormalChargeMatch = _multipleFormalChargeMatch;
		}

		/**
		 * Creates a descriptive string
		 * 
		 * @return
		 */
		String info() {
			String rtn = sln + " [formal charge " + formalCharge + " partial charge "
					+ partialCharge + "]";
			return rtn;
		}

		/**
		 * Searches for this charged group in the current molecule. Sets partial
		 * and formal charges for matches.
		 * 
		 * @param charge
		 * @return
		 */
		void matchMolecule(Charge _charge) {
			charge = _charge;
			Molecule molecule = charge.molecule;

			match(molecule);
		}

		/*
		 * Does the matching work, including atom-by atom matching to determine
		 * if formal charges are set multiple times. (non-Javadoc)
		 * 
		 * @see com.cairn.molecule.PatternMatch#process()
		 */
		@Override
		public void process() {
			int match[] = queryMatches[nMatches - 1];
			logger.debug("Match ");

			Atom atom = target.getAtom(match[0]);

			if (charge.matches[atom.getNo()])
				return;
			charge.matches[atom.getNo()] = true;

			if (formalCharge != 0) {
				atom.setFormalCharge(formalCharge);
			}
			atom.setPartialCharge(partialCharge);

			charge.molecule.getInfoMessageLogger().infoMessageln(3,
					"atom " + atom.info() + " matches " + info());

			if (multipleFormalChargeMatch || formalCharge == 0)
				return;

			logger.debug("checking formal charge single match");

			// Remove formal charge from other atoms in the matching pattern
			for (int i = 1; i < query.getnAtoms(); i++) {
				Atom tAtom = target.getAtom(match[i]);
				if (tAtom.getFormalCharge() != null) {
					logger.debug("Removing formal charge on atom " + tAtom.info());
					tAtom.setFormalCharge(null);
				}
			}
		}

		/**
		 * Loads charge groups from a text resource file (charge.txt)
		 * 
		 * Charge groups are defined by sln, formal and partial charges and a
		 * yes/no indicating if the formal charge can appear more than once in a
		 * pattern.
		 * 
		 * @return
		 */
		private static ChargeGroup[] loadChargeGroups() {
			ArrayList<ChargeGroup> chargeGroups = new ArrayList<ChargeGroup>();

			BufferedReader in = new BufferedReader(new InputStreamReader(
					(new ChargeGroup()).getClass().getResourceAsStream("charge.txt")));

			try {
				String line = in.readLine();

				while (line != null) {
					line = line.trim();
					if (line.equals("") || line.startsWith("#")) {
						line = in.readLine();
						continue;
					}

					logger.debug("Read: " + line + "\n");
					String vals[] = line.split("\\s+");

					String sln = vals[0];
					int formalCharge = Integer.parseInt(vals[1]);
					double partialCharge = Double.parseDouble(vals[2]);
					String str = vals[3];
					boolean multipleFormalChargeMatch = false;
					if (str.equalsIgnoreCase("yes"))
						multipleFormalChargeMatch = true;

					logger.debug("def: sln " + sln + " formalCharge " + formalCharge
							+ " partialCharge " + partialCharge + " single match "
							+ multipleFormalChargeMatch);
					chargeGroups.add(new ChargeGroup(sln, formalCharge, partialCharge,
							multipleFormalChargeMatch));

					line = in.readLine();
				}
			} catch (IOException ex) {
				System.out.println(ex.toString());
				System.exit(0);

			}
			return chargeGroups.toArray(new ChargeGroup[chargeGroups.size()]);

		}
	}
}