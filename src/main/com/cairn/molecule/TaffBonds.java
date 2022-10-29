package com.cairn.molecule;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import org.apache.log4j.Logger;

import com.cairn.gape.utils.InfoMessageLogger;

/**
 * A class to store all bonds in the Tripos forcefield
 * 
 * @author Gareth Jones
 *
 */
class TaffBonds {
	private final List<TaffBond> parameters = new ArrayList<>();

	private static Logger logger = Logger.getLogger(TaffBonds.class);;

	static {
		// logger.setLevel(Level.DEBUG);
	}

	TaffBonds() {
		loadParameterFile();
	}

	void addParameters(Molecule mol) {
		InfoMessageLogger infoMessageLogger = mol.getInfoMessageLogger();

		for (Bond bond : mol.getBonds()) {
			bond.setTaff(null);
			int no = 0;
			for (TaffBond tb : parameters) {
				no++;
				if (bond.sybylBondType(mol) == tb.getBond()
						&& Taff.matchType(bond.getAtom1().getAtomType(), tb.getAtom1())
						&& Taff.matchType(bond.getAtom2().getAtomType(), tb.getAtom2())) {
					logger.debug("matched bond " + bond.getNo() + " to taff entry " + no);
					bond.setTaff(tb);
				}
				if (bond.sybylBondType(mol) == tb.getBond()
						&& Taff.matchType(bond.getAtom1().getAtomType(), tb.getAtom2())
						&& Taff.matchType(bond.getAtom2().getAtomType(), tb.getAtom1())) {
					logger.debug("matched bond " + bond.getNo() + " to taff entry " + no);
					bond.setTaff(tb);
				}
			}
			if (bond.getTaff() == null) {
				infoMessageLogger.infoMessageln(3, "No parameters for bond "
						+ bond.getNo() + " [using l=1.5; k=600]");
				bond.setTaff(new TaffBond(bond.getAtom1().getAtomType(),
						bond.getAtom2().getAtomType(), bond.getBondType(), 1.5, 600));
			}
		}
	}

	private void loadParameterFile() {
		parameters.clear();
		try {
			// URL file = getClass().getResource("TAFF_BOND_STRETCH.txt");
			// BufferedReader in =
			// new BufferedReader
			// (new InputStreamReader(file.openStream()));
			BufferedReader in = new BufferedReader(new InputStreamReader(
					getClass().getResourceAsStream("TAFF_BOND_STRETCH.txt")));
			String line;
			while ((line = in.readLine()) != null) {
				StringTokenizer st = new StringTokenizer(line);
				AtomType.Type atom1 = AtomType.sybType(st.nextToken()).getType();
				AtomType.Type atom2 = AtomType.sybType(st.nextToken()).getType();
				BondType.Type bond = BondType.sybType(st.nextToken()).getType();
				double len = (new Double(st.nextToken())).doubleValue();
				double k = (new Double(st.nextToken())).doubleValue();
				TaffBond tb = new TaffBond(atom1, atom2, bond, len, k);
				parameters.add(tb);
			}
		} catch (IOException ex) {
			System.err.println("Failed to open TAFF_BOND_STRETCH");
		}
	}
}