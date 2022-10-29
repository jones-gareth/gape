package com.cairn.molecule;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import org.apache.log4j.Logger;

import com.cairn.gape.utils.InfoMessageLogger;

class TaffTorsions {

	private static final Logger logger = Logger.getLogger(TaffTorsions.class);

	static {
		// logger.setLevel(Level.DEBUG);
	}

	private final List<TaffTorsion> parameters = new ArrayList<>();

	TaffTorsions() {
		loadParameterFile();
	}

	void addParameters(Molecule mol) {
		InfoMessageLogger infoMessageLogger = mol.getInfoMessageLogger();
		for (Torsion torsion : mol.getTorsions()) {
			torsion.setTaff(null);
			int no = 0;
			for (TaffTorsion tt : parameters) {
				no++;
				if (Taff.matchType(torsion.getAtom1().getAtomType(), tt.getAtom1())
						&& Taff.matchType(torsion.getAtom2().getAtomType(), tt.getAtom2())
						&& Taff.matchType(torsion.getAtom3().getAtomType(), tt.getAtom3())
						&& Taff.matchType(torsion.getAtom4().getAtomType(), tt.getAtom4())
						&& torsion.getBond().sybylBondType(mol) == tt.getBond()) {
					if (logger.isDebugEnabled()) {
						logger.debug("matched torsion to taff entry " + no);
						logger.debug(torsion.info());
						logger.debug(" --> " + tt.info());
					}
					if (torsion.getTaff() != null) {
						if (tt.getWeight() > torsion.getTaff().getWeight()) {
							torsion.setTaff(tt);
							torsion.setReverse(false);
						}
					} else {
						torsion.setTaff(tt);
						torsion.setReverse(false);
					}
				}
				if (Taff.matchType(torsion.getAtom1().getAtomType(), tt.getAtom4())
						&& Taff.matchType(torsion.getAtom2().getAtomType(), tt.getAtom3())
						&& Taff.matchType(torsion.getAtom3().getAtomType(), tt.getAtom2())
						&& Taff.matchType(torsion.getAtom4().getAtomType(), tt.getAtom1())
						&& torsion.getBond().sybylBondType(mol) == tt.getBond()) {
					if (logger.isDebugEnabled()) {
						logger.debug("matched reverse torsion to taff entry " + no);
						logger.debug(torsion.info() + " --> " + tt.info());
					}
					if (torsion.getTaff() != null) {
						if (tt.getWeight() > torsion.getTaff().getWeight()) {
							torsion.setTaff(tt);
							torsion.setReverse(true);
						}
					} else {
						torsion.setTaff(tt);
						torsion.setReverse(true);
					}
				}
			}
			if (torsion.getTaff() == null) {
				infoMessageLogger.infoMessageln(3,
						"No parameters for torsion " + torsion.info());
			} else if (logger.isDebugEnabled()) {
				logger.debug("final match reverse:" + torsion.isReverse()
						+ " torsion to taff entry " + no);
				logger.debug(torsion.info() + " --> " + torsion.getTaff().info());
			}
		}
	}

	private void loadParameterFile() {
		parameters.clear();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(
					getClass().getResourceAsStream("TAFF_TORS.txt")));
			String line;
			while ((line = in.readLine()) != null) {
				logger.debug("Line: " + line);
				StringTokenizer st = new StringTokenizer(line);
				AtomType.Type atom1 = AtomType.sybType(st.nextToken()).getType();
				AtomType.Type atom2 = AtomType.sybType(st.nextToken()).getType();
				AtomType.Type atom3 = AtomType.sybType(st.nextToken()).getType();
				AtomType.Type atom4 = AtomType.sybType(st.nextToken()).getType();
				BondType.Type bond = BondType.sybType(st.nextToken()).getType();
				double k = (new Double(st.nextToken())).doubleValue();
				double p = (new Double(st.nextToken())).doubleValue();
				TaffTorsion tt = new TaffTorsion(atom1, atom2, atom3, atom4, bond, k, p);
				logger.debug("Added entry " + tt.info());
				parameters.add(tt);
			}
		} catch (IOException ex) {
			logger.error("Failed to open TAFF_TORS");
		}
	}
}
