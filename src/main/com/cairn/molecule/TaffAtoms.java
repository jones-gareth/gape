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
 * A class to store all torsions in the Tripos forcefield
 * 
 * @author Gareth Jones
 *
 */
class TaffAtoms {
	private static Logger logger = Logger.getLogger(TaffAtoms.class);

	static {
		// logger.setLevel(Level.DEBUG);
	}

	private final List<TaffAtom> parameters = new ArrayList<>();

	TaffAtoms() {
		loadParameterFile();
	}

	void addParameters(Molecule mol) {
		InfoMessageLogger infoMessageLogger = mol.getInfoMessageLogger();

		for (Atom atom : mol.getAtoms()) {
			atom.setTaff(null);
			int no = 0;
			for (TaffAtom ta : parameters) {
				no++;
				if (Taff.matchType(atom.getAtomType(), ta.getAtom())) {
					logger.debug("matched atom " + atom.getNo() + " to taff entry " + no);
					atom.setTaff(ta);
				}
			}
			if (atom.getTaff() == null)
				infoMessageLogger.infoMessageln(3,
						"No parameters for atom " + atom.getNo());
		}
	}

	private void loadParameterFile() {
		parameters.clear();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(
					getClass().getResourceAsStream("TAFF_VDW.txt")));
			String line;
			while ((line = in.readLine()) != null) {
				StringTokenizer st = new StringTokenizer(line);
				AtomType.Type atom = AtomType.sybType(st.nextToken()).getType();
				double r = Double.parseDouble(st.nextToken());
				double k = Double.parseDouble(st.nextToken());
				TaffAtom ta = new TaffAtom(atom, r, k);
				parameters.add(ta);
			}
		} catch (IOException ex) {
			logger.error("Failed to open TAFF_VDW");
		}
	}
}