package com.cairn.molecule;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import org.apache.log4j.Logger;

/**
 * A class to store all out of plane atoms in the Tripos forcefield
 * 
 * @author Gareth Jones
 *
 */
class TaffOops {
	private final List<TaffOop> parameters = new ArrayList<>();

	static final Logger logger = Logger.getLogger(TaffOops.class);

	TaffOops() {
		loadParameterFile();
	}

	void addParameters(Molecule mol) {
		for (Atom atom : mol.getAtoms()) {
			atom.setTaffOop(null);

			int no = 0;
			for (TaffOop toop : parameters) {
				no++;
				if (Taff.matchType(atom.getAtomType(), toop.getAtom())) {
					logger.debug("matched atom " + atom.getNo() + " to taff entry " + no);
					atom.setTaffOop(toop);
				}
			}
		}
	}

	private void loadParameterFile() {
		parameters.clear();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(
					getClass().getResourceAsStream("TAFF_OOP_BEND.txt")));
			String line;
			while ((line = in.readLine()) != null) {
				logger.debug("Line: " + line);
				StringTokenizer st = new StringTokenizer(line);
				AtomType.Type atom = AtomType.sybType(st.nextToken()).getType();
				double k = (new Double(st.nextToken())).doubleValue();
				TaffOop toop = new TaffOop(atom, k);
				parameters.add(toop);
			}
		} catch (IOException ex) {
			System.err.println("Failed to open TAFF_OOP_BEND");
		}
	}
}