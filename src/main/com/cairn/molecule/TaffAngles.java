package com.cairn.molecule;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import com.cairn.gape.utils.InfoMessageLogger;

/**
 * A class to store all angles in the Tripos forcefield
 * 
 * @author Gareth Jones
 *
 */
class TaffAngles {

	private final List<TaffAngle> parameters = new ArrayList<TaffAngle>();;

	TaffAngles() {
		loadParameterFile();
	}

	void addParameters(Molecule mol) {

		for (Angle angle : mol.getAngles()) {

			angle.setTaff(null);

			InfoMessageLogger infoMessageLogger = mol.getInfoMessageLogger();

			for (TaffAngle ta : parameters) {

				assert angle != null;
				assert angle.getAtom1() != null;
				assert ta != null;

				if (Taff.matchType(angle.getAtom1().getAtomType(), ta.getAtom1())
						&& Taff.matchType(angle.getMid().getAtomType(), ta.getMid())
						&& Taff.matchType(angle.getAtom2().getAtomType(),
								ta.getAtom2())) {
					if (angle.getTaff() != null) {
						if (ta.getWeight() > angle.getTaff().getWeight())
							angle.setTaff(ta);
					} else {
						angle.setTaff(ta);
					}
				}
				if (Taff.matchType(angle.getAtom2().getAtomType(), ta.getAtom1())
						&& Taff.matchType(angle.getMid().getAtomType(), ta.getMid())
						&& Taff.matchType(angle.getAtom1().getAtomType(),
								ta.getAtom2())) {
					if (angle.getTaff() != null) {
						if (ta.getWeight() > angle.getTaff().getWeight())
							angle.setTaff(ta);
					} else {
						angle.setTaff(ta);
					}
				}
			}

			if (angle.getTaff() == null) {
				infoMessageLogger.infoMessage(3,
						"No parameters for angle " + angle.info());
				AtomType.Geometry geom = angle.getMid().getType().getGeometry();
				if (geom == AtomType.Geometry.TET) {
					angle.setTaff(new TaffAngle(angle.getAtom1().getAtomType(),
							angle.getMid().getAtomType(), angle.getAtom2().getAtomType(),
							109.5, 0.02));
					infoMessageLogger.infoMessageln(3, " [using a=109.5; k=0.02]");
				} else if (geom == AtomType.Geometry.LIN) {
					angle.setTaff(new TaffAngle(angle.getAtom1().getAtomType(),
							angle.getMid().getAtomType(), angle.getAtom2().getAtomType(),
							180.0, 0.04));
					infoMessageLogger.infoMessageln(3, " [using a=180; k=0.04]");
				} else if (geom == AtomType.Geometry.TRI) {
					angle.setTaff(new TaffAngle(angle.getAtom1().getAtomType(),
							angle.getMid().getAtomType(), angle.getAtom2().getAtomType(),
							120.0, 0.03));
					infoMessageLogger.infoMessageln(3, " [using a=120; k=0.03]");
				}
			}
		}
	}

	private void loadParameterFile() {
		parameters.clear();
		try {
			// URL file = getClass().getResource("TAFF_ANGLE_BEND.txt");
			// BufferedReader in =
			// new BufferedReader
			// (new InputStreamReader(file.openStream()));
			BufferedReader in = new BufferedReader(new InputStreamReader(
					getClass().getResourceAsStream("TAFF_ANGLE_BEND.txt")));
			String line;
			while ((line = in.readLine()) != null) {
				StringTokenizer st = new StringTokenizer(line);
				AtomType.Type atom1 = AtomType.sybType(st.nextToken()).getType();
				AtomType.Type atom2 = AtomType.sybType(st.nextToken()).getType();
				AtomType.Type atom3 = AtomType.sybType(st.nextToken()).getType();
				double a = (new Double(st.nextToken())).doubleValue();
				double k = (new Double(st.nextToken())).doubleValue();
				TaffAngle ta = new TaffAngle(atom1, atom2, atom3, a, k);
				parameters.add(ta);
			}
		} catch (IOException ex) {
			System.err.println("Failed to open TAFF_ANGLE_BEND");
		}
	}
}