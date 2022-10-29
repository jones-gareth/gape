package com.cairn.molecule.test;

import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.commons.lang3.mutable.MutableBoolean;
import org.apache.commons.math3.util.Precision;
import org.apache.log4j.Logger;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import com.cairn.molecule.Atom;
import com.cairn.molecule.Charge;
import com.cairn.molecule.Molecule;

/**
 * Unit test to check charges added to a compound
 * 
 * @author Gareth Jones
 *
 */
public class TestCharge {
	private static final Logger logger = Logger.getLogger(TestCharge.class);
	static {
		// logger.setLevel(Level.DEBUG);
	}

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
	}

	private class ChargeAtom {

		private String atomName;
		private int formalCharge;
		private double partialCharge;

		public ChargeAtom(String atomName, int formalCharge, double partialCharge) {
			super();
			this.atomName = atomName;
			this.formalCharge = formalCharge;
			this.partialCharge = partialCharge;
		}

	}

	private boolean testAtom(Atom atom, int formalCharge, double partialCharge) {
		logger.debug("Testing atom " + atom.getLabel() + " for formal charge "
				+ formalCharge + " partial charge " + partialCharge);
		boolean pass = true;
		int atomFormalCharge = atom.getFormalCharge() != null ? atom.getFormalCharge()
				: 0;
		if (formalCharge != atomFormalCharge) {
			logger.error("Formal charge mismatch for atom " + atom.getLabel()
					+ " expected " + formalCharge + " got " + atomFormalCharge);
			pass = false;
		}
		if (!sameCharge(partialCharge, atom.getPartialCharge())) {
			logger.error("Partial charge mismatch for atom " + atom.getLabel()
					+ " expected " + partialCharge + " got " + atom.getPartialCharge());
			pass = false;
		}

		return pass;
	}

	private boolean test(String mol2File, ChargeAtom[] chargeAtoms) throws Exception {
		Molecule molecule = TestUtil.getMoleculeFromSybylResource(mol2File);

		Map<String, ChargeAtom> chargeAtomsMap = Arrays.asList(chargeAtoms).stream()
				.collect(Collectors.toMap(x -> x.atomName, x -> x));
		Charge charge = new Charge(molecule);
		charge.charge();

		MutableBoolean pass = new MutableBoolean(true);
		molecule.getAtoms().forEach(atom -> {

			ChargeAtom chargeAtom = chargeAtomsMap.get(atom.getLabel());
			if (chargeAtom == null) {
				if (!testAtom(atom, 0, .0)) {
					pass.setFalse();
				}
			} else {
				if (!testAtom(atom, chargeAtom.formalCharge, chargeAtom.partialCharge)) {
					pass.setFalse();
				}
				chargeAtomsMap.remove(atom.getLabel());
			}

		});

		assert chargeAtomsMap.size() == 0 : "Some atom definitions are untested";
		return pass.booleanValue();
	}

	private boolean sameCharge(double charge1, double charge2) {
		return Precision.equals(charge1, charge2, 0.0001);
	}

	@Test
	public void testOmega1() throws Exception {
		assertTrue(test("/com/cairn/molecule/test/data/omega_1.mol2",
				new ChargeAtom[] {}));
	}

	@Test
	public void testOmega12() throws Exception {
		assertTrue(test("/com/cairn/molecule/test/data/omega_12.mol2",
				new ChargeAtom[] { new ChargeAtom("N5", 1, 1.0) }));
	}

	@Test
	public void testPyrrole() throws Exception {
		assertTrue(test("/com/cairn/molecule/test/data/pyrrole.mol2",
				new ChargeAtom[] {}));
	}

}
