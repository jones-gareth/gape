package com.cairn.molecule.test;

import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.commons.lang3.mutable.MutableBoolean;
import org.apache.log4j.Logger;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import com.cairn.gape.feature.HydrogenBondingType;
import com.cairn.gape.molecule.GaMolecule;

/**
 * Check that the identified donors and acceptors are as expected for a number
 * of compounds.
 * 
 * @author Gareth Jones
 *
 */
public class TestHydrogenBondingType {
	private static final Logger logger = Logger.getLogger(TestHydrogenBondingType.class);
	static {
		// logger.setLevel(Level.DEBUG);
	}

	private class AtomInfo {

		private String atomLabel, hBondTypeName;

		public AtomInfo(String atomLabel, String hBondTypeName) {
			super();
			this.atomLabel = atomLabel;
			this.hBondTypeName = hBondTypeName;
		}

	}

	private boolean testMolecule(String resource, AtomInfo[] donorAtomInfoArray,
			AtomInfo[] acceptorAtomInfoArray) throws Exception {

		logger.info("Checking molecule " + resource);
		GaMolecule molecule = TestUtil.getGaMoleculeFromSybylResource(resource);

		HydrogenBondingType.searchMolecule(molecule);

		logger.info("Checking donors");
		boolean check1 = testMolecule(molecule, donorAtomInfoArray, true);
		logger.info("Checking acceptors");
		boolean check2 = testMolecule(molecule, acceptorAtomInfoArray, false);

		return check1 && check2;
	}

	private boolean testMolecule(GaMolecule molecule, AtomInfo[] atomInfoArray,
			boolean donor) throws Exception {
		Map<String, AtomInfo> atomInfoMap = Arrays.asList(atomInfoArray).stream()
				.collect(Collectors.toMap(x -> x.atomLabel, x -> x));
		assert atomInfoArray.length == atomInfoMap.size();

		MutableBoolean pass = new MutableBoolean(true);
		molecule.getAtoms()
				.stream()
				.forEach(
						atom -> {
							AtomInfo atomInfo = atomInfoMap.get(atom.getLabel());
							HydrogenBondingType hBondType = donor ? atom.getDonorType()
									: atom.getAcceptorType();
							if (hBondType != null) {
								logger.debug("Atom " + atom.getLabel()
										+ " has hydrogen bonding type "
										+ hBondType.getName());
							}
							if (atomInfo != null) {
								atomInfoMap.remove(atom.getLabel());
								if (hBondType == null) {
									logger.warn("Atom " + atom.getLabel()
											+ " has no hydrogen bond type expected "
											+ atomInfo.hBondTypeName);
									pass.setFalse();
								} else if (!atomInfo.hBondTypeName.equals(hBondType
										.getName())) {
									logger.warn("Atom " + atom.getLabel()
											+ " has hydrogen bond type "
											+ hBondType.getName() + " expected "
											+ atomInfo.hBondTypeName
											+ " correct: new AtomInfo(\""
											+ atom.getLabel() + "\", \""
											+ hBondType.getName() + "\"), ");
									pass.setFalse();
								}
							} else {
								if (hBondType != null) {
									logger.warn("Atom " + atom.getLabel()
											+ " has hydrogen bond type "
											+ hBondType.getName() + " expected no type"
											+ " correct: new AtomInfo(\""
											+ atom.getLabel() + "\", \""
											+ hBondType.getName() + "\"), ");
									pass.setFalse();

								}
							}
						});

		assert atomInfoMap.size() == 0 : "unchecked donors or acceptors";

		return pass.booleanValue();
	}

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		HydrogenBondingType.loadParameters();
	}

	@Test
	public void testOmega1() throws Exception {
		AtomInfo[] donors = new AtomInfo[] { new AtomInfo("N3", "Het5 NH") };
		AtomInfo[] acceptors = new AtomInfo[] { new AtomInfo("O1", "Sulphoxide"),
				new AtomInfo("N2", "Symmetric Het5 N"),
				new AtomInfo("N1", "Symmetric Het6 N") };
		assertTrue(testMolecule("omega_1.mol2", donors, acceptors));
	}

	@Test
	public void testOmega12() throws Exception {
		AtomInfo[] donors = new AtomInfo[] { new AtomInfo("N5", "sp3 N+H2") };
		AtomInfo[] acceptors = new AtomInfo[] { new AtomInfo("N1", "Symmetric Het6 N"),
				new AtomInfo("N2", "Asymetric Het6 N"),
				new AtomInfo("N3", "Asymetric Het6 N") };
		assertTrue(testMolecule("omega_12.mol2", donors, acceptors));
	}

	@Test
	public void testPyrrole() throws Exception {
		AtomInfo[] donors = new AtomInfo[] { new AtomInfo("N1", "Het5 NH") };
		AtomInfo[] acceptors = new AtomInfo[] {};
		assertTrue(testMolecule("pyrrole.mol2", donors, acceptors));
	}

	@Test
	public void test1boz() throws Exception {
		AtomInfo[] donors = new AtomInfo[] { new AtomInfo("N1'", "sp2 N+H"),
				new AtomInfo("N4'", "sp2 N+H2"), new AtomInfo("N2'", "Primary Amine NH2") };
		AtomInfo[] acceptors = new AtomInfo[] { new AtomInfo("N8'", "Symmetric Het6 N"),
				new AtomInfo("N3'", "Symmetric Het6 N"), new AtomInfo("O2'", "Ether"),
				new AtomInfo("O5'", "Ether") };
		assertTrue(testMolecule("1boz_aligned.mol2", donors, acceptors));
	}

	@Test
	public void test2dhf() throws Exception {
		AtomInfo[] donors = new AtomInfo[] { new AtomInfo("N1", "sp2 N+H"),
				new AtomInfo("NA2", "sp2 N+H2"),
				new AtomInfo("N3", "Secondary Amide NH"),
				new AtomInfo("N10", "Phenyl NH"), new AtomInfo("N", "Secondary Amide NH") };
		AtomInfo[] acceptors = new AtomInfo[] { new AtomInfo("OH4", "Ketone"),
				new AtomInfo("N8", "Symmetric Het6 N"), new AtomInfo("O", "Ketone"),
				new AtomInfo("OE1", "Carboxylate"), new AtomInfo("OE2", "Carboxylate"),
				new AtomInfo("O1", "Carboxylate"), new AtomInfo("O2", "Carboxylate") };
		assertTrue(testMolecule("2dhf_aligned.mol2", donors, acceptors));
	}

	@Test
	public void test1dlr() throws Exception {
		AtomInfo[] donors = new AtomInfo[] { new AtomInfo("N1", "sp2 N+H"),
				new AtomInfo("N2", "Primary Amine NH2"), new AtomInfo("N4", "sp2 N+H2") };
		AtomInfo[] acceptors = new AtomInfo[] { new AtomInfo("N3", "Symmetric Het6 N"),
				new AtomInfo("N8", "Symmetric Het6 N"), new AtomInfo("O2'", "Ether"),
				new AtomInfo("O5'", "Ether") };
		assertTrue(testMolecule("1dlr_aligned.mol2", donors, acceptors));
	}

}
