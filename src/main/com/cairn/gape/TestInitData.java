package com.cairn.gape;

import com.cairn.gape.molecule.GaMolecule;

/**
 * Test class for Checking that chromosome initialization works.
 * 
 * @author Gareth Jones
 * 
 */
public class TestInitData {

	// Initialization data- not only valid for glaxo.mol in Dean test set

	// corresponds to these chromosomes.

	// 00100010 00101111 00010000 10100000 01100010 01111110 01101101 01110000
	// 10000000 11011100

	// 00100010 00101111 00010000 10100000 01100010 01111110 01101101 11110000
	// 10000000 11011100

	// 00110010 00111111 00010000 11100000 01100010 01111110 01101101 11110000
	// 10000000 11011100

	static boolean initData[][] = new boolean[][] {
			new boolean[] { false, false, true, false, false, false, true,
					false, false, false, true, false, true, true, true, true,
					false, false, false, true, false, false, false, false,
					true, false, true, false, false, false, false, false,
					false, true, true, false, false, false, true, false, false,
					true, true, true, true, true, true, false, false, true,
					true, false, true, true, false, true, false, true, true,
					true, false, false, false, false, true, false, false,
					false, false, false, false, false, true, true, false, true,
					true, true, false, false },
			new boolean[] { false, false, true, false, false, false, true,
					false, false, false, true, false, true, true, true, true,
					false, false, false, true, false, false, false, false,
					true, false, true, false, false, false, false, false,
					false, true, true, false, false, false, true, false, false,
					true, true, true, true, true, true, false, false, true,
					true, false, true, true, false, true, true, true, true,
					true, false, false, false, false, true, false, false,
					false, false, false, false, false, true, true, false, true,
					true, true, false, false },
			new boolean[] { false, false, true, true, false, false, true,
					false, false, false, true, true, true, true, true, true,
					false, false, false, true, false, false, false, false,
					true, true, true, false, false, false, false, false, false,
					true, true, false, false, false, true, false, false, true,
					true, true, true, true, true, false, false, true, true,
					false, true, true, false, true, true, true, true, true,
					false, false, false, false, true, false, false, false,
					false, false, false, false, true, true, false, true, true,
					true, false, false } };

	/**
	 * Runs the GA as an application.
	 * 
	 * @param args
	 */
	public static void main(String args[]) {

		String configFile = "/home/gjones/gape/test_init/minimizer.conf";
		String moleculeFile = "/home/gjones/gape/test_init/glaxo.mol2";

		GaMolecule molecule = new GaMolecule();
		molecule.setFindFeatures(false);
		molecule.loadFile(moleculeFile);
		TorsionalMinimizer tm = new TorsionalMinimizer();
		tm.minimize(configFile, molecule, initData);
	}

}
