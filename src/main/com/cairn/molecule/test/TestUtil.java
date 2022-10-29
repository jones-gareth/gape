package com.cairn.molecule.test;

import java.io.BufferedReader;
import java.io.InputStreamReader;

import com.cairn.gape.molecule.GaMolecule;
import com.cairn.molecule.Molecule;

public class TestUtil {

	public static Molecule getMoleculeFromSybylResource(String resource) throws Exception {
		if (!resource.startsWith("/")) {
			resource = "/com/cairn/molecule/test/data/" + resource;
		}
		BufferedReader in = new BufferedReader(new InputStreamReader(
				TestCharge.class.getResourceAsStream(resource)));
		Molecule molecule = new Molecule();
		molecule.loadSybylMol2(in);
		in.close();
		molecule.build();
		return molecule;
	}

	public static GaMolecule getGaMoleculeFromSybylResource(String resource)
			throws Exception {
		if (!resource.startsWith("/")) {
			resource = "/com/cairn/molecule/test/data/" + resource;
		}
		BufferedReader in = new BufferedReader(new InputStreamReader(
				TestCharge.class.getResourceAsStream(resource)));
		GaMolecule molecule = new GaMolecule();
		molecule.loadSybylMol2(in);
		in.close();
		molecule.init();
		return molecule;
	}
}
