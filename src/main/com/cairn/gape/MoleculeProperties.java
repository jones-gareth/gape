package com.cairn.gape;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;


import com.cairn.gape.feature.Feature.FeatureType;
import com.cairn.gape.feature.HydrogenBondingType;
import com.cairn.gape.molecule.GaMolecule;
import com.cairn.gape.utils.InfoMessageLogger;
import com.cairn.molecule.Atom;
import com.cairn.molecule.Molecule;
import com.opencsv.CSVWriter;

/**
 * Application class to determine simple molecular properties (number of atoms,
 * rotatable bonds, donors etc) for a number of compounds and write them out to
 * a CSV file.
 * 
 * @author Gareth Jones
 *
 */
public class MoleculeProperties {

	private static final InfoMessageLogger infoMessageLogger = new InfoMessageLogger();

	public static void main(String[] args) throws IOException {
		List<GaMolecule> molecules = GaMolecule.loadFiles(args,
				Molecule.FileType.UNKNOWN, Molecule.Source.FILE, infoMessageLogger, true);

		HydrogenBondingType.loadParameters();

		try (CSVWriter csvWriter = new CSVWriter(new FileWriter("mol_properties.csv"),
				',', '"', '"', "\n")) {
			String[] titles = new String[] { "MOLECULE", "N_ATOMS", "N_HYDROGENS",
					"N_HEAVY_ATOMS", "N_ROTATABLE_BONDS", "N_HYDROPHOBIC_ATOMS",
					"N_ACCEPTORS", "N_DONORS", "N_HYDROPHILIC_ATOMS" };
			csvWriter.writeNext(titles);

			for (GaMolecule molecule : molecules) {
				// Remove pharmacophore information
				infoMessageLogger.infoMessageln("\nProcessing " + molecule.getName());
				molecule.removeLonePairs();

				infoMessageLogger.infoMessageln(3,
						"\nFinding Dean and Mills Donors and Acceptors");
				HydrogenBondingType.searchMolecule(molecule, infoMessageLogger);
				infoMessageLogger.infoMessageln(3, "");

				molecule.assignRotatableBonds(false);

				infoMessageLogger.infoMessageln(2, "Finding Features");
				molecule.findFeatures();

				int nHydrogens = molecule.countHydrogens();
				int nAtoms = molecule.getnAtoms();
				int nHeavyAtoms = molecule.getHeavyAtoms().size();
				int nRotatableBonds = molecule.getnRotatableBonds();
				int nHydrophobicAtoms = molecule
						.getFeatureMappings(FeatureType.HYDROPHOBIC_ATOM).getFeatures()
						.size();

				int nAcceptors = 0, nDonors = 0, nDonorHydrogens = 0;
				for (Atom atom : molecule.getAtoms()) {
					if (atom.getDonorType() != null) {
						nDonors++;
					}
					if (atom.getAcceptorType() != null) {
						nAcceptors++;
					}
					if (atom.isDonorHydrogen()) {
						nDonorHydrogens++;
					}
				}

				System.out.println("Molecule " + molecule.getName());
				System.out.println("nHydrogens " + nHydrogens);
				System.out.println("nAtoms " + nAtoms);
				System.out.println("nHeavyAtoms " + nHeavyAtoms);
				System.out.println("nRotatableBonds " + nRotatableBonds);
				System.out.println("nAcceptors " + nAcceptors);
				System.out.println("nDonors " + nDonors);
				System.out.println("nDonorHydrogens " + nDonorHydrogens);
				System.out.println("nHydrophobicAtoms " + nHydrophobicAtoms);
				System.out.println();

				String[] values = new String[] { molecule.getName(),
						String.valueOf(nAtoms), String.valueOf(nHydrogens),
						String.valueOf(nHeavyAtoms), String.valueOf(nRotatableBonds),
						String.valueOf(nHydrophobicAtoms), String.valueOf(nAcceptors),
						String.valueOf(nDonors), String.valueOf(nAcceptors + nDonors) };
				csvWriter.writeNext(values);
			}
		}
	}
}
