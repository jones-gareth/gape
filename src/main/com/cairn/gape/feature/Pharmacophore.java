package com.cairn.gape.feature;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.cairn.gape.Superposition;
import com.cairn.gape.chromosome.SuperpositionChromosome;
import com.cairn.gape.feature.Feature.FeatureType;
import com.cairn.gape.molecule.GaMolecule;
import com.cairn.gape.utils.InfoMessageLogger;
import com.cairn.molecule.Atom;
import com.cairn.molecule.AtomType;
import com.cairn.molecule.Bond;
import com.cairn.molecule.BondType;
import com.cairn.molecule.Molecule;
import com.cairn.molecule.Molecule.FileType;
import com.cairn.molecule.Molecule.Source;

/**
 * Class for describing a pharmacophore.
 * 
 * @author Gareth Jones
 * 
 */
public class Pharmacophore {
	private volatile int nMolecules;

	private volatile List<Feature> allFeatures;

	private volatile Map<Integer, List<Feature>> featuresByMolecule;

	private static final InfoMessageLogger infoMessageLogger = new InfoMessageLogger();

	/**
	 * Empty constructor.
	 */
	public Pharmacophore() {
		if (HydrogenBondingType.getDonorAcceptorTypes() == null)
			HydrogenBondingType.loadParameters();
	}

	/**
	 * Creates Pharmacophore from saved file
	 * 
	 * @param File
	 */
	public Pharmacophore(String file) {
		this();
		try {
			BufferedReader in = new BufferedReader(new FileReader(file));

			readDefinition(in);
			in.close();
		} catch (FileNotFoundException ex) {
			throw new RuntimeException(ex.toString());
		} catch (IOException ex) {
			throw new RuntimeException(ex.toString());
		}
	}

	/**
	 * Create the pharmacophore from a GAPE chromosome
	 * 
	 * @param chrom
	 */
	public Pharmacophore(SuperpositionChromosome chrom) {
		this();

		allFeatures = new ArrayList<Feature>();
		featuresByMolecule = new HashMap<>();

		if (chrom != SuperpositionChromosome.getCurrentChromosome())
			throw new RuntimeException("Pharmacphore: not currently fitted chromosome");
		nMolecules = chrom.getMolecules().size();
		int molNo = 0;
		for (GaMolecule molecule : chrom.getMolecules()) {
			List<Feature> list = new ArrayList<Feature>();
			featuresByMolecule.put(molNo, list);

			for (FeatureMapping features : molecule.getFeatureMappings().values()) {
				for (Feature feature : features.getFeatures()) {
					if (feature.isPharmPoint()) {
						feature.pharmFeatureGeometry = feature.getPharmFeatureGeometry();
						list.add(feature);
						allFeatures.add(feature);
					}
				}
			}
			molNo++;
		}
	}

	/**
	 * Reads a pharmacophore file or string and creates features.
	 * 
	 * @param file
	 * @see Feature#featureFromPharmDescription(String, int);
	 */
	private void readDefinition(BufferedReader in) throws RuntimeException {
		ArrayList<Feature> featureList = new ArrayList<Feature>();
		try {

			String line = in.readLine();
			if (!line.equals("Pharmacophore"))
				throw new RuntimeException("not a pharmacophore file or string");

			String nMols = getValue("no_molecules", in.readLine());
			nMolecules = Integer.parseInt(nMols);
			for (int i = 0; i < nMolecules; i++) {
				String molNo = getValue("molecule_number", in.readLine());
				int moleculeNo = Integer.parseInt(molNo);
				moleculeNo--;
				if (moleculeNo != i)
					throw new RuntimeException(
							"Pharmacophore: readFile: molecule number mismatch");
				while (true) {
					line = in.readLine();
					if (line.startsWith("end_features"))
						break;
					String featureString = getValue("feature", line);
					Feature f = Feature.featureFromPharmDescription(featureString,
							moleculeNo);
					featureList.add(f);
				}
			}
		} catch (IOException ex) {
			throw new RuntimeException("Pharmacophore: readFile " + ex);
		}

		allFeatures = featureList;
		addFeaturesByMolecule();
	}

	/**
	 * Creates the featuresByMolecule hash from the allfeatures list. Adds
	 * feature set numbers. Feature set numbers are unique to each molecule.
	 */
	private void addFeaturesByMolecule() {
		featuresByMolecule = new HashMap<>();

		Map<FeatureType, Integer> featureSetNos = new HashMap<>();
		// this double loop looks pretty inefficient, but I don't suppose it
		// matters.
		for (int i = 0; i < nMolecules; i++) {
			featureSetNos.clear();
			ArrayList<Feature> list = new ArrayList<Feature>();
			featuresByMolecule.put(i, list);
			for (Feature feature : allFeatures) {
				if (feature.moleculeNo != i)
					continue;
				FeatureType featureSet = feature.getFeatureType();
				int featureSetNo = 0;
				if (featureSetNos.containsKey(featureSetNo))
					featureSetNo = featureSetNos.get(featureSetNo) + 1;
				feature.featureSetNo = featureSetNo;
				featureSetNos.put(featureSet, featureSetNo);
				list.add(feature);
			}
		}

	}

	/**
	 * Creates a Pharmacophore from an array of moleucles which have
	 * pharmacophore descriptions. When you read a GaMolecule pharm descriptions
	 * will be read in, but the pharmacophore won't be created unless you use
	 * this routine.
	 * 
	 * @param molecules
	 */
	public void createFromMoleculePharmDescriptions(GaMolecule molecules[]) {
		ArrayList<Feature> featureList = new ArrayList<Feature>();
		nMolecules = molecules.length;

		for (int i = 0; i < nMolecules; i++) {
			for (int j = 0; j < molecules[i].getNPharmPoints(); j++) {
				String pharmDescription = molecules[i].getPharmDescription(j);
				Feature f = Feature.featureFromPharmDescription(pharmDescription, i);
				featureList.add(f);
			}
		}

		allFeatures = featureList;
		addFeaturesByMolecule();
	}

	/**
	 * Parses a "key = value" line and returns value.
	 * 
	 * @param key
	 * @param line
	 * @return
	 */
	private String getValue(String key, String line) {
		line = line.trim();
		if (!line.startsWith(key + " = "))
			throw new RuntimeException("key " + key + " not found in line " + line);
		int equals = line.indexOf('=');
		String value = line.substring(equals + 1).trim();
		return value;
	}

	/**
	 * Creates a molecule of dummy atoms representing feature points.
	 * 
	 * @return
	 */
	public Molecule toMolecule() {

		List<Atom> atoms = new ArrayList<>();
		List<Bond> bonds = new ArrayList<>();
		List<double[]> coords = new ArrayList<>();
		int nAtoms = 0;
		int nBonds = 0;

		for (Feature feature : allFeatures) {
			PharmFeatureGeometry geometry = feature.pharmFeatureGeometry;

			if (geometry.getGeometry() == PharmFeatureGeometry.Geometry.MULTI_VECTOR) {

				MultiVectorPharmFeatureGeometry mvg = (MultiVectorPharmFeatureGeometry) geometry;
				int no = mvg.getNVectors();

				Atom center = new Atom(nAtoms, AtomType.Type.DU);
				atoms.add(center);
				coords.add(mvg.getCenter());
				nAtoms++;
				for (int i = 0; i < no; i++) {
					Atom vertex = new Atom(nAtoms, AtomType.Type.DU);
					coords.add(mvg.getEnd(i));
					atoms.add(vertex);
					nAtoms++;
					Bond bond = new Bond(nBonds, center, vertex, BondType.Type.SINGLE);
					bonds.add(bond);
					nBonds++;
				}

			} else {
				int no = geometry.getNPoints();
				if (no == 1) {
					Atom atom = new Atom(nAtoms, AtomType.Type.DU);
					atoms.add(atom);
					coords.add(geometry.getPoint(0));
					nAtoms++;
				} else if (no == 2) {
					Atom atom1 = new Atom(nAtoms, AtomType.Type.DU);
					Atom atom2 = new Atom(nAtoms + 1, AtomType.Type.DU);
					Bond bond = new Bond(nBonds, atom1, atom2, BondType.Type.SINGLE);
					atoms.add(atom1);
					atoms.add(atom2);
					bonds.add(bond);
					coords.add(geometry.getPoint(0));
					coords.add(geometry.getPoint(1));
					nAtoms += 2;
					nBonds++;
				} else if (no == 3) {
					Atom atom1 = new Atom(nAtoms, AtomType.Type.DU);
					Atom atom2 = new Atom(nAtoms + 1, AtomType.Type.DU);
					Atom atom3 = new Atom(nAtoms + 2, AtomType.Type.DU);
					Bond bond1 = new Bond(nBonds, atom1, atom2, BondType.Type.SINGLE);
					Bond bond2 = new Bond(nBonds + 1, atom1, atom3, BondType.Type.SINGLE);
					Bond bond3 = new Bond(nBonds + 2, atom2, atom3, BondType.Type.SINGLE);
					atoms.add(atom1);
					atoms.add(atom2);
					atoms.add(atom3);
					bonds.add(bond1);
					bonds.add(bond2);
					bonds.add(bond3);
					coords.add(geometry.getPoint(0));
					coords.add(geometry.getPoint(1));
					coords.add(geometry.getPoint(2));
					nAtoms += 3;
					nBonds += 3;
				}
			}
		}

		Molecule mol = new Molecule("Pharmacophore molecule", atoms, bonds, coords);

		return mol;
	}

	/**
	 * Main routine- used for testing and debugging. Alo provides a number of
	 * utilities for the user to use in manipulating pharmacophores and
	 * extracting pharmacophore descriptions from GAPE strucutre files.
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		try {
			HydrogenBondingType.loadParameters();

			if (args.length == 0) {
				System.out.println("Usage: Pharmacophore command [args..]");
				System.out.println("commands:");
				System.out.println("-extract-pharm <molfile> <pharmfile>");
				System.out.println("-extract-base-only-pharm <molfile> <pharmfile>");
				System.out.println("-create <pharm_file>");
				System.out.println("-transform <pharm_file> <new_pharmfile>");
				System.out.println("-to-molecule <pharm_file> <pharm_mol_file>");
				System.out
						.println("-extract-pharm-to-molecule <mol_file> <pharm_mol_file>");

				System.exit(0);
			}

			String command = args[0];

			if (command.equalsIgnoreCase("-extract-pharm")) {
				String molFile = args[1];
				String pharmFile = args[2];
				List<GaMolecule> molecules = loadFiles(molFile);
				moleculesToFile(molecules, pharmFile, false);
			}

			else if (command.equalsIgnoreCase("-extract-base-only-pharm")) {
				String molFile = args[1];
				String pharmFile = args[2];
				List<GaMolecule> molecules = loadFiles(molFile);
				moleculesToFile(molecules, pharmFile, true);
			}

			else if (command.equalsIgnoreCase("-create")) {
				String pharmFile = args[1];
				new Pharmacophore(pharmFile);
			}

			else if (command.equalsIgnoreCase("-transform")) {
				String pharmFile = args[1];
				String outFile = args[2];
				Pharmacophore pharm = new Pharmacophore(pharmFile);
				pharm.pharmacophoreToFile(outFile);
			}

			else if (command.equalsIgnoreCase("-to-molecule")) {
				String pharmFile = args[1];
				String molFile = args[2];

				Pharmacophore pharm = new Pharmacophore(pharmFile);
				Molecule mol = pharm.toMolecule();
				mol.write(molFile, "Pharmacphore points");
			}

			else if (command.equalsIgnoreCase("-extract-pharm-to-molecule")) {
				String molFile = args[1];
				String molPharmFile = args[2];
				List<GaMolecule> molecules = loadFiles(molFile);
				String str = moleculesToString(molecules, false);
				BufferedReader in = new BufferedReader(new StringReader(str));
				Pharmacophore pharm = new Pharmacophore();
				pharm.readDefinition(in);
				in.close();
				Molecule mol = pharm.toMolecule();
				mol.write(molPharmFile, "Pharmacphore points");
			}

			else {
				System.out.println("unknown command " + command);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * Loads and setups input structures
	 * 
	 * @param files
	 * @return
	 */
	private static List<GaMolecule> loadFiles(String file) {
		List<GaMolecule> molecules = GaMolecule.loadFiles(new String[] { file },
				FileType.UNKNOWN, Source.FILE);
		Superposition problem = new Superposition();
		problem.initFromResource();

		for (GaMolecule molecule : molecules) {
			molecule.setProblem(problem);
			setupMolecule(molecule);
		}
		return molecules;
	}

	/**
	 * Sets up a molecule. Stripped down version of {@link GaMolecule#setup}.
	 * 
	 * @param molecule
	 *            Structure to initialize.
	 */
	private static void setupMolecule(GaMolecule molecule) {

		molecule.removeLonePairs();
		HydrogenBondingType.searchMolecule(molecule, infoMessageLogger);

		LonePairAddition.addLonePairs(molecule);
		infoMessageLogger.infoMessageln(2, "Finding Features");
		molecule.findFeatures();

	}

	/**
	 * Extracts pharmacophore information from molecules and dumps it out to the
	 * file. Works just by extracting pharmacophore descriptions, so these need
	 * to be present.
	 * 
	 * @param molecules
	 * @param file
	 * @param baseOnly
	 */
	private static void moleculesToFile(List<GaMolecule> molecules, String file,
			boolean baseOnly) {
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(file));
			// BufferedWriter out = new BufferedWriter(new
			// PrintWriter(System.out));

			out.write(moleculesToString(molecules, baseOnly));
			out.close();

		} catch (IOException ex) {
			System.out.println(ex);
		}
	}

	/**
	 * Extracts pharmacophore information from molecules and return it as a
	 * string Works just by extracting pharmacophore descriptions, so these need
	 * to be present. Set baseOnly to extract only the baseMolecule.
	 * 
	 * @param molecules
	 * @param baseOnly
	 * @return
	 */
	private static String moleculesToString(List<GaMolecule> molecules, boolean baseOnly) {
		StringBuffer buf = new StringBuffer();
		buf.append("Pharmacophore\n");
		if (baseOnly)
			buf.append("no_molecules = 1\n");
		else
			buf.append("no_molecules = " + molecules.size() + "\n");
		for (int i = 0; i < molecules.size(); i++) {
			GaMolecule mol = molecules.get(i);
			if (baseOnly) {
				if (mol.getName().endsWith(" Base Molecule"))
					buf.append("molecule_number = 1\n");
				else
					continue;
			} else
				buf.append("molecule_number = " + String.valueOf(i + 1) + "\n");
			for (int j = 0; j < mol.getNPharmDescriptions(); j++) {
				buf.append("feature = " + mol.getPharmDescription(j) + "\n");
			}
			buf.append("end_features\n");
		}

		return buf.toString();
	}

	/**
	 * Writes out the pharmacophore to a file. In this case the definitions are
	 * created from feature objects. So you can compare the pharmacophore file
	 * extracted from an overlay with that rebuilt from this class.
	 * 
	 * @param file
	 */
	private void pharmacophoreToFile(String file) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(file));
			out.write("Pharmacophore\n");
			out.write("no_molecules = " + nMolecules + "\n");
			for (int i = 0; i < nMolecules; i++) {
				out.write("molecule_number = " + String.valueOf(i + 1) + "\n");
				List<Feature> list = featuresByMolecule.get(i);
				for (Feature feature : list) {
					out.write("feature = " + feature.featureLabel() + "\n");
				}
				out.write("end_features\n");
			}
			out.close();
		} catch (IOException ex) {
			System.out.println(ex);
		}
	}

	/**
	 * @return the allFeatures
	 */
	public List<Feature> getAllFeatures() {
		return allFeatures;
	}

	/**
	 * @return the featuresByMolecule
	 */
	public Map<Integer, List<Feature>> getFeaturesByMolecule() {
		return featuresByMolecule;
	}

}
