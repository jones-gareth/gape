package com.cairn.gape.molecule;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;
import com.cairn.gape.MultiMolSuperposition;
import com.cairn.gape.feature.HydrogenBondingType;
import com.cairn.gape.feature.LonePairAddition;
import com.cairn.gape.utils.InfoMessageLogger;
import com.cairn.molecule.Atom;
import com.cairn.molecule.Bond;
import com.cairn.molecule.Molecule;
import com.cairn.molecule.RotatableBond;

/**
 * A class to store multi conformer molecules for use in GAPE.
 * 
 * @author Gareth Jones
 * 
 */
public class MultiGaMolecule extends GaMolecule {

	private static Logger logger;
	static {
		logger = Logger.getLogger(MultiGaMolecule.class);
		logger.setLevel(Level.DEBUG);
	}

	private final List<Conformer> conformers = new ArrayList<Conformer>();;
	private volatile int currentConformerNo;

	private class Conformer {
		private List<double[]> coords;
		private double energy;
	}

	/**
	 * @return the number of conformers
	 */
	public int nConformers() {
		return conformers.size();
	}

	/**
	 * Copies conformational coordinates to the molecule
	 * 
	 * @param conformerNo
	 */
	public void setConformer(int conformerNo) {
		Conformer conformer = conformers.get(conformerNo);
		List<double[]> conformerCoords = conformer.coords;
		for (int i = 0; i < getnAtoms(); i++) {
			if (conformerCoords.size() == i) {
				conformerCoords.add(new double[4]);
			}
			assert conformerCoords.get(i) != null : "null conformer coordinates";
			assert getCoord(i) != null : "null molecule coords";
			Coord.copy(conformerCoords.get(i), getCoord(i));
		}
		currentConformerNo = conformerNo;
	}

	/**
	 * Copy molecules coordinates back to current conformer- used for adding
	 * lone pairs.
	 */
	private void copyCoordinatesToConformer() {
		Conformer conformer = conformers.get(currentConformerNo);
		List<double[]> conformerCoords = conformer.coords;
		for (int i = 0; i < getnAtoms(); i++) {
			if (conformerCoords.size() == i) {
				conformerCoords.add(new double[4]);
			}
			assert conformerCoords.get(i) != null : "null conformer coordinates";
			assert getCoord(i) != null : "null molecule coords";
			Coord.copy(getCoord(i), conformerCoords.get(i));
		}
	}

	/**
	 * @return any energy associated with the current conformation.
	 */
	public double getConformerEnergy() {
		return conformers.get(currentConformerNo).energy;
	}

	/**
	 * Creates an array of molecules. Specify filename(s), structure format,
	 * source, log level and default output stream. Uses init routine to prepare
	 * molecules (unlike loadFiles in Molecule)- if init is set.
	 * 
	 * @param files
	 * @param type
	 * @param source
	 * @param logLevel
	 * @param out
	 * @return
	 * 
	 * @see Molecule#init()
	 * @see #loadSybylMol2(BufferedReader)
	 * @see #loadSybylMol2(BufferedReader)
	 */
	public static List<MultiGaMolecule> loadMultiMolFiles(String files[], FileType type,
			Source source, InfoMessageLogger infoMessageLogger, boolean init) {

		List<MultiGaMolecule> molecules = new ArrayList<MultiGaMolecule>();
		List<MultiGaMolecule> conformers = new ArrayList<MultiGaMolecule>();
		MultiGaMolecule currentMolecule = null;

		int molNo = 0;
		for (int i = 0; i < files.length; i++) {
			infoMessageLogger.infoMessageln("loading " + files[i]);
			BufferedReader in = openReader(files[i], source);
			if (in == null)
				continue;
			if (type == FileType.UNKNOWN) {
				type = Molecule.getType(files[i]);
			}
			while (true) {
				MultiGaMolecule molecule = new MultiGaMolecule();
				molecule.setInfoMessageLogger(infoMessageLogger);
				try {
					if (type == FileType.MOL2) {
						molecule.loadSybylMol2(in);
					} else if (type == FileType.SDF) {
						molecule.loadSdfMol(in);
					}
				} catch (IOException | MolReadException ex) {
					break;
				}
				if (molecule.getName().endsWith("Pharmacophore"))
					continue;

				if (init)
					molecule.init();
				else
					molecule.build();

				if (currentMolecule == null || !equalMolecules(molecule, currentMolecule)) {
					if (currentMolecule != null) {
						// create multiconformer molecule
						molNo++;
						if (currentMolecule.getName().equals("")
								|| currentMolecule.getName().equals("****"))
							currentMolecule.setName("Structure_" + molNo);
						currentMolecule.addConformers(conformers);
						conformers.clear();
					}
					currentMolecule = molecule;
					molecules.add(currentMolecule);
				}

				conformers.add(molecule);
			}

			// create multiconformer molecule
			molNo++;
			if (currentMolecule.getName().equals("")
					|| currentMolecule.getName().equals("****"))
				currentMolecule.setName("Structure_" + molNo);
			currentMolecule.addConformers(conformers);
			conformers.clear();

			try {
				in.close();
			} catch (IOException ex) {
				System.err.println("IO error closing " + ex);
			}
		}

		return molecules;
	}

	/**
	 * Not this assumes the both conformers have the same numbering- if this is
	 * not true a graph isomorphism/exact match search may be required.
	 * 
	 * @param molecule1
	 * @param molecule2
	 * @return true if two molecules are different conformers of the same
	 *         molecule.
	 */
	public static boolean equalMolecules(GaMolecule molecule1, GaMolecule molecule2) {

		if (molecule1.getnAtoms() != molecule2.getnAtoms())
			return false;
		if (molecule1.getnBonds() != molecule2.getnBonds())
			return false;

		List<Atom> atoms1 = molecule1.getAtoms();
		List<Atom> atoms2 = molecule2.getAtoms();
		for (int i = 0; i < molecule1.getnAtoms(); i++)
			if (atoms1.get(i).getAtomType() != atoms2.get(i).getAtomType())
				return false;

		List<Bond> bonds1 = molecule1.getBonds();
		List<Bond> bonds2 = molecule2.getBonds();
		for (int i = 0; i < molecule1.getnBonds(); i++)
			if (bonds1.get(i).getBondType() != bonds2.get(i).getBondType())
				return false;

		return true;

	}

	/**
	 * Adds other molecules as conformers.
	 * 
	 * @param conformerMolecules
	 */
	public void addConformers(List<MultiGaMolecule> conformerMolecules) {
		if (CollectionUtils.isEmpty(conformerMolecules))
			return;

		for (MultiGaMolecule molecule : conformerMolecules) {
			assert equalMolecules(this, molecule) : "addConformers: conformers do not match";
			Conformer conformer = new Conformer();
			// one of the conformers is this molecule- we need to copy (not just
			// reference) it's coordinates.
			if (molecule == this) {
				conformer.coords = new ArrayList<>(getCoords().size());
				;
				for (int i = 0; i < getnAtoms(); i++) {
					conformer.coords.add(new double[4]);
					Coord.copy(getCoord(i), conformer.coords.get(i));
				}
			} else {
				conformer.coords = molecule.getCoords();
			}
			String energyStr = molecule.getSdfField("mmff94s");
			if (StringUtils.isNotBlank(energyStr))
				conformer.energy = Double.parseDouble(energyStr);
			conformers.add(conformer);
		}

		infoMessageLogger.infoMessageln(0, "Added " + conformerMolecules.size()
				+ " conformation(s) to " + getName());
		return;

	}

	/**
	 * Adds lone pairs to all conformers.
	 * 
	 */
	private void addLonePairs() {
		// add lone pairs to molecule using first conformer coordinates
		setConformer(0);
		double hBondLen = getProblem().getDoubleValue("h_bond_len");
		LonePairAddition.addLonePairs(this, hBondLen);
		// sync coordinates back to conformer
		copyCoordinatesToConformer();

		// update lone pair coordinates in each other conformer
		for (int i = 1; i < conformers.size(); i++) {
			// rebuild the lone pair coordinates
			setConformer(i);
			LonePairAddition.updateLonePairs(this);
			copyCoordinatesToConformer();
		}
	}

	/**
	 * Sets molecule parameters using settings from program instance. Also does
	 * stuff like finding starting energy and randomizing.
	 * 
	 */
	@Override
	public void setup() {
		setUseGray(getProblem().getBooleanValue("graycode"));

		if (isUseGray())
			infoMessageLogger.infoMessageln(2, "Using Gray encoding");
		else
			infoMessageLogger.infoMessageln(2, "Using binary encoding");

		if (countHydrogens() == 0)
			throw new RuntimeException("SD Molecule " + getName() + " has no hydrogens\n"
					+ "GAPE requires that SD files contain all "
					+ "hydrogens or fill_valence is set");

		if (findFeatures) {
			infoMessageLogger.infoMessageln(3,
					"\nFinding Dean and Mills Donors and Acceptors");
			HydrogenBondingType.searchMolecule(this, infoMessageLogger);
			infoMessageLogger.infoMessageln(3, "");

			// We do want lone pairs to rotate if possible so add them
			// before assigning rotatable bonds

			addLonePairs();
		}

		// Find All Rotatable Bonds
		assignRotatableBonds(false);
		// only keep those for terminal donors or acceptors.
		infoMessageLogger
				.infoMessageln(3,
						"Retaining only terminal rotatable bonds to donor hydrogens or lone pairs");
		List<RotatableBond> newRotatableBonds = new ArrayList<RotatableBond>();
		for (RotatableBond rotatableBond : getRotatableBonds()) {
			if (isTerminalRotatable(rotatableBond))
				newRotatableBonds.add(rotatableBond);
		}
		setRotatableBonds(newRotatableBonds);
		infoMessageLogger.infoMessageln(2, "Found " + newRotatableBonds.size()
				+ " rotatable bonds");

		setIgnoreTorsion(true);
		if (findFeatures) {
			infoMessageLogger.infoMessageln(2, "Finding Features");
			findFeatures();
		}

		if (getProblem().hasKey("relax_molecule")) {
			setRelaxMolecule(getProblem().getBooleanValue("relax_molecule"));
			setRelaxMaxDistance(getProblem().getDoubleValue("relax_max_distance"));
			setRelaxMaxAngle(getProblem().getDoubleValue("relax_max_angle"));
		}

		if (isRelaxMolecule()) {
			infoMessageLogger.infoMessageln(2, "Relaxing molecule after fitting");
			infoMessageLogger.infoMessageln(2, "Max distance " + getRelaxMaxDistance()
					+ " Max angle " + getRelaxMaxAngle());
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.molecule.GaMolecule#conformationalEnergy()
	 */
	@Override
	public double conformationalEnergy() {
		// just use the determined conformer energy (if available).
		return getConformerEnergy();
	}

	/**
	 * Writes out a mol2 of sdf file based on file suffix. Creates an entry for
	 * each conformation.
	 * 
	 * @param format
	 * @param out
	 * @param comment
	 * @param incLP
	 * @throws IOException
	 */
	public void writeAllConformations(FileType format, Writer out, String comment,
			boolean incLP) throws IOException {
		String saveName = getName();

		for (int i = 0; i < nConformers(); i++) {
			setConformer(i);
			setName(saveName + "_" + String.valueOf(i + 1));
			if (format == Molecule.FileType.SDF)
				writeSdfMol(out, comment);
			else
				writeSybylMol2(out, comment, incLP);
		}

		setName(saveName);
	}

	/**
	 * Writes out a mol2 of sdf file based on file suffix. Creates an entry for
	 * each conformation.
	 * 
	 * @param molecules
	 * @param file
	 * @param comment
	 */
	static public void writeAllConformations(List<MultiGaMolecule> molecules,
			String file, String comment) {
		Molecule.FileType type = Molecule.getType(file);
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(file));
			for (MultiGaMolecule molecule : molecules) {
				molecule.writeAllConformations(type, out, comment, true);
			}
			out.close();
		} catch (IOException ex) {
			throw new RuntimeException("IOException: " + ex.toString());
		}
	}

	/**
	 * Test routine
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		try {
			if (args.length < 1) {
				System.out.println("Usage: " + MultiGaMolecule.class.getName()
						+ " -version | <conf file> [molecule files..]");
				System.exit(0);
			}

			String confFile = args[0];
			String[] molFiles = ArrayUtils.remove(args, 0);
			MultiMolSuperposition superposition = new MultiMolSuperposition();
			superposition.init(confFile);

			List<MultiGaMolecule> molecules = loadMultiMolFiles(molFiles,
					FileType.UNKNOWN, Source.FILE, superposition.getInfoMessageLogger(),
					true);
			List<GaMolecule> gaMolecules = molecules.stream().map(m -> m)
					.collect(Collectors.toList());
			superposition.setupMolecules(gaMolecules);
			writeAllConformations(molecules, "gape_conformers.mol2", null);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
