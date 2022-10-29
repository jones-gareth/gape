package com.cairn.gape.feature;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import org.apache.log4j.Logger;

import com.cairn.gape.molecule.GaMolecule;
import com.cairn.gape.utils.InfoMessageLogger;
import com.cairn.molecule.Atom;
import com.cairn.molecule.AtomType;
import com.cairn.molecule.MolPattern;
import com.cairn.molecule.Molecule;
import com.cairn.molecule.PatternMatch;

/**
 * Class to represent donor and acceptor types. These types are loaded from a
 * text file.
 * 
 * The default hydrogen bonding types are based on Mill and Dean (JCAMD, v10,
 * 1996, pp607-622). The default hydrogen bond definitions are in a resource
 * file (com.cairn.gape.donor_acceptor_types.txt).
 * 
 * @author Gareth Jones
 * 
 */
public class HydrogenBondingType {

	// array of donor/acceptor types
	private static HydrogenBondingType donorAcceptorTypes[];

	// type is donor/acceptor, weight is the number of atoms in the sln, or
	// assigned weight
	private final int id, type, weight;

	// name and query sln
	private final String name, sln;

	private volatile MolPattern pattern;

	private volatile DonorAcceptorPatternMatch match;

	private final double probability;

	// acceptor gemoetries
	public enum AcceptorGeometry {
		AG_NONE, DIR, PLANE, CONE
	}

	private AcceptorGeometry geometry;

	// type definitions
	public static final int DONOR = 1, ACCEPTOR = 2;

	private static final Logger logger = Logger.getLogger(HydrogenBondingType.class);

	/**
	 * Loads in the default hydrogen bond definitions
	 * 
	 */
	public static void loadParameters() {
		loadParameterFile();
	}

	/**
	 * Loads in a specific file/resource containing hydrogen bonding types
	 * 
	 * @param file
	 */
	public static void loadParameters(String file) {
		loadParameterFile(file);
	}

	/**
	 * Creates a hydrogen bonding type from file definition
	 * 
	 * @param id
	 * @param t
	 * @param n
	 * @param s
	 * @param p
	 * @param g
	 */
	private HydrogenBondingType(int id, int t, String n, String s, double p,
			AcceptorGeometry g, Integer w) {

		this.id = id;
		type = t;
		name = n;
		sln = s;
		probability = p;
		geometry = g;

		pattern = MolPattern.generateMolPattern(sln);

		match = new DonorAcceptorPatternMatch(pattern);
		if (w == null) {
			weight = pattern.getPatternMol().getnAtoms();
		} else {
			weight = w;
		}
		logger.debug(info());
	}

	/**
	 * @return an informative description about the type.
	 */
	public String info() {
		String info = (type == DONOR) ? "Donor" : "Acceptor";
		info += " " + name + " [Prob " + probability + "]";
		if (type == ACCEPTOR) {
			info += " [geometry ";
			if (geometry == AcceptorGeometry.AG_NONE)
				info += "none";
			else if (geometry == AcceptorGeometry.DIR)
				info += "dir";
			else if (geometry == AcceptorGeometry.PLANE)
				info += "plane";
			else if (geometry == AcceptorGeometry.CONE)
				info += "cone";
			info += "]";
		}
		info += " " + sln;
		info += " [wt " + weight + "]";

		return info;
	}

	/**
	 * Searches a molecule and finds and marks all instances of this donor or
	 * acceptor
	 * 
	 * @param mol
	 * */
	public void findType(Molecule mol) {
		logger.debug("searching for " + sln);
		match.match(mol);
	}

	/**
	 * This class handles the matching of the donor/acceptor sln against a
	 * molecule
	 */
	private class DonorAcceptorPatternMatch extends PatternMatch {
		DonorAcceptorPatternMatch(MolPattern p) {
			super(p);
		}

		/*
		 * Sets the hydrogen bonding type of the matching acceptor or donor
		 * hydrogen. In the event of a tie with another hydrogen-bonding type
		 * select the one with the highest weight (number of atoms in sln
		 * pattern) then highest probability. (non-Javadoc)
		 * 
		 * @see com.cairn.molecule.PatternMatch#process()
		 */
		@Override
		public void process() {
			int match[] = queryMatches[nMatches - 1];
			Atom daAtom = target.getAtom(match[0]);
			if (type == DONOR) {
				if (daAtom.getDonorType() == null)
					daAtom.setDonorType(HydrogenBondingType.this);
				else if (daAtom.getDonorType().weight < weight)
					daAtom.setDonorType(HydrogenBondingType.this);
				else if (daAtom.getDonorType().weight == weight
						&& daAtom.getDonorType().probability < probability)
					daAtom.setDonorType(HydrogenBondingType.this);

				for (Atom test : daAtom.getNotDummyNeighbours()) {
					if (test.getAtomType() == AtomType.Type.H)
						test.setDonorHydrogen(true);
				}
			} else {
				if (daAtom.getAtomType() == AtomType.Type.NPL3
						&& daAtom.getnNotDummyNeighbours() == 3)
					return;
				if (daAtom.getAcceptorType() == null)
					daAtom.setAcceptorType(HydrogenBondingType.this);
				else if (daAtom.getAcceptorType().weight < weight)
					daAtom.setAcceptorType(HydrogenBondingType.this);
				else if (daAtom.getAcceptorType().weight == weight
						&& daAtom.getAcceptorType().probability < probability)
					daAtom.setAcceptorType(HydrogenBondingType.this);
			}
		}
	}

	/**
	 * Fills in the array of hydrogen bonding types from a file.
	 * 
	 * @param file
	 */
	private static void loadParameterFile(String file) {
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(file)));
			System.out.println("\nReading hydrogen bonding definitions from " + file);
			loadParameterFile(in);
		} catch (IOException ex) {
			throw new RuntimeException("Failed to read " + file + " " + ex);
		}
	}

	/**
	 * Fills in the array of hydrogen bonding types from the default resource
	 * (donor_acceptor_types.txt)
	 * 
	 */
	private static void loadParameterFile() {
		BufferedReader in = new BufferedReader(
				new InputStreamReader(
						HydrogenBondingType.class
								.getResourceAsStream("donor_acceptor_types.txt")));
		loadParameterFile(in);
	}

	/**
	 * Loads in a parameter file and fills the hydrogen bond array. Defintions
	 * are tab-separated lines. Blank lines, lines that begin with whitespace or
	 * # are ignored.
	 * 
	 * For donors: type name probability sln For acceptors: type name geometry
	 * probability sln
	 * 
	 * @param in
	 */
	private static void loadParameterFile(BufferedReader in) {
		List<HydrogenBondingType> parameters = new ArrayList<HydrogenBondingType>();

		try {
			String line = null;

			int no = 0;

			while ((line = in.readLine()) != null) {
				logger.debug("Line: " + line);
				if (line.startsWith("#") || line.startsWith(" ") || line.startsWith("\t"))
					continue;
				line = line.trim();
				if (line.equals(""))
					continue;

				StringTokenizer st = new StringTokenizer(line, "\t");
				String typeStr = st.nextToken().trim();
				int type = 0;
				if (typeStr.equalsIgnoreCase("donor"))
					type = DONOR;
				else if (typeStr.equalsIgnoreCase("acceptor"))
					type = ACCEPTOR;
				else
					throw new RuntimeException("Unknown Donor/Acceptor type " + typeStr);
				String name = st.nextToken().trim();
				double p = Double.valueOf(st.nextToken().trim()).doubleValue();
				AcceptorGeometry geom = null;
				if (type == ACCEPTOR) {
					String geomStr = st.nextToken();
					if (geomStr.equalsIgnoreCase("none"))
						geom = AcceptorGeometry.AG_NONE;
					else if (geomStr.equalsIgnoreCase("dir"))
						geom = AcceptorGeometry.DIR;
					else if (geomStr.equalsIgnoreCase("plane"))
						geom = AcceptorGeometry.PLANE;
					else if (geomStr.equalsIgnoreCase("cone"))
						geom = AcceptorGeometry.CONE;
					else
						throw new RuntimeException("Unknown geometry " + geomStr);
				}
				String sln = st.nextToken().trim();

				Integer weight = null;
				if (st.hasMoreTokens()) {
					weight = Integer.parseInt(st.nextToken());
				}

				HydrogenBondingType daType = new HydrogenBondingType(no, type, name, sln,
						p, geom, weight);

				parameters.add(daType);
				no++;
			}
		} catch (IOException ex) {
			throw new RuntimeException("Failed to read donor acceptor definitions " + ex);
		}

		donorAcceptorTypes = new HydrogenBondingType[parameters.size()];
		for (int i = 0; i < donorAcceptorTypes.length; i++)
			donorAcceptorTypes[i] = parameters.get(i);
	}

	/**
	 * Searches a molecule for all hydrogen bonding types
	 * 
	 * @see #findType(Molecule)
	 * 
	 * @param mol
	 */
	public static void searchMolecule(GaMolecule mol) {
		searchMolecule(mol, mol.getInfoMessageLogger());
	}

	/**
	 * Searches a molecule for all hydrogen bonding types
	 * 
	 * @see #findType(Molecule)
	 * 
	 * @param mol
	 * @param log
	 * @param out
	 */
	public static void searchMolecule(GaMolecule mol, InfoMessageLogger infoMessageLogger) {

		if (donorAcceptorTypes == null)
			throw new RuntimeException("no donor acceptor types defined");
		for (int i = 0; i < donorAcceptorTypes.length; i++) {
			donorAcceptorTypes[i].findType(mol);
		}
		for (Atom atom : mol.getAtoms()) {
			if (atom.getDonorType() != null)
				infoMessageLogger.infoMessageln(3, atom.info() + ":"
						+ atom.getDonorType().info());
			if (atom.getAcceptorType() != null)
				infoMessageLogger.infoMessageln(3, atom.info() + ":"
						+ atom.getAcceptorType().info());

		}
	}

	/**
	 * Removes donor and acceptor references from all atoms in a molecule
	 * 
	 * @param mol
	 */
	public static void freeMolecule(GaMolecule mol) {
		for (Atom atom : mol.getAtoms()) {
			atom.setDonorType(null);
			atom.setAcceptorType(null);
		}
	}

	/**
	 * For testing and debugging. Identifies all hydrogen bonding types in a
	 * structure file.
	 * 
	 * @param args
	 */
	public static void main(String args[]) {
		if (args.length != 1) {
			System.err.println("Usage: HydrogenBondingType <structure_file>");
			System.exit(0);
		}

		HydrogenBondingType.loadParameters();

		List<GaMolecule> molecules = GaMolecule.loadFiles(new String[] { args[0] },
				Molecule.FileType.UNKNOWN, Molecule.Source.FILE);

		for (GaMolecule mol : molecules) {
			System.out.println("Searching " + mol.getName());
			// DonorAcceptorType.findMoleculeDonorAcceptors(mol);
			mol.init();
			HydrogenBondingType.searchMolecule(mol);
		}

	}

	/**
	 * Indentifies the hydrogen bonding type from a name string. The string is
	 * usually a pharmacophore defintion. Enables recreation of features from
	 * pharmacophore strings
	 * 
	 * @param name
	 * @return
	 */
	public static HydrogenBondingType getTypeFromName(String name) {
		name = name.toUpperCase().replace(' ', '_');
		for (int i = 0; i < donorAcceptorTypes.length; i++) {
			HydrogenBondingType type = donorAcceptorTypes[i];
			String test = type.name.toUpperCase().replace(' ', '_');
			if (test.equals(name))
				return type;
		}
		throw new RuntimeException("can't find Mills and Dean type for " + name);
	}

	public AcceptorGeometry getGeometry() {
		return geometry;
	}

	public int getId() {
		return id;
	}

	public String getName() {
		return name;
	}

	public double getProbability() {
		return probability;
	}

	public static HydrogenBondingType[] getDonorAcceptorTypes() {
		return donorAcceptorTypes;
	}

	/**
	 * @return the type
	 */
	public int getType() {
		return type;
	}

	/**
	 * @return the weight
	 */
	public int getWeight() {
		return weight;
	}

	/**
	 * @return the sln
	 */
	public String getSln() {
		return sln;
	}

}
