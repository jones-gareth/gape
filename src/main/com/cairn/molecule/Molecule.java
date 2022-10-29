package com.cairn.molecule;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.io.Writer;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.StringTokenizer;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;
import com.cairn.gape.utils.InfoMessageLogger;
import com.cairn.molecule.Atom.SdStereo;

/**
 * Class for representing molecule structure
 * 
 * @author Gareth Jones
 * 
 */
public class Molecule {

	private static Logger logger = Logger.getLogger(Molecule.class);

	static {
		// logger.setLevel(Level.DEBUG);
	}

	private volatile String name, chargeType;

	// mol2 molecule type
	private volatile String type;

	private volatile int nSubs, nLonePairs;

	private volatile int nOutputAtoms, nOutputBonds;

	private volatile List<Atom> atoms = new ArrayList<>();

	private volatile List<double[]> coords = new ArrayList<>();

	private volatile List<Bond> bonds = new ArrayList<>();

	private volatile List<Ring> rings = new ArrayList<>();

	private volatile List<Angle> angles = new ArrayList<>();;

	private volatile List<Torsion> torsions = new ArrayList<>();

	private volatile List<RotatableBond> rotatableBonds = new ArrayList<>();

	private volatile Sssr sssr;

	protected volatile InfoMessageLogger infoMessageLogger = new InfoMessageLogger();

	public static class MolReadException extends Exception {

		private MolReadException(String message) {
			super(message);
		}

	}

	/**
	 * Can handle MOL2 and SDF files
	 */
	public enum FileType {
		MOL2, SDF, UNKNOWN;

		/**
		 * @return common suffix
		 */
		public String getSuffix() {
			switch (this) {
			case MOL2:
				return "mol2";
			case SDF:
				return "sdf";
			default:
				return null;
			}
		}
	}

	/**
	 * Can load structures from a local file or a URL.
	 */
	public enum Source {
		FILE, URL
	}

	private boolean loaded, structureChanged;

	protected String fileName, baseName;

	protected HashMap<String, String> sdfFields;

	static protected java.text.NumberFormat nf;

	private boolean sdChiral;

	private boolean _initDone = false;

	// Flags to control how molecules are loaded
	protected static final ThreadLocal<Boolean> assignTypesFlag = ThreadLocal
			.withInitial(() -> true);
	protected static final ThreadLocal<Boolean> addHydrogensFlag = ThreadLocal
			.withInitial(() -> true);
	protected static final ThreadLocal<Boolean> solvateFlag = ThreadLocal
			.withInitial(() -> true);
	protected static final ThreadLocal<Boolean> chargeFlag = ThreadLocal
			.withInitial(() -> true);

	// this is for debugging setting atom types
	private final boolean validateAtomTypes = true;

	static {
		nf = java.text.NumberFormat.getInstance();
		nf.setMaximumFractionDigits(3);
		nf.setGroupingUsed(false);
	}

	static final private ThreadLocal<AddHydrogens> addHydrogens = ThreadLocal
			.withInitial(() -> new AddHydrogens());

	static final private ThreadLocal<Solvate> solvate = ThreadLocal
			.withInitial(() -> new Solvate());

	static final private ThreadLocal<Charge> charge = ThreadLocal
			.withInitial(() -> new Charge());

	/**
	 * Empty Constructor
	 * 
	 */
	public Molecule() {
		;
	}

	/**
	 * Create a new molecule which references an existing molecule.
	 * 
	 * @param m
	 */
	public Molecule(Molecule m) {
		createReference(m);
	}

	public Molecule(String name, List<Atom> atoms, List<double[]> coords) {
		this.name = name;
		this.atoms = atoms;
		this.coords = coords;
		for (Atom atom : atoms) {
			atom.setMolecule(this);
		}
	}

	public Molecule(String name, List<Atom> atoms, List<Bond> bonds,
			List<double[]> coords) {
		this.name = name;
		this.atoms = atoms;
		this.bonds = bonds;
		this.coords = coords;
		for (Atom atom : atoms) {
			atom.setMolecule(this);
		}

	}

	/**
	 * Creates a molecule by reading a file. Assumed to be in MOL2 format
	 * 
	 * @param file
	 */
	public Molecule(String file) {
		this(file, FileType.MOL2, Source.FILE);
	}

	/**
	 * Creates a molecule by reading a file or fetch a URL. Specify format and
	 * source.
	 * 
	 * @param file
	 * @param type
	 * @param source
	 */
	public Molecule(String file, FileType type, Source source) {
		load(file, type, source);
	}

	/**
	 * Creates a new molecule by reading from a buffered reader. Specify format.
	 * 
	 * @param in
	 * @param type
	 */
	public Molecule(BufferedReader in, FileType type) {
		load(in, type);
	}

	/**
	 * Returns true if the molecule was loaded successfully.
	 * 
	 * @return
	 */
	public boolean isLoaded() {
		return loaded;
	}

	/**
	 * Guess the structure file type from a filename.
	 * 
	 * @param file
	 * @return
	 */
	public static FileType getType(String file) {
		if (file.toUpperCase().endsWith(".MOL2"))
			return FileType.MOL2;
		else if (file.toUpperCase().endsWith(".SDF"))
			return FileType.SDF;
		if (file.toUpperCase().endsWith(".MOL2.GZ"))
			return FileType.MOL2;
		else if (file.toUpperCase().endsWith(".SDF.GZ"))
			return FileType.SDF;

		else {
			System.out.println("Don't know filetype of " + file);
			return FileType.UNKNOWN;
		}
	}

	/**
	 * Loads a single molecule from a location.
	 * 
	 * @param location
	 * @param type
	 * @param source
	 */
	public synchronized void load(String location, FileType type, Source source) {
		BufferedReader in = openReader(location, source);
		if (type == FileType.UNKNOWN)
			type = getType(location);
		load(in, type);

		try {
			in.close();
		} catch (IOException ex) {
			System.out.println("IO error closing " + ex);
		}
	}

	/**
	 * Opens a BufferedReader to structural data.
	 * 
	 * @param location
	 * @param source
	 * @return
	 */
	public static BufferedReader openReader(String location, Source source) {
		BufferedReader in = null;
		if (source == Source.FILE) {
			try {
				if (location.toUpperCase().endsWith(".GZ")) {
					in = new BufferedReader(
							new InputStreamReader(new FileInputStream(location)));
				} else
					in = new BufferedReader(new FileReader(location));
			} catch (IOException ex) {
				if (in != null) {
					IOUtils.closeQuietly(in);
				}
				logger.error("File " + location + " not found or IO exception", ex);
			}
		} else if (source == Source.URL) {
			try {
				URL purl = new URL(location);
				URLConnection URLcon = purl.openConnection();
				URLcon.setDoInput(true);
				URLcon.setUseCaches(true);
				URLcon.setAllowUserInteraction(true);
				InputStream urlIn = URLcon.getInputStream();
				if (location.toUpperCase().endsWith(".GZ")) {
					GZIPInputStream zipStream = new GZIPInputStream(urlIn);
					in = new BufferedReader(new InputStreamReader(zipStream));
				} else
					in = new BufferedReader(
							new InputStreamReader(URLcon.getInputStream()));
			} catch (MalformedURLException urlex) {
				logger.error("Malformed url: " + location, urlex);
			} catch (IOException ioex) {
				logger.error("IO Exception", ioex);
			}
		}
		return in;
	}

	/**
	 * Determines the centroid of the molecule
	 * 
	 * @param heavy
	 *            set true to use heavy atoms only.
	 * @return
	 */
	public double[] calculateCentroid(boolean heavy) {
		double x = 0, y = 0, z = 0, cnt = 0;

		int nAtoms = atoms.size();
		for (int no = 0; no < nAtoms; no++) {
			Atom atom = atoms.get(no);
			if (atom.getType().isNotDummy()) {
				if (heavy && !atom.getType().isHeavy()) {
					continue;
				}
			}
			x += coords.get(no)[0];
			y += coords.get(no)[1];
			z += coords.get(no)[2];
			cnt += 1;
		}
		if (cnt > 0) {
			x /= cnt;
			y /= cnt;
			z /= cnt;
		}
		return new double[] { x, y, z };
	}

	/**
	 * Loads a single Molecule from a local MOL2 file
	 * 
	 * @param f
	 */
	public synchronized void loadFile(String f) {
		fileName = f;
		String b = (new File(f)).getName();
		baseName = b.substring(0, b.length() - 5);
		if (!fileName.toUpperCase().endsWith(".MOL2")) {
			System.out.println("Failed to load " + fileName + " not a mol2 file");
			;
		}
		// System.out.println("loading " + fileName);
		BufferedReader in = openReader(fileName, Source.FILE);
		if (in == null) {
			System.out.println("Failed to open " + fileName);
			return;
		}
		load(in, FileType.MOL2);
		if (name.equals("****") || name.equals(""))
			name = baseName;
	}

	/**
	 * Given a reader loads in a MOL2 or SDF molecule. Choose to set atom types.
	 * 
	 * @param in
	 * @param type
	 */
	public synchronized void load(BufferedReader in, FileType type) {
		try {
			if (type == FileType.MOL2) {
				loadSybylMol2(in);
			} else if (type == FileType.SDF) {
				loadSdfMol(in);
			}
		} catch (IOException ex) {
			logger.warn("IO error reading " + ex);
		} catch (MolReadException ex) {
			logger.warn("Format error or end of stream reading " + ex);
		}
		init();
	}

	/**
	 * Finds ring systems SSSR, atom neighbour lists. Sets atom types. Depending
	 * on flags it will also add hydrogens, solvate and determine common charge
	 * groups.
	 * 
	 */
	public synchronized void init() {

		if (_initDone)
			throw new IllegalStateException("Molecule already initialized");
		_initDone = true;

		build();

		// add Hydrogens
		if (addHydrogensFlag.get()) {
			addHydrogens();
		}

		// solvate
		if (solvateFlag.get()) {
			solvate();
		}

		// find charged groups
		if (chargeFlag.get()) {
			charge();
		}

	}

	/**
	 * Solvates molecule by identifying and charging common acids and bases.
	 * 
	 * @see Solvate
	 */
	public synchronized void solvate() {
		if (!assignTypesFlag.get()) {
			assignAtomTypes();
		}

		solvate.get().setMolecule(this);
		solvate.get().solvate();
	}

	/**
	 * Charges molecule by searching for common charge groups.
	 * 
	 * @see Charge
	 */
	public synchronized void charge() {
		charge.get().setMolecule(this);
		charge.get().charge();
	}

	/**
	 * Fills valance by adding hydrogens to this molecule.
	 * 
	 * @see AddHydrogens
	 */
	public synchronized void addHydrogens() {
		addHydrogens.get().setMolecule(this);
		addHydrogens.get().addHydrogens();

		assignAtomTypes();
	}

	/**
	 * Finds ring systems SSSR, atom neighbour lists. Sets atom types. Called
	 * when loading the molecule and also if the strucuture has been altered by
	 * functions like addBond, addAtom, deleteAtom. Sets atom types.
	 * 
	 */
	public synchronized void build() {
		logger.debug("assigning atom neighbours");
		assignAtomNeighbours();
		logger.debug("Finding rings");
		RingFinder rf = new RingFinder(this);
		rf.findRings();
		logger.debug("Doing SSSR");
		sssr = new Sssr(this);
		sssr.doSssr();

		if (assignTypesFlag.get()) {
			assignAtomTypes();
		} else {
			// do aromatic perception and planar ring perception even if we
			// don't check types.
			percieveAromaticity();
		}
		logger.debug("finding groups");
		findGroups();
	}

	/**
	 * Rebuilds neighbour lists, SSSR etc (cals build), if the structure has
	 * been altered by functions like addBond, addAtom, deleteAtom. Sets atom
	 * types.
	 * 
	 * @param setTypes
	 * @see #init(boolean)
	 */
	public synchronized void update() {
		if (structureChanged) {
			build();
			structureChanged = false;
		}
	}

	/**
	 * Frees up data structures- atom neighbour lists and rings and forcefield
	 * structures. If we're keeping lots of molecules about it may make sense to
	 * free up extra structures then re-create them when required/
	 */
	public void free() {

	}

	/**
	 * Loads a molecule from an SD structure
	 * 
	 * @param in
	 * @throws IOException
	 */
	public synchronized void loadSdfMol(BufferedReader in)
			throws IOException, MolReadException {

		String line = readLine(in).trim();
		if (line.equals("$$$$")) {
			line = readLine(in).trim();
		}
		name = line;
		line = readLine(in);
		line = readLine(in);
		line = readLine(in);
		int nAtoms = Integer.valueOf(line.substring(0, 3).trim()).intValue();
		int nBonds = Integer.valueOf(line.substring(3, 6).trim()).intValue();
		if (line.length() > 15) {
			String chiralStr = line.substring(12, 15).trim();
			if (StringUtils.isNotEmpty(chiralStr)) {
				int chiral = Integer.valueOf(chiralStr).intValue();
				if (chiral == 1)
					sdChiral = true;
			}
		}
		atoms.clear();
		coords.clear();
		int no = 0;
		for (int i = 0; i < nAtoms; i++) {
			line = readLine(in);

			double x = Double.valueOf(line.substring(0, 10).trim()).doubleValue();
			double y = Double.valueOf(line.substring(10, 20).trim()).doubleValue();
			double z = Double.valueOf(line.substring(20, 30).trim()).doubleValue();
			String element = line.substring(31, 34).trim();

			Atom.SdStereo sdStereo = Atom.SdStereo.NONE;

			if (line.length() > 42) {
				String atomStereoParityString = line.substring(39, 42).trim();
				if (StringUtils.isNotEmpty(atomStereoParityString)) {
					int atomStereoParity = Integer.valueOf(atomStereoParityString);

					if (atomStereoParity == 1)
						sdStereo = Atom.SdStereo.ODD;
					else if (atomStereoParity == 2)
						sdStereo = Atom.SdStereo.EVEN;
					else if (atomStereoParity == 3)
						sdStereo = Atom.SdStereo.EITHER;
				}
			}

			// ingnore hydrogen count query only
			// int hydrogenCount = Integer.valueOf(line.substring(42,
			// 45).trim());
			// ignore stereo care -query only.
			// int stereoCareBox = Integer.valueOf(line.substring(45,
			// 48).trim());

			int massDifference = 0;
			if (line.length() > 36) {
				String massDifferenceStr = line.substring(34, 36).trim();
				if (StringUtils.isNotEmpty(massDifferenceStr))
					massDifference = Integer.valueOf(massDifferenceStr);
			}

			int formalChargeVal = 0;
			if (line.length() > 39) {
				String formalChargeString = line.substring(36, 39).trim();
				if (StringUtils.isNotEmpty(formalChargeString))
					formalChargeVal = Integer.valueOf(formalChargeString);
			}

			int valence = 0;
			if (line.length() > 51) {
				String valenceString = line.substring(48, 51).trim();
				if (StringUtils.isNotEmpty(valenceString))
					valence = Integer.valueOf(valenceString);
			}

			logger.debug("x " + x + " y " + y + " z " + z + " atom " + element);
			Atom atom = new Atom(this, no, element, element);
			atoms.add(atom);
			coords.add(new double[] { x, y, z, 1.0 });

			atom.setSdMassDifference(massDifference);
			if (formalChargeVal > 0) {
				atom.setFormalCharge(4 - formalChargeVal);
			}
			if (valence > 0) {
				atom.setSdValence(valence == 15 ? 0 : valence);
			}
			atom.setSdStereo(sdStereo);

			no++;
		}
		bonds.clear();
		no = 0;
		for (int i = 0; i < nBonds; i++) {
			line = readLine(in);
			logger.debug(line);
			int a1 = Integer.valueOf(line.substring(0, 3).trim()).intValue() - 1;
			int a2 = Integer.valueOf(line.substring(3, 6).trim()).intValue() - 1;
			int t = Integer.valueOf(line.substring(6, 9).trim()).intValue();
			Atom atom1 = atoms.get(a1);
			Atom atom2 = atoms.get(a2);
			logger.debug("a1 " + a1 + " a2 " + a2 + " t " + t);
			Bond bond = new Bond(no, atom1, atom2, t);
			bonds.add(bond);
			if (line.length() > 9) {
				int end = line.length() - 1;
				if (end > 12)
					end = 12;
				try {
					int s = Integer.valueOf(line.substring(9, end).trim()).intValue();
					if (s == 1)
						bond.setStereo(Bond.Stereo.UP);
					else if (s == 4)
						bond.setStereo(Bond.Stereo.EITHER);
					else if (s == 6)
						bond.setStereo(Bond.Stereo.DOWN);
				} catch (NumberFormatException ex) {
				}
			}
			no++;
		}
		while (true) {
			line = in.readLine();
			if (line == null) {
				return;
			}
			if (line.startsWith("M  END"))
				break;
			if (line.startsWith("$$$$"))
				throw new MolReadException("Can't find <M END> before end of record");
		}
		relabelAtoms();
		loaded = true;
		logger.debug("Loaded molecule");
		// Chemdraw sd files miss terminating $$$$
		try {
			while (true) {
				line = readLine(in);
				if (line.startsWith("$$$$"))
					break;
				if (line.startsWith("> ")) {
					int s = line.indexOf('<');
					int e = line.indexOf('>', s + 1);
					String field = line.substring(s + 1, e);
					String value = null;
					while (line != null) {
						line = readLine(in).trim();
						if (line.equals(""))
							break;
						if (value == null)
							value = line;
						else
							value += "\n" + line;
					}
					logger.debug("SD Field " + field + " SD value " + value);
					if (value != null) {
						addSdfField(field, value);
					}
				}
			}
		} catch (IOException ex) {
			;
		} catch (MolReadException ex) {
			;
		}
	}

	/**
	 * Adds a SDF field.
	 * 
	 * @param field
	 *            Field name
	 * @param value
	 *            Field value
	 */
	public void addSdfField(String field, String value) {
		if (sdfFields == null)
			sdfFields = new HashMap<String, String>();
		sdfFields.put(field, value);
	}

	/**
	 * Returns true if this SDF field is present
	 * 
	 * @param field
	 * @return
	 */
	public boolean hasSdfField(String field) {
		if (sdfFields == null)
			return false;
		return sdfFields.containsKey(field);
	}

	/**
	 * Return's the value of an SDF field (or null if the field does not exist)
	 * 
	 * @param field
	 * @return
	 */
	public String getSdfField(String field) {
		if (sdfFields == null)
			return null;
		if (sdfFields.containsKey(field))
			return sdfFields.get(field);
		else
			return null;
	}

	/**
	 * Removes an SDF field
	 * 
	 * @param field
	 */
	public void removeSdfField(String field) {
		if (sdfFields != null)
			sdfFields.remove(field);
	}

	/**
	 * Remove all sdf fields
	 */
	public void clearSdfFields() {
		if (sdfFields != null)
			sdfFields.clear();
	}

	/**
	 * Reads a line in a structure file
	 * 
	 * @param in
	 * @return
	 * @throws IOException
	 */
	private String readLine(BufferedReader in) throws IOException, MolReadException {
		String val = in.readLine();
		if (val == null)
			throw new MolReadException("Read Error");
		return val;
	}

	/**
	 * Reads in a structure from MOL2 format.
	 * 
	 * @param in
	 * @throws IOException
	 */
	public synchronized void loadSybylMol2(BufferedReader in)
			throws IOException, MolReadException {
		StringTokenizer st;
		String line;
		getMol2Section(in, "@<TRIPOS>MOLECULE");
		name = in.readLine().trim();
		line = in.readLine();
		st = new StringTokenizer(line);
		int nAtoms = Integer.valueOf(st.nextToken()).intValue();
		int nBonds = Integer.valueOf(st.nextToken()).intValue();
		type = in.readLine().trim();
		chargeType = in.readLine().trim();
		getMol2Section(in, "@<TRIPOS>ATOM");
		atoms.clear();
		coords.clear();
		for (int i = 0; i < nAtoms; i++) {
			st = new StringTokenizer(in.readLine());
			int no = Integer.valueOf(st.nextToken()).intValue() - 1;
			assert no == i;
			String l = st.nextToken();
			double x = Double.valueOf(st.nextToken()).doubleValue();
			double y = Double.valueOf(st.nextToken()).doubleValue();
			double z = Double.valueOf(st.nextToken()).doubleValue();
			String n = st.nextToken();
			Atom atom = new Atom(this, no, l, n);
			atoms.add(atom);
			coords.add(new double[] { x, y, z, 1.0 });
			if (st.hasMoreTokens()) {
				int subNo = Integer.valueOf(st.nextToken()).intValue();
				String subName = st.nextToken();
				atom.setSubNo(subNo);
				atom.setSubName(subName);
			}
		}
		getMol2Section(in, "@<TRIPOS>BOND");
		bonds.clear();
		for (int i = 0; i < nBonds; i++) {
			st = new StringTokenizer(in.readLine());
			int no = Integer.valueOf(st.nextToken()).intValue() - 1;
			assert i == no;
			int a1 = Integer.valueOf(st.nextToken()).intValue() - 1;
			int a2 = Integer.valueOf(st.nextToken()).intValue() - 1;
			String n = st.nextToken();
			Atom atom1 = atoms.get(a1);
			Atom atom2 = atoms.get(a2);
			bonds.add(new Bond(no, atom1, atom2, n));
		}

		loaded = true;
	}

	/**
	 * Reads ahead to find a section within a MOL2 file
	 * 
	 * @param in
	 * @param section
	 * @throws IOException
	 */
	private void getMol2Section(BufferedReader in, String section)
			throws IOException, MolReadException {
		while (true) {
			String val = in.readLine();
			if (val == null)
				throw new MolReadException("Can't find section " + section);
			if (val.startsWith(section))
				return;
		}
	}

	/**
	 * Writes out the structure to a file in MOL2 format
	 * 
	 * @param file
	 * @param comment
	 */
	public void writeSybylMol2File(String file, String comment) {
		try {
			file = file.replace(' ', '_');
			FileWriter out = new FileWriter(new File(file));
			writeSybylMol2(out, comment);
			out.close();
		} catch (IOException ex) {
			System.out.println("writeMol2File IO exception: " + ex);
		}
	}

	/**
	 * Returns a String of the structure in MOL2 format.
	 * 
	 * @return
	 */
	public String sybylMol2String() {
		StringWriter writer = new StringWriter();
		try {
			writeSybylMol2(writer, null);
		} catch (IOException ex) {
			;
		}
		return writer.toString();
	}

	/**
	 * Writes out the structure in MOL2 format (Lone pairs are included)
	 * 
	 * @param out
	 * @param info
	 * @throws IOException
	 */
	public void writeSybylMol2(Writer out, String info) throws IOException {
		writeSybylMol2(out, info, true);
	}

	/**
	 * Sets lone pairs for output.
	 * 
	 * @param outputLP
	 *            set to true to include lone pairs in output.
	 */
	private void outputLonePairs(boolean outputLP) {
		for (Atom atom : atoms) {
			if (atom.getAtomType() == AtomType.Type.LP)
				atom.setOutput(outputLP);
		}
	}

	/**
	 * Renumber atoms and bond for writing out. You can have an atom or bond in
	 * the connection table and not write out atoms of bonds by setting output
	 * to false. Sets atom and bond numbers to -1 where output is false and
	 * makes numbers of other atoms and bonds consistent.
	 */
	private synchronized void renumberForOutput() {

		// TODO alter atom and bond so we have an index and a number.
		int cnt = 0;
		for (Atom atom : atoms) {
			if (!atom.isOutput())
				atom.setNo(-1);
			else if (atom.getNo() == -1)
				continue;
			else
				atom.setNo(cnt++);
		}
		nOutputAtoms = cnt;
		cnt = 0;
		for (Bond bond : bonds) {
			if (!bond.getAtom1().isOutput() || !bond.getAtom2().isOutput())
				bond.setNo(-1);
			else if (!bond.isOutput())
				continue;
			else
				bond.setNo(cnt++);

		}
		nOutputBonds = cnt;
	}

	/**
	 * Restores correct numbering of atoms and bonds.
	 * 
	 * @see #renumberForNoLonePairs()
	 */
	private synchronized void renumber() {
		int no = 0;
		for (Atom atom : atoms) {
			atom.setNo(no++);
		}
		no = 0;
		for (Bond bond : bonds) {
			bond.setNo(no++);
		}
	}

	/**
	 * Writes out the structure in MOL2 format
	 * 
	 * @param out
	 * @param info
	 * @param incLP
	 * @throws IOException
	 */
	public void writeSybylMol2(Writer out, String info, boolean incLP)
			throws IOException {

		// renumber for no lone pairs
		outputLonePairs(incLP);
		renumberForOutput();

		// edit to count allowed
		int nA = nOutputAtoms;
		int nB = nOutputBonds;

		// Date not supported by gcj on win32

		String date = java.text.DateFormat.getDateTimeInstance().format(new Date());
		out.write("#       Name:                   " + name + "\n");
		out.write("#       Creating User Name:     <unknown>\n");
		out.write("#       Creation Time:          " + date + "\n");
		out.write("#       \n");
		if (info != null)
			out.write("#       " + info + "\n");
		else
			out.write("#\n");
		out.write("\n\n");

		out.write("@<TRIPOS>MOLECULE\n");
		out.write(name + "\n");
		out.write("   " + nA + "    " + nB + "    0\n");
		out.write("SMALL\nUSER_CHARGES\n\n\n");

		out.write("@<TRIPOS>ATOM\n");
		for (int i = 0; i < atoms.size(); i++) {
			Atom atom = atoms.get(i);
			int no = atom.getNo() + 1;
			if (no == 0)
				continue;
			double[] coord = coords.get(i);
			out.write("    " + no + "   " + atom.getLabel() + "    " + nf.format(coord[0])
					+ "    " + nf.format(coord[1]) + "    " + nf.format(coord[2]) + "    "
					+ atom.getType().getName());
			if (atom.getSubName() != null)
				out.write("    " + atom.getSubNo() + "    " + atom.getSubName());
			else
				out.write("     1     MOLECULE");
			out.write("    " + atom.getPartialCharge());
			out.write("  \n");
		}

		out.write("@<TRIPOS>BOND\n");
		for (int i = 0; i < bonds.size(); i++) {
			Bond bond = bonds.get(i);
			int no = bond.getNo() + 1;
			if (no == 0)
				continue;
			int a1 = bond.getAtom1().getNo() + 1;
			int a2 = bond.getAtom2().getNo() + 1;
			out.write("   " + no + "    " + a1 + "   " + a2 + "   "
					+ bond.getType().getName() + "   \n");
		}

		// restore atom numbering
		renumber();

	}

	/**
	 * Creates the neighbour lists for all atoms
	 * 
	 * @see #getNeighbours(Atom)
	 */
	public void assignAtomNeighbours() {
		for (Atom atom : atoms) {
			atom.getNeighbours();
		}
	}

	/**
	 * Return true if the two atoms are not dummy and bonded
	 * 
	 * @param a1
	 * @param a2
	 * @return
	 */
	public boolean bonded(Atom a1, Atom a2) {
		if (a1 == a2)
			return false;
		return a1.getNotDummyNeighbours().stream().anyMatch(x -> a2 == x);
	}

	/**
	 * Return true if the two atoms are 1-3 bonded
	 * 
	 * @param a1
	 * @param a2
	 * @return
	 */
	public boolean bonded13(Atom a1, Atom a2) {
		for (Atom neighbour1 : a1.getNotDummyNeighbours()) {
			for (Atom neighbour2 : a2.getNotDummyNeighbours()) {
				if (neighbour1 == neighbour2)
					return true;
			}
		}
		return false;
	}

	/**
	 * Return true if the two atoms are 1-4 bonded
	 * 
	 * @param a1
	 * @param a2
	 * @return
	 */
	public boolean bonded14(Atom a1, Atom a2) {
		for (Atom neighbour1 : a1.getNotDummyNeighbours()) {
			for (Atom neighbour2 : a2.getNotDummyNeighbours()) {
				if (bonded(neighbour1, neighbour2))
					return true;
			}
		}
		return false;
	}

	/**
	 * Relabels the atoms to C1, C2, H1, H2 etc. The ordering is the same as the
	 * input file.
	 */
	public void relabelAtoms() {
		Map<AtomType, List<Atom>> elementalMap = atoms.stream().collect(
				Collectors.groupingBy(atom -> atom.getType().getElementalType()));
		for (Entry<AtomType, List<Atom>> entry : elementalMap.entrySet()) {
			int no = 1;
			for (Atom atom : entry.getValue()) {
				atom.setLabel(entry.getKey().getName() + no++);
			}
		}
	}

	/**
	 * Finds all torsions in a molecule
	 * 
	 */
	synchronized void assignTorsions() {
		torsions.clear();
		for (Bond bond : bonds) {
			Atom a2 = bond.getAtom1();
			Atom a3 = bond.getAtom2();
			for (Atom a1 : a2.getNotDummyNeighbours()) {
				if (a1 == a3)
					continue;
				for (Atom a4 : a3.getNotDummyNeighbours()) {
					if (a4 == a2)
						continue;
					logger.debug("Added Torsion");
					Torsion t = new Torsion(this, bond, a1, a2, a3, a4);
					torsions.add(t);
				}
			}
		}

		logger.debug(torsions.size() + " torsions");
	}

	/**
	 * Finds rotatable bonds in a molecule
	 * 
	 */
	public synchronized void assignRotatableBonds() {
		assignRotatableBonds(false);
	}

	/**
	 * Finds rotatable bonds in a molecule.
	 * 
	 * @param flipAmideBonds
	 *            set to allow amide bonds to flip between trans and cis
	 */
	public synchronized void assignRotatableBonds(boolean flipAmideBonds) {
		rotatableBonds.clear();
		// Assign torsions before rotatable bonds
		if (torsions == null)
			assignTorsions();
		for (Bond bond : bonds) {
			Bond.Rotatable rotateType = bond.isRotatable(flipAmideBonds);
			if (rotateType == Bond.Rotatable.NO)
				continue;

			List<Torsion> bondTorsions = torsions.stream()
					.filter(t -> t.getBond() == bond).collect(Collectors.toList());

			boolean flip = rotateType == Bond.Rotatable.FLIP ? true : false;
			RotatableBond rb = new RotatableBond(bond, this, bondTorsions, flip);
			infoMessageLogger.infoMessageln(3, rb.info());

			rotatableBonds.add(rb);
		}

		infoMessageLogger.infoMessageln(2, rotatableBonds.size() + " Rotatable Bonds");

	}

	/**
	 * Loops through the bond table to create an array of angles within the
	 * molecule.
	 */
	synchronized void assignAngles() {
		angles.clear();
		int nBonds = getnBonds();
		for (int i = 0; i < nBonds; i++)
			for (int j = i + 1; j < nBonds; j++) {
				Bond b1 = bonds.get(i);
				Bond b2 = bonds.get(j);
				Atom mid = null, a1 = null, a2 = null;
				if (b1.getAtom1() == b2.getAtom1()) {
					mid = b1.getAtom1();
					a1 = b1.getAtom2();
					a2 = b2.getAtom2();
				}
				if (b1.getAtom2() == b2.getAtom2()) {
					mid = b1.getAtom2();
					a1 = b1.getAtom1();
					a2 = b2.getAtom1();
				}
				if (b1.getAtom1() == b2.getAtom2()) {
					mid = b1.getAtom1();
					a1 = b1.getAtom2();
					a2 = b2.getAtom1();
				}
				if (b1.getAtom2() == b2.getAtom1()) {
					mid = b1.getAtom2();
					a1 = b1.getAtom1();
					a2 = b2.getAtom2();
				}
				if (mid == null)
					continue;
				Angle a = new Angle(a1, mid, a2);
				angles.add(a);
				logger.debug("Adding Angle");
			}
	}

	/**
	 * Looks through the bond table and adds atoms cyclicly bonded to atom into
	 * the neighbours array. Returns count of neighbours found.
	 * 
	 * @param atom
	 * @param neighbours
	 * @return
	 */
	int getRingNeighbours(Atom atom, Atom[] neighbours) {
		int count = 0;
		for (Bond bond : getBonds()) {
			if (!bond.isInRing())
				continue;
			if (atom == bond.getAtom2()) {
				neighbours[count] = bond.getAtom1();
				count++;
			}
			if (atom == bond.getAtom1()) {
				neighbours[count] = bond.getAtom2();
				count++;
			}
		}
		return count;
	}

	/**
	 * Returns the bond between two atoms (or null if no bond is present).
	 * 
	 * @param a1
	 * @param a2
	 * @return
	 */
	public Bond getBond(Atom a1, Atom a2) {
		Optional<Bond> bond = bonds.stream()
				.filter(b -> (b.getAtom1() == a1 && b.getAtom2() == a2)
						|| (b.getAtom1() == a2 && b.getAtom2() == a1))
				.findAny();
		return bond.orElse(null);
	}

	/**
	 * This only removes Lone Pairs we've added- it assumes they're last!
	 */
	public synchronized void removeLonePairs() {
		boolean changed = false;
		int no = 0;
		for (Atom atom : atoms) {
			atom.setnLonePairs(0);
			atom.setLonePairs(null);
			if (atom.getAtomType() == AtomType.Type.LP) {
				changed = true;
				break;
			}
			no++;
		}
		if (no != atoms.size()) {
			atoms = new ArrayList<>(atoms.subList(0, no));
			coords = new ArrayList<>(coords.subList(0, no));
		}

		no = 0;
		for (Bond bond : bonds) {
			if (bond.getAtom1().getAtomType() == AtomType.Type.LP
					|| bond.getAtom2().getAtomType() == AtomType.Type.LP) {
				changed = true;
				break;
			}
			no++;

		}
		if (no != bonds.size()) {
			bonds = new ArrayList<>(bonds.subList(0, no));
		}

		nLonePairs = 0;
		if (changed)
			structureChanged = true;
	}

	/**
	 * Finds all the atoms on one side (defined by atom) of the bond. The path
	 * stops whenever atoms in stopList are encountered.
	 * 
	 * @param atom
	 * @param bond
	 * @param stopList
	 * @return
	 */
	public List<Atom> findNodes(Atom atom, Bond bond, List<Atom> stopList) {
		List<Atom> v = new ArrayList<Atom>();
		Atom otherAtom = (bond.getAtom1() == atom) ? bond.getAtom2() : bond.getAtom1();
		findNodesArrayList(atom, otherAtom, atom, stopList, v);
		return v;
	}

	/**
	 * Finds all the atoms on one side (defined by atom) of the bond.
	 * 
	 * @param atom
	 * @param bond
	 * @return
	 */
	public List<Atom> findNodes(Atom atom, Bond bond) {
		return findNodes(atom, bond, null);
	}

	/**
	 * Recursive function to implement findNodes. Tracks starting, ending,
	 * current atom and arraylist of atoms. If the current atom is not the
	 * starting or ending atom, in the stoplist or not already in the list it is
	 * added to the list and it's neighbours checked.
	 * 
	 * @param start
	 * @param end
	 * @param current
	 * @param stopList
	 * @param v
	 */
	private void findNodesArrayList(Atom start, Atom end, Atom current,
			List<Atom> stopList, List<Atom> v) {
		if (v.contains(current)) {
			return;
		}
		if (stopList != null && stopList.contains(current)) {
			return;
		}

		if (current != start)
			v.add(current);
		for (Bond bond : getBonds()) {
			Atom check = null;
			if (bond.getAtom1() == current)
				check = bond.getAtom2();
			else if (bond.getAtom2() == current)
				check = bond.getAtom1();
			else
				continue;
			if (check == end)
				continue;
			findNodesArrayList(start, end, check, stopList, v);
		}
	}

	/**
	 * Prints out information about which fragments appear in atoms in the
	 * molecule. Fragments are structures such as Benzene, Ribose.
	 * 
	 * @see Fragments
	 */
	public void printFragments() {
		atoms.stream()
				.forEachOrdered(a -> infoMessageLogger.infoMessageln(a.fragmentsInfo()));
	}

	/**
	 * Returns an array of all heavy atoms in this molecule.
	 * 
	 * @return
	 */
	public List<Atom> getHeavyAtoms() {
		return atoms.stream().filter(a -> a.isHeavy()).collect(Collectors.toList());

	}

	/**
	 * Initializes data from another molecule. The structures from the other
	 * molecule are references so changes to the other molecule are reflected in
	 * this one (and vice-versa).
	 * 
	 * @param m
	 */
	public synchronized void createReference(Molecule m) {
		name = m.name;
		type = m.type;
		chargeType = m.chargeType;
		nSubs = m.nSubs;
		atoms = m.atoms;
		coords = m.coords;
		bonds = m.bonds;
		rings = m.rings;
		angles = m.angles;
		torsions = m.torsions;
		rotatableBonds = m.rotatableBonds;
		loaded = true;
		fileName = m.fileName;
		baseName = m.baseName;

	}

	/**
	 * @see getBond
	 * @param a1
	 * @param a2
	 * @return bond between the two atoms.
	 */
	public Bond findBond(Atom a1, Atom a2) {
		return getBond(a1, a2);
	}

	/**
	 * Flattens bonds such as amide bonds and bonds round planar nitrogens. Call
	 * this after finding Rotatable bonds (bonds which can be flattened are
	 * identified then)
	 * 
	 * @see RotatableBond#flattenBond(Molecule, Atom, Atom)
	 */
	public synchronized void flattenBonds() {
		bonds.stream().forEach(bond -> {
			if (bond.isCanFlatten()) {
				infoMessageLogger.infoMessageln(4, bond.info());
				logger.debug("can flatten true for " + bond.info());
				RotatableBond.flattenBond(this, bond.getAtom1(), bond.getAtom2());
			}
		});
	}

	/**
	 * Finds common groups (COOH planar NH2, NO2, Npl3 link) and marks atoms.
	 * Call this after finding rings and setting atom types
	 * 
	 * @see Atom#findGroups()
	 */
	public synchronized void findGroups() {
		atoms.stream().forEach(atom -> atom.findGroups());
	}

	/**
	 * Adds an atom to this molecule. Sets structureChanged- the neighbour lists
	 * for atoms will not be complete until update is called.
	 * 
	 * @param atom
	 * @param coord
	 */
	public synchronized void addAtom(Atom atom, double coord[]) {
		atom.setNo(getnAtoms());
		atoms.add(atom);
		coords.add(coord);
		structureChanged = true;
	}

	/**
	 * Adds a new bond to this molecule. Sets Sets structureChanged- the
	 * neighbour lists for atoms will not be complete until update is called.
	 * 
	 * @param bond
	 */
	public synchronized void addBond(Bond bond) {
		bond.setNo(getnBonds());
		;
		bonds.add(bond);
		structureChanged = true;
	}

	/**
	 * Checks each atom and bond in the molecule and sets it to the correct
	 * Tripos atom or bondtype. Does aromatic ring perception.
	 * 
	 * @see Atom#checkAtomType()
	 * @see Bond#checkBondType()
	 */
	public synchronized void assignAtomTypes() {

		// set atoms
		for (Atom atom : atoms) {
			AtomType.Type check = atom.checkAtomType();
			if (check != atom.getAtomType()) {
				if (check == null) {
					infoMessageLogger.infoMessageln(4,
							"Can't guess type for atom " + atom.info());
					logger.warn("Can't guess type for atom");
				} else {
					AtomType newType = AtomType.sybType(check);
					infoMessageLogger.infoMessageln(4, "Changing atom " + atom.info()
							+ " to type " + newType.getName());
					atom.setType(newType);
				}
			}
		}

		// can use validate atom types to check that atom types are correctly
		// and consistently assigned
		if (validateAtomTypes) {
			for (Atom atom : atoms) {
				AtomType.Type check = atom.checkAtomType();
				if (check != atom.getAtomType()) {
					if (check == null) {
						infoMessageLogger.infoMessageln(2,
								"Can't guess type for atom " + atom.info());
						logger.warn("Can't guess type for atom", new Throwable());
					} else {
						throw new RuntimeException(
								"Failed to validate atom type for atom " + atom.info()
										+ " first pass " + atom.getAtomType()
										+ " second pass " + check);
					}
				}
			}
		}

		// check bonds
		for (Bond bond : bonds) {
			BondType.Type check = bond.checkBondType();
			if (check != bond.getBondType()) {
				if (check == null) {
					infoMessageLogger.infoMessageln(4,
							"Can't guess type for bond " + bond.info());
				} else {
					BondType newType = BondType.sybType(check);
					infoMessageLogger.infoMessageln(4,
							"Changing " + bond.info() + " to type " + newType.getName());
					bond.setType(newType);
					// if a bond type changed we need to rebuild the neighbor
					// lists
					structureChanged = true;
				}
			}
		}

		percieveAromaticity();

		// The aromatic ring atom and bond types are set after perception-
		// setting during perception will disrupt perception of fused aromatic
		// rings.

		rings.stream().forEach(ring -> {
			ring.setAromaticTypes();
		});

		// Try and sort out charged acid groups
		atoms.stream().forEach(atom -> atom.setChargedAcidGroups());

	}

	/**
	 * Perceive aromatic rings (does not set aromatic ring bond types)
	 */
	private void percieveAromaticity() {
		// Find aromatic rings and planar rings
		rings.stream().forEach(ring -> {
			infoMessageLogger.infoMessageln(4, "Checking ring " + ring.info());
			ring.percieveSp2();
		});

		// determine aromaticity
		Huckel.huckel(this);
	}

	/**
	 * Deletes bond no from the molecule. Sets structureChanged. Atom neighbour
	 * lists will be inconsistent until the structure is updated.
	 * 
	 * @param no
	 */
	public synchronized void deleteBond(Bond bond) {
		bonds.remove(bond.getNo());
		int idx = 0;
		for (Bond b : bonds) {
			b.setNo(idx++);
		}
		structureChanged = true;
	}

	/**
	 * Deletes atom number no from a molecule. Rebuilds the bond table and
	 * renumbers atoms. If init is true, finds rings and rebuild atom neighbour
	 * lists. Otherwise structureChanged is set. In this case atom neighbour
	 * lists will be inconsistent until the structure is updated.
	 * 
	 * @param no
	 * @param init
	 */
	public void deleteAtom(Atom atom, boolean rebuild) {
		atoms.remove(atom.getNo());
		coords.remove(atom.getNo());

		List<Bond> checkBonds = new ArrayList<>(bonds);
		checkBonds.stream()
				.filter(bond -> bond.getAtom1() == atom || bond.getAtom2() == atom)
				.forEach(bond -> bonds.remove(bond));

		int idx = 0;
		for (Bond b : bonds) {
			b.setNo(idx++);
		}
		idx = 0;
		for (Atom a : atoms) {
			a.setNo(idx++);
		}

		if (rebuild)
			build();
		else
			structureChanged = true;
	}

	/**
	 * Deletes atom number no from a molecule. Rebuilds the bond table and
	 * renumbers atoms.
	 * 
	 * @param no
	 */
	public void deleteAtom(Atom atom) {
		deleteAtom(atom, true);
	}

	/**
	 * Copies Atoms, Bonds and Coordinates to a new molecule
	 * 
	 * @param m
	 *            New Molecule.
	 */
	public void copyMinimum(Molecule m) {
		m.nLonePairs = nLonePairs;

		for (Atom atom : atoms) {
			m.atoms.add(atom.copy(m));
			double[] coord = new double[4];
			Coord.copy(coords.get(atom.getNo()), coord);
			m.coords.add(coord);
		}

		for (Bond bond : bonds) {
			m.bonds.add(bond.copy(m));
		}

	}

	/**
	 * Merges one molecule into this one.
	 * 
	 * @param other
	 * @param init
	 */
	public void merge(Molecule other, boolean init) {

		int startAtomNo = getnAtoms();
		int bondNo = getnBonds();

		int atomNo = startAtomNo;
		for (Atom atom : other.atoms) {
			Atom newAtom = atom.copy(this);
			newAtom.setNo(atomNo++);
			atoms.add(newAtom);
		}
		for (double[] coord : other.coords) {
			coords.add(coord);
		}

		for (Bond bond : other.bonds) {
			Atom a1 = atoms.get(startAtomNo + bond.getAtom1().getNo());
			Atom a2 = atoms.get(startAtomNo + bond.getAtom2().getNo());
			Bond newBond = new Bond(bondNo++, a1, a2, bond.getName());
			bonds.add(newBond);
		}

		if (init) {
			init();
		} else {
			structureChanged = true;
		}
	}

	/**
	 * Writes out an SDF entry to a file
	 * 
	 * @param file
	 * @param comment
	 */
	public void writeSdfFile(String file, String comment) {
		try {
			file = file.replace(' ', '_');
			FileWriter out = new FileWriter(new File(file));
			writeSdfMol(out, comment);
			out.close();
		} catch (IOException ex) {
			System.err.println("writeSdfFile IO exception: " + ex);
		}
	}

	/**
	 * Returns an SDF entry for the molecule as a String
	 * 
	 * @return
	 */
	public String sdfString() {
		StringWriter writer = new StringWriter();
		try {
			writeSdfMol(writer, null);
		} catch (IOException ex) {
			;
		}
		return writer.toString();
	}

	/**
	 * Writes out a SDF entry for this molecule to a Writer.
	 * 
	 * @param writer
	 * @param info
	 * @throws IOException
	 */
	public void writeSdfMol(Writer writer, String info) throws IOException {
		PrintWriter out = new PrintWriter(writer);

		// Header
		out.write(name + "\n");
		// out.write(" -GAPE---08220413213D\n");
		out.write("\n");
		out.write("\n");

		// Renumber to remove LPs
		outputLonePairs(false);
		renumberForOutput();

		// Counts Line
		out.printf("%3d%3d", nOutputAtoms, nOutputBonds);
		// out.write(" 0 0 0 0 0 0 0999 V2000\n");

		out.write("  0  0");
		if (sdChiral)
			out.write("  1");
		else
			out.write("  0");
		out.write("  0  0  0  0  0999 V2000\n");

		// Atom Block
		for (Atom atom : atoms) {

			// LP
			if (atom.getNo() == -1)
				continue;

			double coord[] = coords.get(atom.getNo());
			out.printf("%10.4f%10.4f%10.4f", coord[0], coord[1], coord[2]);
			out.printf(" %-3s", atom.getType().sdType());

			int massDiffference = atom.getSdMassDifference() != null
					? atom.getSdMassDifference() : 0;
			if (massDiffference < 0)
				out.printf("-%1d", -massDiffference);
			else
				out.printf(" %1d", -massDiffference);
			int charge = 0;
			if (atom.getFormalCharge() != null) {
				charge = 4 - atom.getFormalCharge();
			}
			out.printf("%3d", charge);
			int stereo = 0;
			Atom.SdStereo sdStereo = atom.getSdStereo();
			if (sdStereo == SdStereo.ODD)
				stereo = 1;
			if (sdStereo == SdStereo.EVEN)
				stereo = 2;
			if (sdStereo == SdStereo.EITHER)
				stereo = 3;
			out.printf("%3d", stereo);
			out.write("  0  0");
			int valance = atom.getSdValence() == null ? 0
					: atom.getSdValence() == 0 ? 15 : atom.getSdValence();
			out.printf("%3d", valance);
			out.write("  0  0  0  0  0  0\n");
		}

		// Bond Block
		for (Bond bond : bonds) {

			// Bond to LP
			if (bond.getNo() == -1)
				continue;

			int a1 = bond.getAtom1().getNo() + 1;
			int a2 = bond.getAtom2().getNo() + 1;
			out.printf("%3d%3d%3d", a1, a2, bond.sdType());
			int stereo = 0;
			Bond.Stereo bondStereo = bond.getStereo();
			if (bondStereo == Bond.Stereo.UP)
				stereo = 1;
			else if (bondStereo == Bond.Stereo.DOWN)
				stereo = 4;
			else if (bondStereo == Bond.Stereo.EITHER)
				stereo = 6;
			else if (bondStereo == Bond.Stereo.CIS_TRANS)
				stereo = 3;
			out.printf("%3d", stereo);
			out.write("  0  0  0\n");

		}

		out.write("M  END\n");

		// restore atom numbering
		renumber();

		if (info != null && !info.equals("")) {
			if (sdfFields == null)
				sdfFields = new HashMap<String, String>();
			if (!sdfFields.containsKey("COMMENT"))
				sdfFields.put("COMMENT", info);
		}

		if (sdfFields != null) {
			for (String k : sdfFields.keySet()) {
				String v = sdfFields.get(k);
				out.write(">  <" + k + ">\n" + v + "\n\n");
			}
		}
		out.write("$$$$\n");
	}

	/**
	 * Creates an array of molecules from local files. Specify filename(s).
	 * Checks types but does not perform solvation etc (i.e. calls build, but
	 * not init).
	 * 
	 * @param files
	 * @return
	 */
	public static List<Molecule> loadFiles(String files[]) {
		return loadFiles(files, false);
	}

	/**
	 * Creates an array of molecules from local files. Specify filename(s).
	 * Checks types but does not perform solvation etc unless init is set (i.e.
	 * calls build, but not init).
	 * 
	 * @param files
	 * @see #build()
	 * @see #loadSybylMol2(BufferedReader)
	 * @see #loadSdfMol(BufferedReader)
	 * 
	 * @return
	 */
	public static List<Molecule> loadFiles(String files[], boolean init) {

		FileType type = FileType.UNKNOWN;

		ArrayList<Molecule> molecules = new ArrayList<Molecule>();
		int molNo = 0;
		for (int i = 0; i < files.length; i++) {
			System.out.println("loading " + files[i]);
			BufferedReader in = openReader(files[i], Source.FILE);
			if (in == null)
				continue;
			if (files[i].toUpperCase().endsWith(".MOL2")
					|| files[i].toUpperCase().endsWith(".MOL2.GZ")) {
				type = FileType.MOL2;
			} else {
				type = FileType.SDF;
			}

			while (true) {
				Molecule molecule = new Molecule();
				try {
					if (type == FileType.MOL2) {
						molecule.loadSybylMol2(in);
					} else if (type == FileType.SDF) {
						molecule.loadSdfMol(in);
					}
				} catch (IOException ex) {
					break;
				} catch (MolReadException ex) {
					break;
				}
				if (init) {
					molecule.init();
				} else {
					molecule.build();
				}
				molecules.add(molecule);
				molNo++;
				if (molecule.name.equals("") || molecule.name.equals("****"))
					molecule.name = "Structure_" + molNo;
			}
			try {
				in.close();
			} catch (IOException ex) {
				System.err.println("IO error closing " + ex);
			}
		}

		return molecules;
	}

	/**
	 * Writes out a mol2 of sdf file based on file suffix.
	 * 
	 * @param molecules
	 * @param file
	 * @param comment
	 */
	static public void write(List<? extends Molecule> molecules, String file,
			String comment) {
		Molecule.FileType type = Molecule.getType(file);
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(file));
			for (Molecule molecule : molecules) {
				switch (type) {
				case MOL2:
					molecule.writeSybylMol2(out, comment, true);
					break;
				case SDF:
					molecule.writeSdfMol(out, comment);
					break;
				default:
					throw new IllegalArgumentException("can't write file type " + type);
				}
			}
			out.close();
		} catch (IOException ex) {
			throw new RuntimeException("IOException: " + ex.toString());
		}
	}

	/**
	 * Writes an array of structures to a MOL2 or SDFile depending on file
	 * format
	 * 
	 * @param molecules
	 *            array of molecules
	 * @param file
	 *            filename (does not include suffix)
	 * @param type
	 *            file type
	 * @param incLP
	 *            set to include lone pairs (applicable to MOL2 only- sdf files
	 *            dopn't contain lone pairs
	 */
	public static void writeMols(Molecule molecules[], String file, FileType type,
			boolean incLP) throws IOException {

		if (type == FileType.MOL2)
			file += ".mol2";
		else
			file += ".sdf";
		Writer out = new FileWriter(file);
		for (int i = 0; i < molecules.length; i++) {
			if (type == FileType.MOL2)
				molecules[i].writeSybylMol2(out, null, incLP);
			else
				molecules[i].writeSdfMol(out, null);
		}
		out.close();
	}

	/**
	 * Writes an array of structures to a MOL2 or SDFile. The file format is
	 * deducted from the file suffix.
	 * 
	 * @param molecules
	 * @param file
	 * @throws IOException
	 */
	public static void writeMols(List<Molecule> molecules, String file)
			throws IOException {

		FileType type = Molecule.getType(file);
		Writer out = new FileWriter(file);
		for (Molecule molecule : molecules) {
			if (type == FileType.MOL2)
				molecule.writeSybylMol2(out, null, false);
			else
				molecule.writeSdfMol(out, null);
		}
		out.close();
	}

	/**
	 * Removes all atoms in the molecule that lie outside the sphere of size
	 * distance centered at [x, y, z].
	 * 
	 * @param x
	 * @param y
	 * @param z
	 * @param distance
	 */
	public void inclusionSphere(double x, double y, double z, double distance) {

		double tolerance = distance * distance;
		double[] center = new double[] { x, y, z, 1.0 };
		ArrayList<Atom> removeAtoms = new ArrayList<Atom>();
		for (int i = 0; i < getnAtoms(); i++) {
			double[] coord = getCoord(i);
			double sqrDistance = Coord.sqrDistance(coord, center);
			if (sqrDistance > tolerance)
				removeAtoms.add(atoms.get(i));
		}

		for (Atom atom : removeAtoms) {
			deleteAtom(atom, false);
		}
	}

	/**
	 * Writes out a mol2 of sdf file based on file suffix.
	 * 
	 * @param file
	 * @param comment
	 */
	public void write(String file, String comment) {
		Molecule.FileType type = Molecule.getType(file);
		switch (type) {
		case MOL2:
			writeSybylMol2File(file, comment);
			break;
		case SDF:
			writeSdfFile(file, comment);
			break;
		default:
			throw new IllegalArgumentException("can't write file type " + type);
		}
	}

	public void write(FileType format, Writer out, String info, boolean incLP)
			throws IOException {
		if (format == Molecule.FileType.SDF)
			writeSdfMol(out, info);
		else
			writeSybylMol2(out, info, incLP);
	}

	/**
	 * @return the nBonds
	 */
	public int getnBonds() {
		return bonds.size();
	}

	/**
	 * @return the nRings
	 */
	public int getnRings() {
		return rings.size();
	}

	/**
	 * @return the nTorsions
	 */
	public int getnTorsions() {
		return torsions.size();
	}

	/**
	 * @return the nRotatableBonds
	 */
	public int getnRotatableBonds() {
		return rotatableBonds.size();
	}

	/**
	 * @return the nLonePairs
	 */
	public int getnLonePairs() {
		return nLonePairs;
	}

	public int getnAtoms() {
		return atoms.size();
	}

	public Atom getAtom(int no) {
		return atoms.get(no);
	}

	/**
	 * @return the atoms
	 */
	public List<Atom> getAtoms() {
		return atoms;
	}

	/**
	 * @return the coords
	 */
	public List<double[]> getCoords() {
		return coords;
	}

	public double[] getCoord(int no) {
		return coords.get(no);
	}

	/**
	 * @return the bonds
	 */
	public List<Bond> getBonds() {
		return bonds;
	}

	public Bond getBond(int no) {
		return bonds.get(no);
	}

	/**
	 * @return the rings
	 */
	public List<Ring> getRings() {
		return rings;
	}

	/**
	 * @return the nAngles
	 */
	public int getnAngles() {
		return angles.size();
	}

	/**
	 * @return the angles
	 */
	public List<Angle> getAngles() {
		return angles;
	}

	/**
	 * @return the torsions
	 */
	public List<Torsion> getTorsions() {
		return torsions;
	}

	/**
	 * @return the rotatableBonds
	 */
	public List<RotatableBond> getRotatableBonds() {
		return rotatableBonds;
	}

	/**
	 * @return the sdfFields
	 */
	public HashMap<String, String> getSdfFields() {
		return sdfFields;
	}

	/**
	 * @return the assignTypesFlag
	 */
	public static boolean getAssignTypesFlag() {
		return assignTypesFlag.get();
	}

	/**
	 * @return the addHydrogensFlag
	 */
	public static boolean getAddHydrogensFlag() {
		return addHydrogensFlag.get();
	}

	/**
	 * @return the solvateFlag
	 */
	public static boolean getSolvateFlag() {
		return solvateFlag.get();
	}

	/**
	 * @return the chargeFlag
	 */
	public static boolean getChargeFlag() {
		return chargeFlag.get();
	}

	public static void setAssignTypesFlag(boolean flg) {
		assignTypesFlag.set(flg);
	}

	/**
	 * @return the addHydrogensFlag
	 */
	public static void setAddHydrogensFlag(boolean flg) {
		addHydrogensFlag.set(flg);
	}

	/**
	 * @return the solvateFlag
	 */
	public static void setSolvateFlag(boolean flg) {
		solvateFlag.set(flg);
	}

	/**
	 * @return the chargeFlag
	 */
	public static void setChargeFlag(boolean flg) {
		chargeFlag.set(flg);
	}

	/**
	 * @return the infoMessageLogger
	 */
	public InfoMessageLogger getInfoMessageLogger() {
		return infoMessageLogger;
	}

	/**
	 * @param infoMessageLogger
	 *            the infoMessageLogger to set
	 */
	public void setInfoMessageLogger(InfoMessageLogger infoMessageLogger) {
		this.infoMessageLogger = infoMessageLogger;
	}

	/**
	 * Get molecule name
	 * 
	 * @return
	 */
	public String getName() {
		return name;
	}

	/**
	 * @return the baseName
	 */
	public String getBaseName() {
		return baseName;
	}

	/**
	 * Set molecule name
	 * 
	 * @param name
	 */
	public void setName(String name) {
		this.name = name;
	}

	/**
	 * @param coords
	 *            the coords to set
	 */
	public void setCoords(List<double[]> coords) {
		this.coords = coords;
	}

	/**
	 * @param atoms
	 *            the atoms to set
	 */
	public void setAtoms(List<Atom> atoms) {
		this.atoms = atoms;
	}

	/**
	 * @param nLonePairs
	 *            the nLonePairs to set
	 */
	public void setnLonePairs(int nLonePairs) {
		this.nLonePairs = nLonePairs;
	}

	/**
	 * @param rotatableBonds
	 *            the rotatableBonds to set
	 */
	public void setRotatableBonds(List<RotatableBond> rotatableBonds) {
		this.rotatableBonds = rotatableBonds;
	}

	void addRing(Ring ring) {
		rings.add(ring);
	}

	void clearRings() {
		rings.clear();
	}
}
