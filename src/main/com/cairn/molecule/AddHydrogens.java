package com.cairn.molecule;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;
import com.cairn.gape.feature.LonePairAddition;
import com.cairn.gape.utils.InfoMessageLogger;

/**
 * Class to fill valance in SDF files by adding hydrogens
 * 
 * @author Gareth Jones
 * 
 */
public class AddHydrogens {
	private static Logger logger = Logger.getLogger(AddHydrogens.class);
	static {
		// logger.setLevel(Level.DEBUG);
	}

	private Molecule molecule;

	// Look up table for bond lengths to hydrogens
	static final HashMap<String, Double> hBondLengths;

	static {
		hBondLengths = new HashMap<String, Double>();
		hBondLengths.put("C.1", 1.008);
		hBondLengths.put("C.2", 1.089);
		hBondLengths.put("C.3", 1.1);
		hBondLengths.put("C.ar", 1.084);
		hBondLengths.put("N.3", 1.08);
		hBondLengths.put("N.am", 1.0);
		hBondLengths.put("O.3", 0.95);
		hBondLengths.put("*", 1.08);
	}

	// Pattern for N in Arginine groups
	private static ThreadLocal<MolPattern> argPattern = new ThreadLocal<MolPattern>() {
		@Override
		protected MolPattern initialValue() {
			MolPattern guaPattern = null;
			try {
				guaPattern = MolPattern.generateMolPattern("N:C[f](:N):N");
			} catch (Exception ex) {
				guaPattern = null;
			}
			return guaPattern;
		}
	};

	/**
	 * Constructor. Does not actually add hydrogens
	 * 
	 * @param m
	 */
	public AddHydrogens(Molecule m) {
		molecule = m;
	}

	/**
	 * Constructor. Does not actually add hydrogens
	 */
	public AddHydrogens() {
		;
	}

	/**
	 * Sets molecule for hydrogen addition.
	 * 
	 * @param m
	 */
	public void setMolecule(Molecule m) {
		molecule = m;
	}

	/**
	 * Fills valence in a molecule by adding hydrogens. Uses static methods
	 * developed for adding Lone Pairs to position the hydrogens correctly in
	 * 3-D space
	 * 
	 * 
	 * @see LonePairAddition
	 */
	public void addHydrogens() {

		// Get atom count no as adding hydrogens will increase number of atoms!
		int nAtoms = molecule.getnAtoms();

		// Get Arginine N atoms
		// The solvate class will handle this for traditional bond arrangements
		ArrayList<Atom> argAtoms = AtomPatternMatch.matchPattern(argPattern.get(),
				molecule);

		for (int i = 0; i < nAtoms; i++) {
			Atom atom = molecule.getAtom(i);

			// Arginine
			if (argAtoms.contains(atom)) {
				if (atom.getnNotDummyNeighbours() == 1)
					addTwoHydrogensToTrigonal(atom, hBondLengths.get("*"), 2);
				else if (atom.getnNotDummyNeighbours() == 2)
					addOneHydrogenToTrigonal(atom, hBondLengths.get("*"));
				continue;
			}

			double nbo = atom.getType().getNeutralBondOrder();
			double nbo2 = atom.getType().getNeutralBondOrder2();
			double bo = atom.getBondOrder();
			Integer formalCharge = atom.getFormalCharge();
			if (formalCharge == null) {
				formalCharge = 0;
			}

			double hCnt = nbo - bo + formalCharge;

			if (hCnt < -0.4 && nbo2 > nbo)
				hCnt = nbo2 - bo + formalCharge;

			if (hCnt < 0.9)
				continue;

			if (logger.isDebugEnabled())
				logger.debug(atom.info() + " hCnt " + hCnt + " nbo " + nbo + " bo "
						+ atom.getBondOrder());

			// Add hydrogens to Carbon
			if (atom.getType().isCarbonType()) {
				int nHydrogens = (int) Math.round(hCnt);

				if (atom.getnAromaticNeighbours() == 2
						&& atom.getnNotDummyNeighbours() == 2) {
					if (nHydrogens != 1)
						throw new RuntimeException("Hydrogen Count error for "
								+ atom.info());
					addOneHydrogenToTrigonal(atom, hBondLengths.get("C.ar"));
					continue;
				}

				if (atom.getnDoubleNeighbours() == 1 && atom.getnSingleNeighbours() == 1
						&& atom.getnNotDummyNeighbours() == 2) {
					if (nHydrogens != 1)
						throw new RuntimeException("Hydrogen Count error for "
								+ atom.info());
					addOneHydrogenToTrigonal(atom, hBondLengths.get("C.2"));
					continue;
				}

				if (atom.getnDoubleNeighbours() == 1
						&& atom.getnNotDummyNeighbours() == 1) {
					if (nHydrogens != 2)
						throw new RuntimeException("Hydrogen Count error for "
								+ atom.info());
					addTwoHydrogensToTrigonal(atom, hBondLengths.get("C.2"), 2);
					continue;
				}

				if (atom.getnSingleNeighbours() == 3
						&& atom.getnNotDummyNeighbours() == 3) {
					if (nHydrogens != 1)
						throw new RuntimeException("Hydrogen Count error for "
								+ atom.info());
					addOneHydrogenToTetrahedral(atom, hBondLengths.get("C.3"));
					continue;
				}

				if (atom.getnSingleNeighbours() == 2
						&& atom.getnNotDummyNeighbours() == 2) {
					if (nHydrogens != 2)
						throw new RuntimeException("Hydrogen Count error for "
								+ atom.info());
					addTwoHydrogensToTetrahedral(atom, hBondLengths.get("C.3"), 2);
					continue;
				}

				if (atom.getnSingleNeighbours() == 1
						&& atom.getnNotDummyNeighbours() == 1) {
					if (nHydrogens != 3)
						throw new RuntimeException("Hydrogen Count error for "
								+ atom.info());
					addThreeHydrogensToTetrahedral(atom, hBondLengths.get("C.3"), 3);
					continue;
				}
				if (atom.getnTripleNeighbours() == 1
						&& atom.getnNotDummyNeighbours() == 1) {
					if (nHydrogens != 1)
						throw new RuntimeException("Hydrogen Count error for "
								+ atom.info());
					addOneHydrogenToLinear(atom, hBondLengths.get("C.1"));
					continue;
				}

				if (atom.isArginineCarbon()) {
					continue;
				}
			}

			// Add hydrogens to oxygens
			else if (atom.getType().isOxygenType()) {
				// Don't add H to aromatically bonded O
				if (bo > 1.4 && bo < 1.6)
					continue;
				int nHydrogens = (int) Math.round(hCnt);

				if (atom.getnSingleNeighbours() == 1
						&& atom.getnNotDummyNeighbours() == 1) {
					if (nHydrogens != 1)
						throw new RuntimeException("Hydrogen Count error for "
								+ atom.info());
					addThreeHydrogensToTetrahedral(atom, hBondLengths.get("O.3"), 1);
					continue;
				}
			}

			// Add hydrogens to nitrogens
			else if (atom.getType().isNitrogenType()) {

				int nHydrogens = (int) Math.round(hCnt);
				// Add 2H to single aromatically bonded N
				if (bo > 1.4 && bo < 1.6)
					nHydrogens = 2;

				if (atom.getnDoubleNeighbours() == 0) {

					// NAM
					boolean added = false;
					for (Atom neighbour : atom.getNotDummyNeighbours()) {
						if (neighbour.isAmideCarbon()
								|| neighbour.isSulphonamideSulphur()) {
							if (atom.getnNotDummyNeighbours() == 1) {
								addTwoHydrogensToTrigonal(atom, hBondLengths.get("N.am"),
										nHydrogens);
							} else {
								if (nHydrogens != 1)
									throw new IllegalStateException(
											"Hydrogen Count error for " + atom.info());
								addOneHydrogenToTrigonal(atom, hBondLengths.get("N.am"));
								;
							}
							added = true;
							break;
						}
					}

					if (added)
						continue;

					// NPL3
					added = false;
					for (Atom neighbour : atom.getNotDummyNeighbours()) {
						if (neighbour.getType().isCarbonType()
								&& (neighbour.getnDoubleNeighbours() >= 1 || neighbour
										.getnAromaticNeighbours() >= 2)) {
							if (atom.getnNotDummyNeighbours() == 1) {
								addTwoHydrogensToTrigonal(atom, hBondLengths.get("*"),
										nHydrogens);
							} else {
								if (nHydrogens != 1)
									throw new IllegalStateException(
											"Hydrogen Count error for " + atom.info());
								addOneHydrogenToTrigonal(atom, hBondLengths.get("*"));
								;
							}
							added = true;
							break;
						}
					}

					if (added)
						continue;

				}

				// Aromatic ntirogens- need to have charge to be protonated
				if (atom.getnAromaticNeighbours() == 2 && nHydrogens == 1) {
					addOneHydrogenToTrigonal(atom, hBondLengths.get("*"));
					continue;
				}

				// Sp2 nitrogens, with two heavy neighbours- need to have charge
				// to be protonated
				if (atom.getnDoubleNeighbours() == 1
						&& atom.getnNotDummyNeighbours() == 2 && nHydrogens == 1) {
					addOneHydrogenToTrigonal(atom, hBondLengths.get("*"));
					continue;
				}

				// Sp2 ntirogens
				if (atom.getnDoubleNeighbours() == 1
						&& atom.getnNotDummyNeighbours() == 1) {
					if (nHydrogens != 1)
						throw new IllegalStateException("Hydrogen Count error for "
								+ atom.info());
					addTwoHydrogensToTrigonal(atom, hBondLengths.get("*"), 1);
					continue;
				}

				// Tetrahedral nitrogens
				if (atom.getnSingleNeighbours() == 1
						&& atom.getnNotDummyNeighbours() == 1) {
					if (nHydrogens != 2 && nHydrogens != 3)
						throw new IllegalStateException("Hydrogen Count error for "
								+ atom.info());
					addThreeHydrogensToTetrahedral(atom, hBondLengths.get("N.3"),
							nHydrogens);
					continue;
				}

				if (atom.getnSingleNeighbours() == 2
						&& atom.getnNotDummyNeighbours() == 2) {
					if (nHydrogens != 1 && nHydrogens != 2)
						throw new IllegalStateException("Hydrogen Count error for "
								+ atom.info());
					addTwoHydrogensToTetrahedral(atom, hBondLengths.get("N.3"),
							nHydrogens);
					continue;
				}

				if (atom.getnSingleNeighbours() == 3
						&& atom.getnNotDummyNeighbours() == 3 && nHydrogens == 1) {
					addOneHydrogenToTetrahedral(atom, hBondLengths.get("N.3"));
					continue;
				}

			}

			// Sulphur with valence 4 and three connections
			else if (atom.getType().isSulphurType()) {
				int nHydrogens = (int) Math.round(hCnt);

				if (atom.getnSingleNeighbours() == 3
						&& atom.getnNotDummyNeighbours() == 3) {
					if (nHydrogens != 1)
						throw new IllegalStateException("Hydrogen Count error for "
								+ atom.info());
					addOneHydrogenToTetrahedral(atom, hBondLengths.get("*"));
					continue;
				}
			}

			logger.warn("Don't know how to add hydrogens to " + atom.info());
			// throw new AddHydrogenException("Failed to add hydrogens to "
			// + atom.info());
		}

		// Need to rebuild atom neighbour lists- should also be able to set
		// correct types.
		molecule.update();
	}

	/**
	 * Adds single hydrogen to an sp2 atom which has two heavy-atom neighbours.
	 * 
	 * @param atom
	 * @param distance
	 * 
	 * @see LonePairAddition#addOnePairToTrigonal(double[], double[], double[],
	 *      double[])
	 */
	void addOneHydrogenToTrigonal(Atom atom, double distance) {
		double protonCoord[] = new double[4];
		double[] coord1 = molecule.getCoord(atom.getNo());
		double[] coord2 = molecule.getCoord(atom.getNotDummyNeighbours().get(0).getNo());
		double[] coord3 = molecule.getCoord(atom.getNotDummyNeighbours().get(1).getNo());

		// Hijack the Lone Pair addition code to add the proton.
		if (!LonePairAddition.addOnePairToTrigonal(coord1, coord2, coord3, protonCoord,
				distance))
			throw new RuntimeException("failed to add hydrogen to " + atom.info());

		addHydrogensToAtom(atom, 1, new double[][] { protonCoord });

	}

	/**
	 * Adds a hydrogen to an sp1 atom which has one heavy-atom neighbour.
	 * 
	 * @param atom
	 * @param distance
	 * 
	 * @see LonePairAddition#addOnePairToLinear(double[], double[], double[])
	 */
	void addOneHydrogenToLinear(Atom atom, double distance) {
		double protonCoord[] = new double[4];
		double[] coord1 = molecule.getCoord(atom.getNo());
		double[] coord2 = molecule.getCoord(atom.getNotDummyNeighbours().get(0).getNo());

		// Hijack the Lone Pair addition code to add the proton.
		LonePairAddition.addOnePairToLinear(coord1, coord2, protonCoord, distance);

		addHydrogensToAtom(atom, 1, new double[][] { protonCoord });
	}

	/**
	 * Adds up to two hydrogens to an sp2 atom which has one heavy atom
	 * neighbour.
	 * 
	 * @param atom
	 * @param distance
	 * @param nHydrogens
	 * 
	 * @see LonePairAddition#addTwoPairsToTrigonal(double[], double[], double[],
	 *      double[], double[])
	 */
	void addTwoHydrogensToTrigonal(Atom atom, double distance, int nHydrogens) {
		double coord[] = molecule.getCoord(atom.getNo());
		Atom otherAtom = atom.getNotDummyNeighbours().get(0);
		Atom thirdAtom = null;
		for (Atom test : otherAtom.getNotDummyNeighbours()) {
			if (test != atom) {
				thirdAtom = test;
			}
		}

		double coord2[] = molecule.getCoord(otherAtom.getNo());
		double coord3[] = molecule.getCoord(thirdAtom.getNo());
		double proton1Coord[] = new double[4];
		double proton2Coord[] = new double[4];

		// Hijack the Lone Pair addition code to add the proton.
		if (!LonePairAddition.addTwoPairsToTrigonal(coord, coord2, coord3, proton1Coord,
				proton2Coord, distance))
			throw new RuntimeException("failed to add hydrogen to " + atom.info());
		addHydrogensToAtom(atom, nHydrogens,
				new double[][] { proton1Coord, proton2Coord });
	}

	/**
	 * Adds one hydrogen to a an sp3 atom which has three heavy-atom neighbours.
	 * 
	 * @param atom
	 * @param distance
	 * 
	 * @see LonePairAddition#addOnePairToTetrahedral(double[], double[],
	 *      double[], double[], double[])
	 */
	void addOneHydrogenToTetrahedral(Atom atom, double distance) {

		double coord[] = molecule.getCoord(atom.getNo());
		double coord2[] = molecule.getCoord(atom.getNotDummyNeighbours().get(0).getNo());
		double coord3[] = molecule.getCoord(atom.getNotDummyNeighbours().get(1).getNo());
		double coord4[] = molecule.getCoord(atom.getNotDummyNeighbours().get(2).getNo());
		double protonCoord[] = new double[4];

		// Hijack the Lone Pair addition code to add the proton.
		LonePairAddition.addOnePairToTetrahedral(coord, coord2, coord3, coord4,
				protonCoord, distance);

		addHydrogensToAtom(atom, 1, new double[][] { protonCoord });
	}

	/**
	 * Adds one or two hydrogens to a sp3 atom which has two heavy-atom
	 * neighbours.
	 * 
	 * @param atom
	 * @param distance
	 * @param nHydrogens
	 * 
	 * @see LonePairAddition#addThreePairsToTetrahedral(double[], double[],
	 *      double[], double[], double[])
	 */
	void addTwoHydrogensToTetrahedral(Atom atom, double distance, int nHydrogens) {

		double coord[] = molecule.getCoord(atom.getNo());
		double coord2[] = molecule.getCoord(atom.getNotDummyNeighbours().get(0).getNo());
		double coord3[] = molecule.getCoord(atom.getNotDummyNeighbours().get(1).getNo());
		double proton1Coord[] = new double[4];
		double proton2Coord[] = new double[4];

		// Hijack the Lone Pair addition code to add the proton.
		LonePairAddition.addTwoPairsToTetrahedral(coord, coord2, coord3, proton1Coord,
				proton2Coord, distance);

		addHydrogensToAtom(atom, nHydrogens,
				new double[][] { proton1Coord, proton2Coord });
	}

	/**
	 * Adds up to three hydrogens to an sp3 atom which has only one heavy atom
	 * heignbour.
	 * 
	 * @param atom
	 * @param distance
	 * @param nHydrogens
	 * 
	 * @see LonePairAddition#addThreePairsToTetrahedral(double[], double[],
	 *      double[], double[], double[])
	 */
	void addThreeHydrogensToTetrahedral(Atom atom, double distance, int nHydrogens) {

		double coord[] = molecule.getCoord(atom.getNo());
		double coord2[] = molecule.getCoord(atom.getNotDummyNeighbours().get(0).getNo());
		double proton1Coord[] = new double[4];
		double proton2Coord[] = new double[4];
		double proton3Coord[] = new double[4];

		// Hijack the Lone Pair addition code to add the proton.
		LonePairAddition.addThreePairsToTetrahedral(coord, coord2, proton1Coord,
				proton2Coord, proton3Coord, distance);

		addHydrogensToAtom(atom, nHydrogens, new double[][] { proton1Coord, proton2Coord,
				proton3Coord });
	}

	/**
	 * Updates molecule atom and bond lists with the new hydrogens.
	 * 
	 * @param atom
	 * @param nHydrogens
	 * @param coords
	 */
	void addHydrogensToAtom(Atom atom, int nHydrogens, double[][] coords) {
		InfoMessageLogger infoMessageLogger = molecule.getInfoMessageLogger();

		infoMessageLogger.infoMessageln(3, "Adding " + nHydrogens
				+ " hydrogen(s) to atom " + atom.info());
		for (int i = 0; i < nHydrogens; i++) {
			if (logger.isDebugEnabled())
				logger.debug("adding hydrogen at " + Coord.info(coords[i]));
			if (coords[i][3] < 0.999 || coords[i][3] > 1.00001)
				throw new RuntimeException("addHydrogensToAtom: bad coordinate");
			Atom proton = new Atom(molecule, molecule.getnAtoms(), AtomType.Type.H);
			Bond bond = new Bond(molecule.getnBonds(), atom, proton, BondType.Type.SINGLE);
			molecule.addAtom(proton, coords[i]);
			molecule.addBond(bond);

		}
	}

	/**
	 * Main test routine. Adds hydrogens to an sd file.
	 * 
	 * @param args
	 */
	public static void main(String args[]) {
		String sdFile = args[0];
		try {

			List<Molecule> mols = Molecule.loadFiles(new String[] { sdFile });
			for (Molecule mol : mols) {
				AddHydrogens addHydrogens = new AddHydrogens(mol);
				addHydrogens.addHydrogens();
			}

			Molecule.write(mols, "h_added_" + sdFile, "Solvated");

		} catch (Exception ex) {
			System.err.println("Exception " + ex);
			ex.printStackTrace();
		}
	}
}
