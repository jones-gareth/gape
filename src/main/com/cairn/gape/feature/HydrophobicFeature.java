package com.cairn.gape.feature;

import java.util.ArrayList;
import java.util.List;

import com.cairn.gape.molecule.GaMolecule;
import com.cairn.gape.utils.InfoMessageLogger;
import com.cairn.molecule.Atom;
import com.cairn.molecule.AtomType;

/**
 * Represents a Hydrophobic atom (Carbon). This feature is used in the fitting
 * process, but is not used in scoring.
 * 
 * @author Gareth Jones
 * 
 */
public class HydrophobicFeature extends Feature {

	private HydrophobicFeature() {
		featureType = FeatureType.HYDROPHOBIC_ATOM;
		featureSetName = "HYDROPHOBIC_ATOM";
		atomFeature = true;
	}

	/**
	 * Contructs a hydrophobic feature.
	 * 
	 * @param m
	 * @param a
	 * @param no
	 */
	public HydrophobicFeature(GaMolecule m, Atom a, int no) {
		this();
		molecule = m;
		atom = a;
		featureSetNo = no;
	}

	/**
	 * Searches a molecule for hydrophobic atoms
	 * 
	 * @param m
	 * @return
	 */
	public static List<Feature> findFeatures(GaMolecule m) {
		InfoMessageLogger infoMessageLogger = m.getInfoMessageLogger();
		ArrayList<Feature> v = new ArrayList<Feature>();
		infoMessageLogger.infoMessage(2, "Looking for Hydrophobic atoms: ");
		int no = 0;
		for (Atom atom : m.getAtoms()) {
			if (atom.getAtomType() == AtomType.Type.C3
					|| atom.getAtomType() == AtomType.Type.C2
					|| atom.getAtomType() == AtomType.Type.CAR) {
				int an = atom.getNo() + 1;
				infoMessageLogger.infoMessage(2, atom.getType().getName() + "," + an
						+ " ");
				HydrophobicFeature feature = new HydrophobicFeature(m, atom, no);
				v.add(feature);

				no++;
			}
		}

		infoMessageLogger.infoMessageln(2, "");
		return v;
	}

	@Override
	public void printInfo() {
		System.out.print(info());
	}

	@Override
	public String info() {
		int no = featureSetNo + 1;
		return "HydroP Atom [" + no + " " + atom.getType().getName() + "]";
	}

	@Override
	public double score(Feature f) {
		return .0;
	}

	@Override
	public String pharmLabel() {
		return featureSetName + " ATOM_" + String.valueOf(atom.getNo() + 1);
	}

	public String label() {
		int no = featureSetNo + 1;
		return "HYDROPHOBIC_ATOM_" + no;
	}

	@Override
	public PharmFeatureGeometry getPharmFeatureGeometry() {
		return new SpherePharmFeatureGeometry(coordinate, atom.getType().getRadius());
	}
}
