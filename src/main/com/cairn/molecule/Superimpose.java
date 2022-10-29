package com.cairn.molecule;

import java.util.Map;
import java.util.Map.Entry;

import com.cairn.common.utils.Coord;
import com.cairn.molecule.Ullman.MatchType;

// A class to compare 2 molecular conformations: superimposes all
// isomorphisms of the two conformations and returns an rms for each
// conformation.

class Superimpose implements Ullman.UllmanCallback {
	protected volatile Molecule confX, confY;
	private volatile double refX[][];
	private volatile Ullman ull;

	public Superimpose() {
		;
	}

	Superimpose(Molecule m1, Molecule m2) {
		search(m1, m2);
	}

	public void search(Molecule m1, Molecule m2) {
		confX = m1;
		confY = m2;
		int nAtoms = confX.getnAtoms();
		refX = new double[nAtoms][4];
		for (int i = 0; i < nAtoms; i++)
			Coord.copy(confX.getCoord(i), refX[i]);

		ull = new Ullman(confX, confY, MatchType.MATCH);
		ull.doUllman(this);

	}

	public static void main(String[] args) {
		try {
			if (args.length != 2) {
				System.out
						.println("usage: Superimpose <conformation1.mol2> <conformation2.mol2>");
				System.exit(0);
			}
			Molecule confX = new Molecule();
			Molecule confY = new Molecule();
			confX.loadFile(args[0]);
			confY.loadFile(args[1]);
			new Superimpose(confX, confY);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void processRms(double rms, int no) {
		System.out.println("Isomorphism " + no + " rms " + rms);

		confX.writeSybylMol2File("Fitted Isomorphism " + no + " " + confX.baseName
				+ ".mol2", "Fitted Molecule");
	}

	@Override
	public void callback(Map<Atom, Atom> atomMapping) {
		double mapX[][] = new double[ull.getMolA().getnAtoms()][];
		double mapY[][] = new double[ull.getMolB().getnAtoms()][];
		int mapping[] = new int[confX.getnAtoms()];
		for (int i = 0; i < mapping.length; i++)
			mapping[i] = -1;

		int no = 0;
		for (Entry<Atom, Atom> entry : atomMapping.entrySet()) {
			int no1 = entry.getKey().getNo();
			int no2 = entry.getValue().getNo();
			mapX[no] = refX[no1];
			mapY[no] = confY.getCoord(no2);
			mapping[no1] = no2;
			no++;
		}

		double trans[][] = new double[4][4];
		Coord.leastSquaresFit(mapX, mapY, trans);
		for (int i = 0; i < confX.getnAtoms(); i++)
			Coord.transPoint(trans, refX[i], confX.getCoord(i));

		no = ull.getIsomorphisms();
		double d = .0, cnt = 0;
		for (int i = 0; i < confX.getnAtoms(); i++) {
			if (mapping[i] == -1)
				continue;
			double p1[] = confX.getCoord(i);
			double p2[] = confY.getCoord(mapping[i]);
			for (int j = 0; j < 3; j++) {
				double diff = p1[j] - p2[j];
				d += diff * diff;
			}
			cnt++;
		}
		double rms = Math.sqrt(d / cnt);

		processRms(rms, no);

	}
}
