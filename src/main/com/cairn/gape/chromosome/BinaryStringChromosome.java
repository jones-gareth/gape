package com.cairn.gape.chromosome;

import com.cairn.gape.ga.GaSupervisor;

/**
 * Represents a Binary string Chromosome.
 * 
 * @author Gareth Jones
 * 
 */
public class BinaryStringChromosome implements Chromosome {
	protected int nBits = 0, nBytes = 0;

	// binary string
	protected byte[] bits = null;

	protected GaSupervisor problem;

	// set allowSwitch if crosover simply exchanges strings (i.e child1=parent2,
	// child2=parent1). Useful for multo-type chromosomes.
	private boolean allowSwitch;

	private static final boolean DEBUG = false;

	// 8 bit Graycode
	private static final int BINOFGRAY256[] = { 0, 1, 3, 2, 7, 6, 4, 5, 15, 14, 12, 13,
			8, 9, 11, 10, 31, 30, 28, 29, 24, 25, 27, 26, 16, 17, 19, 18, 23, 22, 20, 21,
			63, 62, 60, 61, 56, 57, 59, 58, 48, 49, 51, 50, 55, 54, 52, 53, 32, 33, 35,
			34, 39, 38, 36, 37, 47, 46, 44, 45, 40, 41, 43, 42, 127, 126, 124, 125, 120,
			121, 123, 122, 112, 113, 115, 114, 119, 118, 116, 117, 96, 97, 99, 98, 103,
			102, 100, 101, 111, 110, 108, 109, 104, 105, 107, 106, 64, 65, 67, 66, 71,
			70, 68, 69, 79, 78, 76, 77, 72, 73, 75, 74, 95, 94, 92, 93, 88, 89, 91, 90,
			80, 81, 83, 82, 87, 86, 84, 85, 255, 254, 252, 253, 248, 249, 251, 250, 240,
			241, 243, 242, 247, 246, 244, 245, 224, 225, 227, 226, 231, 230, 228, 229,
			239, 238, 236, 237, 232, 233, 235, 234, 192, 193, 195, 194, 199, 198, 196,
			197, 207, 206, 204, 205, 200, 201, 203, 202, 223, 222, 220, 221, 216, 217,
			219, 218, 208, 209, 211, 210, 215, 214, 212, 213, 128, 129, 131, 130, 135,
			134, 132, 133, 143, 142, 140, 141, 136, 137, 139, 138, 159, 158, 156, 157,
			152, 153, 155, 154, 144, 145, 147, 146, 151, 150, 148, 149, 191, 190, 188,
			189, 184, 185, 187, 186, 176, 177, 179, 178, 183, 182, 180, 181, 160, 161,
			163, 162, 167, 166, 164, 165, 175, 174, 172, 173, 168, 169, 171, 170 };

	// 4 bit graycode
	// private static final int BINOFGRAY16[] = { 0, 1, 3, 2, 7, 6, 4, 5, 15,
	// 14, 12, 13,
	// 8, 9, 11, 10 };

	protected BinaryStringChromosome() {
	}

	/**
	 * Creates a new binary string cheomosome.
	 * 
	 * @param p
	 *            GA problem
	 * @param n
	 *            number of bits
	 * @throws GaException
	 */
	public BinaryStringChromosome(GaSupervisor p, int n) {
		problem = p;
		nBits = n;
		createEmpty();
	}

	@Override
	public Chromosome create(GaSupervisor problemHandle) {
		this.problem = problemHandle;
		createEmpty();
		initialize();
		return this;
	}

	@Override
	public Chromosome create(GaSupervisor problemHandle, Object initialData) {
		this.problem = problemHandle;
		createEmpty();

		boolean[] data = null;
		try {
			data = (boolean[]) initialData;
		} catch (ClassCastException ex) {
			throw new IllegalStateException(
					"BinaryStringChromosome: create initial data must be of type boolean[]");
		}

		if (data.length != nBits)
			throw new IllegalStateException(
					"BinaryStringChromosome: create initial data must be of length nBits");

		for (int i = 0; i < nBytes; i++) {
			bits[i] = 00;
		}

		for (int i = 0; i < nBits; i++)
			if (data[i])
				swapBit(i);
		return this;
	}

	// Dummy methods so we can use Chromosome interface
	@Override
	public double getFitness() {
		throw new RuntimeException("method not defined");
	}

	@Override
	public String fitnessInfo() {
		return "";
	}

	@Override
	public double distance(Chromosome c) {
		throw new RuntimeException("method not defined");
	}

	@Override
	public boolean sameNiche(Chromosome c) {
		throw new RuntimeException("method not defined");
	}

	@Override
	public boolean ok() {
		throw new RuntimeException("method not defined");
	}

	@Override
	public double rebuild() {
		throw new RuntimeException("method not defined");
	}

	public Chromosome getChromosome(Chromosome c) {
		throw new RuntimeException("method not defined");
	}

	@Override
	public void freeChromosome() {
		throw new RuntimeException("method not defined");
	}

	@Override
	public int getChromosomeId() {
		throw new RuntimeException("method not defined");
	}

	@Override
	public Chromosome getChromosome() {
		throw new RuntimeException("method not defined");
	}

	public void initialize() {
		for (int i = 0; i < nBytes; i++) {
			bits[i] = 00;
		}

		for (int i = 0; i < nBits; i++)
			if (problem.normalRand() < 0.5)
				swapBit(i);
	}

	/**
	 * Creates storage for string.
	 * 
	 * @throws GaException
	 */
	public void createEmpty() {
		if (nBits == 0)
			throw new IllegalStateException("initialize: nBits is zero");
		nBytes = nBits / 8;
		if (nBits % 8 > 0)
			nBytes++;
		bits = new byte[nBytes];
	}

	/*
	 * Compares two binary string chromosomes
	 * 
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.ga.Chromosome#equals(com.cairn.gape.ga.Chromosome
	 * )
	 */
	@Override
	public boolean equals(Chromosome c) {
		BinaryStringChromosome c2 = (BinaryStringChromosome) c;
		byte bits2[] = c2.bits;
		for (int i = 0; i < nBytes; i++) {
			if (bits[i] != bits2[i])
				return false;
		}
		return true;
	}

	/**
	 * Performs mutation. The per-bit mutation probability is 1/length. The
	 * operation is repeated until a mutation occurs.
	 */
	public void mutate() {
		double pMutate = 1.0 / nBits;
		boolean mutated = false;
		for (int i = 0; i < nBits; i++) {
			if (problem.normalRand() < pMutate) {
				swapBit(i);
				mutated = true;
			}
		}
		if (!mutated)
			mutate();
	}

	/**
	 * Toggles a bit.
	 * 
	 * @param pos
	 */
	private void swapBit(int pos) {
		int bytePos = pos / 8;
		int bitPos = pos % 8;

		int b = bits[bytePos];

		int i = b & 1 << bitPos;
		b |= 1 << bitPos;
		b &= ~i;
		bits[bytePos] = (byte) b;
	}

	/**
	 * Performs one point crossover.
	 * 
	 * @param parent2
	 * @param child1
	 * @param child2
	 */
	public void crossover(BinaryStringChromosome parent2, BinaryStringChromosome child1,
			BinaryStringChromosome child2) {

		//
		// crossover
		//
		// performs one-point binary string crossover on parent1 and parent2 to
		// produce child1 and child2.

		int length = nBits / 8;
		int extra = nBits % 8;

		if (extra > 0)
			length++;

		// choose a cross point
		if (problem == null) {
			System.out.println("help");
		}
		int site = problem.randomInt(1, allowSwitch ? nBits : nBits - 1);
		if (DEBUG)
			System.out.println("Cross point " + site);

		byte p1[] = this.bits;
		byte p2[] = parent2.bits;

		byte c1[] = null;
		byte c2[] = null;

		boolean switchFlg = true;
		if (allowSwitch)
			switchFlg = problem.randomBoolean();
		if (switchFlg) {
			c1 = child1.bits;
			c2 = child2.bits;
		} else {
			c1 = child2.bits;
			c2 = child1.bits;
		}

		// x is the number of bytes before the cross point. y is any remaining
		// bits.
		int x = site / 8;
		int y = site % 8;

		// byte portions of child before cross point
		int i;
		for (i = 0; i < x; i++) {
			c1[i] = p1[i];
			c2[i] = p2[i];
		}

		if (x == length)
			return;

		if (y > 0) {
			// now do bits within the byte that the cross point occurs
			int _c1 = 0, _c2 = 0;
			int _p1 = p1[i];
			int _p2 = p2[i];

			// first the bits before the cross point
			int j = 0;
			for (j = 0; j < y; j++) {
				if (bitSet(_p1, j))
					_c1 |= 1 << j;
				if (bitSet(_p2, j))
					_c2 |= 1 << j;
			}

			/* next the bits after the cross point */
			for (j = y; j < 8; j++) {
				if (bitSet(_p1, j))
					_c2 |= 1 << j;
				if (bitSet(_p2, j))
					_c1 |= 1 << j;
			}

			c1[i] = (byte) _c1;
			c2[i] = (byte) _c2;
			i = x + 1;
		} else
			i = x;

		// next the bytes after the cross point
		for (; i < length; i++) {
			c1[i] = p2[i];
			c2[i] = p1[i];
		}

	}

	/**
	 * Checks a bit within a byte. Used by crossover op.
	 * 
	 * @param val
	 * @param pos
	 * @return
	 */
	private boolean bitSet(int val, int pos) {
		if ((((val) >> (pos)) & 01) == 01)
			return true;
		return false;
	}

	/**
	 * @param pos
	 *            position on the string
	 * @return the value of a bit
	 */
	public boolean bitSet(int pos) {
		int val = bits[pos / 8];
		pos = pos % 8;
		return bitSet(val, pos);
	}

	/*
	 * Copies this binary string to another
	 * 
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.ga.Chromosome#copyGene(com.cairn.gape.ga.Chromosome
	 * )
	 */
	@Override
	public void copyGene(Chromosome c2) {
		BinaryStringChromosome c = (BinaryStringChromosome) c2;
		if (c.bits == null || c.nBits != nBits) {
			c.nBits = nBits;
			c.nBytes = nBytes;
			c.bits = new byte[nBytes];
		}
		for (int i = 0; i < nBytes; i++)
			c.bits[i] = bits[i];
	}

	/*
	 * Formats the binary string
	 * 
	 * (non-Javadoc)
	 * 
	 * @see com.cairn.gape.ga.Chromosome#geneInfo()
	 */
	@Override
	public String geneInfo() {
		String rtn = "";
		int no = 0;
		for (int i = 0; i < nBytes; i++) {
			int val = bits[i];
			for (int j = 0; j < 8; j++) {
				if (no == nBits)
					break;
				if (bitSet(val, j))
					rtn += '1';
				else
					rtn += '0';
				no++;
			}
			rtn += ' ';
		}
		return rtn + "\n";
	}

	/**
	 * @param pos
	 *            position on string.
	 * @return byte binary value starting at pos. Reads 8 bits from string.
	 */
	public int posToInt(int pos) {
		int val = 0;
		for (int i = 0; i < 8; i++) {
			if (bitSet(pos + i))
				val += (1 << i);
		}
		return val;
	}

	/**
	 * @param pos
	 * @param size
	 * @return binary value from decoding by reading size bits from position pos
	 */
	public int posToInt(int pos, int size) {
		int val = 0;
		for (int i = 0; i < size; i++) {
			if (bitSet(pos + i))
				val += (1 << i);
		}
		return val;
	}

	/**
	 * @param pos
	 *            position on string.
	 * @return 8-bit Graycode value starting at pos. Reads 8 bits from string.
	 */
	public int posToGray(int pos) {
		int val = posToInt(pos);
		return BINOFGRAY256[val];
	}

	public int getNBits() {
		return nBits;
	}

	public GaSupervisor getProblem() {
		return problem;
	}

	public void setProblem(GaSupervisor problem) {
		this.problem = problem;
	}

	// int byteToGray(int byteNo) {
	// int val = (int) bits[byteNo];
	// val += 128;
	// return BINOFGRAY256[val];
	// }

	// int byteToInt(int byteNo) {
	// int val = (int) bits[byteNo];
	// val += 128;
	// return val;
	// }

}
