package com.cairn.gape.utils;

/**
 * A breath-first search algorithm. The search space is defined by a series of
 * ranges (maximum and minimum values- such that a position can have any value
 * between minimum and maximum-1). The algorithm steps through all available
 * values.
 * 
 * @author Gareth Jones
 * 
 */
public class BreadthSearch {
	private int maximums[], minimums[];
	private int currentMapping[];
	private int stepPosition, mappingLength;
	private long nMappings, mappingCount;

	/**
	 * @param i
	 * @return top value at a position
	 */
	private int getMaxValue(int i) {
		return maximums[i];
	}

	/**
	 * If minimums are not set, we start at 0.
	 * 
	 * @param i
	 * @return minimum value at a position
	 */
	private int getMinValue(int i) {
		if (minimums == null)
			return 0;
		return minimums[i];
	}

	/**
	 * Sets up the start vector
	 */
	public void setStart() {
		currentMapping = new int[mappingLength];
		nMappings = 1;
		for (int i = 0; i < mappingLength; i++) {
			int min = getMinValue(i);
			int max = getMaxValue(i);
			nMappings *= (max - min);
			currentMapping[i] = min;
		}

		stepPosition = 0;
	}

	/**
	 * Finds the next set of values in the search space.
	 * 
	 * @return false is we've finished searching.
	 */
	public boolean createNextMapping() {
		boolean stepPositionChanged = false;

		// find the next position we can increment
		while (currentMapping[stepPosition] + 1 >= getMaxValue(stepPosition)) {
			stepPosition++;
			stepPositionChanged = true;
			if (stepPosition >= mappingLength) {
				assert nMappings - 1 == mappingCount;
				return false;
			}
		}
		mappingCount++;

		// increment at step position
		currentMapping[stepPosition]++;

		// if the step position has changed reset previous positions, and start
		// again at the beginning
		if (stepPositionChanged) {
			for (int i = 0; i < stepPosition; i++)
				currentMapping[i] = getMinValue(i);
			stepPosition = 0;
		}
		return true;
	}

	/**
	 * @return a string describing the current point in the search space
	 */
	public String mappingInfo() {
		StringBuffer buf = new StringBuffer();
		buf.append("No:").append(mappingCount + 1).append("/").append(nMappings)
				.append(" [");
		buf.append(currentMapping[0]);
		for (int i = 1; i < mappingLength; i++)
			buf.append(",").append(currentMapping[i]);
		buf.append("]");
		return buf.toString();
	}

	// getters and setters

	public void setMaximums(int[] maximums) {
		this.maximums = maximums;
		mappingLength = maximums.length;
	}

	public void setMinimums(int[] minimums) {
		this.minimums = minimums;
	}

	public int[] getMaximums() {
		return maximums;
	}

	public int[] getMinimums() {
		return minimums;
	}

	public int[] getCurrentMapping() {
		return currentMapping;
	}

	public int getStepPosition() {
		return stepPosition;
	}

	public long getNMappings() {
		return nMappings;
	}

	/**
	 * Test routine- just loops through the search space printing out values.
	 * 
	 * @param args
	 */
	public static void main(String args[]) {
		int maximums[] = new int[] { 99, 50, 200 };
		int minimums[] = new int[] { -1, 0, -1 };

		BreadthSearch breadthSearch = new BreadthSearch();
		breadthSearch.setMaximums(maximums);
		breadthSearch.setMinimums(minimums);
		breadthSearch.setStart();

		System.out.println("nMappings " + breadthSearch.getNMappings());
		System.out.println(breadthSearch.mappingInfo());

		while (breadthSearch.createNextMapping()) {
			System.out.println(breadthSearch.mappingInfo());
		}
	}

}
