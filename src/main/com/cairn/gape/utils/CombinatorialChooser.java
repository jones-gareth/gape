package com.cairn.gape.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * A class to enumerate all possible selections of n (unordered) items from a
 * list. Enumerates all items from n choose r calling a callback for each
 * selection.
 * 
 * @author Gareth Jones
 *
 * @param <T>
 */
public class CombinatorialChooser<T> {
	private final int nChoose;
	private final List<T> values;
	private int count;
	private final CombinatorialChooserCallback<T> callback;

	/**
	 * Callback for when we find a selection
	 * 
	 * @param <T>
	 */
	@FunctionalInterface
	public interface CombinatorialChooserCallback<T> {
		void callback(List<T> values);
	}

	/**
	 * 
	 * @param nChoose
	 *            number of items to pick
	 * @param values
	 *            list of possible items
	 * @param callback
	 *            calls this for each possible selection
	 */
	public CombinatorialChooser(int nChoose, List<T> values,
			CombinatorialChooserCallback<T> callback) {
		super();
		this.nChoose = nChoose;
		this.values = values;
		this.callback = callback;
	}

	/**
	 * Start the search
	 */
	public void enumerate() {
		// Create array for storing chosen values
		List<T> chosen = new ArrayList<>(nChoose);
		IntStream.range(0, nChoose).forEach(it -> chosen.add(null));

		// do the biz
		count = 0;
		choose(values, nChoose, chosen);
	}

	/**
	 * Recursive function that does the work
	 * 
	 * @param currentValues
	 * @param length
	 * @param chosen
	 */
	private void choose(List<T> currentValues, int length, List<T> chosen) {
		for (int i = 0; i < currentValues.size(); i++) {
			// add current value to chosen values
			chosen.set(chosen.size() - length, currentValues.get(i));
			if (length == 1) {
				// got to last item- call callback
				callback.callback(chosen);
				count++;
			} else {
				// pick remaining items from next positions
				choose(currentValues.subList(i + 1, currentValues.size()), length - 1,
						chosen);
			}
		}
	}

	/**
	 * @return the count
	 */
	public int getCount() {
		return count;
	}

	/**
	 * A combinatorial chooser that processes item indexes rather than lists of
	 * items
	 * 
	 * @param nChoose
	 * @param nAvailable
	 * @param callback
	 * @return
	 */
	public static CombinatorialChooser<Integer> combinatorialIndexChooser(int nChoose,
			int nAvailable, CombinatorialChooserCallback<Integer> callback) {
		List<Integer> indices = IntStream.range(0, nAvailable).boxed()
				.collect(Collectors.toList());
		return new CombinatorialChooser<Integer>(nChoose, indices, callback);
	}

	/**
	 * Test routine
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		CombinatorialChooser<String> chooser = new CombinatorialChooser<String>(3,
				Arrays.asList("a", "b", "c", "d", "e", "f", "g", "h"),
				(vals) -> System.out
						.println(vals.stream().collect(Collectors.joining(", "))));
		chooser.enumerate();
		System.out.println("Found " + chooser.getCount() + " combinations");

		CombinatorialChooser<Integer> indexChooser = CombinatorialChooser
				.combinatorialIndexChooser(3, 8,
						(vals) -> System.out
								.println(vals.stream().map(i -> String.valueOf(i))
										.collect(Collectors.joining(", "))));
		indexChooser.enumerate();
		System.out.println("Found " + indexChooser.getCount() + " combinations");
	}
}
