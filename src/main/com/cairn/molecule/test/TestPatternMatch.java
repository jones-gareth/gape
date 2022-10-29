package com.cairn.molecule.test;

import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import com.cairn.molecule.Molecule;
import com.cairn.molecule.PatternMatch;

/**
 * Simple unit test to evaluate pattern matches
 * 
 * @author Gareth Jones
 *
 */
public class TestPatternMatch {
	private static final Logger logger = Logger.getLogger(TestPatternMatch.class);

	/**
	 * @throws java.lang.Exception
	 */
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	/**
	 * @throws java.lang.Exception
	 */
	@Before
	public void setUp() throws Exception {
	}

	private static class TestPatternMatcher extends PatternMatch {
		public TestPatternMatcher(String sln) {
			super(sln);
		}

		private Set<String> matchingLabels = new HashSet<>();

		@Override
		public void process() {
			int match[] = queryMatches[nMatches - 1];
			String targetLabel = target.getAtom(match[0]).getLabel();
			matchingLabels.add(targetLabel);
			logger.debug("Query matches " + targetLabel);
		}

	}

	private boolean test(String fileName, String pattern, String[] matchingLabels)
			throws Exception {
		logger.info("Matching " + pattern + " against " + fileName);
		Molecule molecule = TestUtil.getMoleculeFromSybylResource(fileName);
		TestPatternMatcher patternMatch = new TestPatternMatcher(pattern);
		patternMatch.match(molecule);

		Set<String> matches = patternMatch.matchingLabels;
		Set<String> matchingLabelsSet = new HashSet<>(Arrays.asList(matchingLabels));
		String matchesStr = matches.stream().collect(Collectors.joining(", "));
		String matchingLabelsStr = matchingLabelsSet.stream().collect(
				Collectors.joining(", "));
		logger.info("Expected " + matchingLabelsStr + " got " + matchesStr);
		return matchingLabelsSet.equals(matches);

	}

	@Test
	public void testOmega1() throws Exception {
		assertTrue(test("omega_1.mol2", "N[1:f]H:-=Hev:-=Hev:-=Hev:-=Hev:-=@1",
				new String[] { "N3" }));
	}

	@Test
	public void testAmideN() throws Exception {
		assertTrue(test("2dhf_aligned.mol2", "N[f]H(-Hev)-C=O",
				new String[] { "N", "N3" }));
	}

	@Test
	public void testAmideO() throws Exception {
		assertTrue(test("2dhf_aligned.mol2", "O=C-N[f]H-Hev", new String[] { "O", "OH4" }));
	}

}