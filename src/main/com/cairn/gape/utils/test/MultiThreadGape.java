package com.cairn.gape.utils.test;

import com.cairn.common.utils.CommonUtils;
import com.cairn.gape.Superposition;
import com.cairn.gape.ga.TestProgressListener;

/**
 * A class for testing multi thread execution of GAPE
 * 
 * @author gjones
 * 
 */
public class MultiThreadGape {

	public static void main(String[] args) {

		if (args.length < 1) {
			System.out.println("Usage: " + MultiThreadGape.class.getName()
					+ " -version | <conf files > ");
			System.exit(0);
		}

		for (final String configFile : args) {
			createGapeThread(configFile);
			System.out.println("Started thread for config file " + configFile);
		}
	}

	/**
	 * @param configFile
	 * @return a new running GAPE thread
	 */
	static Thread createGapeThread(final String configFile) {

		Thread t = new Thread() {
			@Override
			public void run() {
				runGape(configFile);
			}

		};

		t.start();
		return t;

	}

	/**
	 * Run a GAPE thread.
	 * 
	 * @param configFile
	 */
	static void runGape(String configFile) {
		Superposition sp = null;
		try {
			sp = new Superposition();

			// attach a progress listener to the gape job
			sp.addProgressListener(new TestProgressListener());
			sp.fitMolecules(configFile, null);
			sp.getInfoMessageLogger().finish();
		} catch (Exception ex) {
			String stackTrace = CommonUtils.getStackTrace(ex);
			try {
				sp.infoMessageln("Exception:");
				sp.infoMessageln(stackTrace);
			} catch (Exception e) {
				;
			}
			System.err.println(MultiThreadGape.class.getName()
					+ ": Exception: " + stackTrace);
		}
	}
}
