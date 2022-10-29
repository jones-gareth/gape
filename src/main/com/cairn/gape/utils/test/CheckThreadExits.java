package com.cairn.gape.utils.test;

import com.cairn.common.utils.CommonUtils;
import com.cairn.gape.Superposition;
import com.cairn.gape.ga.TestProgressListener;

/**
 * Checks to see if gape can handle sequential jobs with the same thread id.
 * Unable to create two threads with the same threadid on Linux.
 * 
 * @author gjones
 * 
 */
public class CheckThreadExits {

	public static void main(String[] args) {

		try {
			String configFiles[] = new String[] {
					"/home/gjones/gape/test/test_multi_thread/trypsin/superposition.conf",
					"/home/gjones/gape/test/test_multi_thread/CDK2/superposition.conf",
					"/home/gjones/gape/test/test_multi_thread/P38/superposition.conf" };

			GapeJob job = new GapeJob(configFiles[0]);
			Thread t = job.createGapeThread();
			Thread.sleep(1000L);
			System.out.println("Thread id for Ist gape job is " + t.getId()
					+ " name " + t.getName());
			t = null;
			Thread.sleep(30000L);
			job.getSuperposition().cancel();
			System.out.println("1st job cancelled");

			Thread.sleep(5000L);
			System.gc();
			System.runFinalization();
			Thread.sleep(5000L);

			job = new GapeJob(configFiles[1]);
			t = job.createGapeThread();
			Thread.sleep(1000L);
			System.out.println("Thread id for 2nd gape job is " + t.getId()
					+ " name " + t.getName());
			t = null;

			Thread.sleep(30000L);
			job.getSuperposition().cancel();
			System.out.println("2nd job cancelled");

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}

/**
 * A class to represent a gape job in a thread
 * 
 * @author gjones
 * 
 */
class GapeJob {
	Superposition superposition;
	String configFile;
	Thread thread;

	public GapeJob(String configFile) {
		super();
		this.configFile = configFile;
	}

	public Thread createGapeThread() {

		thread = new Thread() {
			@Override
			public void run() {
				runGape(configFile);
			}

		};

		thread.start();
		return thread;

	}

	private void runGape(String configFile) {

		try {
			superposition = new Superposition();

			superposition.addProgressListener(new TestProgressListener());
			superposition.fitMolecules(configFile, null);
			superposition.getInfoMessageLogger().finish();
		} catch (Exception ex) {
			String stackTrace = CommonUtils.getStackTrace(ex);
			try {
				superposition.infoMessageln("Exception:");
				superposition.infoMessageln(stackTrace);
			} catch (Exception e) {
				;
			}
			System.err.println(MultiThreadGape.class.getName()
					+ ": Exception: " + stackTrace);
		}
	}

	/**
	 * @return the sp
	 */
	public Superposition getSuperposition() {
		return superposition;
	}

	/**
	 * @return the configFile
	 */
	public String getConfigFile() {
		return configFile;
	}

	/**
	 * @return the thread
	 */
	public Thread getThread() {
		return thread;
	}

}