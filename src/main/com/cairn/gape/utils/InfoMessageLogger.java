package com.cairn.gape.utils;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;

/**
 * A class for logging informational messages
 * 
 * @author gjones
 * 
 */
public class InfoMessageLogger {

	protected int logLevel = 4;
	// protected PrintStream infoOutput = System.out;
	protected PrintStream infoOutput;

	/**
	 * Sets the info logger to a file
	 */
	public boolean setupInfoOutput(String logFile) {

		try {
			FileOutputStream fileOut = new FileOutputStream(logFile);
			infoOutput = new PrintStream(fileOut);
			return true;
		} catch (FileNotFoundException ex) {
			;
		}

		return false;
	}

	public void infoMessage(int minLevel, String message) {
		if (logLevel > minLevel)
			infoMessage(message);
	}

	public void infoMessageln(int minLevel, String message) {
		if (logLevel > minLevel)
			infoMessageln(message);
	}

	public void infoMessage(String msg) {
		if (infoOutput == null)
			System.out.print(msg);
		else
			infoOutput.print(msg);
	}

	public void infoMessageln(String msg) {
		infoMessage(msg + "\n");
	}

	public void infoMessageln() {
		infoMessage("\n");
	}

	/**
	 * @return the logLevel
	 */
	public int getLogLevel() {
		return logLevel;
	}

	/**
	 * @param logLevel
	 *            the logLevel to set
	 */
	public void setLogLevel(int logLevel) {
		this.logLevel = logLevel;
	}

	/**
	 * @param infoOutput
	 *            the infoOutput to set
	 */
	public void setInfoOutput(PrintStream infoOutput) {
		this.infoOutput = infoOutput;
	}

	public void finish() {
		try {
			if (infoOutput != null)
				infoOutput.close();
		} catch (Exception e) {
			;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#finalize()
	 */
	@Override
	protected void finalize() throws Throwable {
		super.finalize();
		finish();
	}

}
