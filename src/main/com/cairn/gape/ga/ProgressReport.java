package com.cairn.gape.ga;

/**
 * A class for containing progress information
 * 
 */
public class ProgressReport {
	private final double percentComplete;
	private final String message;

	public ProgressReport(double percentComplete, String message) {
		super();
		this.percentComplete = percentComplete;
		this.message = message;
	}

	/**
	 * @return the percentComplete
	 */
	public double getPercentComplete() {
		return percentComplete;
	}

	/**
	 * @return the message
	 */
	public String getMessage() {
		return message;
	}

}