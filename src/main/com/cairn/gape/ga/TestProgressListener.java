package com.cairn.gape.ga;

import com.cairn.gape.ga.BaseSupervisor.ProgressListener;

/**
 * Test progress listener
 * 
 * @author gjones
 * 
 */
public class TestProgressListener implements ProgressListener {

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * com.cairn.gape.ga.BaseSupervisor.ProgressListener#reportProgress
	 * (com.cairn.gape.ga.ProgressReport)
	 */
	@Override
	public void reportProgress(ProgressReport progressReport) {
		long threadId = Thread.currentThread().getId();
		System.out.println("Thread " + threadId + " Percent complete: "
				+ progressReport.getPercentComplete());
		System.out.println("Message         : " + progressReport.getMessage());

	}

}
