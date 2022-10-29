package com.cairn.gape.utils.config;

import java.io.Writer;
import java.util.Map;

public class ActivitySection extends ConfigSection {
	protected boolean useActivities;
	protected Map<String, Double> activities;

	public ActivitySection(String n, String d) {
		super(n, d);
	}

	@Override
	public void printEntries(Writer out, boolean showDescription) {
		super.printEntries(out, showDescription);

		try {
			out.write("\n");
			out.write(" Activity lines are of the form: "
					+ "activity <pKi> <molecule_name>\n");
			out.write("\n");
			if (useActivities) {
				for (String name : activities.keySet()) {
					Double pki = activities.get(name);
					if (pki != null)
						out.write("activity " + pki + " " + name + "\n");
				}
			}
		} catch (Exception e) {
			logger.error("IO exception writing info ", e);
		}
	}

	/**
	 * @return the useActivities
	 */
	public boolean isUseActivities() {
		return useActivities;
	}

	/**
	 * @param useActivities
	 *            the useActivities to set
	 */
	public void setUseActivities(boolean useActivities) {
		this.useActivities = useActivities;
	}

	/**
	 * @return the activities
	 */
	public Map<String, Double> getActivities() {
		return activities;
	}

	/**
	 * @param activities
	 *            the activities to set
	 */
	public void setActivities(Map<String, Double> activities) {
		this.activities = activities;
	}

}
