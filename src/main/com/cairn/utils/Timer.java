package com.cairn.common.utils;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.text.NumberFormat;

/**
 * Thread timer that uses the thread system managed bean.
 * 
 * @author Gareth Jones
 * 
 */
public class Timer {
	private double start;
	private double startUser;
	private double end;
	private double endUser;
	private double total;
	private double user;
	private double system;
	static final NumberFormat nf = NumberFormat.getInstance();;

	static {
		nf.setMaximumFractionDigits(2);
		nf.setGroupingUsed(false);
	}

	public Timer() {
		reset();
	}

	private double nanoTimeToSeconds(long nano) {
		return nano * 1e-9;
	}

	public synchronized void reset() {
		ThreadMXBean threadMXBean = ManagementFactory.getThreadMXBean();
		if (!threadMXBean.isCurrentThreadCpuTimeSupported()) {
			start = startUser = .0;
		} else {
			start = nanoTimeToSeconds(threadMXBean.getCurrentThreadCpuTime());
			startUser = nanoTimeToSeconds(threadMXBean.getCurrentThreadUserTime());
		}

	}

	public synchronized void interval() {
		ThreadMXBean threadMXBean = ManagementFactory.getThreadMXBean();
		if (!threadMXBean.isCurrentThreadCpuTimeSupported()) {
			end = endUser = .0;
		} else {
			end = nanoTimeToSeconds(threadMXBean.getCurrentThreadCpuTime());
			endUser = nanoTimeToSeconds(threadMXBean.getCurrentThreadUserTime());
		}
		total = end - start;
		user = endUser - startUser;
		system = total - user;
	}

	public synchronized String info() {
		if (user > .0)
			return "Total: " + nf.format(total) + "s [user: " + nf.format(user)
					+ " system: " + nf.format(system) + "]";
		else
			return "Total: " + nf.format(total);
	}

	public String total() {
		return nf.format(total);
	}

}
