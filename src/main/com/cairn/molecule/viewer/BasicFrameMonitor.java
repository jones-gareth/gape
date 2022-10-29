package com.cairn.molecule.viewer;

import java.awt.Window;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

public class BasicFrameMonitor extends WindowAdapter {
	@Override
	public void windowClosing(WindowEvent e) {
		Window w = e.getWindow();
		w.setVisible(false);
		w.dispose();
	}
}
