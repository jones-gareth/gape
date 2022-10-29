package com.cairn.molecule.viewer;

import java.awt.Component;

import javax.swing.UIManager;
import javax.swing.UIManager.LookAndFeelInfo;

public class Utils {

	/**
	 * Sets the default look and feel (currently System)
	 */
	public static final void setLookAndFeel(Component c) {

		// On Linux I've had stack overflows in GTK look and feel
		// and also lots of initially blank applets.
		if (System.getProperty("os.name").equals("Linux")) {
			if (false) { // forget about nimbus for now- it's horrible
				try {

					for (LookAndFeelInfo info : UIManager.getInstalledLookAndFeels()) {
						if ("Nimbus".equals(info.getName())) {
							UIManager.setLookAndFeel(info.getClassName());
							return;
						}
					}
				} catch (Exception e) {
					// If Nimbus is not available, you can set the GUI to
					// another
					// look and feel.
				}
			} else
				return;
		}

		try {
			// Set System L&F
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		} catch (Exception e) {
			System.err.println("Unable to set look and feel");
		}

	}

}
