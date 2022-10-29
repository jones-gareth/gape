package com.cairn.gape.utils.config;

import java.awt.BorderLayout;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JTextField;

import com.cairn.gape.utils.config.StructureFileEntryComponent.StructureFileListener;

/**
 * Special section for assigning activities to molecules.
 * 
 */
class ActivitySectionComponent extends ConfigSectionComponent implements
		StructureFileListener {

	private JComponent info;

	private JCheckBox useActivitiesCb;

	private JTextField activityFields[];

	ActivitySectionComponent(ConfigEditor configEditor,
			ConfigSection configSection) {
		super(configEditor, configSection);
	}

	private void enableFields() {
		boolean v = useActivitiesCb.isSelected();
		for (int i = 0; i < activityFields.length; i++)
			activityFields[i].setEnabled(v);
	}

	private void init() {
		for (ConfigEntryComponent<?> component : getConfigEntryComponents()) {

			if (component.getConfigEntry().getName().equals("USE_ACTIVITIES")) {
				useActivitiesCb = (JCheckBox) ((BooleanConfigEntryComponent) component)
						.getEditor();
			}
		}

		if (useActivitiesCb == null)
			System.err.println("no use Activities entry");

		ItemListener listeners[] = useActivitiesCb.getItemListeners();
		for (int i = 0; i < listeners.length; i++)
			useActivitiesCb.removeItemListener(listeners[i]);

		useActivitiesCb.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {
				enableFields();
			}
		});
	}

	@Override
	public void display() {
		if (this.configEditor.structureFileEntryComponent == null) {
			ConfigEditor.logger.warn("ActivityPane: no structure file menu\n");
			return;
		}

		this.configEditor.structureFileEntryComponent.addListener(this);
		ActivitySection activitySection = (ActivitySection) configSection;

		String mols[] = this.configEditor.structureFileEntryComponent.molNames;

		if (mols == null) {
			setLayout(new BorderLayout());
			ConfigEditor.logger.debug("No mols defined");
			info = this.configEditor.infoText("You need to have picked a "
					+ "structure data file before " + "assigning activities");
			add(BorderLayout.NORTH, info);
			return;
		}

		super.display();
		init();

		activityFields = new JTextField[mols.length];

		layout.add(new JLabel("MOLECULE"), 1, y);
		layout.add(new JLabel("ACTIVITY"), 2, y++);
		for (int i = 0; i < mols.length; i++) {
			layout.add(new JLabel(mols[i]), 1, y);
			activityFields[i] = new JTextField("", 20);
			if (activitySection.getActivities().containsKey(mols[i])) {
				Double act = activitySection.getActivities().get(mols[i]);
				activityFields[i].setText(act.toString());
			}
			layout.add(activityFields[i], 2, y++);
		}

		enableFields();
	}

	public void tellListenerNames(String molNames[]) {
		removeAll();
		display();
		revalidate();
	}
}