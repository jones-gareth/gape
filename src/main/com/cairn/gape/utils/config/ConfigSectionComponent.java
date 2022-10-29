package com.cairn.gape.utils.config;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.BevelBorder;

import com.cairn.common.utils.UtilLayout;

/**
 * This class contains information about a section together with code to display
 * in a panel.
 */
public class ConfigSectionComponent extends JPanel {
	protected UtilLayout layout;
	protected int y = 0;
	protected final ConfigEditor configEditor;
	protected final ConfigSection configSection;
	protected final List<ConfigEntryComponent<?>> configEntryComponents = new ArrayList<ConfigEntryComponent<?>>();

	public ConfigSectionComponent(ConfigEditor editor,
			ConfigSection configSection) {
		super();
		this.configEditor = editor;
		this.configSection = configSection;
	}

	/**
	 * Draw editor for this section
	 * 
	 */
	public void display() {
		setLayout(new BorderLayout());
		JComponent ta = configEditor.infoText(configSection.getDescription()
				.replaceAll("\n", ""));
		JPanel p1 = new JPanel(new FlowLayout(FlowLayout.LEFT));
		ta.setBorder(new BevelBorder(BevelBorder.RAISED));
		p1.add(ta);
		add(BorderLayout.NORTH, p1);

		JPanel settingsPanel = new JPanel();
		layout = new UtilLayout(settingsPanel);

		y = 0;
		for (int i = 0; i < configEntryComponents.size(); i++) {

			final ConfigEntryComponent<?> entryComponent = configEntryComponents
					.get(i);
			JLabel l = new JLabel(entryComponent.getConfigEntry().getName()
					.replace('_', ' '));
			layout.add(l, 1, y);
			entryComponent.setDefault();
			layout.add(entryComponent.getEditor(), 2, y);

			if (entryComponent.isOptional()) {
				// handle optional variables
				final JCheckBox optCb = new JCheckBox();
				boolean use = entryComponent.getConfigEntry().isUse();
				optCb.setSelected(use);
				entryComponent.setChecked(use);
				optCb.addItemListener(new ItemListener() {
					public void itemStateChanged(ItemEvent e) {
						entryComponent.setChecked(optCb.isSelected());
					}
				});
				layout.add(optCb, 0, y);
			}

			y++;
			if (entryComponent.getHelp() != null) {
				// print description
				JComponent helpTa = configEditor.infoText(entryComponent
						.getConfigEntry().getHelp().replaceAll("\n", ""));
				layout.c.fill = GridBagConstraints.BOTH;
				layout.add(helpTa, 1, y++, 3, 1);
				layout.c.fill = GridBagConstraints.NONE;
			}
		}

		JPanel p = new JPanel(new FlowLayout(FlowLayout.LEFT));
		p.add(settingsPanel);
		// JScrollPane sp = new JScrollPane(p);
		add(BorderLayout.CENTER, p);
	}

	public void addConfigEntryComponent(ConfigEntryComponent<?> component) {
		configEntryComponents.add(component);
	}

	/**
	 * @return the configEntryComponents
	 */
	public List<ConfigEntryComponent<?>> getConfigEntryComponents() {
		return configEntryComponents;
	}

	/**
	 * @return the configSection
	 */
	public ConfigSection getConfigSection() {
		return configSection;
	}

}
