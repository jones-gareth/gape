package com.cairn.molecule;

import java.util.List;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.log4j.Logger;

/**
 * A class to read molecules from one file and write them to another. Command
 * line arguments are used to control adding hydrogens, solvation etc.
 * 
 * @author Gareth Jones
 *
 */
public class ConvertFiles {
	private static final Logger logger = Logger.getLogger(ConvertFiles.class);

	public static void main(String[] args) throws Exception {

		Option inFileOption = new Option("i", "inFile", true, "Input file");
		Option outFileOption = new Option("o", "outFile", true, "Output file");
		Option buildOption = new Option("b", "build", false,
				"Build only (no init performed)");
		Option assignTypesOption = new Option("t", "assignTypes", false, "Assign Types");
		Option addHydrogensOption = new Option("y", "addHydrogens", false,
				"Fill valance with hydrogens");
		Option solvateOption = new Option("s", "solvate", false, "Solvate");
		Option chargeOption = new Option("c", "charge", false, "Set formal charges");
		Option numericNameSuffixOption = new Option("ns", "numericNameSuffix", true,
				"Add a suffix to all molecule names that look like conformer names (e.g compound_1)");
		Option allOption = new Option("a", "all", false,
				"Same as assign types, charge, addHydrogens and solvate");
		Option help = new Option("h", "help", false, "Print this message");

		Options options = new Options();
		options.addOption(inFileOption);
		options.addOption(outFileOption);
		options.addOption(buildOption);
		options.addOption(assignTypesOption);
		options.addOption(addHydrogensOption);
		options.addOption(solvateOption);
		options.addOption(chargeOption);
		options.addOption(help);
		options.addOption(allOption);
		options.addOption(numericNameSuffixOption);

		CommandLineParser parser = new BasicParser();
		CommandLine line = null;
		try {
			// parse the command line arguments
			line = parser.parse(options, args);
		} catch (ParseException exp) {
			System.err.println("Parsing failed.  Reason: " + exp.getMessage());
			return;
		}

		if (line.hasOption("help")) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(ConvertFiles.class.getName(), options);
			return;
		}

		String inFile = line.getOptionValue("i");
		if (inFile == null) {
			logger.error("Input file is required");
			return;
		}
		String outFile = line.getOptionValue("o");
		if (outFile == null) {
			logger.error("Output file is required");
			return;
		}
		boolean addHydrogens = line.hasOption("y");
		boolean assignTypes = line.hasOption("t");
		boolean solvate = line.hasOption("s");
		boolean charge = line.hasOption("c");
		boolean build = line.hasOption("b");
		boolean all = line.hasOption("a");
		if (all) {
			addHydrogens = assignTypes = solvate = charge = true;
		}
		if (!build) {
			logger.info("addHydrogens = " + addHydrogens);
			logger.info("assignTypes = " + assignTypes);
			logger.info("solvate = " + solvate);
			logger.info("charge = " + charge);
		} else {
			logger.info("build only- no init");
		}

		String numericNameSuffix = line.getOptionValue("numericNameSuffix");

		Molecule.setAddHydrogensFlag(addHydrogens);
		Molecule.setAssignTypesFlag(assignTypes);
		Molecule.setChargeFlag(charge);
		Molecule.setSolvateFlag(solvate);

		List<Molecule> molecules = Molecule.loadFiles(new String[] { inFile }, !build);

		for (Molecule molecule : molecules) {
			if (numericNameSuffix != null) {
				String name = molecule.getName();
				if (name.matches(".*_\\d+?$")) {
					name = name + numericNameSuffix;
					molecule.setName(name);
				}
			}
		}

		Molecule.writeMols(molecules, outFile);

	}
}
