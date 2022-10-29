package com.cairn.gape.molecule;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import org.apache.commons.collections.CollectionUtils;
import org.apache.log4j.Logger;

import com.cairn.common.utils.Coord;
import com.cairn.gape.chromosome.BinaryStringChromosome;
import com.cairn.gape.chromosome.IntegerStringChromosome;
import com.cairn.gape.chromosome.SuperpositionChromosome;
import com.cairn.gape.feature.AcceptorAtomFeature;
import com.cairn.gape.feature.AromaticRingFeature;
import com.cairn.gape.feature.DonorHydrogenFeature;
import com.cairn.gape.feature.Feature;
import com.cairn.gape.feature.Feature.FeatureType;
import com.cairn.gape.feature.FeatureMapping;
import com.cairn.gape.feature.HydrogenBondingType;
import com.cairn.gape.feature.HydrophobicFeature;
import com.cairn.gape.feature.LonePairAddition;
import com.cairn.gape.feature.UserFeatureSet;
import com.cairn.gape.ga.GaSupervisor;
import com.cairn.gape.utils.GaussianList;
import com.cairn.gape.utils.InfoMessageLogger;
import com.cairn.molecule.Atom;
import com.cairn.molecule.AtomType;
import com.cairn.molecule.Bond;
import com.cairn.molecule.Molecule;
import com.cairn.molecule.RotatableBond;
import com.cairn.molecule.Taff;
import com.cairn.molecule.Torsion;

/**
 * This class is represents a molecule. Extends Molecule class to accommodate
 * operations needed in GA computations.
 *
 * @author Gareth Jones
 * @see com.cairn.intranet.molcule.Molecule
 */
public class GaMolecule extends Molecule {

    public static final double PASS1_DISTANCE = 2.0, MAX_PASS_DISTANCE = 5.0;
    private static final Logger logger = Logger.getLogger(GaMolecule.class);
    private static final boolean CHECK_TORSION_MATCHES_PREDICTED = true;
    private static final double PROX_DISTANCE = 1.5;

    static {
        // logger.setLevel(Level.DEBUG);
    }

    /**
     * Creates an array of molecules. Specify filename(s), structure format and
     * source. Does not prepare molecules.
     *
     * @param files
     * @param type
     * @param source
     * @return
     */
    public static List<GaMolecule> loadFiles(String files[], FileType type, Source source) {
        return loadFiles(files, type, source, new InfoMessageLogger(), false);
    }

    /**
     * Creates an array of molecules. Specify filename(s), structure format,
     * source, log level and default output stream. Uses init routine to prepare
     * molecules (unlike loadFiles in Molecule)- if init is set.
     *
     * @param files
     * @param type
     * @param source
     * @param logLevel
     * @param out
     * @return
     * @see Molecule#init()
     * @see #loadSybylMol2(BufferedReader)
     * @see #loadSybylMol2(BufferedReader)
     */
    public static List<GaMolecule> loadFiles(String files[], FileType type,
                                             Source source, InfoMessageLogger infoMessageLogger, boolean init) {

        List<GaMolecule> molecules = new ArrayList<GaMolecule>();
        int molNo = 0;
        for (int i = 0; i < files.length; i++) {
            infoMessageLogger.infoMessageln("loading " + files[i]);
            BufferedReader in = openReader(files[i], source);
            if (in == null)
                continue;
            if (type == FileType.UNKNOWN) {
                type = Molecule.getType(files[i]);
            }
            while (true) {
                GaMolecule molecule = new GaMolecule();
                molecule.infoMessageLogger = infoMessageLogger;
                try {
                    if (type == FileType.MOL2) {
                        molecule.loadSybylMol2(in);
                    } else if (type == FileType.SDF) {
                        molecule.loadSdfMol(in);
                    }
                } catch (IOException | MolReadException ex) {
                    break;
                }
                if (molecule.getName().endsWith("Pharmacophore"))
                    continue;
                if (init)
                    molecule.init();
                else
                    molecule.build();
                molecules.add(molecule);
                molNo++;
                if (molecule.getName().equals("") || molecule.getName().equals("****"))
                    molecule.setName("Structure_" + molNo);
            }
            try {
                in.close();
            } catch (IOException ex) {
                System.err.println("IO error closing " + ex);
            }
        }

        return molecules;
    }

    private static void debug(Supplier<String> message) {
        if (logger.isDebugEnabled()) {
            logger.debug(message.get());
        }
    }
    private final List<double[]> reference = new ArrayList<>();

    private final List<String> pharmDescriptions = new ArrayList<>();

    private final List<Feature> pharmFeatures = new ArrayList<>();

    private final List<FreeCorner> freeCorners = new ArrayList<>();

    private final List<Torsion> freeCornerTorsions = new ArrayList<>();

    private final List<LargeCycle> largeCycles = new ArrayList<>();

    private final Map<FeatureType, FeatureMapping> featureMappings = new LinkedHashMap<>();
    private final ThreadLocal<double[][]> in = ThreadLocal
            .withInitial(() -> new double[4][4]);
    private final ThreadLocal<double[][]> out = ThreadLocal
            .withInitial(() -> new double[4][4]);
    private final ThreadLocal<double[][]> xrot = ThreadLocal
            .withInitial(() -> new double[4][4]);
    private final ThreadLocal<double[][]> yrot = ThreadLocal
            .withInitial(() -> new double[4][4]);
    private final ThreadLocal<double[][]> zrot = ThreadLocal
            .withInitial(() -> new double[4][4]);
    private final ThreadLocal<double[][]> product1 = ThreadLocal
            .withInitial(() -> new double[4][4]);
    private final ThreadLocal<double[][]> product2 = ThreadLocal
            .withInitial(() -> new double[4][4]);
    private final ThreadLocal<double[][]> transform = ThreadLocal
            .withInitial(() -> new double[4][4]);
    // Gape will need to process hydrogen bonding types, but other apps may not.
    protected boolean findFeatures = true;
    private volatile Taff taff;
    private volatile GaSupervisor problem;
    private volatile boolean pairsToCheck[][], useGray, flipFreeCorners, rigid = false,
            fixed = false, ignoreVdwAttractive = false, ignoreTorsion = false,
            randomize = true, outputRandomized = false, flipAmideBonds = false,
            relaxMolecule = false, breakLargeCycles = true;
    private volatile GaussianList atomicGaussians;
    private volatile double activity, weight, relaxMaxDistance = 0.5,
            relaxMaxAngle = 5.0, startingEnergy, randomizedEnergy, baseTorsionEnergy,
            baseVdwEnergy, conformationalEnergy;
    private volatile List<Feature> allFeatures;
    private volatile List<Feature> allMappingFeatures;

    /**
     * Empty constructor
     */
    public GaMolecule() {
    }

    public GaMolecule(String name, List<Atom> atoms, List<Bond> bonds,
                      List<double[]> coords) {
        super(name, atoms, bonds, coords);
    }

    /**
     * Creates molecule for us in GA algorithm
     *
     * @param p algorithm handle
     */
    GaMolecule(GaSupervisor p) {
        this();
        problem = p;
        infoMessageLogger = p.getInfoMessageLogger();
    }

    /**
     * Pairwise acceptor atom score
     *
     * @param mol
     * @return
     * @see #featureOverlay(int, GaMolecule)
     */
    public double acceptorAtomScore(GaMolecule mol) {
        return featureOverlay(FeatureType.ACCEPTOR_ATOM, mol);
    }

    public void addPharmDescription(String pharmDescription) {
        pharmDescriptions.add(pharmDescription);
    }

    /**
     * Return full torsion energy for molecule, includes torsions for rigid
     * portions of the molecule.
     *
     * @return
     */
    public double allTorsionEnergy() {
        if (ignoreTorsion)
            return .0;

        double eTor = getTorsions().stream().map(t -> taff.torsionEnergy(t))
                .reduce(0.0, (sum, energy) -> sum + energy);
        logger.debug(" Torsion energy " + eTor);
        return eTor;
    }

    /**
     * Pairwise aromatic ring score.
     *
     * @param mol
     * @return
     * @see #featureOverlay(int, GaMolecule)
     */
    public double aromaticRingScore(GaMolecule mol) {
        return featureOverlay(FeatureType.AROMATIC_RING, mol);
    }

    /**
     * Returns conformational energy (Torsional + VDW) for a molecule. Only
     * considers energy for flexible part of molecule.
     *
     * @return
     */
    public double conformationalEnergy() {

        double eTor = baseTorsionEnergy + torsionEnergy();
        double eVdw = baseVdwEnergy + vdwEnergy();
        // double eTor = torsionEnergy();
        // double eVdw = vdwEnergy();

        logger.debug("eTor " + eTor + " eVdw " + eVdw);
        // printTorsionalEnergy();
        conformationalEnergy = eTor + eVdw;
        return conformationalEnergy;
    }

    /**
     * Determines (squared) distances outside of constraint. Currently,
     * constraints are only implemented for large cycles.
     *
     * @return
     */
    public double constraintDistance() {
        double score = largeCycles.stream().map(largeCycle -> largeCycle.penalty())
                .reduce(0.0, (sum, penalty) -> sum + penalty);
        return score;
    }

    /**
     * Reset to reference conformation
     */
    public void copyReferenceCorordinates() {
        for (int i = 0; i < getnAtoms(); i++) {
            Coord.copy(reference.get(i), getCoord(i));
        }
    }

    /**
     * Counts the number of Hydrogens in a molecule. Used to check that and SD
     * file contains all hydrogens.
     */
    public int countHydrogens() {
        int hcnt = 0;
        for (Atom atom : getAtoms()) {
            if (atom.getAtomType() == AtomType.Type.H) {
                hcnt++;
            }
        }
        return hcnt;
    }

    /**
     * Pairwise donor hydrogen score
     *
     * @param mol
     * @return
     * @see #featureOverlay(int, GaMolecule)
     */
    public double donorHydrogenScore(GaMolecule mol) {
        return featureOverlay(FeatureType.DONOR_INTERACTION_POINT, mol);
    }

    /**
     * Pairwise feature score between two molecules. Note this routine does not
     * include pharmacophore scaling as required by GAPE. For this more
     * complicated overlay score see SuperpositionChromosome.
     *
     * @param featureNo
     * @param mol
     * @return
     * @see SuperpositionChromosome#featureOverlay(int)
     */
    public double featureOverlay(FeatureType featureNo, GaMolecule mol) {
        double score = .0;
        for (Feature feature : featureMappings.get(featureNo).getFeatures()) {
            for (Feature otherFeature : mol.featureMappings.get(featureNo).getFeatures()) {
                score += feature.score(otherFeature);
            }
        }

        return score;
    }

    /**
     * Finds all the features in a molecule. Fills the featureMappings,
     * allFeatures and allMappingFeatures.
     */
    public void findFeatures() {

        // Fill feature mapping classes
        featureMappings.clear();

        List<Feature> v = HydrophobicFeature.findFeatures(this);
        featureMappings.put(FeatureType.HYDROPHOBIC_ATOM, new FeatureMapping(
                FeatureType.HYDROPHOBIC_ATOM, v));

        v = DonorHydrogenFeature.findFeatures(this);
        featureMappings.put(FeatureType.DONOR_INTERACTION_POINT, new FeatureMapping(
                FeatureType.DONOR_INTERACTION_POINT, v));

        v = AcceptorAtomFeature.findFeatures(this);
        featureMappings.put(FeatureType.ACCEPTOR_ATOM, new FeatureMapping(
                FeatureType.ACCEPTOR_ATOM, v));

        v = AromaticRingFeature.findFeatures(this);
        featureMappings.put(FeatureType.AROMATIC_RING, new FeatureMapping(
                FeatureType.AROMATIC_RING, v));

        // user features
        if (problem != null) {
            for (int i = 0; i < problem.getnUserFeatureTypes(); i++) {
                FeatureType featureType = Feature.userFeatureType(i);
                UserFeatureSet userFeatureSet = problem.getUserFeatureSet(featureType);
                v = userFeatureSet.findUserFeatures(this);
                FeatureType type = FeatureType.valueOf("USER_FEATURES"
                        + String.valueOf(i + 1));
                // int id = userFeatureSet.getFeatureSetNo();
                featureMappings.put(type, new FeatureMapping(type, v));
            }
        }

        // all features
        allFeatures = featureMappings.values().stream().map(FeatureMapping::getFeatures)
                .reduce(new ArrayList<Feature>(), (list, fs) -> {
                    list.addAll(fs);
                    return list;
                });

        allMappingFeatures = allFeatures.stream().filter(Feature::isMappingFeature)
                .collect(Collectors.toList());

    }

    /**
     * Finds all the free corners in this molecule. Sets up the freeCorners
     * array
     */
    public void findFreeCorners() {
        infoMessageLogger.infoMessageln(3, "Finding Free Corners");
        freeCorners.clear();
        freeCornerTorsions.clear();
        for (Atom atom : getAtoms()) {
            logger.debug("checking atom " + atom.info());
            FreeCorner f = new FreeCorner(atom, this);
            if (f.isFreeCorner()) {
                freeCorners.add(f);
                infoMessageLogger.infoMessageln(3, atom.info());
                f.addTorsions(freeCornerTorsions, this);
            }
        }
    }

    /**
     * Matches features to pharmacophore descriptions.
     */
    public void findPharmFeatures() {
        if (pharmDescriptions == null) {
            pharmFeatures.clear();
            return;
        }

        pharmFeatures.clear();

        for (Feature feature : allFeatures) {
            String label = feature.pharmLabel();
            for (String pharmDescription : pharmDescriptions) {
                if (pharmDescription.startsWith(label)) {
                    pharmFeatures.add(feature);
                    feature.setPharmPoint(true);
                }
            }
        }

        if (pharmFeatures.size() != pharmDescriptions.size()) {
            throw new RuntimeException("Unable to find all pharmacophore points");
        }
    }

    /**
     * Decodes an integer string chromosome and builds an overlay. Set remap for
     * a new chromosome. Duplicate mappings mean that this routine may not be
     * reproducible for existing chromsomes. Returns true if we can build an
     * overlay.
     *
     * @param c
     * @param start
     * @param baseMolecule
     * @param remap
     * @return
     */
    public boolean fitMolecule(IntegerStringChromosome c, int start,
                               GaMolecule baseMolecule, boolean remap) {
        int nPoints = baseMolecule.allMappingFeatures.size();
        double xPoints[][] = new double[nPoints][];
        double yPoints[][] = new double[nPoints][];

        List<Feature> baseFeatures = baseMolecule.allMappingFeatures;
        Map<FeatureType, FeatureMapping> baseFeatureMappings = baseMolecule.featureMappings;

        boolean debugLog = logger.isDebugEnabled();
        if (debugLog)
            logger.debug("fitting chromosome:");
        int pos = 0;
        int nMapped = 0;

        for (FeatureMapping baseFeatureMapping : baseFeatureMappings.values()) {
            if (!baseFeatureMapping.isMapping())
                continue;
            for (int k = 0; k < baseFeatureMapping.getNFeatures(); k++) {
                int val = c.getValues(start + pos);

                if (debugLog)
                    logger.debug(val + " ");

                Feature baseFeature = baseFeatures.get(pos);
                if (val != -1) {
                    Feature otherFeature = featureMappings.get(
                            baseFeatureMapping.getFeatureType()).getFeatures(val);
                    xPoints[pos] = otherFeature.calculateCoordinate();
                    yPoints[pos] = baseFeature.calculateCoordinate();
                    baseFeature.setOtherFeature(otherFeature);
                    baseFeature.setMapped(true);
                    nMapped++;
                } else
                    baseFeature.setMapped(false);

                pos++;
            }
        }

        if (debugLog) {
            logger.debug("Pre First Pass Mapped features");
            baseFeatures.stream().filter(f -> f != null && f.isMapped())
                    .forEachOrdered((f) -> {
                        f.calculateSqrDist();
                        logger.debug("Mapped feature " + f.mappingInfo());
                    });
        }

        if (nMapped < 3) {
            if (debugLog)
                logger.debug("First Pass less than 3 fitting points\n");
            return false;
        }

        double trans[][] = new double[4][4];
        Coord.leastSquaresFit(xPoints, yPoints, trans);
        // In some cases the SVD fails- if there are only a few points that
        // are linear.

        if (remap) {
            // 2nd pass fitting

            nMapped = 0;
            BestThree bestThree = new BestThree();

            // first look for up to three points within PASS1_DISTANCE
            for (int i = 0; i < xPoints.length; i++) {
                if (xPoints[i] == null)
                    continue;

                Coord.transPointInPlace(trans, xPoints[i]);

                double sqrDist = Coord.sqrDistance(xPoints[i], yPoints[i]);
                Feature baseFeature = baseFeatures.get(i);
                baseFeature.setSqrDist(sqrDist);
                bestThree.checkBest(i, sqrDist);
                if (sqrDist < PASS1_DISTANCE * PASS1_DISTANCE) {
                    nMapped++;
                    baseFeature.setMapped(true);
                } else {
                    baseFeature.setMapped(false);
                }
            }

            // If we can't find three points then use the best three points
            // within MAX_PASS_DISTANCE
            if (nMapped < 3) {
                if (debugLog)
                    logger.debug("Best Three " + bestThree.firstNo + " "
                            + bestThree.secondNo + " " + bestThree.thirdNo);
                // If the third point is too far away then the mapping has
                // failed
                if (bestThree.thirdSqrDist > MAX_PASS_DISTANCE * MAX_PASS_DISTANCE) {
                    if (debugLog)
                        logger.debug("fitMolecule 2nd pass failed");
                    return false;
                }
                baseFeatures.get(bestThree.firstNo).setMapped(true);
                baseFeatures.get(bestThree.secondNo).setMapped(true);
                baseFeatures.get(bestThree.thirdNo).setMapped(true);
            }

            if (debugLog) {
                logger.debug("Mapped features");
                baseFeatures.stream().filter(f -> f != null && f.isMapped())
                        .forEachOrdered((f) -> {
                            f.calculateSqrDist();
                            logger.debug("Mapped feature " + f.mappingInfo());
                        });
            }

            nMapped = 0;
            for (int i = 0; i < nPoints; i++) {
                if (xPoints[i] == null)
                    continue;
                Feature baseFeature = baseFeatures.get(i);
                if (baseFeature.isMapped()) {
                    nMapped++;
                    xPoints[i] = baseFeature.getOtherFeature().calculateCoordinate();
                } else {
                    xPoints[i] = null;
                    yPoints[i] = null;
                }
            }

            // second fitting
            Coord.leastSquaresFit(xPoints, yPoints, trans);

            for (int i = 0; i < xPoints.length; i++) {
                if (xPoints[i] == null)
                    continue;

                Coord.transPointInPlace(trans, xPoints[i]);

                double sqrDist = Coord.sqrDistance(xPoints[i], yPoints[i]);
                if (sqrDist > MAX_PASS_DISTANCE * MAX_PASS_DISTANCE) {
                    if (debugLog)
                        logger.debug("fitMolecule 2nd pass refitting failed");
                    return false;
                }
                if (debugLog)
                    logger.debug("Mapped point " + i + " sqrDist " + sqrDist);
            }

        }

        for (double[] coord : getCoords()) {
            Coord.transPointInPlace(trans, coord);
        }

        // Remap the chromosome
        if (remap) {
            for (int i = 0; i < baseFeatures.size(); i++)
                if (!baseFeatures.get(i).isMapped())
                    c.setValues(start + i, -1);
        }

        if (debugLog) {
            logger.debug("Remapped Chromosome");
            for (int i = 0; i < baseFeatures.size(); i++) {
                int val = c.getValues(start + i);
                logger.debug("V" + val + " ");
                Feature baseFeature = baseFeatures.get(i);
                if (baseFeature != null && baseFeature.isMapped()) {
                    baseFeature.calculateCoordinate();
                    baseFeature.getOtherFeature().calculateCoordinate();
                    baseFeature.calculateSqrDist();
                    logger.debug(baseFeature.mappingInfo());
                }
            }
        }

        return true;
    }

    /**
     * Frees up extra space used by algorithms
     */
    @Override
    public void free() {
        taff = null;
        // if (!useMillsAndDeanTypes) {
        // DonorAcceptorType.freeMoleculeDonorAcceptors(this);
        // } else {
        HydrogenBondingType.freeMolecule(this);
        // }
        freeFreeCorners();
        pairsToCheck = null;
        freeFeatures();
        reference.clear();
        atomicGaussians = null;
        super.free();
    }

    /**
     * Frees up space associated with free corners
     */
    public void freeFreeCorners() {
        freeCornerTorsions.clear();
        ;
    }

    /**
     * Returns common molecular volume using atomic gaussians.
     *
     * @param mol
     * @return
     */
    public double gaussianIntegral(GaMolecule mol) {

        // atomicGaussians must be current
        GaussianList atoms = atomicGaussians;
        GaussianList otherAtoms = mol.atomicGaussians;

        double volume = atoms.overlapVolume(otherAtoms);
        // if (debugLog)
        // logger.debug("1st order Gaussian Integral " + volume);

        if (true)
            return volume;

        // haven't figured out how to get all the higher order contribs yet!!

        boolean debugLog = logger.isDebugEnabled();
        GaussianList overlay = atoms.intersection(otherAtoms);
        volume = overlay.totalVolume();
        // Subtract 2 order intersections
        overlay = overlay.intersection();
        double diff = overlay.totalVolume();
        volume -= diff;
        if (debugLog)
            logger.debug("2nd order Gaussian Integral " + diff + " new vol " + volume);

        // Add 3 order intersections
        overlay = overlay.intersection();
        diff = overlay.totalVolume();
        volume += diff;
        if (debugLog)
            logger.debug("3rd order Gaussian Integral " + diff + " new vol " + volume);

        if (debugLog)
            logger.debug("Gaussian Integral " + volume);
        return volume;
    }

    /**
     * Given a binary chromosome encoding torsions, decodes to produce a
     * molecular conformation. In this case the chromosome encodes only one
     * mnolecule.
     *
     * @param c
     * @see BinaryStringChromosome
     */
    public void generateConformation(BinaryStringChromosome c) {
        generateConformation(c, 0);
    }

    /**
     * Given a binary chromosome encoding torsions, decodes to produce a
     * molecular conformation.
     *
     * @param c
     * @param start Offset to torsions in the chromosome (used if the chromosome
     *              encode for a number of molecules.
     * @see BinaryStringChromosome
     */
    public void generateConformation(BinaryStringChromosome c, int start) {

        // generate new conformation
        int no = 0;
        for (RotatableBond rb : getRotatableBonds()) {
            int bVal = 0;
            if (useGray)
                bVal = c.posToGray(no * 8 + start);
            else
                bVal = c.posToInt(no * 8 + start);
            double angle = rb.rotateBond(bVal);
            if (CHECK_TORSION_MATCHES_PREDICTED) {
                double a = angle * 180 / Math.PI;
                logger.debug("Angle " + no + " " + a);
                if (rb.getTordistInfo() != null) {
                    double predict = rb.getTordistInfo().getTorAngle() * 180 / Math.PI;
                    double real = rb.getTordistInfo().getTorAngle() * 180 / Math.PI;
                    double ref = rb.getTordistInfo().getRefAngle() * 180 / Math.PI;
                    logger.trace("Reverse " + rb.getTordistInfo().isReverse());
                    logger.trace(" Predict " + predict + " real " + real + " reference "
                            + ref);
                    double diff = real - predict;
                    if (diff > 180)
                        diff -= 360;
                    if (diff < -180)
                        diff += 360;
                    if (Math.abs(diff) > 1.0e-5)
                        throw new IllegalStateException(
                                "Predicted Torsion does not match real");
                }
            }

            no++;
        }

        int pos = start + getnRotatableBonds() * 8;
        no = 0;
        for (FreeCorner freeCorner : freeCorners) {
            if (c.bitSet(pos + no)) {
                freeCorner.flipCorner();
            }
            no++;
        }
    }

    public double getActivity() {
        return activity;
    }

    public void setActivity(double activity) {
        this.activity = activity;
    }

    /**
     * @param i
     * @return a feature from this molecule
     */
    public Feature getAllFeature(int i) {
        return allFeatures.get(i);
    }

    /**
     * @return all features in this molecule
     */
    public List<Feature> getAllFeatures() {
        return allFeatures;
    }

    /**
     * @param i
     * @return A feature that that can be used in chromosome encoding
     */
    public Feature getAllMappingFeature(int i) {
        return allMappingFeatures.get(i);
    }

    /**
     * @return all features from this molecule that can be used in a chromsome
     * encoding
     */
    public List<Feature> getAllMappingFeatures() {
        return allMappingFeatures;
    }

    /**
     * Generates atomic gaussian for use in gaussian volume score. You need to
     * recompute them for each overlay. The gaussian overlay is performed only
     * on hydrophobic atoms.
     *
     * @return
     */
    public GaussianList getAtomicGaussians() {
        // atomicGaussians = GaussianList.atomicGaussians(this);
        atomicGaussians = GaussianList.hydrophobicAtomicGaussians(this);
        return atomicGaussians;
    }

    public double getBaseTorsionEnergy() {
        return baseTorsionEnergy;
    }

    public double getBaseVdwEnergy() {
        return baseVdwEnergy;
    }

    /**
     * Returns the number of bits that are required to encode a moleuclar
     * conformation
     *
     * @return
     */
    public int getConformationalBitLength() {
        int nBonds = getRotatableBonds().size();
        int nCorners = rigid || fixed ? 0 : freeCorners.size();
        int nBits = nBonds * 8 + nCorners;
        return nBits;
    }

    /**
     * @return the conformationalEnergy
     */
    public double getConformationalEnergy() {
        return conformationalEnergy;
    }

    /**
     * @param conformationalEnergy the conformationalEnergy to set
     */
    public void setConformationalEnergy(double conformationalEnergy) {
        this.conformationalEnergy = conformationalEnergy;
    }

    /**
     * @return the featureMappings
     */
    public Map<FeatureType, FeatureMapping> getFeatureMappings() {
        return featureMappings;
    }

    /**
     * @return the featureMappings
     */
    public FeatureMapping getFeatureMappings(FeatureType type) {
        return featureMappings.get(type);
    }

    public FreeCorner getFreeCorner(int i) {
        return freeCorners.get(i);
    }

    public List<FreeCorner> getFreeCorners() {
        return freeCorners;
    }

    public List<LargeCycle> getLargeCycles() {
        return largeCycles;
    }

    public void setLargeCycles(List<LargeCycle> largeCycles) {
        this.largeCycles.clear();
        this.largeCycles.addAll(largeCycles);
    }

    public int getNFeatures() {
        return allFeatures.size();
    }

    public int getNFreeCorners() {
        return freeCorners.size();
    }

    /**
     * @return the number of features that can be used in chromsome encoding
     */
    public int getNMappingFeatures() {
        return allMappingFeatures.size();
    }

    public int getNPharmDescriptions() {
        if (pharmDescriptions == null)
            return 0;
        return pharmDescriptions.size();
    }

    /**
     * @return number of features found from reading in pharm descriptions
     * @see #findPharmFeatures()
     */
    public int getNPharmFeatures() {
        return pharmFeatures.size();
    }

    public int getNPharmPoints() {
        int nPharmPoints = 0;
        for (Feature feature : allFeatures)
            if (feature.isPharmPoint())
                nPharmPoints++;
        return nPharmPoints;
    }

    public String getPharmDescription(int i) {
        return pharmDescriptions.get(i);
    }

    public List<String> getPharmDescriptions() {
        return pharmDescriptions;
    }

    public void setPharmDescriptions(List<String> pharmDescriptions) {
        this.pharmDescriptions.clear();
        if (pharmDescriptions != null)
            this.pharmDescriptions.addAll(pharmDescriptions);
    }

    /**
     * @param i
     * @return pharmacophore feature found from reading in a molecule with
     * pharmacophore descriptions
     * @see #findPharmFeatures()
     */
    public Feature getPharmFeature(int i) {
        return pharmFeatures.get(i);
    }

    /**
     * @return pharmacophore features found from reating in a molecule with
     * pharmacophore descriptions
     * @see #findPharmFeatures()
     */
    public List<Feature> getPharmFeatures() {
        return pharmFeatures;
    }

    public GaSupervisor getProblem() {
        return problem;
    }

    /**
     * Sets GA problem instance/algorithm handle
     *
     * @param p
     */
    public void setProblem(GaSupervisor p) {
        problem = p;
        infoMessageLogger = p.getInfoMessageLogger();
    }

    public List<double[]> getReference() {
        return reference;
    }

    public void setReference(List<double[]> reference) {
        this.reference.clear();
        this.reference.addAll(reference);
    }

    public double[] getReference(int i) {
        return reference.get(i);
    }

    public double getWeight() {
        return weight;
    }

    public void setWeight(double weight) {
        this.weight = weight;
    }

    /**
     * @return ture if the setup routines will try and set up hydrogen bonding
     * features.
     */
    public boolean isFindFeatures() {
        return findFeatures;
    }

    /**
     * Call before setup to determine if features will be found for the
     * molecule.
     *
     * @param findFeatures
     */
    public void setFindFeatures(boolean findFeatures) {
        this.findFeatures = findFeatures;
    }

    public boolean isFixed() {
        return fixed;
    }

    public void setFixed(boolean fixed) {
        this.fixed = fixed;
    }

    /**
     * @return the ignoreTorsion
     */
    public boolean isIgnoreTorsion() {
        return ignoreTorsion;
    }

    /**
     * @param ignoreTorsion the ignoreTorsion to set
     */
    protected void setIgnoreTorsion(boolean ignoreTorsion) {
        this.ignoreTorsion = ignoreTorsion;
    }

    public boolean isRandomize() {
        return randomize;
    }

    public void setRandomize(boolean randomize) {
        this.randomize = randomize;
    }

    public boolean isRelaxMolecule() {
        return relaxMolecule;
    }

    /**
     * @param relaxMolecule the relaxMolecule to set
     */
    protected void setRelaxMolecule(boolean relaxMolecule) {
        this.relaxMolecule = relaxMolecule;
    }

    public boolean isRigid() {
        return rigid;
    }

    public void setRigid(boolean rigid) {
        this.rigid = rigid;
    }

    /**
     * @return the useGray
     */
    public boolean isUseGray() {
        return useGray;
    }

    /**
     * @param useGray the useGray to set
     */
    protected void setUseGray(boolean useGray) {
        this.useGray = useGray;
    }

    /**
     * Overrides {@link com.cairn.intranet.molecue.Molecule#loadSdfMol}.
     * Calls that method then looks for an sdfield with pharmacophore
     * information.
     *
     * @param in Reader to MOL2 file.
     */
    @Override
    public void loadSdfMol(BufferedReader in) throws IOException, MolReadException {
        super.loadSdfMol(in);

        if (hasSdfField("PHARMACOPHORE")) {
            String pharm = getSdfField("PHARMACOPHORE");
            pharmDescriptions.clear();
            pharmDescriptions.addAll(Arrays.asList(pharm.split("\n")));
            // Remove from sdfields
            removeSdfField("PHARMACOPHORE");
        }
    }

    /**
     * Overrides {@link com.cairn.intranet.molecue.Molecule#loadSybylMol2}.
     * Calls that method then looks for a tripos comments section with
     * pharmacophore information.
     *
     * @param in Reader to MOL2 file.
     */
    @Override
    public void loadSybylMol2(BufferedReader in) throws IOException, MolReadException {
        super.loadSybylMol2(in);

        boolean inComments = false, inPharm = false;

        in.mark(16384);
        while (true) {
            String val = in.readLine();
            if (val == null)
                return;
            if (val.startsWith("#"))
                return;
            if (val.startsWith("@<TRIPOS>COMMENT")) {
                inComments = true;
                break;
            }
            if (val.startsWith("@<TRIPOS>MOLECULE")) {
                in.reset();
                return;
            }
        }

        if (!inComments)
            return;
        while (true) {
            String val = in.readLine();
            if (val == null)
                return;
            if (val.startsWith("#"))
                return;
            if (val.startsWith("START_PHARMACOPHORE")) {
                inPharm = true;
                break;
            }
        }

        if (!inPharm)
            return;
        pharmDescriptions.clear();
        while (true) {
            String val = in.readLine();
            if (val == null)
                return;
            if (val.startsWith("#"))
                return;
            if (val.startsWith("END_PHARMACOPHORE"))
                break;
            pharmDescriptions.add(val);
        }

    }

    /**
     * Lists torsions and energy
     */
    public void printTorsionalEnergy() {
        getTorsions().stream().forEachOrdered(
                t -> {
                    infoMessageLogger.infoMessageln("Torsion " + t.info() + " Energy "
                            + t.getEnergy());
                });
    }

    /**
     * Randomizes molecule by applying random rotations to each rotatable bond
     */
    public void randomizeMol() {
        getRotatableBonds().stream().forEachOrdered(rb -> {
            int bVal = problem.randomInt(0, 255);
            rb.rotateBond(bVal);
        });
    }

    /**
     * Purturbs the position of a molecule based on information in a binary
     * chromosome.
     *
     * @param c
     * @param start positon of the purturbation encoding on the chromosome.
     */
    public void relax(BinaryStringChromosome c, int start) {
        if (!relaxMolecule)
            return;

        double r = byteToNormal(c, start) * relaxMaxDistance;
        double polar = byteToNormal(c, start + 8) * 2 * Math.PI;
        double azimuth = byteToNormal(c, start + 16) * 2 * Math.PI;
        double scale = relaxMaxAngle * Math.PI / 180.0;
        double rot1 = byteToNormal(c, start + 24) * scale;
        double rot2 = byteToNormal(c, start + 32) * scale;
        double rot3 = byteToNormal(c, start + 40) * scale;

        // Polar co-ordinates for ligand
        double rz = r * Math.cos(polar);
        double ry = r * Math.sin(polar) * Math.sin(azimuth);
        double rx = r * Math.sin(polar) * Math.cos(azimuth);

        int nAtoms = getnAtoms();
        double cx = 0, cy = 0, cz = 0;
        for (int i = 0; i < nAtoms; i++) {
            double[] coord = getCoord(i);
            cx += coord[0];
            cy += coord[1];
            cz += coord[2];
        }
        cx = cx / nAtoms;
        cy = cx / nAtoms;
        cz = cx / nAtoms;

        double[][] _in = in.get();
        double[][] _out = out.get();
        double[][] _xrot = xrot.get();
        double[][] _yrot = yrot.get();
        double[][] _zrot = zrot.get();

        Coord.identity(_in);
        Coord.identity(_out);
        Coord.identity(_xrot);
        Coord.identity(_yrot);
        Coord.identity(_zrot);

        _in[3][0] = -cx;
        _in[3][1] = -cy;
        _in[3][2] = -cz;
        _out[3][0] = cx + rx;
        _out[3][1] = cy + ry;
        _out[3][2] = cz + rz;

        double cosrot1 = Math.cos(rot1);
        double cosrot2 = Math.cos(rot2);
        double cosrot3 = Math.cos(rot3);
        double sinrot1 = Math.sin(rot1);
        double sinrot2 = Math.sin(rot2);
        double sinrot3 = Math.sin(rot3);

        _xrot[1][1] = _xrot[2][2] = cosrot1;
        _xrot[1][2] = -sinrot1;
        _xrot[2][1] = sinrot1;

        _yrot[0][0] = _yrot[2][2] = cosrot2;
        _yrot[0][2] = -sinrot2;
        _yrot[2][0] = sinrot2;

        _zrot[0][0] = _zrot[1][1] = cosrot3;
        _zrot[0][1] = sinrot3;
        _zrot[1][0] = -sinrot3;

        double[][] _product1 = product1.get();
        double[][] _product2 = product2.get();
        double[][] _transform = transform.get();

        Coord.product(_in, _xrot, _product1);
        Coord.product(_product1, _yrot, _product2);
        Coord.product(_product2, _zrot, _product1);
        Coord.product(_product1, _out, _transform);

        // writeSybylMol2File(name+"_pre_relax.mol2", null);
        for (double[] coord : getCoords()) {
            Coord.transPointInPlace(_transform, coord);
        }
        // writeSybylMol2File(name+"_post_relax.mol2", null);
    }

    public void setReference(int i, double[] reference) {
        this.reference.set(i, reference);
    }

    /**
     * Sets molecule parameters using settings from program instance. Also does
     * stuff like finding starting energy and randomizing.
     */
    public void setup() {
        useGray = problem.getBooleanValue("graycode");
        flipFreeCorners = problem.getBooleanValue("flip_free_corners");
        if (problem.hasKey("break_large_cycles"))
            breakLargeCycles = problem.getBooleanValue("break_large_cycles");

        if (infoMessageLogger.getLogLevel() > 3) {
            if (rigid)
                infoMessageLogger.infoMessageln("Rigid molecule");
            if (fixed)
                infoMessageLogger.infoMessageln("Fixed molecule");
            if (useGray)
                infoMessageLogger.infoMessageln("Using Gray encoding");
            else
                infoMessageLogger.infoMessageln("Using binary encoding");
            if (flipFreeCorners)
                infoMessageLogger.infoMessageln("Flipping ring corners");
            else
                infoMessageLogger.infoMessageln("Not flipping ring corners");
            if (breakLargeCycles)
                infoMessageLogger.infoMessageln("Breaking large cycles");
        }

        if (countHydrogens() == 0)
            throw new RuntimeException("SD Molecule " + getName() + " has no hydrogens\n"
                    + "GAPE requires that SD files contain all "
                    + "hydrogens or fill_valence is set");

        // handle large cycles
        if (!rigid && breakLargeCycles) {
            LargeCycle.findlargeCycles(this);
            update();
        }

        // Create forcefield before adding lone pairs, as we don't
        // want them to contribute
        taff = new Taff(this);

        if (findFeatures) {
            infoMessageLogger.infoMessageln(3,
                    "\nFinding Dean and Mills Donors and Acceptors");
            HydrogenBondingType.searchMolecule(this, infoMessageLogger);
            infoMessageLogger.infoMessageln(3, "");

            // We do want lone pairs to rotate if possible so add them
            // before assigning rotatable bonds
            double hBondLen = problem.getDoubleValue("h_bond_len");
            LonePairAddition.addLonePairs(this, hBondLen);
        }

        if (problem.hasKey("flip_amide_bonds")) {
            flipAmideBonds = problem.getBooleanValue("flip_amide_bonds");
        }
        if (flipAmideBonds)
            infoMessageLogger.infoMessageln(3, "Flipping Amide Bonds");

        // Find Rotatable Bonds
        if (!rigid && !fixed) {
            assignRotatableBonds(flipAmideBonds);
        } else {
            // starting with all rotatable bonds
            assignRotatableBonds(false);
            // only keep those for terminal donors or acceptors.
            infoMessageLogger
                    .infoMessageln(3,
                            "Retaining only terminal rotatable bonds to donor hydrogens or lone pairs");
            List<RotatableBond> terminalRotatableBonds = getRotatableBonds().stream()
                    .filter(rb -> isTerminalRotatable(rb)).collect(Collectors.toList());
            setRotatableBonds(terminalRotatableBonds);
            infoMessageLogger.infoMessageln(2, "Found " + terminalRotatableBonds.size()
                    + " rotatable bonds");
        }

        if (problem.hasKey("ignore_vdw_attractive")) {
            ignoreVdwAttractive = problem.getBooleanValue("ignore_vdw_attractive");
            taff.setIgnoreVdwAttractive(ignoreVdwAttractive);
            if (ignoreVdwAttractive)
                infoMessageLogger.infoMessageln(3, "Ignoring attractive VDW");
        }

        if (problem.hasKey("ignore_torsion")) {
            ignoreTorsion = problem.getBooleanValue("ignore_torsion");
            if (ignoreTorsion)
                infoMessageLogger.infoMessageln(3, "Ignoring Torsion energies");
        }

        if (!rigid && !fixed && problem.hasKey("flatten_bonds")
                && problem.getBooleanValue("flatten_bonds")) {
            infoMessageLogger.infoMessageln(3, "Flattening Bonds");
            flattenBonds();
        }

        if (findFeatures) {
            infoMessageLogger.infoMessageln(2, "Finding Features");
            findFeatures();
        }

        if (!rigid && !fixed && flipFreeCorners)
            findFreeCorners();

        startingEnergy = taff.molEnergy();

        if (logger.isDebugEnabled()) {
            writeSybylMol2File("Starting " + getName() + ".mol2", "Starting GA Molecule");
        }

        infoMessageLogger.infoMessageln(3, "Starting Energy " + startingEnergy);

        if (randomize) {
            randomizeMol();
            randomizedEnergy = taff.molEnergy();
            int no = 1;
            while (!rigid && randomizedEnergy > 1e7) {
                if (no > 100)
                    throw new RuntimeException(
                            "Failed to Randomize molecule to reasonable energy conformer");
                randomizeMol();
                randomizedEnergy = taff.molEnergy();
                no++;
            }
            infoMessageLogger.infoMessageln(3, "Randomized Energy " + randomizedEnergy);
            if (problem.hasKey("output_randomized")
                    && problem.getBooleanValue("output_randomized"))
                outputRandomized = true;
            else
                outputRandomized = false;
            if (outputRandomized)
                writeSybylMol2File("Randomized " + getName() + ".mol2",
                        "Randomized GA Molecule");
            else
                infoMessageLogger.infoMessageln(3, "Not outputing randomized structure");
        }

        reference.clear();
        int nAtoms = getnAtoms();
        for (int i = 0; i < nAtoms; i++) {
            double[] ref = new double[4];
            Coord.copy(getCoord(i), ref);
            reference.add(ref);
        }

        pairsToCheck = new boolean[nAtoms][nAtoms];
        getRotatableBonds().stream().forEach(rb -> setupPairs(rb));
        freeCorners.stream().forEach(freeCorner -> {
            setupPairs(freeCorner.getRotationAB());
            setupPairs(freeCorner.getRotationCD());
            setupPairs(freeCorner.getRotationXC());
        });
        largeCycles.stream().forEach(largeCycle -> {
            // Don't check pairs in LargeCycles breaks
            int no1 = largeCycle.getAtom1().getNo();
            int no2 = largeCycle.getAtom2().getNo();
            pairsToCheck[no1][no2] = pairsToCheck[no2][no1] = false;
        });

        baseTorsionEnergy = allTorsionEnergy() - torsionEnergy();
        baseVdwEnergy = taff.moleculeVdwEnergy() - vdwEnergy();

        if (problem.hasKey("relax_molecule")) {
            relaxMolecule = problem.getBooleanValue("relax_molecule");
            relaxMaxDistance = problem.getDoubleValue("relax_max_distance");
            relaxMaxAngle = problem.getDoubleValue("relax_max_angle");
        }

        if (relaxMolecule) {
            infoMessageLogger.infoMessageln(3, "Relaxing molecule after fitting");
            infoMessageLogger.infoMessageln(3, "Max distance " + relaxMaxDistance
                    + " Max angle " + relaxMaxAngle);
        }
    }

    /**
     * Determines the torsional energy for this molecule. Only determines energy
     * for torsions affected by a rotatable bond or free corner.
     *
     * @return
     */
    public double torsionEnergy() {
        if (ignoreTorsion)
            return .0;

        double rotatableTorsionEnergy = getRotatableBonds().stream()
                .map(rb -> rb.torsionalEnergy(taff))
                .reduce(0.0, (sum, energy) -> sum + energy);
        double freeCornerTorsionEnergy = freeCornerTorsions.stream()
                .map(t -> taff.torsionEnergy(t))
                .reduce(0.0, (sum, energy) -> sum + energy);

        double eTorsion = rotatableTorsionEnergy + freeCornerTorsionEnergy;
        logger.debug(" Torsion energy " + eTorsion);
        return eTorsion;
    }

    /**
     * Returns conformational energy (Torsional + VDW) for a molecule. Considers
     * whole molecule- flexible & rigid fractions.
     *
     * @return
     */
    public double totalConformationalEnergy() {
        double eTor = allTorsionEnergy();
        double eVdw = taff.moleculeVdwEnergy();
        // printTorsionalEnergy();
        logger.debug("eTor " + eTor + " eVdw " + eVdw);
        return eTor + eVdw;
    }

    /**
     * Returns common molecular volume based on hard atomic spheres as in GASP.
     *
     * @param mol
     * @return
     */
    public double volumeIntegral(GaMolecule mol) {

        double integral = .0;

        for (Atom atomA : getAtoms()) {
            for (Atom atomB : mol.getAtoms()) {

                if (!atomA.isNotDummy() || !atomB.isNotDummy())
                    continue;

                double[] coordsA = getCoord(atomA.getNo());
                double[] coordsB = mol.getCoord(atomB.getNo());

                double radiusA = atomA.getType().getRadius();
                double radiusB = atomB.getType().getRadius();

                double srqDistanceAB = Coord.sqrDistance(coordsA, coordsB);
                double test = (radiusA + radiusB) * (radiusA + radiusB);
                if (srqDistanceAB > test)
                    continue;

                double distanceAB = Math.sqrt(srqDistanceAB);

                if ((radiusA - distanceAB) >= radiusB)
                    integral += (4.0 / 3.0 * Math.PI * radiusB * radiusB * radiusB);

                else if ((radiusB - distanceAB) > radiusA)
                    integral += (4.0 / 3.0 * Math.PI * radiusA * radiusA * radiusA);

                else {
                    double xStar = ((radiusA * radiusA) - (radiusB * radiusB) + srqDistanceAB)
                            / (2.0 * distanceAB);

                    integral += (Math.PI
                            / 3.0
                            * ((2.0 * radiusA * radiusA * radiusA)
                            + (2.0 * radiusB * radiusB * radiusB) + (distanceAB * srqDistanceAB)) - Math.PI
                            * ((radiusA * radiusA * xStar) + (distanceAB - xStar)
                            * ((radiusB * radiusB) + (distanceAB * xStar))));
                }
            }
        }

        logger.debug("Volume Integral " + integral);
        return integral;
    }

    /*
     * Calls parent method and, if applicable, adds an sd field with
     * pharmacophore information. Also rebuilds molecule, if required.
     *
     * (non-Javadoc)
     *
     * @see com.cairn.molecule.Molecule#writeSdfMol(java.io.Writer,
     * java.lang.String)
     */
    @Override
    public void writeSdfMol(Writer writer, String info) throws IOException {
        removeSdfField("PHARMACOPHORE");
        if (CollectionUtils.isNotEmpty(pharmDescriptions)) {
            String pharmInfo = pharmDescriptions.stream().collect(
                    Collectors.joining("\n"));
            addSdfField("PHARMACOPHORE", pharmInfo.toString());
        }

        rebuildForWriting();
        super.writeSdfMol(writer, info);
        undoRebuildForWriting();
    }

    /**
     * Overrides {@link com.cairn.intranet.molecue.Molecule#writeSybylMol2}
     * . Calls that method and, if applicable, adds a tripos comments section
     * with pharmacophore information. Also rebuilds molecule, if required.
     *
     * @param out   Writer to MOL2 file.
     * @param info  Comment string.
     * @param incLP Set to include lone pairs.
     */
    @Override
    public void writeSybylMol2(Writer out, String info, boolean incLP) throws IOException {

        // TODO expect molecule name in argument
        rebuildForWriting();
        super.writeSybylMol2(out, info, incLP);
        undoRebuildForWriting();

        if (CollectionUtils.isNotEmpty(pharmDescriptions)) {
            out.write("@<TRIPOS>COMMENT\n");
            out.write("START_PHARMACOPHORE\n");
            for (String pharmDescription : pharmDescriptions) {
                out.write(pharmDescription + "\n");
            }
            out.write("END_PHARMACOPHORE\n");
        }
    }

    /**
     * @return the relaxMaxAngle
     */
    protected double getRelaxMaxAngle() {
        return relaxMaxAngle;
    }

    /**
     * @param relaxMaxAngle the relaxMaxAngle to set
     */
    protected void setRelaxMaxAngle(double relaxMaxAngle) {
        this.relaxMaxAngle = relaxMaxAngle;
    }

    /**
     * @return the relaxMaxDistance
     */
    protected double getRelaxMaxDistance() {
        return relaxMaxDistance;
    }

    /**
     * @param relaxMaxDistance the relaxMaxDistance to set
     */
    protected void setRelaxMaxDistance(double relaxMaxDistance) {
        this.relaxMaxDistance = relaxMaxDistance;
    }

    /**
     * In the multimol version of gape, or for rigid compounds the only bonds we
     * need to rotate are those to allow donor-hydrogens or lone pairs to rotate
     * into position.
     *
     * @param rotatableBond
     * @return true if this
     */
    protected boolean isTerminalRotatable(RotatableBond rotatableBond) {
        Atom atom = rotatableBond.getAtom1();
        boolean acceptor = atom.getAcceptorType() != null;
        boolean donor = atom.getDonorType() != null;

        if (!donor && !acceptor)
            return false;
        for (Atom otherAtom : rotatableBond.getAtom1List()) {
            Bond bond = getBond(atom, otherAtom);
            if (!bond.isTerminal())
                return false;
        }

        infoMessageLogger.infoMessageln(3,
                "Retained rotating bond " + rotatableBond.info());
        return true;
    }

    /**
     * Frees up memory from feature classes and arrays
     */
    void freeFeatures() {
        featureMappings.clear();
        allFeatures = null;
    }

    /**
     * Simple hydrophobic overlay score based on proximity of atoms.
     *
     * @param mol
     * @return
     */
    double hydrophobicProximity(GaMolecule mol) {
        List<Feature> otherFeatures = mol.allFeatures;
        double score = .0;
        for (Feature feature : allFeatures) {
            // for (int i = 0; i < allFeatures.length; i++) {
            if (feature.getFeatureType() != FeatureType.HYDROPHOBIC_ATOM)
                continue;
            for (Feature otherFeature : otherFeatures) {
                // for (int j = 0; j < otherFeatures.length; j++) {
                if (otherFeature.getFeatureType() != FeatureType.HYDROPHOBIC_ATOM)
                    continue;
                double featureCoord[] = feature.getSavedCoordinate();
                double otherFeatureCoord[] = otherFeature.getSavedCoordinate();
                double sqrDist = Coord.sqrDistance(featureCoord, otherFeatureCoord);
                if (sqrDist < PROX_DISTANCE * PROX_DISTANCE) {
                    if (logger.isDebugEnabled()) {
                        int no = feature.getAtom().getNo() + 1;
                        int otherNo = otherFeature.getAtom().getNo() + 1;
                        logger.debug("Proximal hydrophobe " + otherNo + " " + no);
                    }
                    score += 1;
                }
            }
        }
        return score;
    }

    /**
     * Determine which pairs of atoms to check for VDW clashes when rotating a
     * bond.
     *
     * @param rb
     */
    void setupPairs(RotatableBond rb) {
        for (Atom a1 : rb.getAtom1List()) {
            if (a1.getAtomType() == AtomType.Type.LP)
                continue;

            for (Atom a2 : rb.getAtom2List()) {
                // Allow bumps between donor hydrogen and acceptor

                if (a1.isDonorHydrogen() && a2.getAcceptorType() != null)
                    continue;

                if (a2.isDonorHydrogen() && a1.getAcceptorType() != null)
                    continue;

                if (a2.getAtomType() == AtomType.Type.LP)
                    continue;
                if (taff.isIgnore14())
                    if (bonded14(a1, a2))
                        continue;

                if (bonded13(a1, a2))
                    continue;
                pairsToCheck[a1.getNo()][a2.getNo()] = pairsToCheck[a2.getNo()][a1
                        .getNo()] = true;
            }
        }
    }

    /**
     * Determine the VDW energy for the molecule. Does not check all pairs- just
     * those in the pairsToCheck array, so that we only check those atoms
     * separated by a rotatable bond or a free corner.
     *
     * @return
     */
    double vdwEnergy() {
        double eVdw = 0;
        for (int i = 0; i < getnAtoms(); i++) {
            Atom a1 = getAtom(i);
            for (int j = 0; j < i; j++) {
                if (!pairsToCheck[i][j])
                    continue;
                Atom a2 = getAtom(j);
                if (a1.getTaff() == null || a2.getTaff() == null)
                    continue;
                double sqrD = Coord.sqrDistance(getCoord(i), getCoord(j));

                // allow donors and acceptors to approach closer

                if (a1.getDonorType() != null && a2.getAcceptorType() != null)
                    sqrD *= .5;
                if (a2.getDonorType() != null && a1.getAcceptorType() != null)
                    sqrD *= .5;

                double vdw = a1.getTaff().getR() + a2.getTaff().getR();
                sqrD = sqrD / (vdw * vdw);
                double pow6 = sqrD * sqrD * sqrD;
                double pow12 = pow6 * pow6;
                double kMean = Math.sqrt(a1.getTaff().getK() * a2.getTaff().getK());
                double atomVdw = kMean * (1.0 / pow12 - 2.0 / pow6);
                if (ignoreVdwAttractive && atomVdw < .0)
                    atomVdw = 0;
                eVdw += atomVdw;
            }
        }
        return eVdw;
    }

    /**
     * Decodes a byte on a binary string chromosome to normalized double [0, 1).
     *
     * @param c
     * @param pos position of the byte on the binary chromsome
     * @return
     */
    private double byteToNormal(BinaryStringChromosome c, int pos) {
        int bVal = 0;
        if (useGray)
            bVal = c.posToGray(pos);
        else
            bVal = c.posToInt(pos);
        return bVal / 256.0;
    }

    /**
     * Changes the molecule so that we can write out a file for the user to
     * view. At the moment this just reforms any large cycles broken.
     */
    private void rebuildForWriting() {
        for (LargeCycle largeCycle : largeCycles) {
            largeCycle.reformBond();
        }
        update();

    }

    /**
     * Undoes the changes from rebuildForWriting. Currently this just break any
     * defined large cycles.
     */
    private void undoRebuildForWriting() {
        for (LargeCycle largeCycle : largeCycles) {
            largeCycle.restoreBond();
        }
        update();

    }

    /**
     * This class is used to hold the closest three mappings when doing fitting.
     */
    class BestThree {
        double firstSqrDist = 1e100, secondSqrDist = 1e100, thirdSqrDist = 1e100;

        int firstNo, secondNo, thirdNo;

        BestThree() {
        }

        void checkBest(int no, double sqrDist) {
            debug(() -> "checking " + no + " sqrDist " + sqrDist);
            if (sqrDist < firstSqrDist) {
                thirdNo = secondNo;
                thirdSqrDist = secondSqrDist;
                secondNo = firstNo;
                secondSqrDist = firstSqrDist;
                firstNo = no;
                firstSqrDist = sqrDist;
            } else if (sqrDist < secondSqrDist) {
                thirdNo = secondNo;
                thirdSqrDist = secondSqrDist;
                secondNo = no;
                secondSqrDist = sqrDist;
            } else if (sqrDist < thirdSqrDist) {
                thirdNo = no;
                thirdSqrDist = sqrDist;
            }
        }

    }

}