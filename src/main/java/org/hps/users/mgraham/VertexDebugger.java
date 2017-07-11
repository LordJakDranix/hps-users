package org.hps.users.mgraham;

import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;
import hep.physics.matrix.SymmetricMatrix;
import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;

import java.io.IOException;
import static java.lang.Math.sqrt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.hps.analysis.dataquality.DataQualityMonitor;

import org.hps.recon.tracking.TrackUtils;
import org.hps.recon.tracking.TrackerHitUtils;
import org.hps.recon.utils.TrackClusterMatcher;
import org.hps.recon.vertexing.BilliorTrack;
import org.hps.recon.vertexing.BilliorVertex;
import org.hps.recon.vertexing.BilliorVertexer;
import org.lcsim.detector.converter.compact.subdetector.SvtStereoLayer;
import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.event.EventHeader;
import org.lcsim.event.MCParticle;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.Track;
import org.lcsim.event.TrackState;
import org.lcsim.event.base.BaseTrackState;
import org.lcsim.fit.helicaltrack.HelicalTrackFit;
import org.lcsim.fit.helicaltrack.HelixUtils;
import org.lcsim.geometry.Detector;
import org.lcsim.geometry.FieldMap;
import org.lcsim.lcio.LCIOWriter;

/**
 *
 */
public class VertexDebugger extends DataQualityMonitor {

    private static Logger LOGGER = Logger.getLogger(VertexDebugger.class.getPackage().getName());

    String finalStateParticlesColName = "FinalStateParticles";
    String mcParticlesColName = "MCParticle";
    String readoutHitCollectionName = "EcalReadoutHits";//these are in ADC counts
    String calibratedHitCollectionName = "EcalCalHits";//these are in energy
    String clusterCollectionName = "EcalClustersCorr";
    private String notrackFile;
    private String helicalTrackHitCollectionName = "HelicalTrackHits";
    private String rotatedTrackHitCollectionName = "RotatedHelicalTrackHits";
    String[] fpQuantNames = {"nEle_per_Event", "nPos_per_Event", "nPhoton_per_Event", "nUnAssociatedTracks_per_Event", "avg_delX_at_ECal", "avg_delY_at_ECal", "avg_E_Over_P", "avg_mom_beam_elec", "sig_mom_beam_elec"};
    private String outputFile;
    private LCIOWriter writer;
    private LCIOWriter notrackwriter;
    private TrackClusterMatcher matcher = new TrackClusterMatcher();
    //private Map<SiSensor, Map<Integer, Hep3Vector>> stripPositions = new HashMap<SiSensor, Map<Integer, Hep3Vector>>(); 
    private List<HpsSiSensor> sensors = null;
    private Map<Integer, List<SvtStereoLayer>> topStereoLayers = new HashMap<Integer, List<SvtStereoLayer>>();
    private Map<Integer, List<SvtStereoLayer>> bottomStereoLayers = new HashMap<Integer, List<SvtStereoLayer>>();
    // Constants
    public static final double SENSOR_LENGTH = 98.33; // mm
    public static final double SENSOR_WIDTH = 38.3399; // mm
    private static final String SUBDETECTOR_NAME = "Tracker";
    boolean doSkim = false;
    TrackerHitUtils trackerHitUtils = new TrackerHitUtils();
    //some counters
    int nRecoEvents = 0;
    int nTotEle = 0;
    int nTotPos = 0;
    int nTotPhotons = 0;
    int nTotUnAss = 0;
    int nTotAss = 0;
    //some summers
    double sumdelX = 0.0;
    double sumdelY = 0.0;
    double sumEoverP = 0.0;
    private final String plotDir = "FinalStateParticles/";
    // double beamEnergy = 1.05; //GeV
    double clTimeMin = 30;//ns
    double clTimeMax = 50;//ns
    double deltaTimeMax = 4;//ns
    double coplanMean = 180;//degrees
    double coplanWidth = 10;//degrees
    double esumMin = 0.4;
    double esumMax = 1.2;
    double phot_nom_x = 42.52;//nominal photon position (px=0)
    boolean requirePositron = false;
    double maxPairs = 5;
    boolean requireSuperFiducial = false;
    int nbins = 50;
    double B_FIELD = 0.23;//Tesla

    protected double[] beamSize = {0.001, 0.130, 0.050}; //rough estimate from harp scans during engineering run production running
    // Beam position variables.
    // The beamPosition array is in the tracking frame
    /* TODO get the beam position from the conditions db */
    protected double[] beamPosition = {-5.0, 0.0, 0.0}; //

    double minPhi = -0.25;
    double maxPhi = 0.25;

    IHistogram1D elePx;
    IHistogram1D elePy;
    IHistogram1D elePz;
    IHistogram1D elePzBeam;
    IHistogram1D elePzBeamTop;
    IHistogram1D elePzBeamBottom;
    IHistogram1D elePTop;
    IHistogram1D elePBottom;
    IHistogram1D eleClEne;

    IHistogram1D posPx;
    IHistogram1D posPy;
    IHistogram1D posPz;
    IHistogram1D posPTop;
    IHistogram1D posPBottom;
    IHistogram1D posClEne;

    IHistogram2D delElePx;
    IHistogram2D delElePy;
    IHistogram2D delElePz;

    IHistogram2D delPosPx;
    IHistogram2D delPosPy;
    IHistogram2D delPosPz;

    IHistogram2D delZvsVtxZ;
    IHistogram2D delMvsVtxZ;
    IHistogram2D delPvsVtxZ;
    IHistogram2D delMvsVtxZPatch;
    IHistogram2D delMUncvsPatchvsZ;

    IHistogram2D delElePxReco;
    IHistogram2D delElePyReco;
    IHistogram2D delElePzReco;

    IHistogram2D delPosPxReco;
    IHistogram2D delPosPyReco;
    IHistogram2D delPosPzReco;
    IHistogram2D chiSqVtxZ;

    IHistogram2D delElePxShift;
    IHistogram2D delElePyShift;
    IHistogram2D delElePzShift;

    IHistogram2D delPosPxShift;
    IHistogram2D delPosPyShift;
    IHistogram2D delPosPzShift;

    IHistogram2D delZvsVtxZShift;
    IHistogram2D delMvsVtxZShift;
    IHistogram2D delPvsVtxZShift;
    IHistogram2D chiSqVtxZShift;

    IHistogram2D delElePxBSC;
    IHistogram2D delElePyBSC;
    IHistogram2D delElePzBSC;

    IHistogram2D delPosPxBSC;
    IHistogram2D delPosPyBSC;
    IHistogram2D delPosPzBSC;

    IHistogram2D delZvsVtxZBSC;
    IHistogram2D delMvsVtxZBSC;
    IHistogram2D delPvsVtxZBSC;
    IHistogram2D chiSqVtxZBSC;

    IHistogram2D delElePxBSCShift;
    IHistogram2D delElePyBSCShift;
    IHistogram2D delElePzBSCShift;
    IHistogram2D chiSqVtxZBSCShift;

    IHistogram2D delPosPxBSCShift;
    IHistogram2D delPosPyBSCShift;
    IHistogram2D delPosPzBSCShift;

    IHistogram2D delZvsVtxZBSCShift;
    IHistogram2D delMvsVtxZBSCShift;
    IHistogram2D delPvsVtxZBSCShift;

    IHistogram2D delZShiftMinusNoShiftvsVtxZBSC;
    IHistogram2D delZShiftMinusNoShiftvsVtxZUnc;

    IHistogram2D delMShiftMinusNoShiftvsVtxZBSC;
    IHistogram2D delMShiftMinusNoShiftvsVtxZUnc;

    private enum Constraint {

        /**
         * Represents a fit with no constraints.
         */
        UNCONSTRAINED,
        /**
         * Represents a fit with beam spot constraints.
         */
        BS_CONSTRAINED,
        /**
         * Represents a fit with target constraints.
         */
        TARGET_CONSTRAINED
    }

    String detName;

    /**
     * The B field map
     */
    FieldMap bFieldMap = null;
    /**
     * Position of the Ecal face
     */
    private double ecalPosition = 0; // mm

    /**
     * Z position to start extrapolation from
     */
    double extStartPos = 700; // mm

    /**
     * The extrapolation step size
     */
    double stepSize = 5.0; // mm
    /**
     * Name of the constant denoting the position of the Ecal face in the
     * compact description.
     */
    private static final String ECAL_POSITION_CONSTANT_NAME = "ecal_dface";

    public void setFinalStateParticlesColName(String fsp) {
        this.finalStateParticlesColName = fsp;
    }

    public void setOutputFilePath(String output) {
        this.outputFile = output;
    }

    public void setNoTracksFilePath(String output) {
        this.notrackFile = output;
    }

    public void setDoSkim(boolean doit) {
        this.doSkim = doit;
    }

    protected void detectorChanged(Detector detector) {

        super.detectorChanged(detector);
        detName = detector.getName();
        double maxFactor = 1.5;
        double feeMomentumCut = 0.75; //this number, multiplied by the beam energy, is the actual cut
        B_FIELD = detector.getFieldMap().getField(new BasicHep3Vector(0, 0, 500)).y();
        // Get the field map from the detector object
        bFieldMap = detector.getFieldMap();
        // Get the position of the Ecal from the compact description
        ecalPosition = detector.getConstants().get(ECAL_POSITION_CONSTANT_NAME).getValue();
        // Get the position of the Ecal from the compact description
        ecalPosition = detector.getConstants().get(ECAL_POSITION_CONSTANT_NAME).getValue();
        LOGGER.setLevel(Level.ALL);
        LOGGER.info("B_FIELD=" + B_FIELD);
        LOGGER.info("Setting up the plotter");
        aida.tree().cd("/");


        /*  Final State Particle Quantities   */
        /*  plot electron & positron momentum separately  */
        elePx = aida.histogram1D("Electron Px (GeV)", nbins, -0.1 * beamEnergy, 0.200 * beamEnergy);
        elePy = aida.histogram1D("Electron Py (GeV)", nbins, -0.1 * beamEnergy, 0.1 * beamEnergy);
        elePz = aida.histogram1D("Electron Pz (GeV)", nbins, 0, beamEnergy * maxFactor);

        posPx = aida.histogram1D("Positron Px (GeV)", nbins, -0.1 * beamEnergy, 0.200 * beamEnergy);
        posPy = aida.histogram1D("Positron Py (GeV)", nbins, -0.1 * beamEnergy, 0.1 * beamEnergy);
        posPz = aida.histogram1D("Positron Pz (GeV)", nbins, 0, beamEnergy * maxFactor);

        delElePxReco = aida.histogram2D("(MC-Reco) Electron Px (GeV) vs MC Z Vertex", nbins, -0.02, 0.02, nbins, -10, 100);
        delElePyReco = aida.histogram2D("(MC-Reco) Electron Py (GeV) vs MC Z Vertex", nbins, -0.02, 0.02, nbins, -10, 100);
        delElePzReco = aida.histogram2D("(MC-Reco) Electron Pz (GeV) vs MC Z Vertex", nbins, -0.1, 0.1, nbins, -10, 100);

        delPosPxReco = aida.histogram2D("(MC-Reco) Positron Px (GeV) vs MC Z Vertex", nbins, -0.02, 0.02, nbins, -10, 100);
        delPosPyReco = aida.histogram2D("(MC-Reco) Positron Py (GeV) vs MC Z Vertex", nbins, -0.02, 0.02, nbins, -10, 100);
        delPosPzReco = aida.histogram2D("(MC-Reco) Positron Pz (GeV) vs MC Z Vertex", nbins, -0.1, 0.1, nbins, -10, 100);

        delElePx = aida.histogram2D("(MC-Fitted) Electron Px (GeV) vs MC Z Vertex", nbins, -0.02, 0.02, nbins, -10, 100);
        delElePy = aida.histogram2D("(MC-Fitted) Electron Py (GeV) vs MC Z Vertex", nbins, -0.02, 0.02, nbins, -10, 100);
        delElePz = aida.histogram2D("(MC-Fitted) Electron Pz (GeV) vs MC Z Vertex", nbins, -0.1, 0.1, nbins, -10, 100);

        delPosPx = aida.histogram2D("(MC-Fitted) Positron Px (GeV) vs MC Z Vertex", nbins, -0.02, 0.02, nbins, -10, 100);
        delPosPy = aida.histogram2D("(MC-Fitted) Positron Py (GeV) vs MC Z Vertex", nbins, -0.02, 0.02, nbins, -10, 100);
        delPosPz = aida.histogram2D("(MC-Fitted) Positron Pz (GeV) vs MC Z Vertex", nbins, -0.1, 0.1, nbins, -10, 100);

        delZvsVtxZ = aida.histogram2D("(MC-Fitted) Z vs MC Z Vertex", nbins, -50, 50, nbins, -10, 100);
        delMvsVtxZ = aida.histogram2D("(MC-Fitted) Mass vs MC Z Vertex", nbins, -0.04, 0.04, nbins, -10, 100);
        delPvsVtxZ = aida.histogram2D("(MC-Fitted) P vs MC Z Vertex", nbins, -1, 1, nbins, -10, 100);
        chiSqVtxZ = aida.histogram2D("Fitted Chi Sq vs MC Z Vertex", nbins, 0, 10, nbins, -10, 100);

        delElePxShift = aida.histogram2D("Shifted (MC-Fitted) Electron Px (GeV) vs MC Z Vertex", nbins, -0.02, 0.02, nbins, -10, 100);
        delElePyShift = aida.histogram2D("Shifted (MC-Fitted) Electron Py (GeV) vs MC Z Vertex", nbins, -0.02, 0.02, nbins, -10, 100);
        delElePzShift = aida.histogram2D("Shifted (MC-Fitted) Electron Pz (GeV) vs MC Z Vertex", nbins, -0.1, 0.1, nbins, -10, 100);

        delPosPxShift = aida.histogram2D("Shifted (MC-Fitted) Positron Px (GeV) vs MC Z Vertex", nbins, -0.02, 0.02, nbins, -10, 100);
        delPosPyShift = aida.histogram2D("Shifted (MC-Fitted) Positron Py (GeV) vs MC Z Vertex", nbins, -0.02, 0.02, nbins, -10, 100);
        delPosPzShift = aida.histogram2D("Shifted (MC-Fitted) Positron Pz (GeV) vs MC Z Vertex", nbins, -0.1, 0.1, nbins, -10, 100);

        delZvsVtxZShift = aida.histogram2D("Shifted (MC-Fitted) Z vs MC Z Vertex", nbins, -50, 50, nbins, -10, 100);
        delMvsVtxZShift = aida.histogram2D("Shifted (MC-Fitted) Mass vs MC Z Vertex", nbins, -0.04, 0.04, nbins, -10, 100);
        delPvsVtxZShift = aida.histogram2D("Shifted (MC-Fitted) P vs MC Z Vertex", nbins, -1, 1, nbins, -10, 100);
        chiSqVtxZShift = aida.histogram2D("Shifted Chi Sq vs MC Z Vertex", nbins, 0, 10, nbins, -10, 100);

        delElePxBSC = aida.histogram2D("BSC (MC-Fitted) Electron Px (GeV) vs MC Z Vertex", nbins, -0.02, 0.02, nbins, -10, 100);
        delElePyBSC = aida.histogram2D("BSC (MC-Fitted) Electron Py (GeV) vs MC Z Vertex", nbins, -0.02, 0.02, nbins, -10, 100);
        delElePzBSC = aida.histogram2D("BSC (MC-Fitted) Electron Pz (GeV) vs MC Z Vertex", nbins, -0.1, 0.1, nbins, -10, 100);

        delPosPxBSC = aida.histogram2D("BSC (MC-Fitted) Positron Px (GeV) vs MC Z Vertex", nbins, -0.02, 0.02, nbins, -10, 100);
        delPosPyBSC = aida.histogram2D("BSC (MC-Fitted) Positron Py (GeV) vs MC Z Vertex", nbins, -0.02, 0.02, nbins, -10, 100);
        delPosPzBSC = aida.histogram2D("BSC (MC-Fitted) Positron Pz (GeV) vs MC Z Vertex", nbins, -0.1, 0.1, nbins, -10, 100);

        delZvsVtxZBSC = aida.histogram2D("BSC (MC-Fitted) Z vs MC Z Vertex", nbins, -50, 50, nbins, -10, 100);
        delMvsVtxZBSC = aida.histogram2D("BSC (MC-Fitted) Mass vs MC Z Vertex", nbins, -0.04, 0.04, nbins, -10, 100);
        delPvsVtxZBSC = aida.histogram2D("BSC (MC-Fitted) P vs MC Z Vertex", nbins, -1, 1, nbins, -10, 100);
        chiSqVtxZBSC = aida.histogram2D("BSC Chi Sq vs MC Z Vertex", nbins, 0, 25, nbins, -10, 100);

        delElePxBSCShift = aida.histogram2D("BSCShift (MC-Fitted) Electron Px (GeV) vs MC Z Vertex", nbins, -0.02, 0.02, nbins, -10, 100);
        delElePyBSCShift = aida.histogram2D("BSCShift (MC-Fitted) Electron Py (GeV) vs MC Z Vertex", nbins, -0.02, 0.02, nbins, -10, 100);
        delElePzBSCShift = aida.histogram2D("BSCShift (MC-Fitted) Electron Pz (GeV) vs MC Z Vertex", nbins, -0.1, 0.1, nbins, -10, 100);

        delPosPxBSCShift = aida.histogram2D("BSCShift (MC-Fitted) Positron Px (GeV) vs MC Z Vertex", nbins, -0.02, 0.02, nbins, -10, 100);
        delPosPyBSCShift = aida.histogram2D("BSCShift (MC-Fitted) Positron Py (GeV) vs MC Z Vertex", nbins, -0.02, 0.02, nbins, -10, 100);
        delPosPzBSCShift = aida.histogram2D("BSCShift (MC-Fitted) Positron Pz (GeV) vs MC Z Vertex", nbins, -0.1, 0.1, nbins, -10, 100);

        delZvsVtxZBSCShift = aida.histogram2D("BSCShift (MC-Fitted) Z vs MC Z Vertex", nbins, -50, 50, nbins, -10, 100);
        delMvsVtxZBSCShift = aida.histogram2D("BSCShift (MC-Fitted) Mass vs MC Z Vertex", nbins, -0.04, 0.04, nbins, -10, 100);
        delPvsVtxZBSCShift = aida.histogram2D("BSCShift (MC-Fitted) P vs MC Z Vertex", nbins, -1, 1, nbins, -10, 100);
        chiSqVtxZBSCShift = aida.histogram2D("BSCShift Chi Sq vs MC Z Vertex", nbins, 0, 25, nbins, -10, 100);

        delMvsVtxZPatch = aida.histogram2D("(MC-Patch) Mass vs MC Z Vertex", nbins, -0.04, 0.04, nbins, -10, 100);
        delMUncvsPatchvsZ = aida.histogram2D("(Patch-Reco) Mass vs MC Z Vertex", nbins, -0.04, 0.04, nbins, -10, 100);

        delZShiftMinusNoShiftvsVtxZBSC = aida.histogram2D("BSC (Shift-No Shift) Z vs MC Z Vertex", nbins, -50, 50, nbins, -10, 100);
        delZShiftMinusNoShiftvsVtxZUnc = aida.histogram2D("Unconstrained (Shift-No Shift) Z vs MC Z Vertex", nbins, -50, 50, nbins, -10, 100);
        delMShiftMinusNoShiftvsVtxZBSC = aida.histogram2D("BSC (Shift-No Shift) Mass vs MC Z Vertex", nbins, -0.04, 0.04, nbins, -10, 100);
        delMShiftMinusNoShiftvsVtxZUnc = aida.histogram2D("Unconstrained (Shift-No Shift) Mass vs MC Z Vertex", nbins, -0.04, 0.04, nbins, -10, 100);

    }

    @Override
    public void process(EventHeader event) {
        /*  make sure everything is there */

//        List<CalorimeterHit> hits;
//        if (event.hasCollection(CalorimeterHit.class, calibratedHitCollectionName))
//            hits = event.get(CalorimeterHit.class, calibratedHitCollectionName);
//        else
//            return; //this might be a non-data event        
        if (!event.hasCollection(ReconstructedParticle.class, finalStateParticlesColName)) {
            if (debug)
                LOGGER.info(finalStateParticlesColName + " collection not found???");
            return;
        }

        if (!event.hasCollection(MCParticle.class, mcParticlesColName)) {
            if (debug)
                LOGGER.info(mcParticlesColName + " collection not found???");
            return;
        }

        //check to see if this event is from the correct trigger (or "all");
        if (!matchTrigger(event))
            return;

        nRecoEvents++;
        if (debug)
            LOGGER.info("##########  Start of VertexDebugger   ##############");
        List<ReconstructedParticle> finalStateParticles = event.get(ReconstructedParticle.class, finalStateParticlesColName);
        if (debug)
            LOGGER.info("This events has " + finalStateParticles.size() + " final state particles");
//        for (ReconstructedParticle fsPart1 : finalStateParticles) {
        List<MCParticle> MCParticleList = event.get(MCParticle.class, mcParticlesColName);

        boolean saveEvent = false;
        boolean saveNoTrack = false;
        //find electron and positron MCParticles
        Hep3Vector vertexPositionMC = null;
        double apMassMC = -9999;
        MCParticle eleMC = null;
        MCParticle posMC = null;

        for (MCParticle mcp : MCParticleList)
//             if (debug)
            //               LOGGER.info("MC PDGID = "+mcp.getPDGID()+"; # daughters = "+mcp.getDaughters().size());
            if (mcp.getPDGID() == 622 && mcp.getDaughters().size() == 2) {
                vertexPositionMC = mcp.getEndPoint();
                apMassMC = mcp.getMass();
                List<MCParticle> daugList = mcp.getDaughters();
                for (MCParticle daug : daugList)
                    if (daug.getPDGID() == 11)
                        eleMC = daug;
                    else if (daug.getPDGID() == -11)
                        posMC = daug;
            }
        if (eleMC == null || posMC == null) {
            if (debug)
                LOGGER.info("Couldn't find the MC e+e- from A'?????  Quitting.");
            return;
        }

        ReconstructedParticle electron = null;
        ReconstructedParticle positron = null;
        double bestMom = 99999;
        for (ReconstructedParticle fsp : finalStateParticles)
            if (fsp.getCharge() > 0) {//found a positron
                if (fsp.getMomentum().y() * posMC.getPY() > 0) //dumb matching..make sure in same half
                    positron = fsp;
            } else if (fsp.getCharge() < 0)
                if (fsp.getMomentum().y() * eleMC.getPY() > 0 && Math.abs(eleMC.getMomentum().magnitude() - fsp.getMomentum().magnitude()) < bestMom) //dumb matching..make sure in same half
                    electron = fsp;
        if (electron == null || positron == null) {
            if (debug)
                LOGGER.info("Couldn't find MCP matched reconed e+ or e- ?????  Quitting.");
            return;
        }
        //ok..made it this far...now lets do the vertexing.  
        if (debug) {
            LOGGER.info("Found A' MC with vertex at :" + vertexPositionMC.x() + "; " + vertexPositionMC.y() + "; " + vertexPositionMC.z());
            LOGGER.info("Found A' MC with mass = " + apMassMC);
            LOGGER.info("Found A' MC electron momentum = " + eleMC.getMomentum().x() + "; " + eleMC.getMomentum().y() + "; " + eleMC.getMomentum().z());
            LOGGER.info("Found A' MC positron momentum = " + posMC.getMomentum().x() + "; " + posMC.getMomentum().y() + "; " + posMC.getMomentum().z());
        }

        if (debug) {
            LOGGER.info("electron momentum = " + electron.getMomentum().x() + "; " + electron.getMomentum().y() + "; " + electron.getMomentum().z());
            LOGGER.info("positron momentum = " + positron.getMomentum().x() + "; " + positron.getMomentum().y() + "; " + positron.getMomentum().z());
        }

        Hep3Vector delPEleMC = VecOp.sub(eleMC.getMomentum(), electron.getMomentum());
        Hep3Vector delPPosMC = VecOp.sub(posMC.getMomentum(), positron.getMomentum());

        // Covert the tracks to BilliorTracks.
        BilliorTrack electronBTrack = toBilliorTrack(electron.getTracks().get(0));
        BilliorTrack positronBTrack = toBilliorTrack(positron.getTracks().get(0));
        // Generate a candidate vertex and particle.
        if (debug)
            LOGGER.info("Unconstrained R=(0,0,0)  ##############");

        BilliorVertex vtxFit = fitVertex(Constraint.UNCONSTRAINED, electronBTrack, positronBTrack);
        Hep3Vector vtxPos = vtxFit.getPosition();
        Hep3Vector pEleFit = vtxFit.getFittedMomentum(0);
        Hep3Vector pPosFit = vtxFit.getFittedMomentum(1);

        double delz = vertexPositionMC.z() - vtxPos.z();
        LOGGER.info("vertexMC z=" + vertexPositionMC.z() + " vtxPos z = " + vtxPos.z());
        double delm = apMassMC - vtxFit.getInvMass();
        double mUnc = vtxFit.getInvMass();
        double chisq = vtxFit.getChi2();
//        double delpEle = eleMC.getMomentum().magnitude() - electron.getMomentum().magnitude();
//        double delpPos = posMC.getMomentum().magnitude() - positron.getMomentum().magnitude();
        Hep3Vector delPEleFit = new BasicHep3Vector(eleMC.getMomentum().x() - pEleFit.y(),
                eleMC.getMomentum().y() - pEleFit.z(), eleMC.getMomentum().z() - pEleFit.x());
        Hep3Vector delPPosFit = new BasicHep3Vector(posMC.getMomentum().x() - pPosFit.y(),
                eleMC.getMomentum().y() - pPosFit.z(), posMC.getMomentum().z() - pPosFit.x());

        delZvsVtxZ.fill(delz, vertexPositionMC.z());
        if (delz < 10) {//only fill deltaM  and delta Ps if the V0 looks like it got the track from A'
            delPvsVtxZ.fill(delPEleFit.magnitude(), vertexPositionMC.z());
            delPvsVtxZ.fill(delPPosFit.magnitude(), vertexPositionMC.z());
            delMvsVtxZ.fill(delm, vertexPositionMC.z());
            chiSqVtxZ.fill(chisq, vertexPositionMC.z());

            delElePxReco.fill(delPEleMC.x(), vertexPositionMC.z());
            delElePyReco.fill(delPEleMC.y(), vertexPositionMC.z());
            delElePzReco.fill(delPEleMC.z(), vertexPositionMC.z());
            delPosPxReco.fill(delPPosMC.x(), vertexPositionMC.z());
            delPosPyReco.fill(delPPosMC.y(), vertexPositionMC.z());
            delPosPzReco.fill(delPPosMC.z(), vertexPositionMC.z());

            delElePx.fill(eleMC.getMomentum().x() - pEleFit.y(), vertexPositionMC.z());
            delElePy.fill(eleMC.getMomentum().y() - pEleFit.z(), vertexPositionMC.z());
            delElePz.fill(eleMC.getMomentum().z() - pEleFit.x(), vertexPositionMC.z());
            delPosPx.fill(posMC.getMomentum().x() - pPosFit.y(), vertexPositionMC.z());
            delPosPy.fill(posMC.getMomentum().y() - pPosFit.z(), vertexPositionMC.z());
            delPosPz.fill(posMC.getMomentum().z() - pPosFit.x(), vertexPositionMC.z());

            // patch the track parameters at the found vertex
            if (debug)
                LOGGER.info("Patch R=(0,0,0)  ##############");
            patchVertex(electron, positron, vtxFit);
            double mPatch = vtxFit.getInvMass();
            delMvsVtxZPatch.fill(apMassMC - mPatch, vertexPositionMC.z());
            delMUncvsPatchvsZ.fill(mPatch - mUnc, vertexPositionMC.z());
//////////////////////////////              
// ok, take the initial tracks and move the reference point to Vz (==Vx in tracking frame)
            // first, get the x,y,z of the track at the perigee
//            double[] newRef = {vtxPos.z(), vtxPos.x(), vtxPos.y()};
            double[] newRef = {vtxPos.z(), vtxPos.x(), 0.0};//the  TrackUtils.getParametersAtNewRefPoint method only shifts in xy...?
            BaseTrackState eleOldTS = (BaseTrackState) electron.getTracks().get(0).getTrackStates().get(0);
            BaseTrackState posOldTS = (BaseTrackState) positron.getTracks().get(0).getTrackStates().get(0);
            double[] eleShiftPars = TrackUtils.getParametersAtNewRefPoint(newRef, eleOldTS);
            double[] posShiftPars = TrackUtils.getParametersAtNewRefPoint(newRef, posOldTS);
            SymmetricMatrix eleShiftCov = TrackUtils.getCovarianceAtNewRefPoint(newRef, eleOldTS.getReferencePoint(), eleOldTS.getParameters(), new SymmetricMatrix(5, eleOldTS.getCovMatrix(), true));
            SymmetricMatrix posShiftCov = TrackUtils.getCovarianceAtNewRefPoint(newRef, posOldTS.getReferencePoint(), posOldTS.getParameters(), new SymmetricMatrix(5, posOldTS.getCovMatrix(), true));
            BaseTrackState eleShiftTS = new BaseTrackState(eleShiftPars, newRef, eleShiftCov.asPackedArray(true), TrackState.AtIP, B_FIELD);
            BaseTrackState posShiftTS = new BaseTrackState(posShiftPars, newRef, posShiftCov.asPackedArray(true), TrackState.AtIP, B_FIELD);
//            BaseTrackState eleShiftTS = new BaseTrackState(eleShiftPars, newRef, eleOldTS.getCovMatrix(), TrackState.AtIP, B_FIELD);
//            BaseTrackState posShiftTS = new BaseTrackState(posShiftPars, newRef, posOldTS.getCovMatrix(), TrackState.AtIP, B_FIELD);
            BilliorTrack electronBTrackShift = toBilliorTrack(eleShiftTS);
            BilliorTrack positronBTrackShift = toBilliorTrack(posShiftTS);
            ////////////////////////////////            
//           get the new fitter
            if (debug)
                LOGGER.info("Unconstrained R=(" + newRef[0] + "," + newRef[1] + "," + newRef[2] + ") ##############");
            BilliorVertex vtxFitShift = fitVertex(Constraint.UNCONSTRAINED, electronBTrackShift, positronBTrackShift);
            Hep3Vector vtxPosShift = vtxFitShift.getPosition();
            Hep3Vector pEleFitShift = vtxFitShift.getFittedMomentum(0);
            Hep3Vector pPosFitShift = vtxFitShift.getFittedMomentum(1);
//            Hep3Vector delPEleShift = VecOp.sub(eleMC.getMomentum(), pEleFitShift);
//            Hep3Vector delPPosShift = VecOp.sub(posMC.getMomentum(), pPosFitShift);
            Hep3Vector delPEleFitShift = new BasicHep3Vector(eleMC.getMomentum().x() - pEleFitShift.y(),
                    eleMC.getMomentum().y() - pEleFitShift.z(), eleMC.getMomentum().z() - pEleFitShift.x());
            Hep3Vector delPPosFitShift = new BasicHep3Vector(posMC.getMomentum().x() - pPosFitShift.y(),
                    eleMC.getMomentum().y() - pPosFitShift.z(), posMC.getMomentum().z() - pPosFitShift.x());
            double delzShift = vertexPositionMC.z() - (vtxPosShift.z() + vtxPos.z());//add on old z-value
            LOGGER.info("UnConstrained shifted vertexMC z=" + vertexPositionMC.z() + " re-fit vtxPos z = " + (vtxPosShift.z() + vtxPos.z()));
            LOGGER.info("delzShift = " + delzShift);
            double delmShift = apMassMC - vtxFitShift.getInvMass();
            double mUncShift = vtxFitShift.getInvMass();
            double chisqShift = vtxFitShift.getChi2();
            delZvsVtxZShift.fill(delzShift, vertexPositionMC.z());
            delPvsVtxZShift.fill(delPEleFitShift.magnitude(), vertexPositionMC.z());
            delPvsVtxZShift.fill(delPPosFitShift.magnitude(), vertexPositionMC.z());
            delMvsVtxZShift.fill(delmShift, vertexPositionMC.z());
            chiSqVtxZShift.fill(chisqShift, vertexPositionMC.z());
            delElePxShift.fill(eleMC.getMomentum().x() - pEleFitShift.y(), vertexPositionMC.z());
            delElePyShift.fill(eleMC.getMomentum().y() - pEleFitShift.z(), vertexPositionMC.z());
            delElePzShift.fill(eleMC.getMomentum().z() - pEleFitShift.x(), vertexPositionMC.z());
            delPosPxShift.fill(posMC.getMomentum().x() - pPosFitShift.y(), vertexPositionMC.z());
            delPosPyShift.fill(posMC.getMomentum().y() - pPosFitShift.z(), vertexPositionMC.z());
            delPosPzShift.fill(posMC.getMomentum().z() - pPosFitShift.x(), vertexPositionMC.z());
            //    Ok...same thing with beam-spot constrained
            if (debug)
                LOGGER.info("Constrained R=(0,0,0) ##############");
            BilliorVertex vtxFitBSC = fitVertex(Constraint.BS_CONSTRAINED, electronBTrack, positronBTrack, new BasicHep3Vector(beamPosition));
            Hep3Vector vtxPosBSC = vtxFitBSC.getPosition();
            Hep3Vector pEleFitBSC = vtxFitBSC.getFittedMomentum(0);
            Hep3Vector pPosFitBSC = vtxFitBSC.getFittedMomentum(1);
//            Hep3Vector delPEleBSC = VecOp.sub(eleMC.getMomentum(), pEleFitBSC);
//            Hep3Vector delPPosBSC = VecOp.sub(posMC.getMomentum(), pPosFitBSC);
            Hep3Vector delPEleFitBSC = new BasicHep3Vector(eleMC.getMomentum().x() - pEleFitBSC.y(),
                    eleMC.getMomentum().y() - pEleFitBSC.z(), eleMC.getMomentum().z() - pEleFitBSC.x());
            Hep3Vector delPPosFitBSC = new BasicHep3Vector(posMC.getMomentum().x() - pPosFitBSC.y(),
                    eleMC.getMomentum().y() - pPosFitBSC.z(), posMC.getMomentum().z() - pPosFitBSC.x());
            double delzBSC = vertexPositionMC.z() - (vtxPosBSC.z());
            LOGGER.info("Constrained R=0  vertexMC z=" + vertexPositionMC.z() + " re-fit vtxPos z = " + (vtxPosBSC.z()));
            LOGGER.info("delzBSC = " + delzBSC);
            double delmBSC = apMassMC - vtxFitBSC.getInvMass();
            double mUncBSC = vtxFitBSC.getInvMass();
            double chisqBSC = vtxFitBSC.getChi2();
            delZvsVtxZBSC.fill(delzBSC, vertexPositionMC.z());
            delPvsVtxZBSC.fill(delPEleFitBSC.magnitude(), vertexPositionMC.z());
            delPvsVtxZBSC.fill(delPPosFitBSC.magnitude(), vertexPositionMC.z());
            delMvsVtxZBSC.fill(delmBSC, vertexPositionMC.z());
            chiSqVtxZBSC.fill(chisqBSC, vertexPositionMC.z());
            delElePxBSC.fill(eleMC.getMomentum().x() - pEleFitBSC.y(), vertexPositionMC.z());
            delElePyBSC.fill(eleMC.getMomentum().y() - pEleFitBSC.z(), vertexPositionMC.z());
            delElePzBSC.fill(eleMC.getMomentum().z() - pEleFitBSC.x(), vertexPositionMC.z());
            delPosPxBSC.fill(posMC.getMomentum().x() - pPosFitBSC.y(), vertexPositionMC.z());
            delPosPyBSC.fill(posMC.getMomentum().y() - pPosFitBSC.z(), vertexPositionMC.z());
            delPosPzBSC.fill(posMC.getMomentum().z() - pPosFitBSC.x(), vertexPositionMC.z());

            //    Ok...same thing with beam-spot constrained
            Hep3Vector beamRelToVtx = new BasicHep3Vector(-vtxPos.z() + beamPosition[0], -vtxPos.x() + beamPosition[1], -vtxPos.y() + beamPosition[2]);
            if (debug)
                LOGGER.info("Constrained R=(" + newRef[0] + "," + newRef[1] + "," + newRef[2] + ") ##############");
            BilliorVertex vtxFitBSCShift = fitVertex(Constraint.BS_CONSTRAINED, electronBTrackShift, positronBTrackShift, beamRelToVtx);
            Hep3Vector vtxPosBSCShift = vtxFitBSCShift.getPosition();
            Hep3Vector pEleFitBSCShift = vtxFitBSCShift.getFittedMomentum(0);
            Hep3Vector pPosFitBSCShift = vtxFitBSCShift.getFittedMomentum(1);
//            Hep3Vector delPEleBSCShift = VecOp.sub(eleMC.getMomentum(), pEleFitBSCShift);
//            Hep3Vector delPPosBSCShift = VecOp.sub(posMC.getMomentum(), pPosFitBSCShift);
            Hep3Vector delPEleFitBSCShift = new BasicHep3Vector(eleMC.getMomentum().x() - pEleFitBSCShift.y(),
                    eleMC.getMomentum().y() - pEleFitBSCShift.z(), eleMC.getMomentum().z() - pEleFitBSCShift.x());
            Hep3Vector delPPosFitBSCShift = new BasicHep3Vector(posMC.getMomentum().x() - pPosFitBSCShift.y(),
                    eleMC.getMomentum().y() - pPosFitBSCShift.z(), posMC.getMomentum().z() - pPosFitBSCShift.x());
            double delzBSCShift = vertexPositionMC.z() - (vtxPosBSCShift.z() + vtxPos.z());//add on old z-value
            LOGGER.info("Constrained shifted vertexMC z=" + vertexPositionMC.z() + " re-fit vtxPos z = " + (vtxPosBSCShift.z() + vtxPos.z()));
            LOGGER.info("delzBSCShif = " + delzBSCShift);
            double delmBSCShift = apMassMC - vtxFitBSCShift.getInvMass();
            double mUncBSCShift = vtxFitBSCShift.getInvMass();
            double chisqBSCShift = vtxFitBSCShift.getChi2();
            delZvsVtxZBSCShift.fill(delzBSCShift, vertexPositionMC.z());
            delPvsVtxZBSCShift.fill(delPEleFitBSCShift.magnitude(), vertexPositionMC.z());
            delPvsVtxZBSCShift.fill(delPPosFitBSCShift.magnitude(), vertexPositionMC.z());
            delMvsVtxZBSCShift.fill(delmBSCShift, vertexPositionMC.z());
            chiSqVtxZBSCShift.fill(chisqBSCShift, vertexPositionMC.z());
            delElePxBSCShift.fill(eleMC.getMomentum().x() - pEleFitBSCShift.y(), vertexPositionMC.z());
            delElePyBSCShift.fill(eleMC.getMomentum().y() - pEleFitBSCShift.z(), vertexPositionMC.z());
            delElePzBSCShift.fill(eleMC.getMomentum().z() - pEleFitBSCShift.x(), vertexPositionMC.z());
            delPosPxBSCShift.fill(posMC.getMomentum().x() - pPosFitBSCShift.y(), vertexPositionMC.z());
            delPosPyBSCShift.fill(posMC.getMomentum().y() - pPosFitBSCShift.z(), vertexPositionMC.z());
            delPosPzBSCShift.fill(posMC.getMomentum().z() - pPosFitBSCShift.x(), vertexPositionMC.z());

            delZShiftMinusNoShiftvsVtxZBSC.fill(delzBSC-delzBSCShift, vertexPositionMC.z());
            delZShiftMinusNoShiftvsVtxZUnc.fill(delz-delzShift, vertexPositionMC.z());
            
            delMShiftMinusNoShiftvsVtxZBSC.fill(delmBSC-delmBSCShift, vertexPositionMC.z());
            delMShiftMinusNoShiftvsVtxZUnc.fill(delm-delmShift, vertexPositionMC.z());
            
        }
    }

    private void setupWriter() {
        // Cleanup existing writer.
        if (writer != null)
            try {
                writer.flush();
                writer.close();
                writer = null;
            } catch (IOException x) {
                System.err.println(x.getMessage());
            }

        // Setup new writer.
        try {
            writer = new LCIOWriter(outputFile);
        } catch (IOException x) {
            throw new RuntimeException("Error creating writer", x);
        }

        try {
            writer.reOpen();
        } catch (IOException x) {
            throw new RuntimeException("Error rewinding LCIO file", x);
        }
    }

    private void setupNoTrackWriter() {
        // Cleanup existing writer.
        if (notrackwriter != null)
            try {
                notrackwriter.flush();
                notrackwriter.close();
                notrackwriter = null;
            } catch (IOException x) {
                System.err.println(x.getMessage());
            }

        // Setup new writer.
        try {
            notrackwriter = new LCIOWriter(notrackFile);
        } catch (IOException x) {
            throw new RuntimeException("Error creating writer", x);
        }

        try {
            notrackwriter.reOpen();
        } catch (IOException x) {
            throw new RuntimeException("Error rewinding LCIO file", x);
        }
    }

    protected void startOfData() {
//        setupWriter();
//        setupNoTrackWriter();
    }

    @Override
    public void endOfData() {
//        try {
//            writer.close();
//        } catch (IOException x) {
//            throw new RuntimeException("Error rewinding LCIO file", x);
//        }
//        try {
//            notrackwriter.close();
//        } catch (IOException x) {
//            throw new RuntimeException("Error rewinding LCIO file", x);
//        }
    }

    private BilliorVertex fitVertex(Constraint constraint, BilliorTrack electron, BilliorTrack positron) {
        return fitVertex(constraint, electron, positron, null);
    }

    private BilliorVertex fitVertex(Constraint constraint, BilliorTrack electron, BilliorTrack positron, Hep3Vector v0) {
        // Create a vertex fitter from the magnetic field.

        BilliorVertexer vtxFitter = new BilliorVertexer(B_FIELD);
        // TODO: The beam size should come from the conditions database.
        vtxFitter.setBeamSize(beamSize);
//        vtxFitter.setBeamPosition(beamPosition);

        vtxFitter.setDebug(true);

//            vtxFitter.setV0(v0.v());
        if (v0 != null)
            vtxFitter.setBeamPosition(v0.v());

        // Perform the vertexing based on the specified constraint.
        switch (constraint) {
            case UNCONSTRAINED:
                vtxFitter.doBeamSpotConstraint(false);
                break;
            case BS_CONSTRAINED:
                vtxFitter.doBeamSpotConstraint(true);
                break;
            case TARGET_CONSTRAINED:
                vtxFitter.doTargetConstraint(true);
                break;
        }

        // Add the electron and positron tracks to a track list for
        // the vertex fitter.
        List<BilliorTrack> billiorTracks = new ArrayList<BilliorTrack>();
        billiorTracks.add(electron);
        billiorTracks.add(positron);

        // Find and return a vertex based on the tracks.
        BilliorVertex vtx = vtxFitter.fitVertex(billiorTracks);

        return vtx;
    }

    private BilliorTrack toBilliorTrack(Track track) {
        // Generate and return the billior track.
        return new BilliorTrack(track);
    }

    private BilliorTrack toBilliorTrack(HelicalTrackFit htf) {
        // Generate and return the billior track.
        return new BilliorTrack(htf);
    }

    private BilliorTrack toBilliorTrack(TrackState track) {
        // Generate and return the billior track.
        return new BilliorTrack(track, 0, 0); // track state doesn't store chi^2 info (stored in the Track object)
    }

    // patchVertex written by Norman Graf...I plucked this from HpsReconParticleDriver
    private void patchVertex(ReconstructedParticle electron, ReconstructedParticle positron, BilliorVertex v) {
//        ReconstructedParticle rp = v.getAssociatedParticle();
//        List<ReconstructedParticle> parts = rp.getParticles();
//        ReconstructedParticle electron = null;
//        ReconstructedParticle positron = null;
//        for (ReconstructedParticle part : parts) {
//            if (part.getCharge() < 0)
//                electron = part;
//            if (part.getCharge() > 0)
//                positron = part;
//        }
        //electron
        Track et = electron.getTracks().get(0);
        double etrackMom = electron.getMomentum().magnitude();
        HelicalTrackFit ehtf = TrackUtils.getHTF(et);
        // propagate this to the vertex z position...
        // Note that HPS y is lcsim z  //mg...I don't think this is correct!!!  First, HPS x is lcsim y!    
        //  And I think the vertex is already in the HPS frame!
        double es = HelixUtils.PathToZPlane(ehtf, v.getPosition().y());
//        double es = HelixUtils.PathToXPlane(ehtf, v.getPosition().z(), 10000., 100).get(0);
        LOGGER.info("vertex z=" + v.getPosition().z() + " vertex y = " + v.getPosition().y());
        Hep3Vector epointOnTrackAtVtx = HelixUtils.PointOnHelix(ehtf, es);
        Hep3Vector edirOfTrackAtVtx = HelixUtils.Direction(ehtf, es);
        Hep3Vector emomAtVtx = VecOp.mult(etrackMom, VecOp.unit(edirOfTrackAtVtx));
        //positron
        Track pt = positron.getTracks().get(0);
        double ptrackMom = positron.getMomentum().magnitude();
        HelicalTrackFit phtf = TrackUtils.getHTF(pt);
        // propagate this to the vertex z position...
        // Note that HPS y is lcsim z
        double ps = HelixUtils.PathToZPlane(phtf, v.getPosition().y());
//        double ps = HelixUtils.PathToXPlane(phtf, v.getPosition().z(), 10000., 100).get(0);
        Hep3Vector ppointOnTrackAtVtx = HelixUtils.PointOnHelix(phtf, ps);
        Hep3Vector pdirOfTrackAtVtx = HelixUtils.Direction(phtf, ps);
        Hep3Vector pmomAtVtx = VecOp.mult(ptrackMom, VecOp.unit(pdirOfTrackAtVtx));

        double mass = invMass(emomAtVtx, pmomAtVtx);
        v.setVertexTrackParameters(emomAtVtx, pmomAtVtx, mass);
    }

    // invMass  is probably defined in the code a hundred times...here it is again. 
    private double invMass(Hep3Vector p1, Hep3Vector p2) {
        double me2 = 0.000511 * 0.000511;
        double esum = sqrt(p1.magnitudeSquared() + me2) + sqrt(p2.magnitudeSquared() + me2);
        double pxsum = p1.x() + p2.x();
        double pysum = p1.y() + p2.y();
        double pzsum = p1.z() + p2.z();

        double psum = Math.sqrt(pxsum * pxsum + pysum * pysum + pzsum * pzsum);
        double evtmass = esum * esum - psum * psum;

        if (evtmass > 0)
            return Math.sqrt(evtmass);
        else
            return -99;
    }

}
