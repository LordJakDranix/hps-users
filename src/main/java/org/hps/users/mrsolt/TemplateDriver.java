/**
 * Template Driver
 */
/**
 * @author mrsolt
 *
 */
package org.hps.users.mrsolt;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import hep.aida.IAnalysisFactory;
import hep.aida.IHistogramFactory;
import hep.aida.IPlotterFactory;
import hep.aida.IPlotter;
import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;
import hep.aida.ITree;
import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;

import org.lcsim.detector.tracker.silicon.HpsSiSensor;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.EventHeader;
import org.lcsim.event.LCRelation;
import org.lcsim.event.RawTrackerHit;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;
import org.hps.conditions.beam.BeamEnergy.BeamEnergyCollection;
import org.hps.recon.tracking.FittedRawTrackerHit;

//Change "TemplateDriver" to your dirver name
public class TemplateDriver extends Driver {


    // Use JFreeChart as the default plotting backend
    static { 
        hep.aida.jfree.AnalysisFactory.register();
    }

    // Plotting
    protected AIDA aida = AIDA.defaultInstance();
    ITree tree; 
    IHistogramFactory histogramFactory; 

    //List of Sensors
    private List<HpsSiSensor> sensors = null;
    
    //Sample histograms
    IHistogram1D Histo_1D;
    IHistogram2D Histo_2D;  
    
    Map<String, IHistogram1D> HistoMap_1D = new HashMap<String,IHistogram1D>();
    Map<Integer, IHistogram2D> HistoMap_2D = new HashMap<Integer,IHistogram2D>();
    
    //Sample Variables
    int sampleInt = 6;
    double sampleDouble;
    boolean sampleBoolean;
    
    int[] sampleIntArr = new int[sampleInt];
    double[] sampleDoubleArr = new double[sampleInt];
    boolean[] sampleBooleanArr = new boolean[sampleInt];
    
    Hep3Vector hep3Vector = new BasicHep3Vector(0,0,0);

    Map<Integer, int[]> numberOfTopTracksMomentumLay = new HashMap<Integer,int[]>();
    Map<Integer, double[]> hitEfficiencyMomentumLayTop = new HashMap<Integer,double[]>();
    
    //Histogram Settings
    double minX = 0;
    double maxX = 10;
    double minY = 0;
    double maxY = 10;
    int nBins = 50;
    
    //Collection Strings
    private String fittedHitsCollectionName = "SVTFittedRawTrackerHits";
   
    //Constants
    public static final double CONSTANT = 0;
    private static final String SUBDETECTOR_NAME = "Tracker";
    
    //Configurable Variables
    boolean debug = false;
    double configVar = 0;

    public void setConfigVar(double configVar) { 
        this.configVar = configVar;
    }
    
    public void  setDebug(boolean debug){
        this.debug = debug;
    }
    
    //Beam Energy
    double ebeam;
    
    public void detectorChanged(Detector detector){
    	
    	aida.tree().cd("/");
    	tree = aida.tree();
        histogramFactory = IAnalysisFactory.create().createHistogramFactory(tree);
    
        //Set Beam Energy
        BeamEnergyCollection beamEnergyCollection = 
                this.getConditionsManager().getCachedConditions(BeamEnergyCollection.class, "beam_energies").getCachedData();        
        ebeam = beamEnergyCollection.get(0).getBeamEnergy();      
        
        // Get the HpsSiSensor objects from the tracker detector element
        sensors = detector.getSubdetector(SUBDETECTOR_NAME)
                          .getDetectorElement().findDescendants(HpsSiSensor.class);
   
        // If the detector element had no sensors associated with it, throw
        // an exception
        if (sensors.size() == 0) {
            throw new RuntimeException("No sensors were found in this detector.");
        }
        
        //Setup the Histograms
        Histo_1D = aida.histogram1D("Histo 1D Name", nBins, minX, maxX);
        Histo_2D = aida.histogram2D("Histo 2D Name", nBins, minX, maxX, nBins, minY, maxY);
        
        for(int i = 0; i < sampleInt; i++){
        	HistoMap_2D.put((i+1),histogramFactory.createHistogram2D("Histo Name" + (i+1), nBins, minX, maxX, nBins, minY, maxY));
        }

    }

    public void process(EventHeader event){
		aida.tree().cd("/");	
        
        // Get the list of fitted hits from the event
        List<LCRelation> fittedHits = event.get(LCRelation.class, fittedHitsCollectionName);
         
        // Map the fitted hits to their corresponding raw hits
        Map<RawTrackerHit, LCRelation> fittedRawTrackerHitMap = new HashMap<RawTrackerHit, LCRelation>();
        
        String rawTrackerHitCollectionName = "SVTRawTrackerHits";
        List<RawTrackerHit> rawHits = event.get(RawTrackerHit.class, rawTrackerHitCollectionName);
        	
        for (RawTrackerHit rawHit : rawHits) {
             
            // Access the sensor associated with the raw hit
            HpsSiSensor sensor = (HpsSiSensor) rawHit.getDetectorElement();
             
            // Retrieve the channel ID of the raw hit
            int channel = rawHit.getIdentifierFieldValue("strip");
        } 
        
                 
        for (LCRelation fittedHit : fittedHits) {
            fittedRawTrackerHitMap.put(FittedRawTrackerHit.getRawTrackerHit(fittedHit), fittedHit);
        }
         
        for (RawTrackerHit rawHit : rawHits) {
             
            // Get the hit amplitude
            double amplitude = FittedRawTrackerHit.getAmp(fittedRawTrackerHitMap.get(rawHit));
             
            // Get the t0 of the hit
            double t0 = FittedRawTrackerHit.getT0(fittedRawTrackerHitMap.get(rawHit));
        }
        
        List<CalorimeterHit> hits = event.get(CalorimeterHit.class, "EcalHits");
        for (CalorimeterHit hit : hits) {
            System.out.println("calorimeter hit has energy " + hit.getCorrectedEnergy() + " GeV.");
        }
              
        List<List<CalorimeterHit>> collections = event.get(CalorimeterHit.class);
        for (List<CalorimeterHit> calhits : collections) {    
            for (CalorimeterHit hit : calhits) {
                System.out.println("calorimeter hit has energy " + hit.getCorrectedEnergy() + " GeV.");
            }
        
        }
    }

    private void function() {
        //Insert function here
    }

    protected void startOfData() { 
    	System.out.println("PrintStart of Data Variables");
    }
    
    public void endOfData(){
        System.out.println("Print End of Data Variables");
    }
}