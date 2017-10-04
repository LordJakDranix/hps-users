package org.hps.users.tvm;

import hep.aida.IHistogram1D;

import java.util.List;

import org.hps.analysis.examples.MCTruthExampleDriver;
import org.lcsim.event.EventHeader;
import org.lcsim.event.MCParticle;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;


public class MCAnalysisDriver extends Driver {
    
    private static AIDA PLOT = AIDA.defaultInstance();
    
    private IHistogram1D particleEnergyH1D;
    
    public MCAnalysisDriver() {
        // add MC truth print as sub-driver
        this.add(new MCTruthExampleDriver());
    }
    
    public void detectorChanged(Detector detector) {
        // detector changed hook
    }
    
    public void startOfData() {    
        particleEnergyH1D = PLOT.histogram1D(
                "MCParticle Energy", /* plot name */ 
                200 /* n bins */, 
                0.0 /* lower edge */, 
                2.0 /* upper edge */);    
    }
    
    public void process(EventHeader event) {
        // call child driver
        super.process(event);     
        
        List<MCParticle> particles = event.get(
                MCParticle.class /* type of object */, 
                "MCParticle" /* collection name */);
        
        for (MCParticle p : particles) {
            particleEnergyH1D.fill(p.getEnergy());
        }
    }

}
