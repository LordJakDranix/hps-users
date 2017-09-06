package org.hps.users.tvm;

import org.hps.analysis.examples.MCTruthExampleDriver;
import org.lcsim.event.EventHeader;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;


public class MCAnalysisDriver extends Driver {
    
    public MCAnalysisDriver() {
        // add MC truth print as sub-driver
        this.add(new MCTruthExampleDriver());
    }
    
    public void detectorChanged(Detector detector) {
        // detector changed hook
    }
    
    public void startOfData() {
        // start of data hook
    }
    
    public void process(EventHeader event) {
        // process an event
        super.process(event);        
    }

}
