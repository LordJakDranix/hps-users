package org.hps.users.mrsolt;

import hep.physics.vec.BasicHep3Matrix;
import hep.physics.vec.Hep3Vector;
import hep.physics.vec.VecOp;

import java.util.List;

import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.util.Driver;

/**
 * Driver to skim selected events from LCIO files based on z vertex and mass
 *
 * @author Matt Solt
 *
 * @version $Id:
 */
public class GeneralLcioEventSkimmer extends Driver{

    private boolean skipEvent = true;
    private int _numberOfEventsWritten;
    private final String unconstrainedV0CandidatesColName = "UnconstrainedV0Candidates";
    protected final BasicHep3Matrix beamAxisRotation = BasicHep3Matrix.identity();
    private double zMax = -50;
    private double massMin = 0.0;
    private double massMax = 1.0;
    
    public void setZMax(double zMax) {
        this.zMax = zMax;
    }
    
    public void setMassMin(double massMin) {
        this.massMin = massMin;
    }
    
    public void setMassMax(double massMax) {
        this.massMax = massMax;
    }

    @Override
    protected void startOfData(){
    	beamAxisRotation.setActiveEuler(Math.PI / 2, -0.0305, -Math.PI / 2);
    }

    @Override
    protected void process(EventHeader event){
        skipEvent = true;
        
        List<ReconstructedParticle> unConstrainedV0List = event.get(ReconstructedParticle.class, unconstrainedV0CandidatesColName);

        for (ReconstructedParticle uncV0 : unConstrainedV0List) {
        	Hep3Vector uncVtx = VecOp.mult(beamAxisRotation, uncV0.getStartVertex().getPosition());
        	double z = uncVtx.z();
        	double m = uncV0.getMass();
        	if(z > zMax && m > massMin && m < massMax){
        		skipEvent = false;
        		break;
        	}
        }
        
        if (skipEvent) {
            throw new Driver.NextEventException();
        } else {
            _numberOfEventsWritten++;
        }
    }

    @Override
    protected void endOfData(){
        System.out.println("Selected " + _numberOfEventsWritten + " events");
    }

}