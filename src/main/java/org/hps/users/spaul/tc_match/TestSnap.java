package org.hps.users.spaul.tc_match;

import java.io.IOException;
import java.util.List;
import java.util.Random;

import org.hps.conditions.beam.BeamEnergy.BeamEnergyCollection;
import org.hps.conditions.database.DatabaseConditionsManager;
import org.hps.conditions.ecal.EcalChannel.EcalChannelCollection;
import org.hps.conditions.ecal.EcalChannel;
import org.hps.conditions.ecal.EcalConditionsConverter;
import org.hps.conditions.ecal.EcalCrystalPosition;
import org.hps.conditions.ecal.EcalCrystalPosition.EcalCrystalPositionCollection;
import org.hps.conditions.svt.SvtMotorPosition;
import org.hps.recon.utils.SnapToEdge;
import org.hps.recon.utils.TrackClusterMatcher;
import org.lcsim.conditions.ConditionsManager;
import org.lcsim.conditions.ConditionsManager.ConditionsNotFoundException;
import org.lcsim.geometry.subdetector.HPSEcal3;
import org.lcsim.util.aida.AIDA;

import hep.aida.IHistogram2D;
import hep.physics.vec.BasicHep3Vector;
import hep.physics.vec.Hep3Vector;


public class TestSnap {
    public static void main(String arg[]) throws ConditionsNotFoundException, IOException{
        SnapToEdge snapper = new SnapToEdge();
        DatabaseConditionsManager manager = DatabaseConditionsManager.getInstance();
        //manager.setDetector("HPS-PhysicsRun2016-v5-3-fieldmap_globalAlign", 7796);

        manager.setDetector("HPS-PhysicsRun2016-Nominal-v5-0-fieldmap", Integer.parseInt(arg[0]));
        //manager.registerConditionsConverter(new EcalConditionsConverter());
        //manager.getEcalSubdetector();
        
        
        EcalCrystalPositionCollection positions = manager.getCachedConditions(EcalCrystalPositionCollection.class, "ecal_crystal_positions").getCachedData();
        snapper.loadEdges(positions);
        
        IHistogram2D snap = AIDA.defaultInstance().histogram2D("snap", 100, -500, 500, 100, -100, 100);
        IHistogram2D unsnap = AIDA.defaultInstance().histogram2D("unsnap", 100, -500, 500, 100, -100, 100);
        IHistogram2D old_vs_new = AIDA.defaultInstance().histogram2D("y vs new y", 100, -100, 100, 100, -100, 100);
        
        IHistogram2D crystalPos = AIDA.defaultInstance().histogram2D("crystal pos", 100, -500, 500, 100, -100, 100);

        IHistogram2D crystalIndex = AIDA.defaultInstance().histogram2D("crystal index", 50, -25, 25, 14, -7, 7);
        
        IHistogram2D crystalIndexVsX[] =new IHistogram2D[5];
        for(int i = 0; i<5; i++) 
            crystalIndexVsX[i]= AIDA.defaultInstance().histogram2D("crystal ix vs x (row" + (i+1)+")", 50, -25, 25, 100, -500, 500);
        IHistogram2D crystalIndexVsY = AIDA.defaultInstance().histogram2D("crystal iy vs y", 14, -7, 7, 100, -100, 100);
        
        
        IHistogram2D loadedEdges = AIDA.defaultInstance().histogram2D("loaded edges", 100, -500, 500, 100, -100, 100);
        
        
        Random r = new Random();
        
        double xmax = 500;
        double xmin = -500;
        double ymax = 100;
        double ymin = -100;
        for(int i =0; i<10000; i++){
            double x = (xmax-xmin)*r.nextDouble()+xmin;
            double y = (ymax-ymin)*r.nextDouble()+ymin;
            
            Hep3Vector X = new BasicHep3Vector(x,y,0);
            Hep3Vector Xnew = snapper.snapToEdge(X);
            double xnew = Xnew.x();
            double ynew = Xnew.y();
            
            snap.fill(xnew, ynew);
            unsnap.fill(x, y);
            old_vs_new.fill(y, ynew);
            
        }
        
        /*for(EcalCrystalPosition pos : positions){
            EcalChannel chan = channels.findChannel(pos.getChannelId());
            int ix = chan.getX();
            int iy = chan.getY();
            crystalIndex.fill(ix, iy);
            
            double cx = pos.getFrontX();
            double cy = pos.getFrontY();
            crystalPos.fill(cx, cy);
            if(iy>0)
            crystalIndexVsX[iy-1].fill(ix, cx);
            crystalIndexVsY.fill(iy, cy);
        }*/
        
        /*for(int i = 0; i<46; i++){
            loadedEdges.fill(matcher.x_inner_top[i], matcher.y_inner_top[i]);
        }*/
        
        AIDA.defaultInstance().saveAs("out.root");
    }
}
