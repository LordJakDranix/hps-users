package org.hps.users.byale;

import hep.aida.IAnalysisFactory;
import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;
import hep.aida.IHistogramFactory;
import hep.aida.IPlotter;
import hep.aida.IPlotterFactory;
//import hep.aida.ITree;
//import hep.aida.ref.rootwriter.RootFileStore;

//import java.io.IOException;
import java.util.List;

import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.MCParticle;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.SimCalorimeterHit;
import org.lcsim.event.SimTrackerHit;
import org.lcsim.event.Vertex;
//import org.lcsim.event.ReconstructedParticle;
//import org.lcsim.event.base.BaseMCParticle;
import org.lcsim.geometry.Detector;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

/**
 *
 * @author byale
 */
public class MollerAnalysis extends Driver {

	//ITree tree;
	IHistogramFactory histogramFactory;
    private AIDA aida = AIDA.defaultInstance();
    private String collectionName = "MCParticle";
    IPlotter plotterMC;
    IHistogram1D px;
    IHistogram1D py;
    IHistogram1D pz;
    IHistogram1D SeedE;
    IHistogram1D EDiff;
    
    IHistogram2D EposEele;
    IHistogram1D ePos;
    IHistogram1D eEle1;
    IHistogram1D eEle2;
    IHistogram1D EepSum;
    IHistogram1D EepSumEndpoint;
    IHistogram1D EeepSum;
    IHistogram1D EepDiff;
    IHistogram2D ECal;
    
    IHistogram2D ECalPosition;
    IHistogram1D ECalPositionX;
    IHistogram1D ECalPositionY;
    IHistogram1D ECalPositionZ;
    IHistogram1D ECalPositionZDiff;
    IHistogram1D trackMCDiff;
    IHistogram1D nMC;    
    IHistogram1D nClusters;
    IHistogram1D massMC;
    IHistogram1D zvtx_uc;
    
@Override
    protected void detectorChanged(Detector detector) {
        System.out.println("MCParicleAnalysis::detectorChanged  Setting up the plotter");
        
     //   tree = IAnalysisFactory.create().createTreeFactory().create();
     //   histogramFactory = IAnalysisFactory.create().createHistogramFactory(tree);
        
        IAnalysisFactory fac = aida.analysisFactory();
        IPlotterFactory pfac = fac.createPlotterFactory("MC Particles");

        aida.tree().cd("/");
//        resetOccupancyMap(); // this is for calculatin
        plotterMC = pfac.create("MC Particle Momentum");
        plotterMC.createRegions(3, 3);
        px = aida.histogram1D("px (GeV)", 100, -0.05, 0.14);
        py = aida.histogram1D("py (GeV)", 100, -0.04, 0.04);
       // pz = aida.histogram1D("pz (GeV)", 200, -0.25, 2.5);
        pz = aida.histogram1D("pz (GeV)", 100, 1, 1.3);
        SeedE = aida.histogram1D("ECalHitCorrEnergy", 200, 0.05, 0.5);
        EDiff = aida.histogram1D("MCMomentum - hitEnergy", 200, 0, 1);
    //    ECalPosition = aida.histogram2D("ECal position", 49, -24.5, 24.5, 13, -6.5, 6.5);
        ECalPosition = aida.histogram2D("ECal position", 50, -200, 100, 50, -100, 100);
        ECalPositionX = aida.histogram1D("ECal positionX MC-hit", 200, -200, 100);
        ECalPositionY = aida.histogram1D("ECal positionY MC-hit", 200, -100, 100);
        ECalPositionZ = aida.histogram1D("MC endpoint z", 200, 0, 2000);
        ECalPositionZDiff = aida.histogram1D("ECal positionZ MC-hit", 500, -1000, 500);
        nMC = aida.histogram1D("nMC", 6, 0, 5);
        nClusters = aida.histogram1D("nClust", 6, 0, 5);
        massMC = aida.histogram1D("massMC", 100, 0.04, 0.06);
        zvtx_uc = aida.histogram1D("zvtx_MC", 100, -15, 5);
        
        plotterMC.region(0).plot(px);
        plotterMC.region(1).plot(py);
        plotterMC.region(2).plot(pz);
        plotterMC.region(3).plot(ECalPositionX);
        plotterMC.region(4).plot(ECalPositionY);
        plotterMC.region(5).plot(ECalPosition);
    //    plotterMC.region(6).plot(SeedE);
    //    plotterMC.region(7).plot(EDiff);
      //  plotterMC.region(7).plot(EDiff);
        plotterMC.region(6).plot(ECalPositionZ);
       // plotterMC.region(7).plot(ECalPositionZDiff);
       // plotterMC.region(7).plot(nMC);
       // plotterMC.region(8).plot(nClusters);
        plotterMC.region(7).plot(massMC);
        plotterMC.region(8).plot(zvtx_uc);
        
        plotterMC.show();

        EposEele = aida.histogram2D("Epos vs Eele", 50, 0, 2.0, 50, 0, 2.0);

    //    ECalPosition = aida.histogram2D("ECalPosition", 50, -22, 22, 50, -5.5, 5.5);

        ePos = aida.histogram1D("Epos", 50, 0, 2.0);
        eEle1 = aida.histogram1D("Eele1", 50, 0, 2.0);
        eEle2 = aida.histogram1D("Eele2", 50, 0, 2.0);
        EepSumEndpoint = aida.histogram1D("Pair energy:  Endpoint", 50, 1.5, 2.0);
        EepSum = aida.histogram1D("Pair energy", 50, 0, 2.0);
        EeepSum = aida.histogram1D("Trident energy", 50, 0, 2.0);

        // pairs1
        EepDiff = aida.histogram1D("Energy Diff", 50, 0, 2.0);


    }

    @Override
    public void process(EventHeader event) {
        /*  make sure everything is there */
        if (!event.hasCollection(MCParticle.class, collectionName))
            return;

        List<SimCalorimeterHit> ecalhits = event.get(SimCalorimeterHit.class, "EcalHits");
        List<SimTrackerHit> trackerhits = event.get(SimTrackerHit.class, "TrackerHits");
        List<MCParticle> mcpList = event.get(MCParticle.class, collectionName);
  //      List<Cluster> clusters = event.get(Cluster.class, "EcalClustersCorr");
     //   List<TargetConstrainedMollerCandidates> = event.get(Cluster.class, "EcalClustersCorr");

        //List<ReconstructedParticle> uc_mollers = event.get(ReconstructedParticle.class,"UnconstrainedMollerCandidates");
        List<ReconstructedParticle> uc_mollers = event.get(ReconstructedParticle.class,"BeamspotConstrainedMollerCandidates");
        //List<ReconstructedParticle> uc_mollers = event.get(ReconstructedParticle.class,"TargetConstrainedMollerCandidates");
        
        List<Vertex> uc_moller_verts = event.get(Vertex.class,"UnconstrainedMollerVertices");
   //     List<Vertex> uc_moller_verts = event.get(Vertex.class,"TargetConstrainedMollerVertices");
        
        //      nClusters.fill(clusters.size());
        nMC.fill(mcpList.size());
        
        MCParticle pos = null;
        MCParticle ele1 = null;
        MCParticle ele2 = null;

   //     for (Cluster cluster : clusters){
   //     	cluster.getCalorimeterHits().get(0).getCorrectedEnergy();
   //     }
        
        
    //    for (SimCalorimeterHit ecalhit : ecalhits){
      
       // 	MCParticle mcp = ((SimCalorimeterHit)cluster.getCalorimeterHits().get(0)).getMCParticle(0);
      //  	MCParticle mcp = ecalhit.getMCParticle(0);
       
       //   for (SimTrackerHit trackerhit : trackerhits){
        	
       // 	MCParticle mcp = trackerhit.getMCParticle();
        
        for (MCParticle mcp : mcpList) {
     //   for(ReconstructedParticle moller : uc_mollers){
        //	for(Vertex uc_moller_vert : uc_moller_verts){
        	        	
       // 	if (mcp.getCharge() < 0 && mcp.getMomentum().magnitude() >= 1 && mcp.getMomentum().magnitude() <= 1.3){ 
       //     	if(mcp.getEndPoint().z()>1200 && (mcp.getCharge() < 0) && (((mcp.getEndPoint().y() <= -20.5)||(mcp.getEndPoint().y() >= 20.5))&&((mcp.getEndPoint().x()<=-93.29)||(mcp.getEndPoint().x()>=28.92)) || ((mcp.getEndPoint().y() <= -33.91)||(mcp.getEndPoint().y() >= 33.91))) ){

            		//if(mcp.getEndPoint().z()>1200){
      
        // ECal Hits
               // px.fill(mcp.getPX());
               // py.fill(mcp.getPY());
               // pz.fill(mcp.getPZ());
                
              //  ECalPosition.fill(mcp.getEndPoint().x(),mcp.getEndPoint().y());
           //     ECalPositionX.fill(mcp.getEndPoint().x() - ecalhit.getPosition()[0]);
           //     ECalPositionY.fill(mcp.getEndPoint().y() - ecalhit.getPosition()[1]);
              //  ECalPositionZ.fill(mcp.getEndPoint().z());
              //  massMC.fill(mcp.getMass());
        
   //     	px.fill(moller.getMomentum().x());
   //         py.fill(moller.getMomentum().y());
   //         pz.fill(moller.getMomentum().z());
             
           //  ECalPosition.fill(mcp.getEndPoint().x(),mcp.getEndPoint().y());
        //     ECalPositionX.fill(mcp.getEndPoint().x() - ecalhit.getPosition()[0]);
        //     ECalPositionY.fill(mcp.getEndPoint().y() - ecalhit.getPosition()[1]);
           //  ECalPositionZ.fill(mcp.getEndPoint().z());
 //            massMC.fill(moller.getMass());
           //  zvtx_uc.fill(moller.getStartVertex().getPosition().z());
             zvtx_uc.fill(mcp.getOriginZ());     
                
                
                
          /*  if (!mcp.getParents().isEmpty())
                if (mcp.getParents().get(0).getPDGID() == 622)
                    if (mcp.getCharge() > 0)
                        pos = mcp;
                    else if (ele1 == null)
                        ele1 = mcp;
                    else
                        ele2 = mcp;   
            */
            
        }
        
        
      /*  if (ele1 != null && ele2 != null && pos != null) {
            EposEele.fill(ele1.getEnergy(), pos.getEnergy());
            EposEele.fill(ele2.getEnergy(), pos.getEnergy());
            //EposEele.fill(Math.max(ele2.getEnergy(),ele1.getEnergy()), pos.getEnergy());

            eEle1.fill(ele1.getEnergy());
            eEle2.fill(ele2.getEnergy());
            ePos.fill(pos.getEnergy());

            EepSum.fill(ele1.getEnergy() + pos.getEnergy());
            EepSum.fill(ele2.getEnergy() + pos.getEnergy());
            EeepSum.fill(ele1.getEnergy() + pos.getEnergy() + ele2.getEnergy());

            EepSumEndpoint.fill(ele1.getEnergy() + pos.getEnergy());
            EepSumEndpoint.fill(ele2.getEnergy() + pos.getEnergy());
            
            EepDiff.fill(ele1.getEnergy() - pos.getEnergy());
            EepDiff.fill(ele2.getEnergy() - pos.getEnergy());

        } */
  
        
    }
  //      public void endOfData() { 
        
  //      String rootFile = "_MollerSLIC_analysis.root";
  //      RootFileStore store = new RootFileStore(rootFile);
  //      try {
  //          store.open();
  //          store.add(tree);
  //          store.close(); 
  //      } catch (IOException e) {
  //          e.printStackTrace();
  //      }
  //  }
}
