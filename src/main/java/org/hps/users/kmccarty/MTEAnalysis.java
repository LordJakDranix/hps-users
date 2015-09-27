package org.hps.users.kmccarty;

import hep.aida.IHistogram1D;
import hep.aida.IHistogram2D;
import hep.physics.vec.VecOp;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import org.hps.record.triggerbank.AbstractIntData;
import org.hps.record.triggerbank.TIData;
import org.hps.record.triggerbank.TriggerModule;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.GenericObject;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.util.Driver;
import org.lcsim.util.aida.AIDA;

public class MTEAnalysis extends Driver {
	// Define track LCIO information.
	private boolean skipBadSVT = true;
	private String bankCollectionName = "TriggerBank";
	private String particleCollectionName = "FinalStateParticles";
	private static final AIDA aida = AIDA.defaultInstance();
	private IHistogram1D[] chargedTracksPlot = {
			aida.histogram1D("MTE Analysis/Møller Event Tracks", 10, -0.5, 9.5),
			aida.histogram1D("MTE Analysis/Trident Event Tracks", 10, -0.5, 9.5),
			aida.histogram1D("MTE Analysis/Elastic Event Tracks", 10, -0.5, 9.5)
	};
	private IHistogram1D[] clusterCountPlot = {
			aida.histogram1D("MTE Analysis/Møller Event Clusters", 10, -0.5, 9.5),
			aida.histogram1D("MTE Analysis/Trident Event Clusters", 10, -0.5, 9.5),
			aida.histogram1D("MTE Analysis/Elastic Event Clusters", 10, -0.5, 9.5)
	};
	private IHistogram1D[] energyPlot = {
			aida.histogram1D("MTE Analysis/Møller Energy Sum Distribution", 220, 0, 2.2),
			aida.histogram1D("MTE Analysis/Trident Energy Sum Distribution", 220, 0, 2.2),
			aida.histogram1D("MTE Analysis/Elastic Energy Distribution", 110, 0, 1.5)
	};
	private IHistogram1D[] electronPlot = {
			aida.histogram1D("MTE Analysis/Møller Electron Energy Distribution", 220, 0, 2.2),
			aida.histogram1D("MTE Analysis/Trident Electron Energy Distribution", 220, 0, 2.2),
	};
	private IHistogram1D positronPlot = aida.histogram1D("MTE Analysis/Trident Positron Energy Distribution", 220, 0, 2.2);
	private IHistogram2D[] energy2DPlot = {
			aida.histogram2D("MTE Analysis/Møller 2D Energy Distribution", 55, 0, 1.1, 55, 0, 1.1),
			aida.histogram2D("MTE Analysis/Trident 2D Energy Distribution", 55, 0, 1.1, 55, 0, 1.1),
	};
	private IHistogram1D timePlot = aida.histogram1D("MTE Analysis/Track Cluster Time Distribution", 4000, 0, 400);
	private IHistogram1D timeCoincidencePlot = aida.histogram1D("MTE Analysis/Møller Time Coincidence Distribution", 1000, 0, 100);
	private IHistogram1D timeCoincidenceAllCutsPlot = aida.histogram1D("MTE Analysis/Møller Time Coincidence Distribution (All Møller Cuts)", 1000, 0, 100);
	private IHistogram1D negTrackCount = aida.histogram1D("MTE Analysis/All Negative Tracks", 10, -0.5, 9.5);
	private IHistogram1D posTrackCount = aida.histogram1D("MTE Analysis/All Positive Event Tracks", 10, -0.5, 9.5);
	private IHistogram1D chargedTrackCount = aida.histogram1D("MTE Analysis/All Event Event Tracks", 10, -0.5, 9.5);
	
	private IHistogram1D trTimeCoincidenceAll          = aida.histogram1D("Trident/Time Coincidence",                   150, 0.0, 15.0);
	private IHistogram1D trTimeCoincidenceFiducial     = aida.histogram1D("Trident/Time Coincidence (Fiducial Region)", 150, 0.0, 15.0);
	private IHistogram1D trEnergySumAll                = aida.histogram1D("Trident/Energy Sum",                         220, 0.0,  1.1);
	private IHistogram1D trEnergySumFiducial           = aida.histogram1D("Trident/Energy Sum (Fiducial Region)",       220, 0.0,  1.1);
	private IHistogram2D trEnergySum2DAll              = aida.histogram2D("Trident/First Cluster Energy vs. Second Cluster Energy",                   220, 0, 1.1, 220,   0, 1.1);
	private IHistogram2D trEnergySum2DFiducial         = aida.histogram2D("Trident/First Cluster Energy vs. Second Cluster Energy (Fiducial Region)", 220, 0, 1.1, 220,   0, 1.1);
	private IHistogram2D trSumCoplanarityAll           = aida.histogram2D("Trident/Hardware Coplanarity vs. Energy Sum",                              220, 0, 1.1, 115,   0, 180);
	private IHistogram2D trSumCoplanarityFiducial      = aida.histogram2D("Trident/Hardware Coplanarity vs. Energy Sum (Fiducial Region)",            220, 0, 1.1, 115,   0, 180);
	private IHistogram2D trSumCoplanarityCalcAll       = aida.histogram2D("Trident/Calculated Coplanarity vs. Energy Sum",                            220, 0, 1.1, 115, 130, 230);
	private IHistogram2D trSumCoplanarityCalcFiducial  = aida.histogram2D("Trident/Calculated Coplanarity vs. Energy Sum (Fiducial Region)",          220, 0, 1.1, 115, 130, 230);
	private IHistogram2D trTimeEnergyAll               = aida.histogram2D("Trident/Cluster Time vs. Cluster Energy",                                  220, 0, 1.1, 100,   0, 100);
	private IHistogram2D trTimeEnergyFiducial          = aida.histogram2D("Trident/Cluster Time vs. Cluster Energy (Fiducial Region)",                220, 0, 1.1, 100,   0, 100);
	
	private TriggerPlotsModule allPlots = new TriggerPlotsModule("All");
	private TriggerPlotsModule møllerPlots = new TriggerPlotsModule("Møller");
	private TriggerPlotsModule tridentPlots = new TriggerPlotsModule("Trident");
	private TriggerPlotsModule elasticPlots = new TriggerPlotsModule("Elastic");
	private static final int MØLLER  = 0;
	private static final int TRIDENT = 1;
	private static final int ELASTIC = 2;
	private boolean verbose = false;
	private boolean excludeNoTrackEvents = false;
	private double timeCoincidenceCut = Double.MAX_VALUE;
	private Map<String, Integer> møllerBitMap = new HashMap<String, Integer>();
	private Map<String, Integer> tridentBitMap = new HashMap<String, Integer>();
	private Map<String, Integer> elasticBitMap = new HashMap<String, Integer>();
	private int møllerEvents = 0;
	private int tridentEvents = 0;
	private int elasticEvents = 0;
	private int totalEvents = 0;
	private int pair1Events = 0;
	private int pair0Events = 0;
	private int singles1Events = 0;
	private int singles0Events = 0;
	private int pulserEvents = 0;
	
	@Override
	public void startOfData() {
		for(int s0 = 0; s0 <= 1; s0++) {
			for(int s1 = 0; s1 <= 1; s1++) {
				for(int p0 = 0; p0 <= 1; p0++) {
					for(int p1 = 0; p1 <= 1; p1++) {
						for(int pulser = 0; pulser <=1; pulser++) {
							// Set each "trigger bit."
							boolean s0bit = (s0 == 1);
							boolean s1bit = (s1 == 1);
							boolean p0bit = (p0 == 1);
							boolean p1bit = (p1 == 1);
							boolean pulserBit = (p1 == 1);
							
							// Generate the bit string.
							String bitString = getBitString(s0bit, s1bit, p0bit, p1bit, pulserBit);
							
							// Set a default value of zero for this bit combination.
							møllerBitMap.put(bitString, 1);
							tridentBitMap.put(bitString, 1);
							elasticBitMap.put(bitString, 1);
						}
					}
				}
			}
		}
	}
	
	@Override
	public void endOfData() {
		System.out.println("Møller  Events   :: " + møllerEvents);
		System.out.println("Trident Events   :: " + tridentEvents);
		System.out.println("Elastic Events   :: " + elasticEvents);
		System.out.println("Total Events     :: " + totalEvents);
		System.out.println("Pair 1 Events    :: " + pair1Events);
		System.out.println("Pair 0 Events    :: " + pair0Events);
		System.out.println("Singles 1 Events :: " + singles1Events);
		System.out.println("Singles 0 Events :: " + singles0Events);
		System.out.println("Pulser Events    :: " + pulserEvents);
		
		System.out.println("Plsr\tS0\tS1\tP0\tP1\tMøller");
		for(Entry<String, Integer> entry : møllerBitMap.entrySet()) {
			System.out.println(entry.getKey() + "\t" + entry.getValue());
		}
		
		System.out.println("Plsr\tS0\tS1\tP0\tP1\tTrident");
		for(Entry<String, Integer> entry : tridentBitMap.entrySet()) {
			System.out.println(entry.getKey() + "\t" + entry.getValue());
		}
		
		System.out.println("Plsr\tS0\tS1\tP0\tP1\tElastic");
		for(Entry<String, Integer> entry : elasticBitMap.entrySet()) {
			System.out.println(entry.getKey() + "\t" + entry.getValue());
		}
	}
	
	private static final String getBitString(boolean s0, boolean s1, boolean p0, boolean p1, boolean pulser) {
		return String.format("%d\t%d\t%d\t%d\t%d", (pulser ? 1 : 0), (s0 ? 1 : 0), (s1 ? 1 : 0), (p0 ? 1 : 0), (p1 ? 1 : 0));
	}
	
	@Override
	public void process(EventHeader event) {
		// Check whether the SVT was active in this event.
		final String[] flagNames = { "svt_bias_good", "svt_burstmode_noise_good", "svt_position_good" };
		boolean svtGood = true;
        for(int i = 0; i < flagNames.length; i++) {
            int[] flag = event.getIntegerParameters().get(flagNames[i]);
            if(flag == null || flag[0] == 0) {
                svtGood = false;
            }
        }
        
        // If the SVT was bad, then skip the event.
        if(!svtGood && skipBadSVT) {
        	return;
        }
        
		if(event.hasCollection(ReconstructedParticle.class, particleCollectionName)) {
			// Get the list of tracks.
			List<ReconstructedParticle> trackList = event.get(ReconstructedParticle.class, particleCollectionName);
			
			// Plot the time stamps of all tracks.
			for(ReconstructedParticle track : trackList) {
				if(track.getClusters().size() != 0) {
					Cluster cluster = track.getClusters().get(0);
					timePlot.fill(cluster.getCalorimeterHits().get(0).getTime());
				}
			}
			
			if(verbose) {
				System.out.println(trackList.size() + " tracks found.");
				for(ReconstructedParticle track : trackList) {
					System.out.printf("Track :: Q = %4.1f; E = %6.3f%n",
							track.getCharge(), track.getEnergy());
				}
			}
			
			// Populate the all cluster plots.
			List<Cluster> topClusters = new ArrayList<Cluster>();
			List<Cluster> botClusters = new ArrayList<Cluster>();
			List<Cluster> clusters = event.get(Cluster.class, "EcalClusters");
			for(Cluster cluster : clusters) {
				allPlots.addCluster(cluster);
				if(cluster.getCalorimeterHits().get(0).getIdentifierFieldValue("iy") > 0) { topClusters.add(cluster); }
				else { botClusters.add(cluster); }
			}
			
			// Make cluster pairs.
			List<Cluster[]> clusterPairs = new ArrayList<Cluster[]>();
			for(Cluster topCluster : topClusters) {
				for(Cluster botCluster : botClusters) {
					clusterPairs.add(new Cluster[] { topCluster, botCluster });
				}
			}
			
			// Populate the all cluster pair plots.
			for(Cluster[] pair : clusterPairs) {
				allPlots.addClusterPair(pair);
			}
			
			// Check each of the event-type conditions.
			boolean isMøller = false;
			boolean isTrident = false;
			boolean isElastic = false;
			
			// Produce all possible pairs of tracks.
			List<ReconstructedParticle[]> pairList = getTrackPairs(trackList);
			
			// Check the Møller condition. A Møller event is expected
			// to have two tracks, both negative, with a net energy
			// within a certain band of the beam energy.
			møllerTrackLoop:
			for(ReconstructedParticle[] pair : pairList) {
				// If trackless events are to be excluded, then require
				// that each "track" have a real track.
				if(excludeNoTrackEvents && (pair[0].getTracks().isEmpty() || pair[1].getTracks().isEmpty())) {
					continue møllerTrackLoop;
				}
				
				// Both tracks are required to be negatively charged.
				if(pair[0].getCharge() >= 0 || pair[1].getCharge() >= 0) {
					continue møllerTrackLoop;
				}
				
				// Both tracks must have clusters associated with them.
				Cluster[] trackClusters = new Cluster[2];
				for(int i = 0; i < 2; i++) {
					// Disallow tracks with no associated clusters.
					if(pair[i].getClusters().size() == 0) {
						continue møllerTrackLoop;
					}
					
					// Store the first cluster associated with the track.
					trackClusters[i] = pair[i].getClusters().get(0);
				}
				
				// Require that the track clusters be within a certain
				// time window of one another.
				CalorimeterHit[] seeds = new CalorimeterHit[2];
				seeds[0] = trackClusters[0].getCalorimeterHits().get(0);
				seeds[1] = trackClusters[1].getCalorimeterHits().get(0);
				timeCoincidencePlot.fill(Math.abs(seeds[0].getTime() - seeds[1].getTime()));
				if(Math.abs(trackClusters[0].getCalorimeterHits().get(0).getTime() - trackClusters[1].getCalorimeterHits().get(0).getTime()) > timeCoincidenceCut) {
					continue møllerTrackLoop;
				}
				
				// Require both tracks to occur within the range of
				// 36.5 and 49 ns.
				if(seeds[0].getTime() < 36.5 || seeds[0].getTime() > 49) {
					continue møllerTrackLoop;
				} if(seeds[1].getTime() < 36.5 || seeds[1].getTime() > 49) {
					continue møllerTrackLoop;
				}
				
				// No track may have an energy that exceeds 900 MeV.
				if(pair[0].getMomentum().magnitude() >= 0.900 || pair[1].getMomentum().magnitude() >= 0.900) {
					continue møllerTrackLoop;
				}
				
				// Get the energy sum.
				double sum = VecOp.add(pair[0].getMomentum(), pair[1].getMomentum()).magnitude();
				
				// "Møller-like" track pairs must have energies within
				// an allowed energy range.
				if(sum < 0.800 || sum > 1.500) {
					continue møllerTrackLoop;
				}
				
				//timeCoincidenceAllCutsPlot.fill(Math.abs(seeds[0].getTime() - seeds[1].getTime()));
				
				// Note that this is a Møller event.
				isMøller = true;
				
				// Populate the Møller plots.
				energyPlot[MØLLER].fill(sum);
				møllerPlots.addClusterPair(trackClusters);
				electronPlot[MØLLER].fill(pair[0].getMomentum().magnitude());
				electronPlot[MØLLER].fill(pair[1].getMomentum().magnitude());
				energy2DPlot[MØLLER].fill(pair[0].getMomentum().magnitude(), pair[1].getMomentum().magnitude());
			}
			
			// Check the elastic condition. Elastic events should be
			// negatively and have an energy approximately equal to
			// the beam energy.
			elasticTrackLoop:
			for(ReconstructedParticle track : trackList) {
				// If trackless events are to be excluded, then require
				// that the "track" has a real track.
				if(excludeNoTrackEvents && track.getTracks().isEmpty()) {
					continue elasticTrackLoop;
				}
				
				// Check the elastic condition.
				if(track.getCharge() < 0 && track.getMomentum().magnitude() >= 0.900) {
					isElastic = true;
					energyPlot[ELASTIC].fill(track.getMomentum().magnitude());
					if(!track.getClusters().isEmpty()) {
						elasticPlots.addCluster(track.getClusters().get(0));
					}
				}
			}
			
			// Check the trident condition. Tridents are events that
			// contain both one positive and one negative track.
			tridentTrackLoop:
			for(ReconstructedParticle[] pair : pairList) {
				// If trackless events are to be excluded, then require
				// that each "track" have a real track.
				if(excludeNoTrackEvents && (pair[0].getTracks().isEmpty() || pair[1].getTracks().isEmpty())) {
					continue tridentTrackLoop;
				}
				
				// Check the trident condition.
				boolean isPosNeg = (pair[0].getCharge() < 0 && pair[1].getCharge() > 0) || (pair[0].getCharge() > 0 && pair[1].getCharge() < 0);
				if(!isPosNeg) { continue tridentTrackLoop; }
				
				// Both tracks must have clusters associated with them.
				Cluster[] trackClusters = new Cluster[pair.length];
				for(int i = 0; i < pair.length; i++) {
					// Disallow tracks with no associated clusters.
					if(pair[i].getClusters().size() == 0) {
						continue tridentTrackLoop;
					}
					
					// Store the first cluster associated with the track.
					trackClusters[i] = pair[i].getClusters().get(0);
				}
				
				// Make sure that the clusters are not the same.
				if(trackClusters[0] == trackClusters[1]) {
					continue tridentTrackLoop;
				}
				
				// Require that the track clusters be within a certain
				// time window of one another.
				CalorimeterHit[] seeds = new CalorimeterHit[2];
				seeds[0] = trackClusters[0].getCalorimeterHits().get(0);
				seeds[1] = trackClusters[1].getCalorimeterHits().get(0);
				timeCoincidencePlot.fill(Math.abs(seeds[0].getTime() - seeds[1].getTime()));
				if(Math.abs(trackClusters[0].getCalorimeterHits().get(0).getTime() - trackClusters[1].getCalorimeterHits().get(0).getTime()) > timeCoincidenceCut) {
					continue tridentTrackLoop;
				}
				
				// Require that the energy of the electron is below
				// 900 MeV.
				boolean electronNotElastic = (pair[0].getCharge() < 0 && pair[0].getMomentum().magnitude() < 0.900)
						|| (pair[1].getCharge() < 0 && pair[1].getMomentum().magnitude() < 0.900);
				if(!electronNotElastic) {
					continue tridentTrackLoop;
				}
				
				// If all tests are passed, this is a trident. Note
				// this and populate the trident plots.
				isTrident = true;
				tridentPlots.addClusterPair(trackClusters);
				if(pair[0].getCharge() > 0) {
					positronPlot.fill(pair[1].getMomentum().magnitude());
					electronPlot[TRIDENT].fill(pair[0].getMomentum().magnitude());
				} else {
					positronPlot.fill(pair[0].getMomentum().magnitude());
					electronPlot[TRIDENT].fill(pair[1].getMomentum().magnitude());
				}
				energyPlot[TRIDENT].fill(VecOp.add(pair[0].getMomentum(), pair[1].getMomentum()).magnitude());
				energy2DPlot[TRIDENT].fill(pair[0].getMomentum().magnitude(), pair[1].getMomentum().magnitude());
				
				// Track which clusters have already been added to the
				// singles plot so that there are no repeats.
				Set<Cluster> plotSet = new HashSet<Cluster>();
				Set<Cluster> plotFiducial = new HashSet<Cluster>();
				
				// Fill the all pairs plots.
				double pairEnergy = trackClusters[0].getEnergy() + trackClusters[1].getEnergy();
				trEnergySumAll.fill(pairEnergy);
				trEnergySum2DAll.fill(trackClusters[1].getEnergy(), trackClusters[0].getEnergy());
				trTimeCoincidenceAll.fill(TriggerModule.getValueTimeCoincidence(trackClusters));
				trSumCoplanarityCalcAll.fill(pairEnergy, getCalculatedCoplanarity(trackClusters));
				trSumCoplanarityAll.fill(pairEnergy, TriggerModule.getValueCoplanarity(trackClusters));
				
				// Fill the singles plots.
				if(!plotSet.contains(trackClusters[0])) {
					plotSet.add(trackClusters[0]);
					trTimeEnergyAll.fill(trackClusters[0].getEnergy(), TriggerModule.getClusterTime(trackClusters[0]));
				} if(!plotSet.contains(trackClusters[1])) {
					plotSet.add(trackClusters[1]);
					trTimeEnergyAll.fill(trackClusters[1].getEnergy(), TriggerModule.getClusterTime(trackClusters[1]));
				}
				
				// Fill the fiducial plots if appropriate.
				if(inFiducialRegion(trackClusters[0]) && inFiducialRegion(trackClusters[1])) {
					trEnergySumFiducial.fill(pairEnergy);
					trEnergySum2DFiducial.fill(trackClusters[1].getEnergy(), trackClusters[0].getEnergy());
					trTimeCoincidenceFiducial.fill(TriggerModule.getValueTimeCoincidence(trackClusters));
					trSumCoplanarityCalcFiducial.fill(pairEnergy, getCalculatedCoplanarity(trackClusters));
					trSumCoplanarityFiducial.fill(pairEnergy, TriggerModule.getValueCoplanarity(trackClusters));
				}
				
				// Fill the singles fiducial plots if appropriate.
				if(!plotFiducial.contains(trackClusters[0]) && inFiducialRegion(trackClusters[0])) {
					plotFiducial.add(trackClusters[0]);
					trTimeEnergyFiducial.fill(trackClusters[0].getEnergy(), TriggerModule.getClusterTime(trackClusters[0]));
				} if(!plotFiducial.contains(trackClusters[1]) && inFiducialRegion(trackClusters[1])) {
					plotFiducial.add(trackClusters[1]);
					trTimeEnergyFiducial.fill(trackClusters[1].getEnergy(), TriggerModule.getClusterTime(trackClusters[1]));
				}
			}
			
			if(verbose) {
				System.out.printf("\tMøller  :: %b%n", isMøller);
				System.out.printf("\tTrident :: %b%n", isTrident);
				System.out.printf("\tElastic :: %b%n", isElastic);
				System.out.println();
			}
			
			// Get the TI bits.
			String bitString = null;
			TIData tiBank = null;
			List<GenericObject> bankList = event.get(GenericObject.class, bankCollectionName);
			for(GenericObject obj : bankList) {
				if(AbstractIntData.getTag(obj) == TIData.BANK_TAG) {
					tiBank = new TIData(obj);
					bitString = getBitString(tiBank.isPulserTrigger(), tiBank.isSingle0Trigger(),
							tiBank.isSingle1Trigger(), tiBank.isPair0Trigger(), tiBank.isPair1Trigger());
					
					if(tiBank.isPair1Trigger()) {
						pair1Events++;
					} else if(tiBank.isPair0Trigger()) {
						pair0Events++;
					} else if(tiBank.isSingle1Trigger()) {
						singles1Events++;
					} else if(tiBank.isSingle0Trigger()) {
						singles0Events++;
					} else if(tiBank.isPulserTrigger()) {
						pulserEvents++;
					}
				}
			}
			if(bitString == null) {
				System.out.println("No TI data found!!");
			}
			
			// Get the number of charged tracks in the event.
			int tracks = 0;
			int posTracks = 0;
			int negTracks = 0;
			for(ReconstructedParticle track : trackList) {
				if(track.getCharge() != 0 && tiBank.isPulserTrigger()) {
					if(excludeNoTrackEvents && !track.getTracks().isEmpty()) {
						tracks++;
						if(track.getCharge() > 0) { posTracks++; }
						else { negTracks++; }
					} else {
						tracks++;
						if(track.getCharge() > 0) { posTracks++; }
						else { negTracks++; }
					}
				}
			}
			
			// Populate the "all tracks" plots.
			posTrackCount.fill(posTracks);
			negTrackCount.fill(negTracks);
			chargedTrackCount.fill(tracks);
			
			// Add the result to the appropriate plots and increment
			// the appropriate trigger bit combination.
			if(isMøller) {
				møllerEvents++;
				chargedTracksPlot[MØLLER].fill(tracks);
				clusterCountPlot[MØLLER].fill(clusters.size());
				
				Integer val = møllerBitMap.get(bitString);
				if(val == null) { møllerBitMap.put(bitString, 1); }
				else { møllerBitMap.put(bitString, val + 1); }
			} else if(isTrident) {
				tridentEvents++;
				chargedTracksPlot[TRIDENT].fill(tracks);
				clusterCountPlot[TRIDENT].fill(clusters.size());
				
				Integer val = tridentBitMap.get(bitString);
				if(val == null) { tridentBitMap.put(bitString, 1); }
				else { tridentBitMap.put(bitString, val + 1); }
			} else if(isElastic) {
				elasticEvents++;
				chargedTracksPlot[ELASTIC].fill(tracks);
				clusterCountPlot[ELASTIC].fill(clusters.size());
				
				Integer val = elasticBitMap.get(bitString);
				if(val == null) { elasticBitMap.put(bitString, 1); }
				else { elasticBitMap.put(bitString, val + 1); }
			}
			totalEvents++;
		}
	}
	
	private static final double getCalculatedCoplanarity(Cluster[] pair) {
		// Define the x- and y-coordinates of the clusters as well as
		// calorimeter center.
		final double ORIGIN_X = 42.52;
		double x[] = { pair[0].getPosition()[0], pair[1].getPosition()[0] };
		double y[] = { pair[0].getPosition()[1], pair[1].getPosition()[1] };
		
        // Get the cluster angles.
        double[] clusterAngle = new double[2];
        for(int i = 0; i < 2; i++) {
        	clusterAngle[i] = Math.atan2(y[i], x[i] - ORIGIN_X) * 180 / Math.PI;
        	if(clusterAngle[i] <= 0) { clusterAngle[i] += 360; }
        }
        
        // Calculate the coplanarity cut value.
        double clusterDiff = clusterAngle[0] - clusterAngle[1];
        return clusterDiff > 0 ? clusterDiff : clusterDiff + 360;
	}
	
	public void setTimeCoincidenceCut(double value) {
		timeCoincidenceCut = value;
	}
	
	public void setExcludeNoTrackEvents(boolean state) {
		excludeNoTrackEvents = state;
	}
	
	public void setSkipBadSVT(boolean state) {
		skipBadSVT = state;
	}
	
	private static final boolean inFiducialRegion(Cluster cluster) {
		// Get the x and y indices for the cluster.
		int ix   = TriggerModule.getClusterXIndex(cluster);
		int absx = Math.abs(TriggerModule.getClusterXIndex(cluster));
		int absy = Math.abs(TriggerModule.getClusterYIndex(cluster));
		
		// Check if the cluster is on the top or the bottom of the
		// calorimeter, as defined by |y| == 5. This is an edge cluster
		// and is not in the fiducial region.
		if(absy == 5) {
			return false;
		}
		
		// Check if the cluster is on the extreme left or right side
		// of the calorimeter, as defined by |x| == 23. This is also
		// and edge cluster is not in the fiducial region.
		if(absx == 23) {
			return false;
		}
		
		// Check if the cluster is along the beam gap, as defined by
		// |y| == 1. This is an internal edge cluster and is not in the
		// fiducial region.
		if(absy == 1) {
			return false;
		}
		
		// Lastly, check if the cluster falls along the beam hole, as
		// defined by clusters with -11 <= x <= -1 and |y| == 2. This
		// is not the fiducial region.
		if(absy == 2 && ix <= -1 && ix >= -11) {
			return false;
		}
		
		// If all checks fail, the cluster is in the fiducial region.
		return true;
	}
	
	private static final List<ReconstructedParticle[]> getTrackPairs(List<ReconstructedParticle> trackList) {
		// Create an empty list for the pairs.
		List<ReconstructedParticle[]> pairs = new ArrayList<ReconstructedParticle[]>();
		
		// Add all possible pairs of tracks.
		for(int i = 0; i < trackList.size(); i++) {
			for(int j = i + 1; j < trackList.size(); j++) {
				pairs.add(new ReconstructedParticle[] { trackList.get(i), trackList.get(j) });
			}
		}
		
		// Return the list of tracks.
		return pairs;
	}
}
