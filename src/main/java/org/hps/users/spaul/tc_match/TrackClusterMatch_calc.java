package org.hps.users.spaul.tc_match;

public class TrackClusterMatch_calc {
    
    double dyMeanBotElecGBL_noL6_2016[] = {-43.1459, 246.16, -536.838, 552.722, -267.945, 49.101, };
    double dySigmBotElecGBL_noL6_2016[] = {-19.1892, 186.524, -469.633, 522.885, -269.012, 52.0043, };
    double dyMeanBotPosiGBL_noL6_2016[] = {126.285, -865.441, 2254.7, -2829.59, 1728, -413.071, };
    double dySigmBotPosiGBL_noL6_2016[] = {-65.0583, 554.976, -1618.47, 2260.79, -1540.78, 411.589, };
    double dyMeanTopElecGBL_noL6_2016[] = {-7.9396, 38.3225, -78.0421, 78.1102, -37.151, 6.62897, };
    double dySigmTopElecGBL_noL6_2016[] = {11.4978, 0.512801, -40.9207, 53.0446, -24.9509, 3.98128, };
    double dyMeanTopPosiGBL_noL6_2016[] = {-46.6203, 314.252, -805.758, 990.56, -588.474, 135.657, };
    double dySigmTopPosiGBL_noL6_2016[] = {13.5897, 5.5759, -92.0985, 155.231, -106.064, 26.6547, };
    double dxMeanBotElecGBL_noL6_2016[] = {89.5302, -503.241, 1047.28, -1052.85, 509.941, -94.733, };
    double dxSigmBotElecGBL_noL6_2016[] = {48.2918, -175.963, 293.211, -242.921, 93.7292, -12.1049, };
    double dxMeanBotPosiGBL_noL6_2016[] = {18.9481, -370.259, 1577.52, -2690.48, 2049.28, -582.213, };
    double dxSigmBotPosiGBL_noL6_2016[] = {-60.611, 518.878, -1414.32, 1771.5, -1041.03, 230.035, };
    double dxMeanTopElecGBL_noL6_2016[] = {128.968, -742.489, 1590.81, -1631.63, 802.959, -151.725, };
    double dxSigmTopElecGBL_noL6_2016[] = {26.2209, -53.5331, 37.5842, 9.29055, -24.4114, 8.82913, };
    double dxMeanTopPosiGBL_noL6_2016[] = {-139.385, 778.484, -1679.36, 1798.11, -950.765, 197.982, };
    double dxSigmTopPosiGBL_noL6_2016[] = {91.8373, -442.036, 950.795, -1053.36, 584.799, -128.003, };
    double dyMeanBotElecGBL_hasL6_2016[] = {4.55212, -32.3026, 82.8215, -102.604, 59.5914, -12.9573, };
    double dySigmBotElecGBL_hasL6_2016[] = {24.7646, -93.0104, 168.237, -156.506, 72.498, -13.1847, };
    double dyMeanBotPosiGBL_hasL6_2016[] = {4.72926, -25.8425, 48.4637, -44.0838, 19.2657, -3.23875, };
    double dySigmBotPosiGBL_hasL6_2016[] = {-0.101683, 43.5468, -112.209, 117.355, -55.5621, 9.82688, };
    double dyMeanTopElecGBL_hasL6_2016[] = {18.5206, -106.714, 220.818, -209.39, 92.8465, -15.6116, };
    double dySigmTopElecGBL_hasL6_2016[] = {4.18864, 28.1707, -101.256, 129.418, -72.7543, 15.129, };
    double dyMeanTopPosiGBL_hasL6_2016[] = {-10.5208, 53.7065, -96.7171, 80.5609, -31.1419, 4.47035, };
    double dySigmTopPosiGBL_hasL6_2016[] = {7.10017, -2.14813, -8.23257, 7.40281, -0.97583, -0.419878, };
    double dxMeanBotElecGBL_hasL6_2016[] = {-26.9005, 109.663, -204.201, 185.67, -82.6174, 14.4039, };
    double dxSigmBotElecGBL_hasL6_2016[] = {17.8116, -60.9044, 102.479, -87.6493, 37.1673, -6.18107, };
    double dxMeanBotPosiGBL_hasL6_2016[] = {-3.91156, 41.5949, -76.9402, 58.3309, -19.6054, 2.44684, };
    double dxSigmBotPosiGBL_hasL6_2016[] = {10.9692, -27.6715, 43.3315, -38.0932, 17.3595, -3.17578, };
    double dxMeanTopElecGBL_hasL6_2016[] = {-22.2302, 81.627, -141.594, 123.058, -52.1465, 8.60953, };
    double dxSigmTopElecGBL_hasL6_2016[] = {-2.89464, 53.9567, -136.678, 147.799, -72.7372, 13.3487, };
    double dxMeanTopPosiGBL_hasL6_2016[] = {13.8746, -53.8677, 112.881, -114.995, 56.2456, -10.4624, };
    double dxSigmTopPosiGBL_hasL6_2016[] = {18.4828, -61.0319, 98.3593, -80.9737, 33.1667, -5.37544, };
    
    int dxMeanlength = 6;
    int dyMeanlength = 6;
    int dxSigmlength = 6;
    int dySigmlength = 6;

    double n_sigma_x(double p, double cx, double tx, boolean hasL6, int charge, boolean isTopTrack){
      double[] dxSigm = null;
      double[] dxMean = null;
      
      if (charge>0) {
        if (isTopTrack) {
          dxMean = !hasL6 ? dxMeanTopPosiGBL_noL6_2016 : dxMeanTopPosiGBL_hasL6_2016;
          dxSigm = !hasL6 ? dxSigmTopPosiGBL_noL6_2016 : dxSigmTopPosiGBL_hasL6_2016;
        }
        else {
          dxMean = !hasL6 ? dxMeanBotPosiGBL_noL6_2016 : dxMeanBotPosiGBL_hasL6_2016;
          dxSigm = !hasL6 ? dxSigmBotPosiGBL_noL6_2016 : dxSigmBotPosiGBL_hasL6_2016;
        }
      }
      else if (charge<0) {
        if (isTopTrack) {
          dxMean = !hasL6 ? dxMeanTopElecGBL_noL6_2016 : dxMeanTopElecGBL_hasL6_2016;
          dxSigm = !hasL6 ? dxSigmTopElecGBL_noL6_2016 : dxSigmTopElecGBL_hasL6_2016;
        }
        else {
          dxMean = !hasL6 ? dxMeanBotElecGBL_noL6_2016 : dxMeanBotElecGBL_hasL6_2016;
          dxSigm = !hasL6 ? dxSigmBotElecGBL_noL6_2016 : dxSigmBotElecGBL_hasL6_2016;
        }
      }

      if (p > 1.7 && hasL6) 
        p=1.7;

      if (p > 1.0 && !hasL6) 
        p=1.0;
                
      
      
      // calculate measured mean and sigma of deltaX for this energy:
      double aDxMean=0,aDxSigm=0;
      for (int ii=dxMeanlength-1; ii>=0; ii--) aDxMean = dxMean[ii] + p*aDxMean;
      for (int ii=dxSigmlength-1; ii>=0; ii--) aDxSigm = dxSigm[ii] + p*aDxSigm;
       
      // calculate nSigma between track and cluster:
      return  (cx - tx - aDxMean) / aDxSigm;
    }

    double ecal_y_min(double x){



      // if(x>-92 && x<-82)
      //return 30+(x - -92);
      //if(x>14 && x<30)
        //return 27.2-(x - 14);
      //else
       if(x>=-92 && x<=30){
       return 42.9;
     }
     else return 29.4;
    }

    boolean do_y_corr = true;

    double n_sigma_y(double p, double cy, double ty,double tx, boolean hasL6, int charge, boolean isTopTrack){
      double[] dySigm = null;
      double[] dyMean = null;

      if (charge>0) {
        if (isTopTrack) {
          dyMean = !hasL6 ? dyMeanTopPosiGBL_noL6_2016 : dyMeanTopPosiGBL_hasL6_2016;
          dySigm = !hasL6 ? dySigmTopPosiGBL_noL6_2016 : dySigmTopPosiGBL_hasL6_2016;
        }
        else {
          dyMean = !hasL6 ? dyMeanBotPosiGBL_noL6_2016 : dyMeanBotPosiGBL_hasL6_2016;
          dySigm = !hasL6 ? dySigmBotPosiGBL_noL6_2016 : dySigmBotPosiGBL_hasL6_2016;
        }
      }
      else if (charge<0) {
        if (isTopTrack) {
          dyMean = !hasL6 ? dyMeanTopElecGBL_noL6_2016 : dyMeanTopElecGBL_hasL6_2016;
          dySigm = !hasL6 ? dySigmTopElecGBL_noL6_2016 : dySigmTopElecGBL_hasL6_2016;
        }
        else {
          dyMean = !hasL6 ? dyMeanBotElecGBL_noL6_2016 : dyMeanBotElecGBL_hasL6_2016;
          dySigm = !hasL6 ? dySigmBotElecGBL_noL6_2016 : dySigmBotElecGBL_hasL6_2016;
        }
      }

      if (p > 1.7 && hasL6) 
        p=1.7;

      if (p > 1.0 && !hasL6) 
        p=1.0;

      //snap to the nearest point where a cluster could be reconstructed at
      if(do_y_corr){
        if(Math.abs(ty)<ecal_y_min(tx)){
          ty = ty/Math.abs(ty)*ecal_y_min(tx);
        }
        else if(ty>85)
          ty = 85;
        else if(ty<-85)
          ty = -85;
      }

      
       // calculate measured mean and sigma of deltaX for this energy:
      double aDyMean=0,aDySigm=0;
      for (int ii=dyMeanlength-1; ii>=0; ii--) aDyMean = dyMean[ii] + p*aDyMean;
      for (int ii=dySigmlength-1; ii>=0; ii--) aDySigm = dySigm[ii] + p*aDySigm;

       // calculate nSigma between track and cluster:
      return  (cy - ty - aDyMean) / aDySigm;
    }
}
