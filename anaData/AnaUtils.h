#ifndef _ANAUTILS_H_
#define _ANAUTILS_H_

#include "AnaIO.h"

using namespace std;

namespace AnaUtils
{

enum{
  //1-2
  gkProton = 1,
  gkNeutron,

  //3-5
  gkPiPlus,
  gkPiZero,
  gkPiMinus,
  
  //6-9
  gkElectron,
  gkPositron,
  gkMuPlus,
  gkMuMinus,

  //10
  gkKaon,

  //11
  gkGamma,

  //12
  gkNeutrino,

  //13
  gkHyperon,

  //14
  gkNucleus,

  //15
  gkNotCaredType
};


int GetParticleType(const int pdg)
{
  int type = -999;
  if(pdg==2212){
    type = gkProton;
  }
  else if(pdg==211){
    type = gkPiPlus;
  }
  else if(pdg==-211){
    type = gkPiMinus;
  }
  else if(pdg==111){
    type = gkPiZero;
  }
  else if(pdg==-11){
    type = gkPositron;
  }
  else if(pdg==11){
    type = gkElectron;
  }
  else if(pdg==-13){
    type = gkMuPlus;
  }
  else if(pdg==13){
    type = gkMuMinus;
  }
  else if(pdg==321||pdg==-321||pdg==310||pdg==130){
    type = gkKaon;
  }
  else if(pdg==2112){
    type = gkNeutron;
  }
  else if(pdg==22){
    type = gkGamma;
  }
  else if(pdg==14 || pdg==-14 || pdg==12 || pdg==-12){
    type = gkNeutrino;
  }
  else if(pdg==3122||pdg==3212||pdg==3222){
    type = gkHyperon;
  }
  else if(pdg>9999){
    type = gkNucleus;
  }
  else{
    cout<<"getParticleType unknown pdg "<<pdg<<endl; exit(1);
    //type = gkNotCaredType;
  }

  return type;
}



double GetTruncatedMean(vector<double> array, const int nsample)
{
  const double fracTrun = 0.6;
  const double nterm = nsample*fracTrun;
  //either nsample<0, or total nterm is too small

  std::sort(array.begin(), array.begin()+nsample);

  double sum =0.0;

  for(unsigned int ii=0; ii< nterm; ii++){
    sum += array[ii];
  }
  return sum / (nterm+1E-10);
}

double GetTruncatedMean(const vector<double> tmparr, const unsigned int nsample0, const unsigned int nsample1, const double lowerFrac, const double upperFrac)
{
  //for proton Bragg peak use 0.4-0.95. Seen by CDF of startE using signal proton samples in drawTracking

  if(nsample1>=tmparr.size()){
    return -999;
  }

  vector<double> array;
  for(unsigned int ii=nsample0; ii<=nsample1; ii++){
    array.push_back(tmparr[ii]);
  }

  std::sort(array.begin(), array.end());

  double sum =0.0;

  const int iter0 = array.size()*lowerFrac;
  const int iter1 = array.size()*upperFrac;

  for(int ii=iter0; ii< iter1; ii++){
    sum += array[ii];
  }
  return sum / ( (iter1-iter0)+1E-10);
}


void GetdEdx(const vector<double> arraydEdx, vector<double> &startE, vector<double> &endE, const unsigned int padding=0)
{
  const unsigned int ncls = arraydEdx.size();

  if(ncls<=padding){
    return;
  }

  //start from [2] because [0] and [1] in both start and last are weird
  for(unsigned int kk=padding; kk<ncls; kk++){
    startE.push_back(arraydEdx[kk]);

    const double endpe = arraydEdx[ncls-1-kk];
    endE.push_back(endpe);
  }
}


bool Manual_beamPos_mc(const double beam_startX, const double beam_startY, const double beam_startZ, const double beam_dirX, const double beam_dirY,   const double beam_dirZ, const double true_dirX,   const double true_dirY, const double true_dirZ,   const double true_startX, const double true_startY, const double true_startZ) 
{
  //https://github.com/calcuttj/PionStudies/blob/master/rDataFrame/eventSelection.h

  //For MC from Owen Goodwins studies
  const double xlow = -3.,  xhigh = 7.,  ylow = -8.,  yhigh = 7.;
  const double zlow = 27.5,  zhigh = 32.5,  coslow = 0.93;

  const double projectX = (true_startX + -1*true_startZ*(true_dirX/true_dirZ) );
  const double projectY = (true_startY + -1*true_startZ*(true_dirY/true_dirZ) );
  const double cos = true_dirX*beam_dirX + true_dirY*beam_dirY + true_dirZ*beam_dirZ;

  if ( (beam_startX - projectX) < xlow )
    return false;
  
  if ( (beam_startX - projectX) > xhigh )
    return false;

  if ( (beam_startY - projectY) < ylow )
    return false;

  if ( (beam_startY - projectY) > yhigh )
    return false;
  
  if (beam_startZ < zlow || zhigh < beam_startZ)
    return false;
  
  if ( cos < coslow)
    return false;

  return true;
}

bool GetBeamPosPass()
{
  //https://github.com/calcuttj/PionStudies/blob/master/rDataFrame/eventSelection.C
  /*
.Define("passBeamCut", manual_beamPos_mc, 
            {"reco_beam_startX", "reco_beam_startY", "reco_beam_startZ",
             "reco_beam_trackDirX", "reco_beam_trackDirY", "reco_beam_trackDirZ",
             "true_beam_startDirX", "true_beam_startDirY", "true_beam_startDirZ",
             "true_beam_startX", "true_beam_startY", "true_beam_startZ"})
   */

  return Manual_beamPos_mc(AnaIO::reco_beam_startX, AnaIO::reco_beam_startY, AnaIO::reco_beam_startZ,
                           AnaIO::reco_beam_trackDirX, AnaIO::reco_beam_trackDirY, AnaIO::reco_beam_trackDirZ,
                           AnaIO::true_beam_startDirX, AnaIO::true_beam_startDirY, AnaIO::true_beam_startDirZ,
                           AnaIO::true_beam_startX, AnaIO::true_beam_startY, AnaIO::true_beam_startZ);
}


bool CutBeamdEdx(const double varSignal)
{
  vector<double> startE, lastE;
  GetdEdx( *AnaIO::reco_beam_calibrated_dEdX, startE, lastE, 2);
  AnaIO::nBeamdEdxCls = startE.size();
  AnaIO::hnBeamdEdxCls->Fill(AnaIO::nBeamdEdxCls);
  if(AnaIO::nBeamdEdxCls<6){
    return false;
  }

  //no bragg peak
  AnaIO::beamStartE0 = startE[0];
  AnaIO::beamStartE1 = startE[1];
  AnaIO::beamStartE2 = startE[2];
  AnaIO::beamStartE3 = startE[3];
  AnaIO::beamStartE4 = startE[4];
  AnaIO::beamStartE5 = startE[5];
  
  AnaIO::beamTMeanStart = GetTruncatedMean(startE, AnaIO::nBeamdEdxCls-6);
  
  AnaIO::hSignalVsStartE0->Fill(AnaIO::beamStartE0, varSignal);
  AnaIO::hSignalVsStartE1->Fill(AnaIO::beamStartE1, varSignal);
  AnaIO::hSignalVsStartE2->Fill(AnaIO::beamStartE2, varSignal);
  AnaIO::hSignalVsStartE3->Fill(AnaIO::beamStartE3, varSignal);
  AnaIO::hSignalVsStartE4->Fill(AnaIO::beamStartE4, varSignal);
  AnaIO::hSignalVsStartE5->Fill(AnaIO::beamStartE5, varSignal);
  
  AnaIO::hSignalVsTMeanStart->Fill(AnaIO::beamTMeanStart, varSignal);
  
  //has Bragg Peak
  AnaIO::beamLastE0 = lastE[0];
  AnaIO::beamLastE1 = lastE[1];
  AnaIO::beamLastE2 = lastE[2];
  AnaIO::beamLastE3 = lastE[3];
  AnaIO::beamLastE4 = lastE[4];
  AnaIO::beamLastE5 = lastE[5];
  
  AnaIO::beamTMeanLast = GetTruncatedMean(lastE, 6);
  
  AnaIO::hSignalVsLastE0->Fill(AnaIO::beamLastE0, varSignal);
  AnaIO::hSignalVsLastE1->Fill(AnaIO::beamLastE1, varSignal);
  AnaIO::hSignalVsLastE2->Fill(AnaIO::beamLastE2, varSignal);
  AnaIO::hSignalVsLastE3->Fill(AnaIO::beamLastE3, varSignal);
  AnaIO::hSignalVsLastE4->Fill(AnaIO::beamLastE4, varSignal);
  AnaIO::hSignalVsLastE5->Fill(AnaIO::beamLastE5, varSignal);
  
  //both need to tune for different energy
  //it is so clean that no need to cut on last since there is no Bragg peak form proton any more
  if(AnaIO::beamTMeanStart>2.8){
    return false;
  }
  
  AnaIO::hSignalVsTMeanLast->Fill(AnaIO::beamTMeanLast, varSignal);
  
  AnaIO::hSigAfterVsStartE0->Fill(AnaIO::beamStartE0, varSignal);
  AnaIO::hSigAfterVsStartE1->Fill(AnaIO::beamStartE1, varSignal);
  AnaIO::hSigAfterVsStartE2->Fill(AnaIO::beamStartE2, varSignal);
  
  AnaIO::hSigAfterVsLastE0->Fill(AnaIO::beamLastE0, varSignal);
  AnaIO::hSigAfterVsLastE1->Fill(AnaIO::beamLastE1, varSignal);
  AnaIO::hSigAfterVsLastE2->Fill(AnaIO::beamLastE2, varSignal);
  
  AnaIO::hSigAfterVsStartE3->Fill(AnaIO::beamStartE3, varSignal);
  AnaIO::hSigAfterVsLastE3->Fill(AnaIO::beamLastE3, varSignal);
  AnaIO::hSigAfterVsStartE4->Fill(AnaIO::beamStartE4, varSignal);
  AnaIO::hSigAfterVsLastE4->Fill(AnaIO::beamLastE4, varSignal);
  AnaIO::hSigAfterVsStartE5->Fill(AnaIO::beamStartE5, varSignal);
  AnaIO::hSigAfterVsLastE5->Fill(AnaIO::beamLastE5, varSignal);
  
  return true;
}

//=== end
}


#endif

