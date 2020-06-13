#include "stdio.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TTree.h"

#include "style.h"
#include "AnaFunctions.h"
#include "AnaUtils.h"
#include "AnaIO.h"

#include <algorithm>    // std::sort
#include <vector>       // std::vector

using namespace AnaFunctions;
using namespace AnaUtils;
//using namespace AnaIO;


bool manual_beamPos_mc(const double beam_startX, const double beam_startY, const double beam_startZ, const double beam_dirX, const double beam_dirY,   const double beam_dirZ, const double true_dirX,   const double true_dirY, const double true_dirZ,   const double true_startX, const double true_startY, const double true_startZ) 
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

bool getBeamPosPass()
{
  //https://github.com/calcuttj/PionStudies/blob/master/rDataFrame/eventSelection.C
  /*
.Define("passBeamCut", manual_beamPos_mc, 
            {"reco_beam_startX", "reco_beam_startY", "reco_beam_startZ",
             "reco_beam_trackDirX", "reco_beam_trackDirY", "reco_beam_trackDirZ",
             "true_beam_startDirX", "true_beam_startDirY", "true_beam_startDirZ",
             "true_beam_startX", "true_beam_startY", "true_beam_startZ"})
   */

  return manual_beamPos_mc(AnaIO::reco_beam_startX, AnaIO::reco_beam_startY, AnaIO::reco_beam_startZ,
                           AnaIO::reco_beam_trackDirX, AnaIO::reco_beam_trackDirY, AnaIO::reco_beam_trackDirZ,
                           AnaIO::true_beam_startDirX, AnaIO::true_beam_startDirY, AnaIO::true_beam_startDirZ,
                           AnaIO::true_beam_startX, AnaIO::true_beam_startY, AnaIO::true_beam_startZ);
}

bool cutBeamdEdx(const double varSignal)
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

vector<TLorentzVector> getFSTruth(const bool kPiZero, int & protonIdx, int & piplusIdx, bool & tmpksig)
{
  //if kPiZero false, return 1piplus (1st and 2nd protons)
  //if kPiZero true, return (1st PiZero) (1st and 2nd protons)
  //only use kPiZero in the end at filling

  const vector<int> * pdg = AnaIO::true_beam_daughter_PDG;
  const vector<double> * px = AnaIO::true_beam_daughter_startPx;
  const vector<double> * py = AnaIO::true_beam_daughter_startPy;
  const vector<double> * pz = AnaIO::true_beam_daughter_startPz;
  TH1I * htype = AnaIO::hselectedDaughterType;
  //===

  const int np = pdg->size();
  AnaIO::nproton = 0;
  AnaIO::nneutron = 0;
  int npiplus = 0;
  AnaIO::nPiZero = 0;
  AnaIO::ngamma = 0;
  int npartialBkg = 0;
  TLorentzVector pPiplus, pPiZero, pProton, p2Proton;

  protonIdx = -999;
  int protonIndices[np];
  double protonmom[np];
  double PiZeromom[np];
  double gammamom[np];
  vector<TVector3> bufferProtonmom;
  vector<TVector3> bufferPiZeromom;
  vector<int> bufferType;

  for(int ii=0; ii<np; ii++){
    const int itype = GetParticleType((*pdg)[ii]);
    const TVector3 tmpp( (*px)[ii], (*py)[ii], (*pz)[ii] );

    if(itype==gkPiPlus){
      npiplus++;
      pPiplus.SetVectM(tmpp, PionMass());

      piplusIdx = ii;
    }
    else if(itype==gkProton){
      //ii is the location in the original vector
      protonIndices[AnaIO::nproton]=ii;

      protonmom[AnaIO::nproton] = tmpp.Mag();
      bufferProtonmom.push_back(tmpp);
      AnaIO::nproton++;
    }
    else if(itype==gkPiZero){
      PiZeromom[AnaIO::nPiZero] = tmpp.Mag();
      bufferPiZeromom.push_back(tmpp);
      AnaIO::nPiZero++;
    }
    else if(itype==gkNeutron){
      AnaIO::nneutron++;
    }
    else if(itype==gkGamma){
      gammamom[AnaIO::ngamma] = tmpp.Mag();
      AnaIO::ngamma++;
    }
    else if(itype==gkPiMinus||itype==gkKaon){
      npartialBkg++;
    }

    bufferType.push_back(itype);
  }

  //proton=======================
  int leadingProtonID = 0, subldProtonID = -999;
  if(AnaIO::nproton>1){
    int protonsortid[AnaIO::nproton];
    TMath::Sort(AnaIO::nproton, protonmom, protonsortid);
    leadingProtonID = protonsortid[0];
    subldProtonID = protonsortid[1];
    //printf("test %d %f %d %f\n", leadingProtonID, protonmom[leadingProtonID], subldProtonID, protonmom[subldProtonID]);
  }
  if(AnaIO::nproton>0){
    pProton.SetVectM(bufferProtonmom[leadingProtonID], ProtonMass());
    protonIdx = protonIndices[leadingProtonID];
  }
  if(AnaIO::nproton>1){
    p2Proton.SetVectM(bufferProtonmom[subldProtonID], ProtonMass());
  }
  //PiZero==============================
  int leadingPiZeroID = 0;
  if(AnaIO::nPiZero>1){
    int PiZerosortid[AnaIO::nPiZero];
    TMath::Sort(AnaIO::nPiZero, PiZeromom, PiZerosortid);
    leadingPiZeroID = PiZerosortid[0];
  }
  if(AnaIO::nPiZero>0){
    pPiZero.SetVectM(bufferPiZeromom[leadingPiZeroID], PiZeroMass());
  }
  //gamma====================================
  int leadinggammaID = 0;
  if(AnaIO::ngamma>1){
    int gammasortid[AnaIO::ngamma];
    TMath::Sort(AnaIO::ngamma, gammamom, gammasortid);
    leadinggammaID = gammasortid[0];
  }
  AnaIO::maxgammaEnergy = gammamom[leadinggammaID];
  //==========================================

  //fill vec regardless of tmpksig
  vector<TLorentzVector> vec;
  vec.push_back(kPiZero?pPiZero:pPiplus);
  vec.push_back(pProton);
  vec.push_back(p2Proton);
    
  if(npiplus!=1){
    piplusIdx = -999;
  }

  //Now need kPiZero
  tmpksig = false;
  if(AnaIO::nproton>=1 && npartialBkg ==0){
    if(kPiZero){
      if(AnaIO::nPiZero>0 && npiplus==0){
        tmpksig = true;
      }
    }
    else{//PiPlus
      if(npiplus==1 && AnaIO::nPiZero==0){
        tmpksig = true;
      }
    }
  }

  if(tmpksig && htype){
    const int nbt = bufferType.size();
    for(int ii=0; ii<nbt; ii++){
      htype->Fill(bufferType[ii]);
    }
  }

  return vec;
}


int getTruthFromRec(const int recidx, int & pdg, double & momentum)
{
  const int truthID = (*AnaIO::reco_daughter_PFP_true_byHits_ID)[recidx];
  int trueidx = -999;
  pdg = -999;
  momentum = -999;
  int counter = 0;
  for(unsigned int ii = 0; ii<AnaIO::true_beam_daughter_ID->size(); ii++){
    if(  (*AnaIO::true_beam_daughter_ID)[ii] == truthID ){
      if(counter){
        printf("getTruthFromRec truthID found again %d %d\n", recidx, truthID); exit(1);
      }
      trueidx = ii;
      counter++;
    }
  }

  if(trueidx>=0){
    pdg = (*AnaIO::true_beam_daughter_PDG)[trueidx];

    const TVector3 vec((*AnaIO::true_beam_daughter_startPx)[trueidx], (*AnaIO::true_beam_daughter_startPy)[trueidx], (*AnaIO::true_beam_daughter_startPz)[trueidx]);
    momentum = vec.Mag();
  }

  return trueidx;
}

int getNTrack(int & nproton, int & ngamma, int & nmichel, const bool kprint=false, const bool ksig=false)
{
  /*
root [11] beamana->Scan("reco_daughter_PFP_ID:reco_daughter_PFP_true_byHits_ID:reco_daughter_allTrack_ID:reco_daughter_allShower_ID:true_beam_daughter_ID:true_beam_daughter_reco_byHits_PFP_ID:true_beam_daughter_reco_byHits_allTrack_ID:true_beam_daughter_reco_byHits_allShower_ID")
***********************************************************************************************************************
*    Row   * Instance * reco_daug * reco_daug * reco_daug * reco_daug * true_beam * true_beam * true_beam * true_beam *
***********************************************************************************************************************
*       13 *        0 *        86 *     40841 *        86 *        86 *     40800 *           *           *           *
*       13 *        1 *       226 *     40855 *       225 *       225 *     40841 *           *           *           *

So
(reco_daughter_PFP_ID = reco_daughter_allTrack_ID = reco_daughter_allShower_ID)
                             !=
 (reco_daughter_PFP_true_byHits_ID = true_beam_daughter_ID)

All 
true_beam_daughter_reco_byHits_PFP_ID:
true_beam_daughter_reco_byHits_allTrack_ID:
true_beam_daughter_reco_byHits_allShower_ID
not set
   */

  const int recsize = AnaIO::reco_daughter_PFP_ID->size();
  const int mbr     = AnaIO::reco_daughter_allTrack_momByRange_proton->size();
  if(recsize!=mbr){
    printf("recsize != mbr %d %d\n", recsize, mbr); exit(1); 
  }

  /*
est rec 1/2 trueID 15 pdg 22 truemomentum 0.002512 recp -999.000000 resolution -397725.603302 nhits 7 trackScore 0.016739 emScore 0.983240 michelScore 0.289565 sum 1.289545
   */
  nproton = 0;
  ngamma = 0;
  nmichel = 0;
  int ntracks = 0;
 
  for(int ii=0; ii<recsize; ii++){

    const vector<double>  dedxarray = (*AnaIO::reco_daughter_allTrack_calibrated_dEdX_SCE)[ii];
    const int dedxsize = dedxarray.size();
     if(dedxsize<3){
      //do not count it as track
      continue;
    }

     const int nhits          = (*AnaIO::reco_daughter_PFP_nHits)[ii];
    /*//they are indeed different: 236 80
    if(nhits != dedxsize){
      printf("nhits!=dedx size %d %d\n", nhits, dedxsize); exit(1);
    }
    */

    double recp        = -999;
    const double startE2 = dedxarray[2];
    const double startE3 = dedxarray[3];

    const double chi2 = (*AnaIO::reco_daughter_allTrack_Chi2_proton)[ii];
    const double ndof = (*AnaIO::reco_daughter_allTrack_Chi2_ndof)[ii];
    const double Chi2NDF = chi2/(ndof+1E-10);

    //all signal protons have nhit below 260
    //all signal protons have startE3 > 9
    //with Chi2NDF cut  44/191 = 23% purity; without 46/197 = 23% purity, slightly higher efficiency
    //test if(startE2>10 && nhits<260 && startE3>9 && Chi2NDF<50){
    if(startE2>10 && nhits<260 && startE3>9){
      recp = (*AnaIO::reco_daughter_allTrack_momByRange_proton)[ii];
      nproton++;
    }
  
    const double trackScore  = (*AnaIO::reco_daughter_PFP_trackScore)[ii];
    if(trackScore>0.5){
      ntracks++;
    }
    const double emScore     = (*AnaIO::reco_daughter_PFP_emScore)[ii];
    if(emScore>0.5){
      ngamma++;
    }
    const double michelScore = (*AnaIO::reco_daughter_PFP_michelScore)[ii];
    if(michelScore>0.5){
      nmichel++;
    }

    int pdg = -999;
    double truemomentum = -999;
    const int trueID = getTruthFromRec(ii, pdg, truemomentum);

    if(kprint){
      printf("test ksig %d rec %d/%d trueID %d pdg %d truemomentum %f recp %f resolution %f startE2 %f startE3 %f nhits %d trackScore %f emScore %f michelScore %f sum %f chi2 %f ndof %f chi2/ndof %f\n", ksig, ii, recsize, trueID, pdg, truemomentum, recp, recp/truemomentum-1, startE2, startE3,  nhits, trackScore, emScore, michelScore, trackScore+emScore+michelScore, chi2, ndof, Chi2NDF);

      if(!ksig){
        if(pdg==2212 && trueID>=0){
          AnaIO::hProtonnHits->Fill(nhits);
          AnaIO::hProtontrackScore->Fill(trackScore);
          AnaIO::hProtonemScore->Fill(emScore);
          AnaIO::hProtonmichelScore->Fill(michelScore);
          AnaIO::hProtonChi2NDF->Fill(Chi2NDF);
          AnaIO::hProtonMomentumRes->Fill(truemomentum, recp/truemomentum-1);
        }

        if(pdg==211 && trueID>=0){
          AnaIO::hPiPlusnHits->Fill(nhits);
          AnaIO::hPiPlustrackScore->Fill(trackScore);
          AnaIO::hPiPlusemScore->Fill(emScore);
          AnaIO::hPiPlusmichelScore->Fill(michelScore);
          AnaIO::hPiPlusChi2NDF->Fill(Chi2NDF);
          //AnaIO::hPiPlusMomentumRes->Fill(truemomentum, recp/truemomentum-1);
        }
      }
    }
    
  }

  if(kprint){
    printf("\n\n");
  }

  return ntracks;
}

double getRecFromTruth(const int truthID, const vector<double> * mombyrange)
{
  vector<double> lastE, startE;

  static int ievent = 0;
  ievent++;

  double rpm = -999;
  int counter = 0;
  for(unsigned int ii = 0; ii < AnaIO::reco_daughter_PFP_true_byHits_ID->size(); ii++){
    if((*AnaIO::reco_daughter_PFP_true_byHits_ID)[ii] == truthID) {
      if((*AnaIO::reco_daughter_allTrack_ID)[ii] != -1) {

        if(counter){
          printf("rpm already set!! %f %f %d %d\n", rpm, (*mombyrange)[ii], counter, ievent); //exit(1);
        }

        rpm = (*mombyrange)[ii];

        AnaIO::chi2 = (*AnaIO::reco_daughter_allTrack_Chi2_proton)[ii];
        AnaIO::ndof = (*AnaIO::reco_daughter_allTrack_Chi2_ndof)[ii];

        if(!AnaIO::reco_daughter_allTrack_calibrated_dEdX_SCE){
          printf("AnaIO::reco_daughter_allTrack_calibrated_dEdX_SCE null!\n"); exit(1);
        }

        GetdEdx( (*AnaIO::reco_daughter_allTrack_calibrated_dEdX_SCE)[ii], startE, lastE, 0);
        //printf("==================================\n");

        const int iter0 = 2;
        AnaIO::lastTruncatedMeanE = GetTruncatedMean(lastE, iter0, iter0-1+10, 0.0, 1.0);
        AnaIO::startTruncatedMeanE10 = GetTruncatedMean(startE, iter0, iter0-1+10, 0.4, 0.95);
        AnaIO::startTruncatedMeanE20 = GetTruncatedMean(startE, iter0, iter0-1+20, 0.4, 0.95);
        AnaIO::startTruncatedMeanE30 = GetTruncatedMean(startE, iter0, iter0-1+30, 0.4, 0.95);
        AnaIO::startTruncatedMeanE40 = GetTruncatedMean(startE, iter0, iter0-1+40, 0.4, 0.95);
        AnaIO::startTruncatedMeanE50 = GetTruncatedMean(startE, iter0, iter0-1+50, 0.4, 0.95);

        counter++;
      }
    }
  }

  AnaIO::ndEdxCls = lastE.size();
  if(AnaIO::ndEdxCls>=6){
    AnaIO::lastE0 = lastE[0];
    AnaIO::lastE1 = lastE[1];
    AnaIO::lastE2 = lastE[2];
    AnaIO::lastE3 = lastE[3];
    AnaIO::lastE4 = lastE[4];
    AnaIO::lastE5 = lastE[5];
    
    AnaIO::startE0 = startE[0];
    AnaIO::startE1 = startE[1];
    AnaIO::startE2 = startE[2];
    AnaIO::startE3 = startE[3];
    AnaIO::startE4 = startE[4];
    AnaIO::startE5 = startE[5];
  }

  return rpm;
}

void anaRec(TList *lout, const TString tag, const int nEntryToStop = -999)
{
  const bool kPiZero = tag.Contains("MPiZero");
  const bool kTrackingProton = !tag.Contains("TrackingPiPlus");

  cout<<"\n\n                       Running anaRec kPiZero "<<kPiZero<<" TrackingProton "<<kTrackingProton<<endl<<endl;

  TTree * tree = AnaIO::GetInputTree("pionana/beamana");
  TTree * tout = AnaIO::GetOutputTree(lout, tag);  
  AnaIO::IniRecHist(lout, tag);

  //==============================================================================================================
  int ientry = 0;
  while(tree->GetEntry(ientry)){
    if(ientry%100000==0){
      printf("myEntries %d\n", ientry);
    }
    
    if(nEntryToStop>0){
      if(ientry>=nEntryToStop){
        printf("\n\n\n************************  Breaking after %d entries ***********************************************\n\n", nEntryToStop);
        break;
      }
    }

    //do it before the loop continues for any reason
    ientry++;
   
    //===========================================================
    //calculate before any cuts! Only filled after ALL cuts!
    const int TruthBeamType = GetParticleType(AnaIO::true_beam_PDG);
    int  protonIdx = -999, piplusIdx = -999;
    bool tmpkSig = false;
    vector<TLorentzVector> vecPiP = getFSTruth(kPiZero, protonIdx, piplusIdx, tmpkSig);

    AnaIO::finPimomentum = vecPiP[0].P();
    AnaIO::finProtonmomentum = vecPiP[1].P();
    AnaIO::fin2Pmom = vecPiP[2].P();

    //no phase space cut
    AnaIO::kSignal = (AnaIO::true_beam_PDG==211) &&  tmpkSig; 
    //with phase space cut
    AnaIO::kSignal = (AnaIO::true_beam_PDG==211) &&  tmpkSig && (AnaIO::finProtonmomentum>0.45 && AnaIO::fin2Pmom<0.45);
    if(!kPiZero){
      AnaIO::kSignal = AnaIO::kSignal && (AnaIO::finPimomentum>0.15);
    }

    //const bool varTrueBeam = (AnaIO::true_beam_PDG==211);
    const bool varSignal = AnaIO::kSignal;

    //for to see within signal
    if(0){
      if(!AnaIO::kSignal){
        continue;
      }
    }
    else{
      //===========================================================
      //0. true beam particle //"In data the beam Instrumentation is able to filter for these events but is not inside the MC" so read data file has this
      if(AnaIO::true_beam_PDG != 211 && AnaIO::true_beam_PDG != -13){
        continue;
      }
      
      //without 20 sig in 36 events; with 18/32, both 56% purity -> same purity, without this has higher efficiency
      if(0){
        //x. primary beam type 
        AnaIO::hRecoBeamType->Fill(AnaIO::reco_beam_type);
        //no effect, shadowed by TMeanStart cut
        //    if(AnaIO::reco_beam_type!=13){//13: Pandora "track like"
        //      continue;
        //    }
        
        //1. beam position MC cut, need MC truth, how is it possible in analysis?
        const bool kBeamPosPass = getBeamPosPass();
        AnaIO::hBeamPosPass->Fill(kBeamPosPass);
        if(!kBeamPosPass){
          continue;
        }
        //-> now signal purity 138/3537 = 3.9%, 2283 pi+ bea, 801 e+ beam
        
        //2. APA3 
        AnaIO::hBeamEndZ->Fill(AnaIO::reco_beam_endZ);
        AnaIO::hBeamEndZPass->Fill(!(AnaIO::reco_beam_endZ>=226));
        if(AnaIO::reco_beam_endZ>=226){
          AnaIO::hBeamEndZPass->Fill(false);
          continue;
        }
        //-> now signal purity 135/3143 = 4.3%, 2102 pi+ beam, 801 e+ beam
      }
    }

    //3. n track daughter
    int cutnproton = 0;
    int cutngamma = 0;
    int cutnmichel = 0;
    AnaIO::nTrack = getNTrack(cutnproton, cutngamma, cutnmichel);

    AnaIO::hBeamNTrack->Fill(AnaIO::nTrack, cutnproton);
    AnaIO::hSignalVsBeamNTrack->Fill(AnaIO::nTrack, varSignal);

    /*
      PiPlus: ngamma = 0, nproton=1 -> 75/218=34% in signal sample; 63/467 = 13% purity, 63/218 = 29% efficiency, p*e = 63*63/218/467 = 3.9%
      PiPlus: ngamma = 0, nproton=1,  noly nhit<260 and startE3>9 and Chi2NDF<50 is proton, no pion-analysis pre-cuts! -> 61/380 = 16% 
      PiPlus: ngamma = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 and Chi2NDF<50 is proton, no pion-analysis pre-cuts! -> 44/191 = 23% purity, e*p = 44.*44./191./218. = 4.6%
  (*) PiPlus: ngamma = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9                is proton, no pion-analysis pre-cuts! -> 46/197 = 23% purity, slightly higher efficiency

      PiZero: ngamma = 2, nmichel=0, nproton=1 -> 32/260=12% in signal sample; 27/68 = 40% purity, 27/260 = 10% efficiency, p*e = 27*27/68/260 = 4.1%
      PiZero: ngamma = 2, nmichel=0, nproton=1, ntrack=1 -> ; 19/40 = 48% purity
      PiZero: ngamma = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 is proton -> 18/37=49% purity
      PiZero: ngamma = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 is proton -> 18/32=56% purity, p*e = 18.*18./32./260. = 3.9%
      PiZero: ngamma = 2, nmichel=0, nproton=1,           noly nhit<260 and startE3>9 is proton -> 21/47=45% purity -> so it is important to have ntrack=1

      PiZero: ngamma = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 and Chi2NDF<50 is proton, no pion-analysis pre-cuts! -> 18/34 = 53% purity
  (*) PiZero: ngamma = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 is proton, no pion-analysis pre-cuts! -> 20/36=56% purity, eff*purity = 20.*20./36./260. = 4.3%
     */

    if(!kPiZero){
      if(cutngamma>0){
        continue;
      }
      //not found in signal, need Michel to tell 2-proton events from p-pi
      /*
      if(cutnmichel!=1){
        continue;
      }
      */
      if(AnaIO::nTrack!=2){
        continue;
      }
    }
    else{
      if(cutngamma!=2){
        continue;
      }

      if(cutnmichel>0){
        continue;
      }

      //with this 56% purity, without 45%
      if(AnaIO::nTrack!=1){
        continue;
      }
    }

    if(cutnproton!=1){
      continue;
    }

    //just for printing
    getNTrack(cutnproton, cutngamma, cutnmichel, true, varSignal);

    //nTrack has Michel object
    /*
    if(kPiZero){
      //the purity profile looks weird for pizeor signal
      //4/260 events signal efficiency
      if(AnaIO::nTrack!=1){
        continue;
      }
    }
    else{
      // 50/218 events signal efficiency
      if(AnaIO::nTrack!=2){
        continue;
      }
    }
    */
    //-> now pi+p signal purity 78/940=8.3%, 611 pi+ beam, 270 e+ beam

    /*   
    //x. Beam dEdx cut shadowed by beam filtering
    AnaIO::hBeamLen->Fill(AnaIO::reco_beam_len);
    AnaIO::hSignalVsLen->Fill(AnaIO::reco_beam_len, varSignal);
    if(!cutBeamdEdx(varSignal)){
      continue;
    }
    */
    //-> now signal purity 167/5274 = 3.2%, 2627 pi+ beam, 1947 e+

   

    //============== Benchmark after ALL cuts !!! =========================
    //benchmark

    AnaIO::hTruthBeamType->Fill(TruthBeamType);
    AnaIO::hTruthSignal->Fill(AnaIO::kSignal);
   
    //======================= NO cuts below ========================

    //--- to test
    //Fill kSignal vs variable; the reason for kSignal but not kBeam is the signal is interacting pion, not all pions. Directly choose matrix to optimize for it
    AnaIO::hSigAfterVsLen->Fill(AnaIO::reco_beam_len, varSignal);

    //============= done loop
    
    tout->Fill();
  }

  cout<<"All entries "<<ientry<<endl;

  GetProfileX(lout);
  DrawHist(lout, "output", tag, true);
}

void anaTruth(TList *lout, const TString tag, const int nEntryToStop = -999)
{
  const bool kPiZero = tag.Contains("MPiZero");
  const bool kTrackingProton = !tag.Contains("TrackingPiPlus");

  cout<<"\n\n                       Running anaTruth kPiZero "<<kPiZero<<" TrackingProton "<<kTrackingProton<<endl<<endl;

  TTree * tree = AnaIO::GetInputTree("pionana/beamana");
  TTree * tout = AnaIO::GetOutputTree(lout, tag);
  AnaIO::IniTruthHist(lout, tag);

  //==============================================================================================================
  int ientry = 0;
  while(tree->GetEntry(ientry)){
    if(ientry%100000==0){
      printf("myEntries %d\n", ientry);
    }
    
    if(nEntryToStop>0){
      if(ientry>=nEntryToStop){
        printf("\n\n\n************************  Breaking after %d entries ***********************************************\n\n", nEntryToStop);
        break;
      }
    }

    //do it before the loop continues for any reason
    ientry++;
   
    //===========================================================

    const int beamtype = GetParticleType(AnaIO::true_beam_PDG);
    AnaIO::hbeamType->Fill(beamtype);

    if(beamtype!=gkPiPlus){
      continue;
    }

    const int nd = AnaIO::true_beam_daughter_PDG->size();
    AnaIO::hndaughter->Fill(nd);
    if(nd==0){
      continue;
    }

    const TLorentzVector iniPiPlus(AnaIO::true_beam_endPx, AnaIO::true_beam_endPy, AnaIO::true_beam_endPz, PionMass());

    AnaIO::iniPimomentum = iniPiPlus.P();
    AnaIO::iniPitheta = iniPiPlus.Theta()*TMath::RadToDeg();

    int  protonIdx = -999, piplusIdx = -999;
    AnaIO::kSignal = false;
    vector<TLorentzVector> vecPiP = getFSTruth(kPiZero, protonIdx, piplusIdx, AnaIO::kSignal);
    AnaIO::hnexcl->Fill(AnaIO::kSignal);

    AnaIO::finPimomentum = vecPiP[0].P();
    AnaIO::finProtonmomentum = vecPiP[1].P();
    AnaIO::fin2Pmom = vecPiP[2].P();

     //signal truth
    if(AnaIO::kSignal){
      AnaIO::hnproton->Fill(AnaIO::nproton);
      AnaIO::hnneutron->Fill(AnaIO::nneutron);
      AnaIO::hnPiZero->Fill(AnaIO::nPiZero);
      
      //re-calculate final pi p theta w.r.t. iniPi
      const int targetA = 40;
      const int targetZ = 18;
      getCommonTKI(targetA, targetZ, &iniPiPlus, &(vecPiP[0]), &(vecPiP[1]), AnaIO::dalphat, AnaIO::dphit, AnaIO::dpt, AnaIO::pn, AnaIO::finPitheta, AnaIO::finProtontheta);
      
      AnaIO::hmomIniPi->Fill(AnaIO::iniPimomentum);
      AnaIO::hmomFinPi->Fill(AnaIO::finPimomentum);
      AnaIO::hmomFinProton->Fill(AnaIO::finProtonmomentum);
      AnaIO::hdalphat->Fill(AnaIO::dalphat);
      AnaIO::hdphit->Fill(AnaIO::dphit);
      AnaIO::hdpt->Fill(AnaIO::dpt);
      AnaIO::hpn->Fill(AnaIO::pn);
    }

    vector<double>* mombyrange = 0x0;
    bool kDoTracking = false;
    double * recmomentum = 0x0;

    if(kTrackingProton){
      mombyrange = AnaIO::reco_daughter_allTrack_momByRange_proton;
      kDoTracking = (protonIdx>=0);
      recmomentum = &AnaIO::recProtonmomentum;
    }
    else{
      mombyrange = AnaIO::reco_daughter_allTrack_momByRange_muon;
      kDoTracking = (piplusIdx>=0);
      recmomentum = &AnaIO::recPiPlusmomentum;
    }

    AnaIO::hksignp->Fill(kDoTracking, AnaIO::kSignal);
 
    //has true proton
    if(kDoTracking){
      const int truthID = (*AnaIO::true_beam_daughter_ID)[protonIdx];
      (*recmomentum) = getRecFromTruth(truthID, mombyrange);
      tout->Fill();
    }
  }

  cout<<"All entries "<<ientry<<endl;

  DrawHist(lout, "output", tag);
}

int main(int argc, char * argv[])
{
  if(argc!=4){
    cout<<"argc!=4 !"<<argc<<endl;
    return 0;
  }

  gSystem->Load("libTree");

  style::SetGlobalStyle();

  const bool kPiZero = atoi(argv[1]);
  const bool kProton = atoi(argv[2]);
  if(kPiZero&&!kProton){
    return 0;
  }
  const bool kTruth = atoi(argv[3]);
  if(!kTruth&&!kProton){
    return 0;
  }

  TString tag = (kPiZero?"MPiZero":"1PiPlus");
  tag+=(kProton?"_TrackingProton":"_TrackingPiPlus");
  tag+=(kTruth?"_anaTruth":"_anaRec");

  TList * lout = new TList;

  TFile * fin = new TFile("input/protoDUNE_reco_flattree.root");
  if(!fin->IsOpen()){
    cout<<"fin not open!"<<endl;
    exit(1);
  }

  if(kTruth){
    anaTruth(lout, tag);
  }
  else{
    anaRec(lout, tag);
  }

  TFile * fout = new TFile(Form("output/outanaData_%s.root", tag.Data()),"recreate");

  lout->Write();
  fout->Save();
  fout->Close();

  return 0;
}
