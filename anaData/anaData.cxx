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

void setFullSignal(const bool kpi0)
{
  int  protonIdx = -999, piplusIdx = -999;
  vector<TLorentzVector> vecPiP = getFSTruth(kpi0, protonIdx, piplusIdx, AnaIO::kSignal);
  
  AnaIO::finPimomentum = vecPiP[0].P();
  AnaIO::finProtonmomentum = vecPiP[1].P();
  AnaIO::fin2Pmom = vecPiP[2].P();
  
  //no phase space cut
  AnaIO::kSignal = AnaIO::kSignal && (AnaIO::true_beam_PDG==211);
  //with phase space cut
  AnaIO::kSignal = AnaIO::kSignal && (AnaIO::finProtonmomentum>0.45 && AnaIO::fin2Pmom<0.45);
  if(!kpi0){
    AnaIO::kSignal = AnaIO::kSignal && (AnaIO::finPimomentum>0.15);
  }
}

int getTruthFromRec(const int recidx, int & pdg, double & momentum)
{
  const int truthID = (*AnaIO::reco_daughter_PFP_true_byHits_ID)[recidx];
  int trueidx = -999;
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

  pdg = -999;
  momentum = -999;

  if(trueidx>=0){
    pdg = (*AnaIO::true_beam_daughter_PDG)[trueidx];

    const TVector3 vec((*AnaIO::true_beam_daughter_startPx)[trueidx], (*AnaIO::true_beam_daughter_startPy)[trueidx], (*AnaIO::true_beam_daughter_startPz)[trueidx]);
    momentum = vec.Mag();
  }

  return trueidx;
}

int getNTrack(const bool kpi0, const bool ksig, int & nproton, int & ngamma, int & nmichel, const bool kprint=false, const bool kfill=false)
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
  int ntracks = 0;
  ngamma = 0;
  nmichel = 0;
  nproton = 0;
 
  for(int ii=0; ii<recsize; ii++){

    //---> need to be done before any counting!!!
    const vector<double>  dedxarray = (*AnaIO::reco_daughter_allTrack_calibrated_dEdX_SCE)[ii];
    const int NdEdx = dedxarray.size();
    const int ndedxcut = kpi0 ? 6 : 16; //need E[3], at least 4 cls
    printf("check cut kpi0 %d ndedxcut %d\n", kpi0, ndedxcut);
    if(NdEdx<ndedxcut){
      //do not count it as anything
      continue;
    }
    //<--- need to be done before any counting!!!

    //track and em scores at 0.5 look reasonable as baseline (after 0.5 the proton fraction look flat)
    //michel is not found (all below 0.5, those above 0.5 are shower and other_type)
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

    //nhit 260 is just clean-up
    const int nhits          = (*AnaIO::reco_daughter_PFP_nHits)[ii];
    /*//they are indeed different: 236 80
    if(nhits != NdEdx){
      printf("nhits!=dedx size %d %d\n", nhits, NdEdx); exit(1);
    }
    */

    //========== proton tagging now!
  
    const double startE2 = dedxarray[2];
    const double startE3 = dedxarray[3];
    const double lastE2 = dedxarray[NdEdx-1-2];
    const double lastE3 = dedxarray[NdEdx-1-3];

    const double chi2 = (*AnaIO::reco_daughter_allTrack_Chi2_proton)[ii];
    const double ndof = (*AnaIO::reco_daughter_allTrack_Chi2_ndof)[ii];
    const double Chi2NDF = chi2/(ndof+1E-10);

    bool isSelProton = false;
    double recp        = -999;
    //all signal protons have nhit below 260
    //all signal protons have startE3 > 9
    //with Chi2NDF cut  44/191 = 23% purity; without 46/197 = 23% purity, slightly higher efficiency
    //test if(startE2>10 && nhits<260 && startE3>9 && Chi2NDF<50){
    if(startE2>10 && nhits<260 && startE3>9){
      isSelProton = true;
      nproton++;

      recp = (*AnaIO::reco_daughter_allTrack_momByRange_proton)[ii];
    }
    else{
      //temporary pi+ momentum
      recp = (*AnaIO::reco_daughter_allTrack_momByRange_muon)[ii];
    }
    //========== proton tagging done!

    int pdg = -999;
    double truemomentum = -999;
    const int trueidx = getTruthFromRec(ii, pdg, truemomentum);

    int fillstktype = -999;
    if(pdg==2212){//proton
      fillstktype = 0;
    }
    else if(pdg==211){//pi+
      fillstktype = 1;
    }
    else if(pdg==22){//gamma
      fillstktype = 2;
    }
    else if(pdg==-999){//shower no true
      fillstktype = 3;
    }
    else{//all others
      fillstktype = 4;
    }

    if(kprint){
      printf("test ksig %d rec %d/%d trueidx %d pdg %d truemomentum %f recp %f resolution %f startE2 %f startE3 %f nhits %d trackScore %f emScore %f michelScore %f sum %f chi2 %f ndof %f chi2/ndof %f\n", ksig, ii, recsize, trueidx, pdg, truemomentum, recp, recp/truemomentum-1, startE2, startE3,  nhits, trackScore, emScore, michelScore, trackScore+emScore+michelScore, chi2, ndof, Chi2NDF);
    }

    if(kfill){
      if(pdg==2212 && isSelProton){
        AnaIO::hProtonMomentumRes->Fill(truemomentum, recp/truemomentum-1);
      }
      else if(pdg==211 && !isSelProton){
        AnaIO::hPiMomentumRes->Fill(truemomentum, recp/truemomentum-1);
      }

      AnaIO::hCutNdEdx->Fill(NdEdx, fillstktype);
      AnaIO::hCutstartE2->Fill(startE2, fillstktype);
      AnaIO::hCutstartE3->Fill(startE3, fillstktype);
      AnaIO::hCutlastE2->Fill(lastE2, fillstktype);
      AnaIO::hCutlastE3->Fill(lastE3, fillstktype);

      AnaIO::hCutnHits->Fill(nhits, fillstktype);
      AnaIO::hCutChi2NDF->Fill(Chi2NDF, fillstktype);
      AnaIO::hCuttrackScore->Fill(trackScore, fillstktype);
      AnaIO::hCutemScore->Fill(emScore, fillstktype);
      AnaIO::hCutmichelScore->Fill(michelScore, fillstktype);
    }
    
  }

  if(kprint){
    printf("\n\n");
  }

  return ntracks;
}

double getRecFromTruth(const int protonIdx, const vector<double> * mombyrange)
{
  const int truthID = (*AnaIO::true_beam_daughter_ID)[protonIdx];

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

  //____________________________________________________________________________________________________

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
   
    //____________________________________________________________________________________________________

    //calculate before any cuts! Only filled after ALL cuts!
    const int TruthBeamType = GetParticleType(AnaIO::true_beam_PDG);
   
    setFullSignal(kPiZero);

    //for to see within signal
    if(0){//====================================================== switch between signal sample and full sample
      if(!AnaIO::kSignal){
        continue;
      }
    }
    else{
      //0. true beam particle //"In data the beam Instrumentation is able to filter for these events but is not inside the MC" so read data file has this
      if(AnaIO::true_beam_PDG != 211 && AnaIO::true_beam_PDG != -13){
        continue;
      }
      
      //without 20 sig in 36 events; with 18/32, both 56% purity -> same purity, without this has higher efficiency
      if(0){//====================================================== switch off pion analysis cuts
        //x. primary beam type 
        AnaIO::hRecoBeamType->Fill(AnaIO::reco_beam_type);
        //no effect, shadowed by TMeanStart cut
        //    if(AnaIO::reco_beam_type!=13){//13: Pandora "track like"
        //      continue;
        //    }
        
        //1. beam position MC cut, need MC truth, how is it possible in analysis?
        const bool kBeamPosPass = GetBeamPosPass();
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
      }//====================================================== switch off pion analysis cuts
    }//====================================================== switch between signal sample and full sample

    //3. n track daughter
    int cutnproton = 0;
    int cutngamma = 0;
    int cutnmichel = 0;

    AnaIO::nTrack = getNTrack(kPiZero, AnaIO::kSignal, cutnproton, cutngamma, cutnmichel, false, true);

    AnaIO::hBeamNTrack->Fill(AnaIO::nTrack, cutnproton);
    AnaIO::hSignalVsBeamNTrack->Fill(AnaIO::nTrack, AnaIO::kSignal);

    /*
      PiPlus: ngamma = 0, nproton=1 -> 75/218=34% in signal sample; 63/467 = 13% purity, 63/218 = 29% efficiency, p*e = 63*63/218/467 = 3.9%
      PiPlus: ngamma = 0, nproton=1,  noly nhit<260 and startE3>9 and Chi2NDF<50 is proton, no pion-analysis pre-cuts! -> 61/380 = 16% 
      PiPlus: ngamma = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 and Chi2NDF<50 is proton, no pion-analysis pre-cuts! -> 44/191 = 23% purity, e*p = 44.*44./191./218. = 4.6%
      PiPlus: ngamma = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9                is proton, no pion-analysis pre-cuts! -> 46/197 = 23% purity, slightly higher efficiency
      PiPlus: ngamma = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=4 instead of 3) is proton, no pion-analysis pre-cuts! -> 49/199 = 25%
      PiPlus: ngamma = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=5 instead of 3) is proton, no pion-analysis pre-cuts! -> 48/201 = 24%
      PiPlus: ngamma = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=6 instead of 3) is proton, no pion-analysis pre-cuts! -> 50/202 = 25%
      PiPlus: ngamma = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=7 instead of 3) is proton, no pion-analysis pre-cuts! -> 50/200 = 25%, e*p = 5.7%
      PiPlus: ngamma = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=8 instead of 3) is proton, no pion-analysis pre-cuts! -> 50/194 = 26%, e*p = 5.9%
      PiPlus: ngamma = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=9 instead of 3) is proton, no pion-analysis pre-cuts! -> 47/176 = 27%, e*p = 5.8%
      PiPlus: ngamma = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=10 instead of 3)is proton, no pion-analysis pre-cuts! -> 46/168 = 27%, e*p = 5.8%
      PiPlus: ngamma = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=13 instead of 3)is proton, no pion-analysis pre-cuts! -> 43/159 = 27%, e*p = 5.3%
      PiPlus: ngamma = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=14 instead of 3)is proton, no pion-analysis pre-cuts! -> 42/149 = 28%, e*p = 5.4%
      PiPlus: ngamma = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=15 instead of 3)is proton, no pion-analysis pre-cuts! -> 40/135 = 30%, e*p = 5.4%
      PiPlus: ngamma = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=17 instead of 3)is proton, no pion-analysis pre-cuts! -> 39/126 = 31%
      PiPlus: ngamma = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=20 instead of 3)is proton, no pion-analysis pre-cuts! -> 34/114 = 30%
      PiPlus: ngamma = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=25 instead of 3)is proton, no pion-analysis pre-cuts! -> 25/82 = 30%
      PiPlus: ngamma = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=30 instead of 3)is proton, no pion-analysis pre-cuts! -> 22/71 = 31%

      (*) PiPlus: ngamma = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=16 instead of 3)is proton, no pion-analysis pre-cuts! -> 40/127 = 32%, e*p = 5.8%


      PiZero: ngamma = 2, nmichel=0, nproton=1 -> 32/260=12% in signal sample; 27/68 = 40% purity, 27/260 = 10% efficiency, p*e = 27*27/68/260 = 4.1%
      PiZero: ngamma = 2, nmichel=0, nproton=1, ntrack=1 -> ; 19/40 = 48% purity
      PiZero: ngamma = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 is proton -> 18/37=49% purity
      PiZero: ngamma = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 is proton -> 18/32=56% purity, p*e = 18.*18./32./260. = 3.9%
      PiZero: ngamma = 2, nmichel=0, nproton=1,           noly nhit<260 and startE3>9 is proton -> 21/47=45% purity -> so it is important to have ntrack=1

      PiZero: ngamma = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 and Chi2NDF<50 is proton, no pion-analysis pre-cuts! -> 18/34 = 53% purity
      PiZero: ngamma = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 is proton, no pion-analysis pre-cuts! -> 20/36=56% purity, eff*purity = 20.*20./36./260. = 4.3%
      PiZero: ngamma = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=4 instead of 3 ) is proton, no pion-analysis pre-cuts! -> 19/33=58%
      PiZero: ngamma = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=5 instead of 3 ) is proton, no pion-analysis pre-cuts! -> 21/34=62%
      PiZero: ngamma = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=7 instead of 3 ) is proton, no pion-analysis pre-cuts! -> 20/33=61%
      PiZero: ngamma = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=8 instead of 3 ) is proton, no pion-analysis pre-cuts! -> 19/31=61%, e*p = 4.5%
      PiZero: ngamma = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=10 instead of 3 )is proton, no pion-analysis pre-cuts! -> 13/20=65%, e*p = 3.3%
      PiZero: ngamma = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=20 instead of 3 )is proton, no pion-analysis pre-cuts! -> 9/10=90%,  e*p = 3.1%

      (*) PiZero: ngamma = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=6 instead of 3 ) is proton, no pion-analysis pre-cuts! -> 22/35=63%, e*p = 5.3%

     */

    if(1){//====================================================== switch of doing selection cuts
      if(cutnproton!=1){
        continue;
      }
      
      if(!kPiZero){
        if(cutngamma>0){
          continue;
        }
        //not found in signal, need Michel to tell 2-proton events from p-pi; ignore michel count for the momentum, things will change after implementing it.
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
    }//====================================================== switch of doing selection cuts

    //just for printing
    getNTrack(kPiZero, AnaIO::kSignal, cutnproton, cutngamma, cutnmichel, true, false);

     /*   
    //x. Beam dEdx cut shadowed by beam filtering
    AnaIO::hBeamLen->Fill(AnaIO::reco_beam_len);
    AnaIO::hSignalVsLen->Fill(AnaIO::reco_beam_len, AnaIO::kSignal);
    if(!cutBeamdEdx(AnaIO::kSignal)){
      continue;
    }
    */

    //____________________________________________________________________________________________________

    //============== Benchmark after ALL cuts !!! =========================
    //benchmark

    AnaIO::hTruthBeamType->Fill(TruthBeamType);
    AnaIO::hTruthSignal->Fill(AnaIO::kSignal);
   
    //======================= NO cuts below ========================

    //--- to test
    //Fill kSignal vs variable; the reason for kSignal but not kBeam is the signal is interacting pion, not all pions. Directly choose matrix to optimize for it
    //AnaIO::hSigAfterVsLen->Fill(AnaIO::reco_beam_len, AnaIO::kSignal);

    //============= done loop
    
    tout->Fill();
  }

  cout<<"All entries "<<ientry<<endl;

  style::Process2DHist(lout);
  style::DrawHist(lout, "output", tag, true);
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
      (*recmomentum) = getRecFromTruth(protonIdx, mombyrange);
      tout->Fill();
    }
  }

  cout<<"All entries "<<ientry<<endl;

  style::DrawHist(lout, "output", tag);
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
