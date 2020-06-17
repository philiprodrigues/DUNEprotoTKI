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
  //AnaIO::kSignal = AnaIO::kSignal && (AnaIO::finProtonmomentum>0.45 && AnaIO::fin2Pmom<0.45);
  AnaIO::kSignal = AnaIO::kSignal && (AnaIO::finProtonmomentum<1 && AnaIO::finProtonmomentum>0.45 && AnaIO::fin2Pmom<0.45);

  //no pi+ phase space cut
  //if(!kpi0){
    //AnaIO::kSignal = AnaIO::kSignal && (AnaIO::finPimomentum>0.15);
  //}

  //p_leading 0.45-1, p_sublieading < 0.45, no cut on pi+
  //p-pi0: 242 events
  //p-pi+: 220 events
}

int getTruthFromRec(const int recidx, int & pdg, double & momentum)
{
  pdg = -999;
  momentum = -999;
  int trueidx = -999;

  if(AnaIO::reco_daughter_PFP_true_byHits_ID && AnaIO::reco_daughter_PFP_true_byHits_ID->size() ){//in data this is size = 0 but not null

    const int truthID = (*AnaIO::reco_daughter_PFP_true_byHits_ID)[recidx];
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
    
  }

  return trueidx;
}

TLorentzVector * getPiZero(const int ksig, const vector<TLorentzVector> & shws,  const bool kprint, const bool kfill)
{
  const int shsize = shws.size();
  if(kfill){
    AnaIO::hRecPi0Nshower->Fill(shsize, !ksig);
  }

  vector<TLorentzVector> piarr;
  int maxEidx = -999;
  double pi0Emax = -999;
  if(shsize>=2){
    for(int ii = 0; ii<shsize; ii++){
      for(int kk = ii+1; kk<shsize; kk++){
        const TLorentzVector tmpgg = shws[ii]+shws[kk];
        if(kprint){
          printf("\n============ %d %d %d \n", shsize, ii, kk);
          tmpgg.Print();        
        }

        if(tmpgg.E() > pi0Emax){
          pi0Emax = tmpgg.E();
          maxEidx = piarr.size();
        }

        piarr.push_back(tmpgg);
      }
    }
  }  

  TLorentzVector * outcopy = 0x0;
  if(maxEidx>=0){
    outcopy = new TLorentzVector(piarr[maxEidx]);
  }

  return outcopy;
}

int getNTrack(const bool kpi0, const bool ksig, int & nproton, int & nshower, int & nmichel, TLorentzVector *  & leadingPi0, const bool kprint=false, const bool kfill=false, const int kdebug = 0)
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
  nshower = 0;
  nmichel = 0;
  nproton = 0;
 
  vector<TLorentzVector> showerArray;

  static bool kShowCut = true;

  for(int ii=0; ii<recsize; ii++){

    int pdg = -999;
    double truemomentum = -999;
    const int trueidx = getTruthFromRec(ii, pdg, truemomentum);

    int fillstktype = -999;
    if(pdg==2212){//proton
      fillstktype = 0;
    }
    else if(pdg==211){//pi+ put it beyond others
      fillstktype = 1;
    }
    else if(pdg==22){//gamma
      fillstktype = 4;
    }
    else if(pdg==-999){//shower no true
      fillstktype = 2;
    }
    else{//all others
      fillstktype = 3;
    }

    //---> need to be done before any counting!!!
    const vector<double>  dedxarray = (*AnaIO::reco_daughter_allTrack_calibrated_dEdX_SCE)[ii];
    const int NdEdx = dedxarray.size();
    if(kfill){
      AnaIO::hCutNdEdx->Fill(NdEdx, fillstktype);
    }
    const int ndedxcut = kpi0 ? 6 : 16; //need E[3], at least 4 cls
    if(kShowCut){
      printf("check cut kpi0 %d ndedxcut %d\n", kpi0, ndedxcut);
    }
    if(NdEdx<ndedxcut){
      //do not count it as anything
      continue;
    }
    //<--- need to be done before any counting!!!

    //track and em scores at 0.5 look reasonable as baseline (after 0.5 the proton fraction look flat)
    //michel is not found (all below 0.5, those above 0.5 are shower and other_type)
    const double trackScore  = (*AnaIO::reco_daughter_PFP_trackScore)[ii];
    if(trackScore>0.5){
      if( (*AnaIO::reco_daughter_allTrack_ID)[ii]==-1 ){
        printf("bad track ID! %d\n", ii); exit(1);
      }
      ntracks++;
    }
    const double emScore     = (*AnaIO::reco_daughter_PFP_emScore)[ii];
    if(emScore>0.5){
      if( (*AnaIO::reco_daughter_allShower_ID)[ii]==-1 ){
        printf("bad shower ID! %d\n", ii); exit(1);
      }
      nshower++;
      
      const TVector3 showerDir( (*AnaIO::reco_daughter_allShower_dirX)[ii], (*AnaIO::reco_daughter_allShower_dirY)[ii], (*AnaIO::reco_daughter_allShower_dirZ)[ii] );
      const TVector3 showerMomentum = showerDir.Unit()*(*AnaIO::reco_daughter_allShower_energy)[ii] * 1E-3; //MeV to GeV
      const TLorentzVector showerLv( showerMomentum, showerMomentum.Mag() );
      showerArray.push_back(showerLv);
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
    const double cutNH = 260;
    if(kpi0){
      //all signal protons have nhit below 260
      //all signal protons have startE3 > 9
      //with Chi2NDF cut  44/191 = 23% purity; without 46/197 = 23% purity, slightly higher efficiency
      //test if(startE2>10 && nhits<260 && startE3>9 && Chi2NDF<50){
      
      const double cutSE2 = 10;
      const double cutSE3 = 9;
      isSelProton = (startE2>cutSE2 && nhits<cutNH && startE3>cutSE3);
      if(kShowCut){
        printf("check cut kpi0 %d proton tag startE2 %.2f nhits %.0f startE3 %.2f\n", kpi0, cutSE2, cutNH, cutSE3); 
      }
    }
    else{
      const double cutCHI = 50;
      isSelProton = (Chi2NDF<cutCHI && nhits<cutNH);
      if(kShowCut){
        printf("check cut kpi0 %d proton tag Chi2NDF %.2f nhits %.0f\n", kpi0, cutCHI, cutNH);
      }
    }

    double recp        = -999;
    if(isSelProton){
      nproton++;
      recp = (*AnaIO::reco_daughter_allTrack_momByRange_proton)[ii];
    }
    else{
      //temporary pi+ momentum
      recp = (*AnaIO::reco_daughter_allTrack_momByRange_muon)[ii];
    }
    //========== proton tagging done!
   

    if(kprint){
      printf("test sig %d rii %d/%d tii %4d pdg %4d sE2 %6.1f sE3 %6.1f lE2 %6.1f lE3 %6.1f nhi %4d tkS %6.1f emS %6.1f miS %6.1f sum %6.1f chf %6.1f\n", ksig, ii, recsize, trueidx, pdg, startE2, startE3, lastE2, lastE3, nhits, trackScore, emScore, michelScore, trackScore+emScore+michelScore, Chi2NDF);
    }

    if(kfill){
      if(pdg==2212 && isSelProton){
        AnaIO::hProtonMomentumRes->Fill(truemomentum, recp/truemomentum-1);
      }
      else if(pdg==211 && !isSelProton){
        AnaIO::hPiMomentumRes->Fill(truemomentum, recp/truemomentum-1);
      }

      if(!kdebug || 
         (kdebug==1 && isSelProton) ||
         (kdebug==2 && !isSelProton)
         ){
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

    kShowCut = false;
  }

  //only fill nshowers when doing fill
  leadingPi0 = getPiZero(ksig, showerArray, false, kfill);

  if(kprint){
    printf("\n\n");
  }

  return ntracks;
}

bool cutTopology(const bool kpi0)
{
    /*
      PiPlus: nshower = 0, nproton=1 -> 75/218=34% in signal sample; 63/467 = 13% purity, 63/218 = 29% efficiency, p*e = 63*63/218/467 = 3.9%
      PiPlus: nshower = 0, nproton=1,  noly nhit<260 and startE3>9 and Chi2NDF<50 is proton, no pion-analysis pre-cuts! -> 61/380 = 16% 
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 and Chi2NDF<50 is proton, no pion-analysis pre-cuts! -> 44/191 = 23% purity, e*p = 44.*44./191./218. = 4.6%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9                is proton, no pion-analysis pre-cuts! -> 46/197 = 23% purity, slightly higher efficiency
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=4 instead of 3) is proton, no pion-analysis pre-cuts! -> 49/199 = 25%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=5 instead of 3) is proton, no pion-analysis pre-cuts! -> 48/201 = 24%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=6 instead of 3) is proton, no pion-analysis pre-cuts! -> 50/202 = 25%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=7 instead of 3) is proton, no pion-analysis pre-cuts! -> 50/200 = 25%, e*p = 5.7%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=8 instead of 3) is proton, no pion-analysis pre-cuts! -> 50/194 = 26%, e*p = 5.9%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=9 instead of 3) is proton, no pion-analysis pre-cuts! -> 47/176 = 27%, e*p = 5.8%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=10 instead of 3)is proton, no pion-analysis pre-cuts! -> 46/168 = 27%, e*p = 5.8%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=13 instead of 3)is proton, no pion-analysis pre-cuts! -> 43/159 = 27%, e*p = 5.3%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=14 instead of 3)is proton, no pion-analysis pre-cuts! -> 42/149 = 28%, e*p = 5.4%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=15 instead of 3)is proton, no pion-analysis pre-cuts! -> 40/135 = 30%, e*p = 5.4%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=17 instead of 3)is proton, no pion-analysis pre-cuts! -> 39/126 = 31%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=20 instead of 3)is proton, no pion-analysis pre-cuts! -> 34/114 = 30%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=25 instead of 3)is proton, no pion-analysis pre-cuts! -> 25/82 = 30%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=30 instead of 3)is proton, no pion-analysis pre-cuts! -> 22/71 = 31%

      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and chi2<20 (NdEdx>=16 instead of 3, no startE2, E3 cuts at all)is proton, no pion-analysis pre-cuts! -> 39/122 = 32%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and chi2<23 (NdEdx>=16 instead of 3, no startE2, E3 cuts at all)is proton, no pion-analysis pre-cuts! -> 41/123 = 33%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and chi2<25 (NdEdx>=16 instead of 3, no startE2, E3 cuts at all)is proton, no pion-analysis pre-cuts! -> 42/125 = 34%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and chi2<30 (NdEdx>=16 instead of 3, no startE2, E3 cuts at all)is proton, no pion-analysis pre-cuts! -> 42/127 = 33%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and chi2<35 (NdEdx>=16 instead of 3, no startE2, E3 cuts at all)is proton, no pion-analysis pre-cuts! -> 43/128 = 34%, e*p = 6.6%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and chi2<40 (NdEdx>=16 instead of 3, no startE2, E3 cuts at all)is proton, no pion-analysis pre-cuts! -> 43/127 = 34%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and chi2<60 (NdEdx>=16 instead of 3, no startE2, E3 cuts at all)is proton, no pion-analysis pre-cuts! -> 44/138 = 32%

      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=16 instead of 3)is proton, no pion-analysis pre-cuts! -> 40/127 = 32%, e*p = 5.8%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and chi2<50 (NdEdx>=16 instead of 3, no startE2, E3 cuts at all)is proton, HAS pion-analysis pre-cuts! -> 39/111 = 35%, e*p = 6.3%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and chi2<50 (NdEdx>=16 instead of 3, no startE2, E3 cuts at all)is proton, no pion-analysis pre-cuts! -> 44/132 = 33%, e*p = 6.7%

      (*) PiPlus New Def:  2) p_lading 0.45-1, p_sl<0.45. PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and chi2<50 (NdEdx>=16 instead of 3, no startE2, E3 cuts at all)is proton, no pion-analysis pre-cuts! -> 43/132 = 33%, e*p = 43.*43./132./220. = 6.4%

      <-- //p_leading 0.45-1, p_sublieading < 0.45, no cut on pi+
      <-- //p-pi+: 220 events

      PiZero: nshower = 2, nmichel=0, nproton=1 -> 32/260=12% in signal sample; 27/68 = 40% purity, 27/260 = 10% efficiency, p*e = 27*27/68/260 = 4.1%
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1 -> ; 19/40 = 48% purity
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 is proton -> 18/37=49% purity
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 is proton -> 18/32=56% purity, p*e = 18.*18./32./260. = 3.9%
      PiZero: nshower = 2, nmichel=0, nproton=1,           noly nhit<260 and startE3>9 is proton -> 21/47=45% purity -> so it is important to have ntrack=1

      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 and Chi2NDF<50 is proton, no pion-analysis pre-cuts! -> 18/34 = 53% purity
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 is proton, no pion-analysis pre-cuts! -> 20/36=56% purity, eff*purity = 20.*20./36./260. = 4.3%
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=4 instead of 3 ) is proton, no pion-analysis pre-cuts! -> 19/33=58%
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=5 instead of 3 ) is proton, no pion-analysis pre-cuts! -> 21/34=62%
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=7 instead of 3 ) is proton, no pion-analysis pre-cuts! -> 20/33=61%
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=8 instead of 3 ) is proton, no pion-analysis pre-cuts! -> 19/31=61%, e*p = 4.5%
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=10 instead of 3 )is proton, no pion-analysis pre-cuts! -> 13/20=65%, e*p = 3.3%
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=20 instead of 3 )is proton, no pion-analysis pre-cuts! -> 9/10=90%,  e*p = 3.1%

      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and chi2<25 ( NdEdx>=6 instead of 3 but no startE2,E3 at all) is proton, no pion-analysis pre-cuts! -> 20/33=61%
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and chi2<30 ( NdEdx>=6 instead of 3 but no startE2,E3 at all) is proton, no pion-analysis pre-cuts! -> 20/33=61%
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and chi2<35 ( NdEdx>=6 instead of 3 but no startE2,E3 at all) is proton, no pion-analysis pre-cuts! -> 20/34=59%

      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=6 instead of 3 ) is proton, HAS pion-analysis pre-cuts! -> 21/32=66%, e*p= 5.3%
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=6 instead of 3 ) is proton, no pion-analysis pre-cuts! -> 22/35=63%, e*p = 5.3%

      (*) PiZero: nshower >= 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=6 instead of 3 ) is proton, no pion-analysis pre-cuts! -> 30/47=64%, e*p = 7.4%, eff = 30/260 = 12%, this is the SAME for both 1) p_leading>0.45, p_sl<0.45, p+>0.15 and 2) p_lading 0.45-1, p_sl<0.45. Eff*pur for new def 2) is 30.*30./47./242 = 7.9%
      
      <-- //p_leading 0.45-1, p_sublieading < 0.45, no cut on pi+
      <-- //p-pi0: 242 events
     */

    int cutnproton = 0;
    int cutnshower = 0;
    int cutnmichel = 0;
    //vector<TLorentzVector> piZeroArray;
    TLorentzVector * leadingPi0 = 0x0;
 
    //debug=1 proton candidate; debug=2 non-proton candidate
    const int kdebug = 0;
    if(kdebug==0){
      //only filling
      AnaIO::nTrack = getNTrack(kpi0, AnaIO::kSignal, cutnproton, cutnshower, cutnmichel, leadingPi0, false, true, 0);
    }
    else{//just no filling, still no debug
      AnaIO::nTrack = getNTrack(kpi0, AnaIO::kSignal, cutnproton, cutnshower, cutnmichel, leadingPi0, false, false, 0);
    }

    int filleventtype = -999;
    if(AnaIO::kSignal){
      filleventtype = 0;
    }
    else if(AnaIO::true_beam_PDG==211){
      filleventtype = 1;
    }
    else{
      filleventtype = 2;
    }

    //----------- do cuts

    AnaIO::hCutnshower->Fill(cutnshower, filleventtype);
    if(kpi0){
      if(cutnshower<2){
        return false;
      }
    }
    else{
      if(cutnshower>0){
        return false;
      }
    }

    AnaIO::hCutnmichel->Fill(cutnmichel, filleventtype);
    if(kpi0){
      if(cutnmichel>0){
        return false;
      }
    }
    else{
      //not found in signal, need Michel to tell 2-proton events from p-pi; ignore michel count for the momentum, things will change after implementing it.
        /*
          if(cutnmichel!=1){
          continue;
          }
        */
    }

    AnaIO::hCutntrack->Fill(AnaIO::nTrack, filleventtype);
    if(kpi0){
      //with this 56% purity, without 45%
      if(AnaIO::nTrack!=1){
        return false;
      }
    }
    else{
      if(AnaIO::nTrack!=2){
        return false;
      }
    }

    AnaIO::hCutnproton->Fill(cutnproton, filleventtype);
    if(cutnproton!=1){
      return false;
    }
    
    //===================================================
    //only checking after-selection distributions
    if(kpi0){
      if(!leadingPi0){
        printf("leadingpi0 null!!\n"); exit(1);
      }
      AnaIO::hCutMpi0->Fill(leadingPi0->M(), filleventtype);
    }

    TLorentzVector *dummypi0 = 0x0;
    //no getting of pizeroarray
    if(kdebug==0){
      //just for printing, no filling
      getNTrack(kpi0, AnaIO::kSignal, cutnproton, cutnshower, cutnmichel, dummypi0, true, false, 0);
    }
    else{
      //debug mode: major background is -999 shower mocking protons, no efficient cut variables found
      //print, fill, and debug
      getNTrack(kpi0, AnaIO::kSignal, cutnproton, cutnshower, cutnmichel, dummypi0, true, true, kdebug);
    }

    return true;
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
      if((*AnaIO::reco_daughter_allTrack_ID)[ii] != -1) {//allTrack force reconstruction is successful assuming track

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
      /*there are indeed -1 due to fail allTrack assumption
      else{
        printf("getRecFromTruth id = -1!! %d %d\n", ii, (*AnaIO::reco_daughter_allTrack_ID)[ii]); exit(1);
      }
      */
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

void anaRec(TString finName, TList *lout, const TString tag, const int nEntryToStop = -999)
{
  bool kMC = true;
  if(!finName.Contains("_mc_")){
    kMC = false;
  }

  TFile * fin = new TFile(finName);
  if(!fin->IsOpen()){
    cout<<"fin not open!"<<endl;
    exit(1);
  }

  const bool kPiZero = tag.Contains("MPiZero");
  const bool kTrackingProton = !tag.Contains("TrackingPiPlus");

  cout<<"\n\n                       Running anaRec kMC "<<kMC<<" kPiZero "<<kPiZero<<" TrackingProton "<<kTrackingProton<<endl<<endl;

  TTree * tree = AnaIO::GetInputTree(fin, "pionana/beamana");
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
    int TruthBeamType = -999;
    if(kMC){
      TruthBeamType = GetParticleType(AnaIO::true_beam_PDG);
      setFullSignal(kPiZero);
    }

    //for to see within signal
    if(0){//====================================================== switch between signal sample and full sample
      if(!AnaIO::kSignal){
        continue;
      }
    }
    else{
      //0. true beam particle //"In data the beam Instrumentation is able to filter for these events but is not inside the MC" so read data file has this
      if(kMC){
        if(AnaIO::true_beam_PDG != 211 && AnaIO::true_beam_PDG != -13){
          continue;
        }
      }
      else{
        const bool data_beamID = GetbeamID((*AnaIO::data_BI_PDG_candidates));
        AnaIO::hCutbeamID->Fill(data_beamID);
        if(!data_beamID){
          continue;
        }
      }
      
      //for pi+ it is worse to include this pion-ana-cut (e*p: 6.7->6.3%), for pi0 it is the same e*p, better purity, worse efficiency.
      if(0){//====================================================== switch off pion analysis cuts
        //x. primary beam type 
        AnaIO::hRecoBeamType->Fill(AnaIO::reco_beam_type);
        //no effect, shadowed by TMeanStart cut
        if(AnaIO::reco_beam_type!=13){//13: Pandora "track like"
          continue;
        }
        
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
    if(!cutTopology(kPiZero)){
      continue;
    }   

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
}

void anaTruth(TString finName, TList *lout, const TString tag, const int nEntryToStop = -999)
{
  if(!finName.Contains("_mc_")){
    printf("anaTruth _mc_ not found in file name! %s\n", finName.Data()); exit(1);
  }

  TFile * fin = new TFile(finName);
  if(!fin->IsOpen()){
    cout<<"fin not open!"<<endl;
    exit(1);
  }

  const bool kPiZero = tag.Contains("MPiZero");
  const bool kTrackingProton = !tag.Contains("TrackingPiPlus");

  cout<<"\n\n                       Running anaTruth kPiZero "<<kPiZero<<" TrackingProton "<<kTrackingProton<<endl<<endl;

  TTree * tree = AnaIO::GetInputTree(fin, "pionana/beamana");
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
    double * recmomentum = 0x0;
    int recidx = -999;
    if(kTrackingProton){
      recidx = protonIdx;
      mombyrange = AnaIO::reco_daughter_allTrack_momByRange_proton;
      recmomentum = &AnaIO::recProtonmomentum;
    }
    else{
      recidx = piplusIdx;
      mombyrange = AnaIO::reco_daughter_allTrack_momByRange_muon;
      recmomentum = &AnaIO::recPiPlusmomentum;
    }

    const bool kDoTracking = (recidx>=0);
    AnaIO::hksignp->Fill(kDoTracking, AnaIO::kSignal);
 
    //has true proton
    if(kDoTracking){
      (*recmomentum) = getRecFromTruth(recidx, mombyrange);
      tout->Fill();
    }
  }

  cout<<"All entries "<<ientry<<endl;

}

void PrintLegend()
{
  TCanvas * c1 = 0x0;

  vector<TString> parType;
  parType.push_back("p");
  parType.push_back("#pi^{+}");
  parType.push_back("EM shower");
  parType.push_back("others");
  c1 = style::DrawLegend(parType, "f");
  c1->Print("output/legend_parType.eps");
  c1->Print("output/legend_parType.pdf");
  c1->Print("output/legend_parType.png");

  vector<TString> evtType;
  evtType.push_back("signal");
  evtType.push_back("background");
  evtType.push_back("#mu^{+} beam");
  c1 = style::DrawLegend(evtType, "f");
  c1->Print("output/legend_evtType.eps");
  c1->Print("output/legend_evtType.pdf");
  c1->Print("output/legend_evtType.png");
}

int main(int argc, char * argv[])
{
  if(argc!=4){
    cout<<"argc!=4 !"<<argc<<endl;
    cout<<"Print legend instead!"<<endl;
    PrintLegend();
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

  TList * mclout = 0x0;
  TList * datalout = 0x0;

  if(1){//switch for test
    mclout = new TList;
  }
  else{
    datalout = new TList;
  }

  //=======================================================================================
  //------------------------- MC
  TString mcfinName = "input/protoDUNE_mc_reco_flattree.root";

  if(kTruth){
    anaTruth(mcfinName, mclout, tag);
  }
  else{
    if(mclout){
      anaRec(mcfinName, mclout, tag);
    }

    if(datalout){
      TString datafinName = "input/protoDUNE_data_reco_flattree.root";
      anaRec(datafinName, datalout, tag);
    }
  }
  //------------------------- Data
  //=======================================================================================

  const double plotscale = 1.0;
  if(mclout){
    if(datalout){

    }
    else{
      style::DrawHist(mclout, plotscale, 0x0, "output", tag, true, false);
    }
  }
  else if(datalout){
    style::DrawHist(datalout, plotscale, 0x0, "output", tag, true, false);
  }

  TFile * fout = new TFile(Form("output/outanaData_%s.root", tag.Data()),"recreate");

  mclout->Write();
  fout->Save();
  fout->Close();

  return 0;
}
