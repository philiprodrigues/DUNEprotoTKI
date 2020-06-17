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
#include "AnaCut.h"

#include <algorithm>    // std::sort
#include <vector>       // std::vector

using namespace AnaFunctions;
using namespace AnaUtils;
using namespace AnaCut;


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

int anaRec(TString finName, TList *lout, const TString tag, const int nEntryToStop = -999)
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
  int selBeamCount = 0;
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
      if(!CutBeamAllInOne(kMC)){
        continue;
      }
    }//====================================================== switch between signal sample and full sample

    //count beam after beam cut before topology cut
    selBeamCount++;

    //3. n track daughter
    if(!CutTopology(kPiZero)){
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

  return selBeamCount; 
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

  double mcBeamCount = -999;
  double dataBeamCount = -999;
  //=======================================================================================
  //------------------------- MC
  TString mcfinName = "input/protoDUNE_mc_reco_flattree.root";

  if(kTruth){
    anaTruth(mcfinName, mclout, tag);
  }
  else{
    if(mclout){
      mcBeamCount = anaRec(mcfinName, mclout, tag);
    }

    if(datalout){
      TString datafinName = "input/protoDUNE_data_reco_flattree.root";
      dataBeamCount = anaRec(datafinName, datalout, tag);
    }
  }
  //------------------------- Data
  //=======================================================================================

  printf("anaRec beamcount data: %.0f mc: %.0f\n", dataBeamCount, mcBeamCount);

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
  if(mclout){
    TDirectory * ld = gDirectory->mkdir("mc");
    ld->cd();
    mclout->Write();
  }
  if(datalout){
    TDirectory * ld = gDirectory->mkdir("data");
    ld->cd();
    datalout->Write();
  }

  fout->Save();
  fout->Close();

  return 0;
}
