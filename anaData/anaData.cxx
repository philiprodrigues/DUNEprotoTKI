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
#include "AnaIO.h"
#include "AnaUtils.h"
#include "AnaCut.h"

#include <algorithm>    // std::sort
#include <vector>       // std::vector

const bool gkOnlySignal = false;
const int gkDataBit = 3;

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
      TruthBeamType = AnaUtils::GetParticleType(AnaIO::true_beam_PDG);
      AnaUtils::SetFullSignal(kPiZero);
    }

    //for to see within signal
    if(gkOnlySignal){//====================================================== switch between signal sample and selected sample
      if(!AnaIO::kSignal){
        continue;
      }
    }
    else{
    
      if(!AnaCut::CutBeamAllInOne(kMC)){
        continue;
      }
      
      //count beam after beam cut before topology cut
      selBeamCount++;
      
      //3. n track daughter
      if(!AnaCut::CutTopology(kPiZero)){
        continue;
      }   
    }//====================================================== switch between signal sample and selected sample

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

    if(kMC){
      AnaIO::hTruthBeamType->Fill(TruthBeamType);
      AnaIO::hTruthSignal->Fill(AnaIO::kSignal);
    }

    //======================= NO cuts below ========================

    //--- to test
    //Fill kSignal vs variable; the reason for kSignal but not kBeam is the signal is interacting pion, not all pions. Directly choose matrix to optimize for it
    //AnaIO::hSigAfterVsLen->Fill(AnaIO::reco_beam_len, AnaIO::kSignal);

    //============= done loop
    
    tout->Fill();
  }

  cout<<"All entries "<<ientry<<endl;

  style::Process2DHist(lout);

  if(kMC){
    /*
seeana010.log:kPiZero 0 signal 224.0 all 224.0 purity 100.0%
seeana110.log:kPiZero 1 signal 246.0 all 246.0 purity 100.0%
     */
    const double nfullsig = kPiZero? 246.0 : 224.0;
    const double nsig = AnaIO::hTruthSignal->GetBinContent(2);
    const double nbk = AnaIO::hTruthSignal->GetBinContent(1);
    const double nall = nsig+nbk;
    const double purity = nsig/nall;
    const double eff = nsig/nfullsig;
    const double ep = eff*purity;
    printf("kPiZero %d fullsig %.1f signal %.1f all %.1f purity %.1f%% eff %.1f%% ep %.1f%%\n", kPiZero, nfullsig, nsig, nall, purity*100, eff*100, ep*100);
  }

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

    const int beamtype = AnaUtils::GetParticleType(AnaIO::true_beam_PDG);
    AnaIO::hbeamType->Fill(beamtype);

    if(beamtype!=AnaUtils::gkPiPlus){
      continue;
    }

    const int nd = AnaIO::true_beam_daughter_PDG->size();
    AnaIO::hndaughter->Fill(nd);
    if(nd==0){
      continue;
    }

    const TLorentzVector iniPiPlus(AnaIO::true_beam_endPx, AnaIO::true_beam_endPy, AnaIO::true_beam_endPz, AnaFunctions::PionMass());

    AnaIO::iniPimomentum = iniPiPlus.P();
    AnaIO::iniPitheta = iniPiPlus.Theta()*TMath::RadToDeg();

    int  protonIdx = -999, piplusIdx = -999;
    AnaIO::kSignal = false;
    vector<TLorentzVector> vecPiP = AnaUtils::GetFSTruth(kPiZero, protonIdx, piplusIdx, AnaIO::kSignal);
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
      AnaFunctions::getCommonTKI(targetA, targetZ, &iniPiPlus, &(vecPiP[0]), &(vecPiP[1]), AnaIO::dalphat, AnaIO::dphit, AnaIO::dpt, AnaIO::pn, AnaIO::finPitheta, AnaIO::finProtontheta);
      
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
      (*recmomentum) = AnaUtils::GetRecFromTruth(recidx, mombyrange);
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
    AnaUtils::PrintLegend();
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

  if(gkDataBit & 1){
    mclout = new TList;
  }
  if(gkDataBit & 2){
    datalout = new TList;
  }

  double mcBeamCount = -999;
  double dataBeamCount = -999;
  //=======================================================================================
  //------------------------- MC
  const TString mcfinName = "input/protoDUNE_mc_reco_flattree.root";

  if(kTruth){
    anaTruth(mcfinName, mclout, tag);
  }
  else{
    if(mclout){
      mcBeamCount = anaRec(mcfinName, mclout, tag);
    }

    if(datalout){
      const TString datafinName = "input/protoDUNE_data_reco_flattree.root";
      dataBeamCount = anaRec(datafinName, datalout, tag);
    }
  }
  //------------------------- Data
  //=======================================================================================
  const double plotscale = dataBeamCount/mcBeamCount;

  printf("anaRec beamcount data: %.0f mc: %.0f plotscale %f\n", dataBeamCount, mcBeamCount, plotscale);

  const bool kfast = false;//true;
  if(mclout){
    if(datalout){
      style::DrawHist(mclout, plotscale, datalout, "output", tag, true, kfast);
    }
    else{
      style::DrawHist(mclout, 1, 0x0, "output", tag, true, kfast);
    }
  }
  else if(datalout){
    style::DrawHist(datalout, 1, 0x0, "output", tag, true, kfast);
  }

  TFile * fout = new TFile(Form("output/outanaData_%s.root", tag.Data()),"recreate");
  if(mclout){
    TDirectory * ld = gDirectory->mkdir("mc");
    ld->cd();
    mclout->Write();
    gDirectory->cd("../");
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
