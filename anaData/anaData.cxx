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
#include "TreeIO.h"

using namespace AnaFunctions;
using namespace AnaUtils;
using namespace TreeIO;

vector<TLorentzVector> getFSTruth(const bool kPiZero, const vector<int> * pdg, const vector<double> * px, const vector<double> * py, const vector<double> * pz, TH1I * htype, int &nproton, int &nneutron, int &nPiZero, int & ngamma, double & maxgammaEnergy, int & protonIdx, int & piplusIdx, bool & kSignal)
{
  //if kPiZero false, return 1piplus (1st and 2nd protons)
  //if kPiZero true, return (1st PiZero) (1st and 2nd protons)
  //only use kPiZero in the end at filling

  const int np = pdg->size();
  nproton = 0;
  nneutron = 0;
  int npiplus = 0;
  nPiZero = 0;
  ngamma = 0;
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
    const int itype = getParticleType((*pdg)[ii]);
    const TVector3 tmpp( (*px)[ii], (*py)[ii], (*pz)[ii] );

    if(itype==gkPiPlus){
      npiplus++;
      pPiplus.SetVectM(tmpp, PionMass());

      piplusIdx = ii;
    }
    else if(itype==gkProton){
      //ii is the location in the original vector
      protonIndices[nproton]=ii;

      protonmom[nproton] = tmpp.Mag();
      bufferProtonmom.push_back(tmpp);
      nproton++;
    }
    else if(itype==gkPiZero){
      PiZeromom[nPiZero] = tmpp.Mag();
      bufferPiZeromom.push_back(tmpp);
      nPiZero++;
    }
    else if(itype==gkNeutron){
      nneutron++;
    }
    else if(itype==gkGamma){
      gammamom[ngamma] = tmpp.Mag();
      ngamma++;
    }
    else if(itype==gkPiMinus||itype==gkKaon){
      npartialBkg++;
    }

    bufferType.push_back(itype);
  }

  //proton=======================
  int leadingProtonID = 0, subldProtonID = -999;
  if(nproton>1){
    int protonsortid[nproton];
    TMath::Sort(nproton, protonmom, protonsortid);
    leadingProtonID = protonsortid[0];
    subldProtonID = protonsortid[1];
    //printf("test %d %f %d %f\n", leadingProtonID, protonmom[leadingProtonID], subldProtonID, protonmom[subldProtonID]);
  }
  if(nproton>0){
    pProton.SetVectM(bufferProtonmom[leadingProtonID], ProtonMass());
    protonIdx = protonIndices[leadingProtonID];
  }
  if(nproton>1){
    p2Proton.SetVectM(bufferProtonmom[subldProtonID], ProtonMass());
  }
  //PiZero==============================
  int leadingPiZeroID = 0;
  if(nPiZero>1){
    int PiZerosortid[nPiZero];
    TMath::Sort(nPiZero, PiZeromom, PiZerosortid);
    leadingPiZeroID = PiZerosortid[0];
  }
  if(nPiZero>0){
    pPiZero.SetVectM(bufferPiZeromom[leadingPiZeroID], PiZeroMass());
  }
  //gamma====================================
  int leadinggammaID = 0;
  if(ngamma>1){
    int gammasortid[ngamma];
    TMath::Sort(ngamma, gammamom, gammasortid);
    leadinggammaID = gammasortid[0];
  }
  maxgammaEnergy = gammamom[leadinggammaID];
  //==========================================

  //fill vec regardless of kSignal
  vector<TLorentzVector> vec;
  vec.push_back(kPiZero?pPiZero:pPiplus);
  vec.push_back(pProton);
  vec.push_back(p2Proton);
    
  if(npiplus!=1){
    piplusIdx = -999;
  }

  //Now need kPiZero
  kSignal = false;
  if(nproton>=1 && npartialBkg ==0){
    if(kPiZero){
      if(nPiZero>0 && npiplus==0){
        kSignal = true;
      }
    }
    else{//PiPlus
      if(npiplus==1 && nPiZero==0){
        kSignal = true;
      }
    }
  }

  if(kSignal){
    const int nbt = bufferType.size();
    for(int ii=0; ii<nbt; ii++){
      htype->Fill(bufferType[ii]);
    }
  }

  return vec;
}

double getRecFromTruth(const int targetid, const vector<int> * reco_daughter_PFP_true_byHits_ID, const vector<int> * reco_daughter_allTrack_ID, const vector<double> * reco_daughter_allTrack_momByRange, const vector<vector<double> >* reco_daughter_allTrack_calibrated_dEdX_SCE, vector<double> &endE, vector<double> &startE, const vector<double>* reco_daughter_allTrack_Chi2_proton, const vector<int>* reco_daughter_allTrack_Chi2_ndof, double & chi2, double & ndof)
{
  double rpm = -999;
  for(unsigned int ii = 0; ii < reco_daughter_PFP_true_byHits_ID->size(); ii++){
    if((*reco_daughter_PFP_true_byHits_ID)[ii] == targetid) {
      if((*reco_daughter_allTrack_ID)[ii] != -1) {
        rpm = (*reco_daughter_allTrack_momByRange)[ii];

        chi2 = (*reco_daughter_allTrack_Chi2_proton)[ii];
        ndof = (*reco_daughter_allTrack_Chi2_ndof)[ii];

        if(!reco_daughter_allTrack_calibrated_dEdX_SCE){
          printf("reco_daughter_allTrack_calibrated_dEdX_SCE null!\n"); exit(1);
        }
        const vector<double> arraydEdx = (*reco_daughter_allTrack_calibrated_dEdX_SCE)[ii];
        const unsigned int ncls = arraydEdx.size();
        for(unsigned int kk=0; kk<ncls; kk++){
          const double endpe = arraydEdx[ncls-1-kk];
          endE.push_back(endpe);

          startE.push_back(arraydEdx[kk]);
          //printf("test %d %d %f\n", ii, kk, endpe);
        }
        //printf("==================================\n");
      }
    }
  }
  return rpm;
}

void anaRec(TList *lout, const TString tag, const int nEntryToStop = -999)
{
  const bool kPiZero = tag.Contains("MPiZero");
  const bool kTrackingProton = !tag.Contains("TrackingPiPlus");

  cout<<"\n\n                       Running anaRec kPiZero "<<kPiZero<<" TrackingProton "<<kTrackingProton<<endl<<endl;

  TTree * tree = GetInputTree("pionana/beamana");
  //test TTree * tout = GetOutputTree(lout, tag);
  IniRecHist(lout, tag);

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

  }

  cout<<"All entries "<<ientry<<endl;

  drawHist(lout, "output", tag);
}

void anaTruth(TList *lout, const TString tag, const int nEntryToStop = -999)
{
  const bool kPiZero = tag.Contains("MPiZero");
  const bool kTrackingProton = !tag.Contains("TrackingPiPlus");

  cout<<"\n\n                       Running anaTruth kPiZero "<<kPiZero<<" TrackingProton "<<kTrackingProton<<endl<<endl;

  TTree * tree = GetInputTree("pionana/beamana");
  TTree * tout = GetOutputTree(lout, tag);
  IniTruthHist(lout, tag);

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

    const int beamtype = getParticleType(true_beam_PDG);
    hbeamType->Fill(beamtype);

    if(beamtype!=gkPiPlus){
      continue;
    }

    const int nd = true_beam_daughter_PDG->size();
    hndaughter->Fill(nd);
    if(nd==0){
      continue;
    }

    const TLorentzVector iniPiPlus(true_beam_endPx, true_beam_endPy, true_beam_endPz, PionMass());

    iniPimomentum = iniPiPlus.P();
    iniPitheta = iniPiPlus.Theta()*TMath::RadToDeg();

    int  protonIdx = -999, piplusIdx = -999;
    kSignal = false;
    vector<TLorentzVector> vecPiP = getFSTruth(kPiZero, true_beam_daughter_PDG, true_beam_daughter_startPx, true_beam_daughter_startPy, true_beam_daughter_startPz, hselectedDaughterType, nproton, nneutron, nPiZero, ngamma, maxgammaEnergy, protonIdx, piplusIdx, kSignal);
    hnexcl->Fill(kSignal);

    finPimomentum = vecPiP[0].P();
    finProtonmomentum = vecPiP[1].P();
    fin2Pmom = vecPiP[2].P();

     //signal truth
    if(kSignal){
      hnproton->Fill(nproton);
      hnneutron->Fill(nneutron);
      hnPiZero->Fill(nPiZero);
      
      //re-calculate final pi p theta w.r.t. iniPi
      const int targetA = 40;
      const int targetZ = 18;
      getCommonTKI(targetA, targetZ, &iniPiPlus, &(vecPiP[0]), &(vecPiP[1]), dalphat, dphit, dpt, pn, finPitheta, finProtontheta);
      
      hmomIniPi->Fill(iniPimomentum);
      hmomFinPi->Fill(finPimomentum);
      hmomFinProton->Fill(finProtonmomentum);
      hdalphat->Fill(dalphat);
      hdphit->Fill(dphit);
      hdpt->Fill(dpt);
      hpn->Fill(pn);
    }

    vector<double>* reco_daughter_allTrack_momByRange = 0x0;
    bool kDoTracking = false;
    double * recmomentum = 0x0;

    if(kTrackingProton){
      reco_daughter_allTrack_momByRange = reco_daughter_allTrack_momByRange_proton;
      kDoTracking = (protonIdx>=0);
      recmomentum = &recProtonmomentum;
    }
    else{
      reco_daughter_allTrack_momByRange = reco_daughter_allTrack_momByRange_muon;
      kDoTracking = (piplusIdx>=0);
      recmomentum = &recPiPlusmomentum;
    }

    hksignp->Fill(kDoTracking, kSignal);
 
    //has true proton
    if(kDoTracking){
      const int targetIdx = (*true_beam_daughter_ID)[protonIdx];
      vector<double> lastE, startE;
      (*recmomentum) = getRecFromTruth(targetIdx, reco_daughter_PFP_true_byHits_ID, reco_daughter_allTrack_ID, reco_daughter_allTrack_momByRange, reco_daughter_allTrack_calibrated_dEdX_SCE, lastE, startE, reco_daughter_allTrack_Chi2_proton, reco_daughter_allTrack_Chi2_ndof, chi2, ndof);

      ndEdxCls = lastE.size();
      if(ndEdxCls>=6){
        lastE0 = lastE[0];
        lastE1 = lastE[1];
        lastE2 = lastE[2];
        lastE3 = lastE[3];
        lastE4 = lastE[4];
        lastE5 = lastE[5];

        startE0 = startE[0];
        startE1 = startE[1];
        startE2 = startE[2];
        startE3 = startE[3];
        startE4 = startE[4];
        startE5 = startE[5];
      }

      tout->Fill();
    }
  }

  cout<<"All entries "<<ientry<<endl;

  drawHist(lout, "output", tag);
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
