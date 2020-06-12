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

#include <algorithm>    // std::sort
#include <vector>       // std::vector

using namespace AnaFunctions;
using namespace AnaUtils;
using namespace TreeIO;

double getTruncatedMean(vector<double> array, const int nsample)
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

void getdEdx(const vector<double> arraydEdx, vector<double> &startE, vector<double> &endE)
{
  const unsigned int ncls = arraydEdx.size();

  if(ncls<3){
    return;
  }

  //start from [2] because [0] and [1] in both start and last are weird
  for(unsigned int kk=2; kk<ncls; kk++){
    startE.push_back(arraydEdx[kk]);

    const double endpe = arraydEdx[ncls-1-kk];
    endE.push_back(endpe);
  }
}

bool cutBeamdEdx(const double varSignal)
{
  vector<double> startE, lastE;
  getdEdx( *reco_beam_calibrated_dEdX, startE, lastE);
  nBeamdEdxCls = startE.size();
  hnBeamdEdxCls->Fill(nBeamdEdxCls);
  if(nBeamdEdxCls<6){
    return false;
  }

  //no bragg peak
  beamStartE0 = startE[0];
  beamStartE1 = startE[1];
  beamStartE2 = startE[2];
  beamStartE3 = startE[3];
  beamStartE4 = startE[4];
  beamStartE5 = startE[5];
  
  beamTMeanStart = getTruncatedMean(startE, nBeamdEdxCls-6);
  
  hSignalVsStartE0->Fill(beamStartE0, varSignal);
  hSignalVsStartE1->Fill(beamStartE1, varSignal);
  hSignalVsStartE2->Fill(beamStartE2, varSignal);
  hSignalVsStartE3->Fill(beamStartE3, varSignal);
  hSignalVsStartE4->Fill(beamStartE4, varSignal);
  hSignalVsStartE5->Fill(beamStartE5, varSignal);
  
  hSignalVsTMeanStart->Fill(beamTMeanStart, varSignal);
  
  //has Bragg Peak
  beamLastE0 = lastE[0];
  beamLastE1 = lastE[1];
  beamLastE2 = lastE[2];
  beamLastE3 = lastE[3];
  beamLastE4 = lastE[4];
  beamLastE5 = lastE[5];
  
  beamTMeanLast = getTruncatedMean(lastE, 6);
  
  hSignalVsLastE0->Fill(beamLastE0, varSignal);
  hSignalVsLastE1->Fill(beamLastE1, varSignal);
  hSignalVsLastE2->Fill(beamLastE2, varSignal);
  hSignalVsLastE3->Fill(beamLastE3, varSignal);
  hSignalVsLastE4->Fill(beamLastE4, varSignal);
  hSignalVsLastE5->Fill(beamLastE5, varSignal);
  
  //both need to tune for different energy
  //it is so clean that no need to cut on last since there is no Bragg peak form proton any more
  if(beamTMeanStart>2.8){
    return false;
  }
  
  hSignalVsTMeanLast->Fill(beamTMeanLast, varSignal);
  
  hSigAfterVsStartE0->Fill(beamStartE0, varSignal);
  hSigAfterVsStartE1->Fill(beamStartE1, varSignal);
  hSigAfterVsStartE2->Fill(beamStartE2, varSignal);
  
  hSigAfterVsLastE0->Fill(beamLastE0, varSignal);
  hSigAfterVsLastE1->Fill(beamLastE1, varSignal);
  hSigAfterVsLastE2->Fill(beamLastE2, varSignal);
  
  hSigAfterVsStartE3->Fill(beamStartE3, varSignal);
  hSigAfterVsLastE3->Fill(beamLastE3, varSignal);
  hSigAfterVsStartE4->Fill(beamStartE4, varSignal);
  hSigAfterVsLastE4->Fill(beamLastE4, varSignal);
  hSigAfterVsStartE5->Fill(beamStartE5, varSignal);
  hSigAfterVsLastE5->Fill(beamLastE5, varSignal);
  
  return true;
}

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

  if(kSignal && htype){
    const int nbt = bufferType.size();
    for(int ii=0; ii<nbt; ii++){
      htype->Fill(bufferType[ii]);
    }
  }

  return vec;
}

int getNTrack(const vector<int> * reco_daughter_allTrack_ID)
{
  int ntrack = 0;
  for(unsigned int ii=0; ii<reco_daughter_allTrack_ID->size(); ii++){
    if((*reco_daughter_allTrack_ID)[ii] != -1) {
      ntrack++;
    }
  }
  return ntrack;
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

        getdEdx( (*reco_daughter_allTrack_calibrated_dEdX_SCE)[ii], startE, endE);
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
  TTree * tout = GetOutputTree(lout, tag);  
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
    //calculate before any cuts! Only filled after ALL cuts!
    const int TruthBeamType = getParticleType(true_beam_PDG);
    int  protonIdx = -999, piplusIdx = -999;
    bool tmpkSig = false;
    vector<TLorentzVector> vecPiP = getFSTruth(kPiZero, true_beam_daughter_PDG, true_beam_daughter_startPx, true_beam_daughter_startPy, true_beam_daughter_startPz, 0x0, nproton, nneutron, nPiZero, ngamma, maxgammaEnergy, protonIdx, piplusIdx, tmpkSig);

    finPimomentum = vecPiP[0].P();
    finProtonmomentum = vecPiP[1].P();
    fin2Pmom = vecPiP[2].P();

    //with phase space cut
    kSignal = (true_beam_PDG==211) &&  tmpkSig && (finProtonmomentum>0.45 && fin2Pmom<0.45);
    if(!kPiZero){
      kSignal = kSignal && (finPimomentum>0.15);
    }

    //const bool varTrueBeam = (true_beam_PDG==211);
    const bool varSignal = kSignal;
    //===========================================================
    //1. Beam dEdx cut
    hBeamLen->Fill(reco_beam_len);
    hSignalVsLen->Fill(reco_beam_len, varSignal);
  
    if(!cutBeamdEdx(varSignal)){
      continue;
    }
    //-> now signal purity 167/5274 = 3.2%, 2627 pi+ beam, 1947 e+

    /*
    //0. true beam particle //only for cross checking previous study
    if(true_beam_PDG != 211 && true_beam_PDG != -13){
      continue;
    }
    */

    //0. primary beam type 
    hRecoBeamType->Fill(reco_beam_type);
    //no effect, shadowed by TMeanStart cut
    /*
    if(reco_beam_type!=13){//13: Pandora "track like"
      continue;
    }
    */

    //1. beam position MC cut, need MC truth, how is it possible in analysis?
    const bool kBeamPosPass = GetBeamPosPass();
    hBeamPosPass->Fill(kBeamPosPass);
    if(!kBeamPosPass){
      continue;
    }
    //-> now signal purity 138/3537 = 3.9%, 2283 pi+ bea, 801 e+ beam

    //2. APA3 
    hBeamEndZ->Fill(reco_beam_endZ);
    hBeamEndZPass->Fill(!(reco_beam_endZ>=226));

    if(reco_beam_endZ>=226){
      hBeamEndZPass->Fill(false);
      continue;
    }
    //-> now signal purity 135/3143 = 4.3%, 2102 pi+ beam, 801 e+ beam
   
    //9. n track daughter
    nTrack = getNTrack(reco_daughter_allTrack_ID);
    hBeamNTrack->Fill(nTrack);
    hSignalVsBeamNTrack->Fill(nTrack, varSignal);
    if(kPiZero){
      //the purity profile looks weird for pizeor signal
      if(nTrack!=1){
        continue;
      }
    }
    else{
      if(nTrack!=2){
        continue;
      }
    }
    //-> now pi+p signal purity 78/940=8.3%, 611 pi+ beam, 270 e+ beam

    //============== Benchmark after ALL cuts !!! =========================
    //benchmark

    hTruthBeamType->Fill(TruthBeamType);
    hTruthSignal->Fill(kSignal);
   
    //======================= NO cuts below ========================

    //--- to test
    //Fill kSignal vs variable; the reason for kSignal but not kBeam is the signal is interacting pion, not all pions. Directly choose matrix to optimize for it
    hSigAfterVsLen->Fill(reco_beam_len, varSignal);

    //============= done loop
    
    tout->Fill();
  }

  cout<<"All entries "<<ientry<<endl;

  getProfileX(lout);
  drawHist(lout, "output", tag, true);
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
