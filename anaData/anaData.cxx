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

using namespace AnaFunctions;

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
  gkNucleus
};

int getParticleType(const int pdg)
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
  else if(pdg==-13){
    type = gkMuPlus;
  }
  else if(pdg==321||pdg==-321||pdg==310){
    type = gkKaon;
  }
  else if(pdg==2112){
    type = gkNeutron;
  }
  else if(pdg==22){
    type = gkGamma;
  }
  else if(pdg==14){
    type = gkNeutrino;
  }
  else if(pdg==3122||pdg==3212||pdg==3222){
    type = gkHyperon;
  }
  else if(pdg>9999){
    type = gkNucleus;
  }
  else{
    cout<<"getParticleType unknown pdg "<<pdg<<endl;
    exit(1);
  }

  return type;
}

vector<TLorentzVector> getPiP1P2(const bool kPiZero, const vector<int> * pdg, const vector<double> * px, const vector<double> * py, const vector<double> * pz, TH1I * htype, int &nproton, int &nneutron, int &nPiZero, int & ngamma, double & maxgammaEnergy, int & protonIdx, int & piplusIdx, bool & kSignal)
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

void drawHist(TList *lout, const TString tag)
{
  TCanvas * c1 = new TCanvas("c1"+tag, "", 1200, 800);
  style::PadSetup(c1);
  gPad->SetTopMargin(0.06);
  gPad->SetRightMargin(0.03);

  for(int ii=0; ii<lout->GetSize(); ii++){
    TH1 * hh = dynamic_cast<TH1*> (lout->At(ii));
    if(!hh){
      continue;
    }

    TH2 * htmp = dynamic_cast<TH2 *>(hh);
    const bool k2d = htmp;
    if(k2d){
      gStyle->SetOptStat(0);
      //gStyle->SetOptStat("enou");
    }
    else{
      gStyle->SetOptStat("enoum");
    }

    const TString tag = hh->GetName();
    cout<<"Printing "<<tag<<endl;

    style::ResetStyle(hh);
    hh->UseCurrentStyle();//can't go after setmarkersize
    hh->SetMaximum(hh->GetBinContent(hh->GetMaximumBin())*1.3);
    hh->SetMinimum(0);
    //hh->Draw(k2d?"box":"text hist");
    hh->UseCurrentStyle();
    hh->SetMarkerSize(3);
    hh->Draw("text hist");

    c1->Print("output/"+tag+".png");
    c1->Print("output/"+tag+".pdf");
    c1->Print("output/"+tag+".eps");
  }

}

double getRecMomentum(const int targetid, const vector<int> * reco_daughter_PFP_true_byHits_ID, const vector<int> * reco_daughter_allTrack_ID, const vector<double> * reco_daughter_allTrack_momByRange, const vector<vector<double> >* reco_daughter_allTrack_calibrated_dEdX_SCE, vector<double> &endE, vector<double> &startE, const vector<double>* reco_daughter_allTrack_Chi2_proton, const vector<int>* reco_daughter_allTrack_Chi2_ndof, double & chi2, double & ndof)
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

TString anaData(TList *lout, const TString tag, const int nEntryToStop = -999)
{
  const bool kPiZero = tag.Contains("MPiZero");
  const bool kTrackingProton = !tag.Contains("TrackingPiPlus");

  cout<<"\n\n                       Running kPiZero "<<kPiZero<<" TrackingProton "<<kTrackingProton<<endl<<endl;

  TTree * tree=(TTree*) gDirectory->Get("pionana/beamana");
  if(!tree){
    cout<<"no tree!"<<endl;
    gDirectory->ls();
    exit(1);
  }

  TTree * tout = new TTree("tree", tag); lout->Add(tout);
  double dalphat, dphit, dpt, pn, iniPimomentum, iniPitheta, finPimomentum, finPitheta, finProtonmomentum, finProtontheta, fin2Pmom, maxgammaEnergy;
  int nproton, nneutron, nPiZero, ngamma;
  bool kSignal;
  double recProtonmomentum, recPiPlusmomentum;
  int ndEdxCls;
  double lastE0, lastE1, lastE2, lastE3, lastE4, lastE5;
  double startE0, startE1, startE2, startE3, startE4, startE5, chi2, ndof;
  tout->Branch("dalphat",&dalphat);
  tout->Branch("dphit",&dphit);
  tout->Branch("dpt",&dpt);
  tout->Branch("pn",&pn);
  tout->Branch("iniPimomentum",&iniPimomentum);
  tout->Branch("iniPitheta",&iniPitheta);
  tout->Branch("finPimomentum",&finPimomentum);
  tout->Branch("finPitheta",&finPitheta);
  tout->Branch("finProtonmomentum",&finProtonmomentum);
  tout->Branch("recProtonmomentum",&recProtonmomentum);
  tout->Branch("recPiPlusmomentum",&recPiPlusmomentum);
  tout->Branch("finProtontheta",&finProtontheta);
  tout->Branch("fin2Pmom",&fin2Pmom);
  tout->Branch("nproton",&nproton);
  tout->Branch("nneutron",&nneutron);
  tout->Branch("nPiZero",&nPiZero);
  tout->Branch("ngamma",&ngamma);
  tout->Branch("kSignal",&kSignal);
  tout->Branch("ndEdxCls",&ndEdxCls);
  tout->Branch("maxgammaEnergy",&maxgammaEnergy);
  tout->Branch("lastE0",&lastE0);
  tout->Branch("lastE1",&lastE1);
  tout->Branch("lastE2",&lastE2);
  tout->Branch("lastE3",&lastE3);
  tout->Branch("lastE4",&lastE4);
  tout->Branch("lastE5",&lastE5);
  tout->Branch("startE0",&startE0);
  tout->Branch("startE1",&startE1);
  tout->Branch("startE2",&startE2);
  tout->Branch("startE3",&startE3);
  tout->Branch("startE4",&startE4);
  tout->Branch("startE5",&startE5);
  tout->Branch("chi2",&chi2);
  tout->Branch("ndof",&ndof);
  //===============================================================
  // variable list
  //===============================================================
  /*
xlu@ubuntu:~/oldHome/workspace/.analysis/protoDUNE_anaPion$ pwd
/home/xlu/oldHome/workspace/.analysis/protoDUNE_anaPion
xlu@ubuntu:~/oldHome/workspace/.analysis/protoDUNE_anaPion$ cat  seestructure.log | sort | grep -v reco | grep true    | grep grand_daughter  -v | grep beam_start -v | grep daughter_end -v | grep beam_slice -v | grep elastic -v | grep beam_Pi0 -v
 new_true_beam_incidentEnergies = (vector<double>*)0x4008030
 new_true_beam_interactingEnergy = -999
 true_beam_daughter_ID = (vector<int>*)0x3b7f5d0
 true_beam_daughter_len = (vector<double>*)0x3b876b0
 true_beam_daughter_nHits = (vector<int>*)0x3bf0210
 true_beam_daughter_PDG = (vector<int>*)0x3b774f0
 true_beam_daughter_Process = (vector<string>*)0x3be0050
 true_beam_daughter_startP = (vector<double>*)0x3bbfcd0
 true_beam_daughter_startPx = (vector<double>*)0x3ba7a30
 true_beam_daughter_startPy = (vector<double>*)0x3bafb10
 true_beam_daughter_startPz = (vector<double>*)0x3bb7bf0
 true_beam_daughter_startX = (vector<double>*)0x3b8f790
 true_beam_daughter_startY = (vector<double>*)0x3b97870
 true_beam_daughter_startZ = (vector<double>*)0x3b9f950
 true_beam_endP  = 0
 true_beam_endProcess = FastScintillation
 true_beam_endPx = 0
 true_beam_endPy = 0
 true_beam_endPz = -0
 true_beam_endX  = -26.7375
 true_beam_endY  = 422.725
 true_beam_endZ  = 8.48351
 true_beam_ID    = 147
 true_beam_IDE_totalDep = 0
 true_beam_incidentEnergies = (vector<double>*)0x3fd7c90
 true_beam_interactingEnergy = -999
 true_beam_nElasticScatters = 0
 true_beam_nHits = 1271
 true_beam_PDG   = -11
 true_beam_process_dSlice = (vector<int>*)0x3dd21c0
 true_beam_processes = (vector<string>*)0x3dc1fb0
 true_beam_process_matched = (vector<int>*)0x3dda2a0
 true_beam_process_slice = (vector<int>*)0x3dca0e0
 true_daughter_nNeutron = 0
 true_daughter_nNucleus = 0
 true_daughter_nPi0 = 0
 true_daughter_nPiMinus = 0
 true_daughter_nPiPlus = 0
 true_daughter_nProton = 0
xlu@ubuntu:~/oldHome/workspace/.analysis/protoDUNE_anaPion$ 
   */

  //to find the rec-truth 
  //true_beam_daughter_ID = (vector<int>*)0x3b7f5d0
  //reco_daughter_PFP_true_byHits_ID = (vector<int>*)0x3798a90
  //reco_daughter_allTrack_ID = (vector<int>*)0x1548bf0
  //reco_daughter_allTrack_momByRange_proton = (vector<double>*)0x3f19020
  /*
true_beam_daughter_ID[0] == 99
for (size_t i = 0; i < reco_daughter_PFP_true_byHits_ID.size(); ++i) {
  if (reco_daughter_PFP_true_byHits_ID[i] == 99) {
    if (reco_daughter_allTrack_ID[i] != -1) {
       reco_momentum = reco_daughter_allTrack_momByRange_proton[i]
    }
  }
}

// reco_daughter_allTrack_calibrated_dEdX_SCE = (vector<vector<double> >*)
   */
  /*
reco_daughter_allTrack_Chi2_proton = (vector<double>*)0x1180370
 reco_daughter_allTrack_Chi2_ndof = (vector<int>*)0x15413a0
   */
  //===============================================================

  //true_beam_daughter_PDG->size()
  //(*true_beam_daughter_PDG)[0]
  vector<int> *true_beam_daughter_PDG=0x0;
  vector<double> *true_beam_daughter_startPx=0x0;
  vector<double> *true_beam_daughter_startPy=0x0;
  vector<double> *true_beam_daughter_startPz=0x0;

  vector<int>* true_beam_daughter_ID = 0x0;
  vector<int>* reco_daughter_PFP_true_byHits_ID = 0x0;
  vector<int>* reco_daughter_allTrack_ID = 0x0;
  vector<double>* reco_daughter_allTrack_momByRange_proton = 0x0;
  vector<double>* reco_daughter_allTrack_momByRange_muon = 0x0;
  vector<vector<double> >* reco_daughter_allTrack_calibrated_dEdX_SCE = 0x0;

  vector<double>* reco_daughter_allTrack_Chi2_proton = 0x0;
  vector<int>* reco_daughter_allTrack_Chi2_ndof = 0x0;

  double true_beam_endPx = -999;
  double true_beam_endPy = -999;
  double true_beam_endPz = -999;

  int true_beam_PDG   = -999;

  int true_daughter_nNeutron = -999;
  int true_daughter_nNucleus = -999;
  int true_daughter_nPi0 = -999;
  int true_daughter_nPiMinus = -999;
  int true_daughter_nPiPlus = -999;
  int true_daughter_nProton = -999;

  tree->SetBranchAddress("true_beam_daughter_PDG", &true_beam_daughter_PDG);
  tree->SetBranchAddress("true_beam_daughter_startPx", &true_beam_daughter_startPx);
  tree->SetBranchAddress("true_beam_daughter_startPy", &true_beam_daughter_startPy);
  tree->SetBranchAddress("true_beam_daughter_startPz", &true_beam_daughter_startPz);

  tree->SetBranchAddress("true_beam_daughter_ID", &true_beam_daughter_ID);
  tree->SetBranchAddress("reco_daughter_PFP_true_byHits_ID", &reco_daughter_PFP_true_byHits_ID);
  tree->SetBranchAddress("reco_daughter_allTrack_ID", &reco_daughter_allTrack_ID);
  tree->SetBranchAddress("reco_daughter_allTrack_momByRange_proton", &reco_daughter_allTrack_momByRange_proton);
  tree->SetBranchAddress("reco_daughter_allTrack_momByRange_muon", &reco_daughter_allTrack_momByRange_muon);
  tree->SetBranchAddress("reco_daughter_allTrack_calibrated_dEdX_SCE", &reco_daughter_allTrack_calibrated_dEdX_SCE);
  tree->SetBranchAddress("reco_daughter_allTrack_Chi2_proton", &reco_daughter_allTrack_Chi2_proton);
  tree->SetBranchAddress("reco_daughter_allTrack_Chi2_ndof", &reco_daughter_allTrack_Chi2_ndof);

  tree->SetBranchAddress("true_beam_endPx", &true_beam_endPx);
  tree->SetBranchAddress("true_beam_endPy", &true_beam_endPy);
  tree->SetBranchAddress("true_beam_endPz", &true_beam_endPz);

  tree->SetBranchAddress("true_beam_PDG", &true_beam_PDG);

  tree->SetBranchAddress("true_daughter_nNeutron", &true_daughter_nNeutron);
  tree->SetBranchAddress("true_daughter_nNucleus", &true_daughter_nNucleus);
  tree->SetBranchAddress("true_daughter_nPi0", &true_daughter_nPi0);
  tree->SetBranchAddress("true_daughter_nPiMinus", &true_daughter_nPiMinus);
  tree->SetBranchAddress("true_daughter_nPiPlus", &true_daughter_nPiPlus);
  tree->SetBranchAddress("true_daughter_nProton", &true_daughter_nProton);

  //==============================================================================================================
  //Histogram definition
  //==============================================================================================================
  TH1I * hbeamType = new TH1I("h0beamType"+tag,  "", 20, -0.5, 19.5); lout->Add(hbeamType);
  TH1I * hndaughter = new TH1I("h1ndaughter"+tag,"", 50, -0.5, 49.5); lout->Add(hndaughter);
  TH1I * hnexcl = new TH1I("h2nexcl"+tag,"",2, -0.5, 1.5); lout->Add(hnexcl);
  //TH2I * hksignp = new TH2I("h2ksignp"+tag,"",10, -0.5, 9.5, 2, -0.5, 1.5); lout->Add(hksignp);
  TH2I * hksignp = new TH2I("h2ksignp"+tag,"",2, -0.5, 1.5, 2, -0.5, 1.5); lout->Add(hksignp);
  TH1I * hselectedDaughterType = new TH1I("h3selectedDaughterType"+tag, "", 21, -0.5, 20.5); lout->Add(hselectedDaughterType);
  TH1I * hnproton = new TH1I("h4nproton"+tag,"",11, -0.5, 10.5); lout->Add(hnproton);
  TH1I * hnneutron = new TH1I("h5nneutron"+tag,"",11, -0.5, 10.5); lout->Add(hnneutron);
  TH1I * hnPiZero = new TH1I("h6nPiZero"+tag,"",11, -0.5, 10.5); lout->Add(hnPiZero);
  TH1D * hmomIniPi = new TH1D("hmomIniPi"+tag,"", 50, 0, 2); lout->Add(hmomIniPi);
  TH1D * hmomFinPi = new TH1D("hmomFinPi"+tag,"", 50, 0, 2); lout->Add(hmomFinPi);
  TH1D * hmomFinProton = new TH1D("hmomFinProton"+tag,"", 50, 0, 2); lout->Add(hmomFinProton);

  const double Ebin[]= {0.000000, 20.000000, 40.000000, 60.000000, 80.000000, 100.000000, 120.000000, 130.000000, 140.000000, 150.000000, 160.000000, 170.000000, 180.000000};
  TH1D *hdalphat = new TH1D("dalphat","", sizeof(Ebin)/sizeof(double)-1, Ebin); lout->Add(hdalphat);

  const double Fbin[]={0.000000, 2.500000, 5.000000, 7.500000, 10.000000, 12.500000, 15.000000, 17.500000, 20.000000, 22.500000, 25.000000, 27.500000, 30.000000, 35.000000, 40.000000, 45.000000, 50.000000, 55.000000, 60.000000, 70.000000, 85.000000, 105.000000, 130.000000, 180.000000};
  TH1D *hdphit = new TH1D("dphit","", sizeof(Fbin)/sizeof(double)-1, Fbin); lout->Add(hdphit);

   const double Gbin[]={0.000000, 0.025000, 0.050000, 0.075000, 0.100000, 0.125000, 0.150000, 0.175000, 0.200000, 0.225000, 0.250000, 0.275000, 0.300000, 0.350000, 0.400000, 0.450000, 0.500000, 0.550000, 0.600000, 0.650000, 0.700000, 0.800000, 1.000000, 1.200000, 2.000000};
   TH1D *hdpt = new TH1D("dpt","", sizeof(Gbin)/sizeof(double)-1, Gbin); lout->Add(hdpt);

   const double Hbin[]={0.000000, 0.025000, 0.050000, 0.075000, 0.100000, 0.125000, 0.150000, 0.175000, 0.200000, 0.225000, 0.250000, 0.275000, 0.300000, 0.350000, 0.400000, 0.450000, 0.500000, 0.550000, 0.600000, 0.650000, 0.700000, 0.800000, 1.000000, 1.200000, 2.000000};
   TH1D *hpn = new TH1D("pn","", sizeof(Hbin)/sizeof(double)-1, Hbin); lout->Add(hpn);

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
    vector<TLorentzVector> vecPiP = getPiP1P2(kPiZero, true_beam_daughter_PDG, true_beam_daughter_startPx, true_beam_daughter_startPy, true_beam_daughter_startPz, hselectedDaughterType, nproton, nneutron, nPiZero, ngamma, maxgammaEnergy, protonIdx, piplusIdx, kSignal);
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
      (*recmomentum) = getRecMomentum(targetIdx, reco_daughter_PFP_true_byHits_ID, reco_daughter_allTrack_ID, reco_daughter_allTrack_momByRange, reco_daughter_allTrack_calibrated_dEdX_SCE, lastE, startE, reco_daughter_allTrack_Chi2_proton, reco_daughter_allTrack_Chi2_ndof, chi2, ndof);

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

  drawHist(lout, tag);

  return tag;
}

int main(int argc, char * argv[])
{
  if(argc!=3){
    cout<<"argc!=3 !"<<argc<<endl;
    return 0;
  }

  gSystem->Load("libTree");

  style::SetGlobalStyle();

  const bool kPiZero = atoi(argv[1]);
  const bool kProton = atoi(argv[2]);
  if(kPiZero&&!kProton){
    return 0;
  }

  TString tag = (kPiZero?"MPiZero":"1PiPlus");
  tag+="_Tracking";
  tag+=(kProton?"Proton":"PiPlus");

  TList * lout = new TList;

  TFile * fin = new TFile("input/protoDUNE_reco_flattree.root");
  if(!fin->IsOpen()){
    cout<<"fin not open!"<<endl;
    exit(1);
  }

  anaData(lout, tag);
 
  TFile * fout = new TFile(Form("output/outanaData_%s.root", tag.Data()),"recreate");

  lout->Write();
  fout->Save();
  fout->Close();

  return 0;
}
