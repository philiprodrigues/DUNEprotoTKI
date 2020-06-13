#ifndef _ANAIO_H_
#define _ANAIO_H_

namespace AnaIO
{
  //======================================= tree out =======================================
  double dalphat; 
  double  dphit;
  double  dpt; 
  double  pn; 
  double  iniPimomentum; 
  double  iniPitheta; 
  double  finPimomentum; 
  double  finPitheta; 
  double  finProtonmomentum; 
  double  finProtontheta; 
  double  fin2Pmom; 
  double  maxgammaEnergy;
  int nproton; 
  int nneutron; 
  int nPiZero; 
  int ngamma;
  bool kSignal;
  double recProtonmomentum;
  double recPiPlusmomentum;
  int ndEdxCls;
  double lastE0; 
  double lastE1; 
  double lastE2; 
  double lastE3; 
  double lastE4; 
  double lastE5;
  double startE0; 
  double startE1; 
  double startE2; 
  double startE3; 
  double startE4; 
  double startE5; 
  double startTruncatedMeanE10;
  double startTruncatedMeanE20;
  double startTruncatedMeanE30;
  double startTruncatedMeanE40;
  double startTruncatedMeanE50;
  double lastTruncatedMeanE;
  double chi2;
  double ndof;

  //--- beam related
  int nBeamdEdxCls;
  double beamLastE0; 
  double beamLastE1; 
  double beamLastE2; 
  double beamLastE3; 
  double beamLastE4; 
  double beamLastE5;
  double beamStartE0; 
  double beamStartE1; 
  double beamStartE2; 
  double beamStartE3; 
  double beamStartE4; 
  double beamStartE5; 
  double beamTMeanLast;
  double beamTMeanStart;

  //======================================= tree in =======================================
  //--------------------------------------- anaTruth ---------------------------------------
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

  vector<int>     *reco_daughter_PFP_ID =0x0;
  vector<int>     *reco_daughter_PFP_nHits =0x0;
  vector<double>  *reco_daughter_PFP_trackScore =0x0;
  vector<double>  *reco_daughter_PFP_emScore =0x0;
  vector<double>  *reco_daughter_PFP_michelScore =0x0;

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

  //--------------------------------------- anaRec ---------------------------------------
  //--- established
  int reco_beam_type = -999;

  Double_t        reco_beam_startX;
  Double_t        reco_beam_startY;
  Double_t        reco_beam_startZ;
  Double_t        reco_beam_trackDirX;
  Double_t        reco_beam_trackDirY;
  Double_t        reco_beam_trackDirZ;
  Double_t        true_beam_startX;
  Double_t        true_beam_startY;
  Double_t        true_beam_startZ;
  Double_t        true_beam_startDirX;
  Double_t        true_beam_startDirY;
  Double_t        true_beam_startDirZ;

  Double_t        reco_beam_endZ;

  vector<double>  *reco_beam_calibrated_dEdX;

  //--- to move above
  Int_t           nTrack; //reco_beam_nTrackDaughters;//not fill in flat tree! Neither is nShowerDaughters

  //--- to test
  Double_t        reco_beam_len;


  //======================================= Truth Hist out =======================================
 TH1I * hbeamType = 0x0;
 TH1I * hndaughter = 0x0;
 TH1I * hnexcl = 0x0;
 TH2I * hksignp = 0x0;
 TH1I * hselectedDaughterType = 0x0;
 TH1I * hnproton = 0x0;
 TH1I * hnneutron = 0x0;
 TH1I * hnPiZero = 0x0;
 TH1D * hmomIniPi = 0x0;
 TH1D * hmomFinPi = 0x0;
 TH1D * hmomFinProton = 0x0;
 TH1D *hdalphat = 0x0;
 TH1D *hdphit = 0x0;
 TH1D *hdpt = 0x0;
 TH1D *hpn = 0x0;

//======================================= Rec Hist out =======================================
 TH1I * hTruthBeamType = 0x0;
 TH1I * hTruthSignal = 0x0;

 //--- established
 TH1I * hnBeamdEdxCls = 0x0;

 TH2D * hSignalVsTMeanStart = 0x0;
 TH2D * hSignalVsTMeanLast = 0x0;

 TH2D * hSignalVsStartE0 = 0x0;
 TH2D * hSignalVsStartE1 = 0x0;
 TH2D * hSignalVsStartE2 = 0x0;
 TH2D * hSignalVsStartE3 = 0x0;
 TH2D * hSignalVsStartE4 = 0x0;
 TH2D * hSignalVsStartE5 = 0x0;
 TH2D * hSignalVsLastE0 = 0x0;
 TH2D * hSignalVsLastE1 = 0x0;
 TH2D * hSignalVsLastE2 = 0x0;
 TH2D * hSignalVsLastE3 = 0x0;
 TH2D * hSignalVsLastE4 = 0x0;
 TH2D * hSignalVsLastE5 = 0x0;

 TH2D * hSigAfterVsStartE0 = 0x0;
 TH2D * hSigAfterVsStartE1 = 0x0;
 TH2D * hSigAfterVsStartE2 = 0x0;
 TH2D * hSigAfterVsStartE3 = 0x0;
 TH2D * hSigAfterVsStartE4 = 0x0;
 TH2D * hSigAfterVsStartE5 = 0x0;
 TH2D * hSigAfterVsLastE0 = 0x0;
 TH2D * hSigAfterVsLastE1 = 0x0;
 TH2D * hSigAfterVsLastE2 = 0x0;
 TH2D * hSigAfterVsLastE3 = 0x0;
 TH2D * hSigAfterVsLastE4 = 0x0;
 TH2D * hSigAfterVsLastE5 = 0x0;

 //--- to move above
 TH1I * hRecoBeamType = 0x0;
 TH1I * hBeamPosPass = 0x0;
 TH1D * hBeamEndZ = 0x0;
 TH1I * hBeamEndZPass = 0x0;

 TH2I * hBeamNTrack = 0x0;
 TH2D * hSignalVsBeamNTrack = 0x0;

 //--- to test
 TH1D * hBeamLen = 0x0;
 TH2D * hSignalVsLen = 0x0;
 TH2D * hSigAfterVsLen = 0x0;

 //--- resolution
 TH2D * hProtonMomentumRes = 0x0;
 TH2D * hPiMomentumRes = 0x0;

 TH2D * hCutnHits = 0x0;
 TH2D * hCutNdEdx = 0x0;
 TH2D * hCutstartE2 = 0x0;
 TH2D * hCutstartE3 = 0x0;
 TH2D * hCutlastE2 = 0x0;
 TH2D * hCutlastE3 = 0x0;
 TH2D * hCuttrackScore = 0x0;
 TH2D * hCutemScore = 0x0;
 TH2D * hCutmichelScore = 0x0;
 TH2D * hCutChi2NDF = 0x0;

//==============================================================================
//==============================================================================

TTree * GetOutputTree(TList * lout, const TString tag)
{
  TTree * tout = new TTree("tree", tag); lout->Add(tout);
 
  //definitely need this to avoid memory-resident Tree
  tout->SetDirectory(0);

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
  tout->Branch("startTruncatedMeanE10",&startTruncatedMeanE10);
  tout->Branch("startTruncatedMeanE20",&startTruncatedMeanE20);
  tout->Branch("startTruncatedMeanE30",&startTruncatedMeanE30);
  tout->Branch("startTruncatedMeanE40",&startTruncatedMeanE40);
  tout->Branch("startTruncatedMeanE50",&startTruncatedMeanE50);
  tout->Branch("lastTruncatedMeanE",&lastTruncatedMeanE);
  tout->Branch("chi2",&chi2);
  tout->Branch("ndof",&ndof);

  tout->Branch("nBeamdEdxCls",&nBeamdEdxCls);
  tout->Branch("beamLastE0",&beamLastE0);
  tout->Branch("beamLastE1",&beamLastE1);
  tout->Branch("beamLastE2",&beamLastE2);
  tout->Branch("beamLastE3",&beamLastE3);
  tout->Branch("beamLastE4",&beamLastE4);
  tout->Branch("beamLastE5",&beamLastE5);
  tout->Branch("beamStartE0",&beamStartE0);
  tout->Branch("beamStartE1",&beamStartE1);
  tout->Branch("beamStartE2",&beamStartE2);
  tout->Branch("beamStartE3",&beamStartE3);
  tout->Branch("beamStartE4",&beamStartE4);
  tout->Branch("beamStartE5",&beamStartE5);

  tout->Branch("beamTMeanStart",&beamTMeanStart);
  tout->Branch("beamTMeanLast",&beamTMeanLast);

  //save input info
  tout->Branch("true_beam_PDG",&true_beam_PDG);
  tout->Branch("reco_beam_len",&reco_beam_len);
  tout->Branch("nTrack",&nTrack);

  return tout;
}

TTree * GetInputTree(const TString tname)
{
  TTree * tree=(TTree*) gDirectory->Get(tname);
  if(!tree){
    cout<<"no tree!"<<endl;
    gDirectory->ls();
    exit(1);
  }

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

  //--------------------------------------- anaTruth ---------------------------------------
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

  tree->SetBranchAddress("reco_daughter_PFP_ID", &reco_daughter_PFP_ID);
  tree->SetBranchAddress("reco_daughter_PFP_nHits", &reco_daughter_PFP_nHits);
  tree->SetBranchAddress("reco_daughter_PFP_trackScore", &reco_daughter_PFP_trackScore);
  tree->SetBranchAddress("reco_daughter_PFP_emScore", &reco_daughter_PFP_emScore);
  tree->SetBranchAddress("reco_daughter_PFP_michelScore", &reco_daughter_PFP_michelScore);

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

  //--------------------------------------- anaRec ---------------------------------------
  tree->SetBranchAddress("reco_beam_type", &reco_beam_type);

  tree->SetBranchAddress("reco_beam_startX", &reco_beam_startX);
  tree->SetBranchAddress("reco_beam_trackDirX", &reco_beam_trackDirX);
  tree->SetBranchAddress("true_beam_startX", &true_beam_startX);
  tree->SetBranchAddress("true_beam_startDirX", &true_beam_startDirX);

  tree->SetBranchAddress("reco_beam_startY", &reco_beam_startY);
  tree->SetBranchAddress("reco_beam_trackDirY", &reco_beam_trackDirY);
  tree->SetBranchAddress("true_beam_startY", &true_beam_startY);
  tree->SetBranchAddress("true_beam_startDirY", &true_beam_startDirY);

  tree->SetBranchAddress("reco_beam_startZ", &reco_beam_startZ);
  tree->SetBranchAddress("reco_beam_trackDirZ", &reco_beam_trackDirZ);
  tree->SetBranchAddress("true_beam_startZ", &true_beam_startZ);
  tree->SetBranchAddress("true_beam_startDirZ", &true_beam_startDirZ);

  tree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ);

  //to check the difference with _allTrack_
  tree->SetBranchAddress("reco_beam_len", &reco_beam_len);
  //tree->SetBranchAddress("reco_beam_allTrack_len", &reco_beam_len);
  tree->SetBranchAddress("reco_beam_calibrated_dEdX", &reco_beam_calibrated_dEdX);
  //tree->SetBranchAddress("reco_beam_allTrack_calibrated_dEdX", &reco_beam_calibrated_dEdX);

  //tree->SetBranchAddress("reco_beam_nTrackDaughters", &reco_beam_nTrackDaughters);

  return tree;
}

void IniRecHist(TList * lout, const TString tag)
{
  const int nPass = 2;
  const double Passmin = -0.5;
  const double Passmax = 1.5;

  hTruthBeamType = new TH1I("a0TruthBeamType"+tag,  "", 20, -0.5, 19.5); lout->Add(hTruthBeamType);
  hTruthSignal = new TH1I("a1TruthSignal"+tag,  "",  nPass, Passmin, Passmax); lout->Add(hTruthSignal);

  //--- established
  //over shadowed by beam pre-filtering
  /*
  hnBeamdEdxCls = new TH1I("c1nBeamdEdxCls"+tag,"", 20, -0.5, 19.5); lout->Add(hnBeamdEdxCls);

  const int nE = 30;
  const double Emin = 0;
  const double lastEmax = 30;
  const double startEmax = 10;

  hSignalVsTMeanStart = new TH2D("c1SignalVsTMeanStartNOHPRFPUR"+tag,"",nE*5, Emin, startEmax, nPass, Passmin, Passmax); lout->Add(hSignalVsTMeanStart);

  hSignalVsStartE0 = new TH2D("z1SignalVsStartE0NOHPRFPUR"+tag,"", nE, Emin, startEmax, nPass, Passmin, Passmax); lout->Add(hSignalVsStartE0);
  hSignalVsStartE1 = new TH2D("z1SignalVsStartE1NOHPRFPUR"+tag,"", nE, Emin, startEmax, nPass, Passmin, Passmax); lout->Add(hSignalVsStartE1);
  hSignalVsStartE2 = new TH2D("z1SignalVsStartE2NOHPRFPUR"+tag,"", nE, Emin, startEmax, nPass, Passmin, Passmax); lout->Add(hSignalVsStartE2);
  hSignalVsStartE3 = new TH2D("z1SignalVsStartE3NOHPRFPUR"+tag,"", nE, Emin, startEmax, nPass, Passmin, Passmax); lout->Add(hSignalVsStartE3);
  hSignalVsStartE4 = new TH2D("z1SignalVsStartE4NOHPRFPUR"+tag,"", nE, Emin, startEmax, nPass, Passmin, Passmax); lout->Add(hSignalVsStartE4);
  hSignalVsStartE5 = new TH2D("z1SignalVsStartE5NOHPRFPUR"+tag,"", nE, Emin, startEmax, nPass, Passmin, Passmax); lout->Add(hSignalVsStartE5);
  hSignalVsLastE0 = new TH2D("z1SignalVsLastE0NOHPRFPUR"+tag,"", nE, Emin, lastEmax, nPass, Passmin, Passmax); lout->Add(hSignalVsLastE0);
  hSignalVsLastE1 = new TH2D("z1SignalVsLastE1NOHPRFPUR"+tag,"", nE, Emin, lastEmax, nPass, Passmin, Passmax); lout->Add(hSignalVsLastE1);
  hSignalVsLastE2 = new TH2D("z1SignalVsLastE2NOHPRFPUR"+tag,"", nE, Emin, lastEmax, nPass, Passmin, Passmax); lout->Add(hSignalVsLastE2);
  hSignalVsLastE3 = new TH2D("z1SignalVsLastE3NOHPRFPUR"+tag,"", nE, Emin, lastEmax, nPass, Passmin, Passmax); lout->Add(hSignalVsLastE3);
  hSignalVsLastE4 = new TH2D("z1SignalVsLastE4NOHPRFPUR"+tag,"", nE, Emin, lastEmax, nPass, Passmin, Passmax); lout->Add(hSignalVsLastE4);
  hSignalVsLastE5 = new TH2D("z1SignalVsLastE5NOHPRFPUR"+tag,"", nE, Emin, lastEmax, nPass, Passmin, Passmax); lout->Add(hSignalVsLastE5);

  hSignalVsTMeanLast = new TH2D("z2SignalVsTMeanLastNOHPRFPUR"+tag,"",nE*5, Emin, lastEmax, nPass, Passmin, Passmax); lout->Add(hSignalVsTMeanLast);
  hSigAfterVsStartE0 = new TH2D("z2SigAfterVsStartE0NOHPRFPUR"+tag,"", nE, Emin, startEmax, nPass, Passmin, Passmax); lout->Add(hSigAfterVsStartE0);
  hSigAfterVsStartE1 = new TH2D("z2SigAfterVsStartE1NOHPRFPUR"+tag,"", nE, Emin, startEmax, nPass, Passmin, Passmax); lout->Add(hSigAfterVsStartE1);
  hSigAfterVsStartE2 = new TH2D("z2SigAfterVsStartE2NOHPRFPUR"+tag,"", nE, Emin, startEmax, nPass, Passmin, Passmax); lout->Add(hSigAfterVsStartE2);
  hSigAfterVsStartE3 = new TH2D("z2SigAfterVsStartE3NOHPRFPUR"+tag,"", nE, Emin, startEmax, nPass, Passmin, Passmax); lout->Add(hSigAfterVsStartE3);
  hSigAfterVsStartE4 = new TH2D("z2SigAfterVsStartE4NOHPRFPUR"+tag,"", nE, Emin, startEmax, nPass, Passmin, Passmax); lout->Add(hSigAfterVsStartE4);
  hSigAfterVsStartE5 = new TH2D("z2SigAfterVsStartE5NOHPRFPUR"+tag,"", nE, Emin, startEmax, nPass, Passmin, Passmax); lout->Add(hSigAfterVsStartE5);
  hSigAfterVsLastE0 = new TH2D("z2SigAfterVsLastE0NOHPRFPUR"+tag,"", nE, Emin, lastEmax, nPass, Passmin, Passmax); lout->Add(hSigAfterVsLastE0);
  hSigAfterVsLastE1 = new TH2D("z2SigAfterVsLastE1NOHPRFPUR"+tag,"", nE, Emin, lastEmax, nPass, Passmin, Passmax); lout->Add(hSigAfterVsLastE1);
  hSigAfterVsLastE2 = new TH2D("z2SigAfterVsLastE2NOHPRFPUR"+tag,"", nE, Emin, lastEmax, nPass, Passmin, Passmax); lout->Add(hSigAfterVsLastE2);
  hSigAfterVsLastE3 = new TH2D("z2SigAfterVsLastE3NOHPRFPUR"+tag,"", nE, Emin, lastEmax, nPass, Passmin, Passmax); lout->Add(hSigAfterVsLastE3);
  hSigAfterVsLastE4 = new TH2D("z2SigAfterVsLastE4NOHPRFPUR"+tag,"", nE, Emin, lastEmax, nPass, Passmin, Passmax); lout->Add(hSigAfterVsLastE4);
  hSigAfterVsLastE5 = new TH2D("z2SigAfterVsLastE5NOHPRFPUR"+tag,"", nE, Emin, lastEmax, nPass, Passmin, Passmax); lout->Add(hSigAfterVsLastE5);
  */
 
  //--- to move above
  hRecoBeamType = new TH1I("h0RecoBeamType"+tag,  "", 30, -4.5, 25.5); lout->Add(hRecoBeamType);

  hBeamPosPass = new TH1I("h1BeamPosPass"+tag, "", 4, -0.5, 3.5); lout->Add(hBeamPosPass);

  hBeamEndZ = new TH1D("h2BeamEndZ"+tag,"",50, 0, 500); lout->Add(hBeamEndZ);
  hBeamEndZPass = new TH1I("h2BeamEndZPass"+tag,"",4, -0.5, 3.5); lout->Add(hBeamEndZPass);

  const int nBeamNTrack = 10;
  const double BeamNTrackmin = -0.5;
  const double BeamNTrackmax = 9.5;
  hBeamNTrack = new TH2I("h6BeamNTrack"+tag,"", nBeamNTrack, BeamNTrackmin, BeamNTrackmax, 5, -0.5, 4.5); lout->Add(hBeamNTrack);
  hSignalVsBeamNTrack = new TH2D("p6SignalVsBeanNTrackNOHPRFPUR"+tag,"", nBeamNTrack, BeamNTrackmin, BeamNTrackmax, nPass, Passmin, Passmax); lout->Add(hSignalVsBeamNTrack);

  //--- to test
  //sensitive to pi+ beam, but not to signal -> sensitive to pi+ beam non-signal, i.e. non-interacting
  /*
  const int nlen = 10;
  const double lenmin = 0;
  const double lenmax = 250;
  hBeamLen = new TH1D("h9BeamLen"+tag,"",nlen, lenmin, lenmax); lout->Add(hBeamLen);
  hSignalVsLen = new TH2D("p9SignalVsLenNOHPRFPUR"+tag,"", nlen, lenmin, lenmax, nPass, Passmin, Passmax); lout->Add(hSignalVsLen);
  hSigAfterVsLen = new TH2D("p9SigAfterVsLenNOHPRFPUR"+tag,"", nlen, lenmin, lenmax, nPass, Passmin, Passmax); lout->Add(hSigAfterVsLen);
  */

  //--- resolution
  hProtonMomentumRes = new TH2D("b0ProtonMomentumResNOH"+tag,"", 20, 0, 2, 20, -1, 1); lout->Add(hProtonMomentumRes);
  hPiMomentumRes = new TH2D("b1PiMomentumResNOH"+tag,"", 20, 0, 2, 20, -1, 1); lout->Add(hPiMomentumRes);

  //--- cut variables
  hCutNdEdx       = new TH2D("c001CutNdEdxSTK"+tag,"", 20, 0, 100, 10, -0.5, 9.5); lout->Add(hCutNdEdx);
  hCuttrackScore  = new TH2D("c002CuttrackScoreSTK"+tag,"", 20, 0, 1, 10, -0.5, 9.5); lout->Add(hCuttrackScore);
  hCutemScore     = new TH2D("c003CutemScoreSTK"+tag,"", 20, 0, 1, 10, -0.5, 9.5); lout->Add(hCutemScore);
  hCutmichelScore = new TH2D("c004CutmichelScoreSTK"+tag,"", 20, 0, 1, 10, -0.5, 9.5); lout->Add(hCutmichelScore);
  hCutnHits       = new TH2D("c005CutnHitsSTK"+tag,"", 20, 0, 400, 10, -0.5, 9.5); lout->Add(hCutnHits);
  hCutChi2NDF     = new TH2D("c006CutChi2NDFSTK"+tag,"", 20, 0, 100, 10, -0.5, 9.5); lout->Add(hCutChi2NDF);
  hCutstartE2     = new TH2D("c007CutstartE2STK"+tag,"", 20, 0, 20, 10, -0.5, 9.5); lout->Add(hCutstartE2);
  hCutstartE3     = new TH2D("c008CutstartE3STK"+tag,"", 20, 0, 20, 10, -0.5, 9.5); lout->Add(hCutstartE3);
  hCutlastE2      = new TH2D("c009CutlastE2STK"+tag,"", 20, 0, 20, 10, -0.5, 9.5); lout->Add(hCutlastE2);
  hCutlastE3      = new TH2D("c010CutlastE3STK"+tag,"", 20, 0, 20, 10, -0.5, 9.5); lout->Add(hCutlastE3);

}

void IniTruthHist(TList * lout, const TString tag)
{
  hbeamType = new TH1I("h0beamType"+tag,  "", 20, -0.5, 19.5); lout->Add(hbeamType);
  hndaughter = new TH1I("h1ndaughter"+tag,"", 50, -0.5, 49.5); lout->Add(hndaughter);
  hnexcl = new TH1I("h2nexcl"+tag,"",2, -0.5, 1.5); lout->Add(hnexcl);
  hksignp = new TH2I("h2ksignp"+tag,"",2, -0.5, 1.5, 2, -0.5, 1.5); lout->Add(hksignp);
  hselectedDaughterType = new TH1I("h3selectedDaughterType"+tag, "", 21, -0.5, 20.5); lout->Add(hselectedDaughterType);
  hnproton = new TH1I("h4nproton"+tag,"",11, -0.5, 10.5); lout->Add(hnproton);
  hnneutron = new TH1I("h5nneutron"+tag,"",11, -0.5, 10.5); lout->Add(hnneutron);
  hnPiZero = new TH1I("h6nPiZero"+tag,"",11, -0.5, 10.5); lout->Add(hnPiZero);
  hmomIniPi = new TH1D("hmomIniPi"+tag,"", 50, 0, 2); lout->Add(hmomIniPi);
  hmomFinPi = new TH1D("hmomFinPi"+tag,"", 50, 0, 2); lout->Add(hmomFinPi);
  hmomFinProton = new TH1D("hmomFinProton"+tag,"", 50, 0, 2); lout->Add(hmomFinProton);

  const double Ebin[]= {0.000000, 20.000000, 40.000000, 60.000000, 80.000000, 100.000000, 120.000000, 130.000000, 140.000000, 150.000000, 160.000000, 170.000000, 180.000000};
  hdalphat = new TH1D("dalphat","", sizeof(Ebin)/sizeof(double)-1, Ebin); lout->Add(hdalphat);

  const double Fbin[]={0.000000, 2.500000, 5.000000, 7.500000, 10.000000, 12.500000, 15.000000, 17.500000, 20.000000, 22.500000, 25.000000, 27.500000, 30.000000, 35.000000, 40.000000, 45.000000, 50.000000, 55.000000, 60.000000, 70.000000, 85.000000, 105.000000, 130.000000, 180.000000};
  hdphit = new TH1D("dphit","", sizeof(Fbin)/sizeof(double)-1, Fbin); lout->Add(hdphit);

   const double Gbin[]={0.000000, 0.025000, 0.050000, 0.075000, 0.100000, 0.125000, 0.150000, 0.175000, 0.200000, 0.225000, 0.250000, 0.275000, 0.300000, 0.350000, 0.400000, 0.450000, 0.500000, 0.550000, 0.600000, 0.650000, 0.700000, 0.800000, 1.000000, 1.200000, 2.000000};
   hdpt = new TH1D("dpt","", sizeof(Gbin)/sizeof(double)-1, Gbin); lout->Add(hdpt);

   const double Hbin[]={0.000000, 0.025000, 0.050000, 0.075000, 0.100000, 0.125000, 0.150000, 0.175000, 0.200000, 0.225000, 0.250000, 0.275000, 0.300000, 0.350000, 0.400000, 0.450000, 0.500000, 0.550000, 0.600000, 0.650000, 0.700000, 0.800000, 1.000000, 1.200000, 2.000000};
   hpn = new TH1D("pn","", sizeof(Hbin)/sizeof(double)-1, Hbin); lout->Add(hpn);
}

}

#endif
