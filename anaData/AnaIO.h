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

  Int_t           nTrack; 

  //not used, keep in case ->
  //--- beam related
  int beamNdEdx;
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
  double beamLastTME;
  double beamStartTME;
  //<--

  //======================================= tree in =======================================
  //--------------------------------------- anaTruth ---------------------------------------
  vector<int>    *true_beam_daughter_PDG=0x0;
  vector<double> *true_beam_daughter_startPx=0x0;
  vector<double> *true_beam_daughter_startPy=0x0;
  vector<double> *true_beam_daughter_startPz=0x0;

  vector<int>* true_beam_daughter_ID = 0x0;
  vector<int>* reco_daughter_PFP_true_byHits_ID = 0x0;
  vector<int>* reco_daughter_allTrack_ID = 0x0;

  vector<vector<double> >* reco_daughter_allTrack_calibrated_dEdX_SCE = 0x0;

  vector<double>* reco_daughter_allTrack_Chi2_proton = 0x0;
  vector<int>*    reco_daughter_allTrack_Chi2_ndof = 0x0;

  double true_beam_endPx = -999;
  double true_beam_endPy = -999;
  double true_beam_endPz = -999;

  int true_beam_PDG   = -999;

  //--------------------------------------- anaRec ---------------------------------------

  vector<int>     *reco_daughter_PFP_ID =0x0;
  vector<int>     *reco_daughter_PFP_nHits =0x0;
  vector<double>  *reco_daughter_PFP_trackScore_collection =0x0;
  vector<double>  *reco_daughter_PFP_emScore_collection =0x0;
  vector<double>  *reco_daughter_PFP_michelScore_collection =0x0;

  vector<int>     *reco_daughter_allShower_ID=0x0;
  vector<double>  *reco_daughter_allShower_dirX=0x0;
  vector<double>  *reco_daughter_allShower_dirY=0x0;
  vector<double>  *reco_daughter_allShower_dirZ=0x0;
  vector<double>  *reco_daughter_allShower_energy=0x0;

  vector<double>* reco_daughter_allTrack_momByRange_proton = 0x0;
  vector<double>* reco_daughter_allTrack_momByRange_muon = 0x0;

  vector<double>* reco_daughter_allTrack_Theta = 0x0;
  vector<double>* reco_daughter_allTrack_Phi = 0x0;

  vector<int>     *data_BI_PDG_candidates=0x0;

  int reco_beam_type = -999;

  Double_t        reco_beam_startX;
  Double_t        reco_beam_startY;
  Double_t        reco_beam_startZ;
  Double_t        reco_beam_trackDirX;
  Double_t        reco_beam_trackDirY;
  Double_t        reco_beam_trackDirZ;
  Double_t        reco_beam_trackEndDirX;
  Double_t        reco_beam_trackEndDirY;
  Double_t        reco_beam_trackEndDirZ;
  Double_t        true_beam_startX;
  Double_t        true_beam_startY;
  Double_t        true_beam_startZ;
  Double_t        true_beam_startDirX;
  Double_t        true_beam_startDirY;
  Double_t        true_beam_startDirZ;

  Double_t        reco_beam_endZ;

  Double_t        data_BI_X;
  Double_t        data_BI_Y;
  Double_t        data_BI_dirX;
  Double_t        data_BI_dirY;
  Double_t        data_BI_dirZ;
  Int_t           data_BI_nMomenta;
  Int_t           data_BI_nTracks;

  //not used, keep in case ->
  vector<double>  *reco_beam_calibrated_dEdX;
  Double_t        reco_beam_len;
  //<-

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

  TH1I * hCutbeamID = 0x0;
  TH2D * hCutBeamType = 0x0;
  TH2D * hCutBeamPosPass = 0x0;
  TH2D * hCutBeamEndZ = 0x0;
  TH2D * hCutBeamEndZPass = 0x0;

  TH2D * hBeamThetaRes = 0x0;

  TH2D * hProtonThetaRes = 0x0;
  TH2D * hPiplusThetaRes = 0x0;

  TH2D * hProtonMomentumRes = 0x0;
  TH2D * hPiplusMomentumRes = 0x0;

  TH2D * hRecBeamTheta = 0x0;

  TH2D * hRecProtonMomentum = 0x0;
  TH2D * hRecProtonTheta = 0x0;

  TH2D * hRecPiplusMomentum = 0x0;
  TH2D * hRecPiplusTheta = 0x0;

  TH2D * hRecProtonStartE2 = 0x0;
  TH2D * hRecProtonStartE3 = 0x0;
  TH2D * hRecPiplusStartE2 = 0x0;
  TH2D * hRecPiplusStartE3 = 0x0;

  TH2D * hRecProtonLastE2 = 0x0;
  TH2D * hRecProtonLastE3 = 0x0;
  TH2D * hRecPiplusLastE2 = 0x0;
  TH2D * hRecPiplusLastE3 = 0x0;

  TH2D * hRecProtonStartTME = 0x0;
  TH2D * hRecProtonLastTME = 0x0;
  TH2D * hRecPiplusStartTME = 0x0;
  TH2D * hRecPiplusLastTME = 0x0;

  TH2D * hRecPi0Nshower = 0x0;

  TH2D * hCutnHits = 0x0;
  TH2D * hCutNdEdx = 0x0;
  TH2D * hCutstartE2 = 0x0;
  TH2D * hCutstartE3 = 0x0;
  TH2D * hCuttrackScore = 0x0;
  TH2D * hCutemScore = 0x0;
  TH2D * hCutmichelScore = 0x0;
  TH2D * hCutChi2NDF = 0x0;

  TH2D * hCutnproton = 0x0;
  TH2D * hCutntrack = 0x0;
  TH2D * hCutnshower = 0x0;
  TH2D * hCutnmichel = 0x0;
  TH2D * hCutMpi0 = 0x0;

  //not used keep in case ->
  TH2D * hBeamNdEdx = 0x0;
  
  TH2D * hBeamStartTME = 0x0;
  TH2D * hBeamLastTME = 0x0;
  
  TH2D * hBeamStartE0 = 0x0;
  TH2D * hBeamStartE1 = 0x0;
  TH2D * hBeamStartE2 = 0x0;
  TH2D * hBeamStartE3 = 0x0;
  TH2D * hBeamStartE4 = 0x0;
  TH2D * hBeamStartE5 = 0x0;
  TH2D * hBeamLastE0 = 0x0;
  TH2D * hBeamLastE1 = 0x0;
  TH2D * hBeamLastE2 = 0x0;
  TH2D * hBeamLastE3 = 0x0;
  TH2D * hBeamLastE4 = 0x0;
  TH2D * hBeamLastE5 = 0x0;
  
  TH2D * hBeamPCStartTME = 0x0;
  TH2D * hBeamPCLastTME = 0x0;
  
  TH2D * hBeamPCStartE0 = 0x0;
  TH2D * hBeamPCStartE1 = 0x0;
  TH2D * hBeamPCStartE2 = 0x0;
  TH2D * hBeamPCStartE3 = 0x0;
  TH2D * hBeamPCStartE4 = 0x0;
  TH2D * hBeamPCStartE5 = 0x0;
  TH2D * hBeamPCLastE0 = 0x0;
  TH2D * hBeamPCLastE1 = 0x0;
  TH2D * hBeamPCLastE2 = 0x0;
  TH2D * hBeamPCLastE3 = 0x0;
  TH2D * hBeamPCLastE4 = 0x0;
  TH2D * hBeamPCLastE5 = 0x0;

  TH2D * hBeamLen = 0x0;
  TH2D * hBeamPCLen = 0x0;
  //<--


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

  tout->Branch("beamNdEdx",&beamNdEdx);
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

  tout->Branch("beamStartTME",&beamStartTME);
  tout->Branch("beamLastTME",&beamLastTME);

  //save input info
  tout->Branch("true_beam_PDG",&true_beam_PDG);
  tout->Branch("reco_beam_len",&reco_beam_len);
  tout->Branch("nTrack",&nTrack);

  return tout;
}

TTree * GetInputTree(TFile * fin, const TString tname)
{
  TTree * tree=(TTree*) fin->Get(tname);
  if(!tree){
    cout<<"no tree!"<<endl;
    gDirectory->ls();
    exit(1);
  }

  //--------------------------------------- anaTruth ---------------------------------------
  tree->SetBranchAddress("true_beam_daughter_PDG", &true_beam_daughter_PDG);
  tree->SetBranchAddress("true_beam_daughter_startPx", &true_beam_daughter_startPx);
  tree->SetBranchAddress("true_beam_daughter_startPy", &true_beam_daughter_startPy);
  tree->SetBranchAddress("true_beam_daughter_startPz", &true_beam_daughter_startPz);

  tree->SetBranchAddress("true_beam_daughter_ID", &true_beam_daughter_ID);
  tree->SetBranchAddress("reco_daughter_PFP_true_byHits_ID", &reco_daughter_PFP_true_byHits_ID);
  tree->SetBranchAddress("reco_daughter_allTrack_ID", &reco_daughter_allTrack_ID);

  tree->SetBranchAddress("reco_daughter_allTrack_calibrated_dEdX_SCE", &reco_daughter_allTrack_calibrated_dEdX_SCE);
  tree->SetBranchAddress("reco_daughter_allTrack_Chi2_proton", &reco_daughter_allTrack_Chi2_proton);
  tree->SetBranchAddress("reco_daughter_allTrack_Chi2_ndof", &reco_daughter_allTrack_Chi2_ndof);

  tree->SetBranchAddress("reco_daughter_PFP_ID", &reco_daughter_PFP_ID);
  tree->SetBranchAddress("reco_daughter_PFP_nHits", &reco_daughter_PFP_nHits);
  tree->SetBranchAddress("reco_daughter_PFP_trackScore_collection", &reco_daughter_PFP_trackScore_collection);
  tree->SetBranchAddress("reco_daughter_PFP_emScore_collection", &reco_daughter_PFP_emScore_collection);
  tree->SetBranchAddress("reco_daughter_PFP_michelScore_collection", &reco_daughter_PFP_michelScore_collection);

  tree->SetBranchAddress("reco_daughter_allShower_ID", &reco_daughter_allShower_ID);
  tree->SetBranchAddress("reco_daughter_allShower_dirX", &reco_daughter_allShower_dirX);
  tree->SetBranchAddress("reco_daughter_allShower_dirY", &reco_daughter_allShower_dirY);
  tree->SetBranchAddress("reco_daughter_allShower_dirZ", &reco_daughter_allShower_dirZ);
  tree->SetBranchAddress("reco_daughter_allShower_energy", &reco_daughter_allShower_energy);

  tree->SetBranchAddress("reco_daughter_allTrack_momByRange_proton", &reco_daughter_allTrack_momByRange_proton);
  tree->SetBranchAddress("reco_daughter_allTrack_momByRange_muon", &reco_daughter_allTrack_momByRange_muon);

  tree->SetBranchAddress("reco_daughter_allTrack_Theta", &reco_daughter_allTrack_Theta);
  tree->SetBranchAddress("reco_daughter_allTrack_Phi", &reco_daughter_allTrack_Phi);

  tree->SetBranchAddress("true_beam_endPx", &true_beam_endPx);
  tree->SetBranchAddress("true_beam_endPy", &true_beam_endPy);
  tree->SetBranchAddress("true_beam_endPz", &true_beam_endPz);

  tree->SetBranchAddress("true_beam_PDG", &true_beam_PDG);

  //--------------------------------------- anaRec ---------------------------------------
  tree->SetBranchAddress("data_BI_PDG_candidates", &data_BI_PDG_candidates);
  tree->SetBranchAddress("reco_beam_type", &reco_beam_type);

  tree->SetBranchAddress("reco_beam_trackDirX", &reco_beam_trackDirX);
  tree->SetBranchAddress("reco_beam_trackDirY", &reco_beam_trackDirY);
  tree->SetBranchAddress("reco_beam_trackDirZ", &reco_beam_trackDirZ);
  tree->SetBranchAddress("reco_beam_trackEndDirX", &reco_beam_trackEndDirX);
  tree->SetBranchAddress("reco_beam_trackEndDirY", &reco_beam_trackEndDirY);
  tree->SetBranchAddress("reco_beam_trackEndDirZ", &reco_beam_trackEndDirZ);

  tree->SetBranchAddress("reco_beam_startX", &reco_beam_startX);
  tree->SetBranchAddress("true_beam_startX", &true_beam_startX);
  tree->SetBranchAddress("true_beam_startDirX", &true_beam_startDirX);

  tree->SetBranchAddress("reco_beam_startY", &reco_beam_startY);
  tree->SetBranchAddress("true_beam_startY", &true_beam_startY);
  tree->SetBranchAddress("true_beam_startDirY", &true_beam_startDirY);

  tree->SetBranchAddress("reco_beam_startZ", &reco_beam_startZ);
  tree->SetBranchAddress("true_beam_startZ", &true_beam_startZ);
  tree->SetBranchAddress("true_beam_startDirZ", &true_beam_startDirZ);

  tree->SetBranchAddress("data_BI_X", &data_BI_X);
  tree->SetBranchAddress("data_BI_Y", &data_BI_Y);
  tree->SetBranchAddress("data_BI_dirX", &data_BI_dirX);
  tree->SetBranchAddress("data_BI_dirY", &data_BI_dirY);
  tree->SetBranchAddress("data_BI_dirZ", &data_BI_dirZ);
  tree->SetBranchAddress("data_BI_nMomenta", &data_BI_nMomenta);
  tree->SetBranchAddress("data_BI_nTracks", &data_BI_nTracks);

  tree->SetBranchAddress("reco_beam_endZ", &reco_beam_endZ);

  //to check the difference with _allTrack_
  tree->SetBranchAddress("reco_beam_len", &reco_beam_len);
  //tree->SetBranchAddress("reco_beam_allTrack_len", &reco_beam_len);
  tree->SetBranchAddress("reco_beam_calibrated_dEdX", &reco_beam_calibrated_dEdX);
  //tree->SetBranchAddress("reco_beam_allTrack_calibrated_dEdX", &reco_beam_calibrated_dEdX);

  return tree;
}

void IniRecHist(TList * lout, const TString tag)
{
  const int nPass = 2;
  const double Passmin = -0.5;
  const double Passmax = 1.5;

  hTruthBeamType = new TH1I("a0TruthBeamType"+tag,  "", 20, -0.5, 19.5); lout->Add(hTruthBeamType);
  hTruthSignal = new TH1I("a1TruthSignal"+tag,  "",  nPass, Passmin, Passmax); lout->Add(hTruthSignal);

  const int nevtType = 10;
  const double evtTypemin = -0.5;
  const double evtTypemax = 9.5;

  const int nE = 50;
  const double Emin = 0;
  const double lastEmax = 20;
  const double startEmax =5;

  hBeamNdEdx   = new TH2D("d000BeamNdEdxSTK"+tag,"", 50, 0, 100, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamNdEdx);

  hBeamStartE0 = new TH2D("d010BeamStartE0STK"+tag,"", nE, Emin, startEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamStartE0);
  hBeamStartE1 = new TH2D("d011BeamStartE1STK"+tag,"", nE, Emin, startEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamStartE1);
  hBeamStartE2 = new TH2D("d012BeamStartE2STK"+tag,"", nE, Emin, startEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamStartE2);
  hBeamStartE3 = new TH2D("d013BeamStartE3STK"+tag,"", nE, Emin, startEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamStartE3);
  hBeamStartE4 = new TH2D("d014BeamStartE4STK"+tag,"", nE, Emin, startEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamStartE4);
  hBeamStartE5 = new TH2D("d015BeamStartE5STK"+tag,"", nE, Emin, startEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamStartE5);
  hBeamStartTME= new TH2D("d016BeamStartTMESTKNTX"+tag,"",nE*5, 1.2, 4.2, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamStartTME);

  hBeamLastE0  = new TH2D("d020BeamLastE0STK"+tag,"", nE, Emin, lastEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamLastE0);
  hBeamLastE1  = new TH2D("d021BeamLastE1STK"+tag,"", nE, Emin, lastEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamLastE1);
  hBeamLastE2  = new TH2D("d022BeamLastE2STK"+tag,"", nE, Emin, lastEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamLastE2);
  hBeamLastE3  = new TH2D("d023BeamLastE3STK"+tag,"", nE, Emin, lastEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamLastE3);
  hBeamLastE4  = new TH2D("d024BeamLastE4STK"+tag,"", nE, Emin, lastEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamLastE4);
  hBeamLastE5  = new TH2D("d025BeamLastE5STK"+tag,"", nE, Emin, lastEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamLastE5);
  hBeamLastTME = new TH2D("d026BeamLastTMESTKNTX"+tag,"",nE*5, 0, 15, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamLastTME);

  hBeamPCStartE0 = new TH2D("d110BeamPCStartE0STK"+tag,"", nE, Emin, startEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamPCStartE0);
  hBeamPCStartE1 = new TH2D("d111BeamPCStartE1STK"+tag,"", nE, Emin, startEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamPCStartE1);
  hBeamPCStartE2 = new TH2D("d112BeamPCStartE2STK"+tag,"", nE, Emin, startEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamPCStartE2);
  hBeamPCStartE3 = new TH2D("d113BeamPCStartE3STK"+tag,"", nE, Emin, startEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamPCStartE3);
  hBeamPCStartE4 = new TH2D("d114BeamPCStartE4STK"+tag,"", nE, Emin, startEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamPCStartE4);
  hBeamPCStartE5 = new TH2D("d115BeamPCStartE5STK"+tag,"", nE, Emin, startEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamPCStartE5);
  hBeamPCStartTME= new TH2D("d116BeamPCStartTMESTKNTX"+tag,"",nE*5, 1.2, 4.2, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamPCStartTME);

  hBeamPCLastE0  = new TH2D("d120BeamPCLastE0STK"+tag,"", nE, Emin, lastEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamPCLastE0);
  hBeamPCLastE1  = new TH2D("d121BeamPCLastE1STK"+tag,"", nE, Emin, lastEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamPCLastE1);
  hBeamPCLastE2  = new TH2D("d122BeamPCLastE2STK"+tag,"", nE, Emin, lastEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamPCLastE2);
  hBeamPCLastE3  = new TH2D("d123BeamPCLastE3STK"+tag,"", nE, Emin, lastEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamPCLastE3);
  hBeamPCLastE4  = new TH2D("d124BeamPCLastE4STK"+tag,"", nE, Emin, lastEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamPCLastE4);
  hBeamPCLastE5  = new TH2D("d125BeamPCLastE5STK"+tag,"", nE, Emin, lastEmax, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamPCLastE5);
  hBeamPCLastTME = new TH2D("d126BeamPCLastTMESTKNTX"+tag,"",nE*5, 0, 15, nevtType, evtTypemin, evtTypemax); lout->Add(hBeamPCLastTME);

  //not used
  /*
  //sensitive to pi+ beam, but not to signal -> sensitive to pi+ beam non-signal, i.e. non-interacting
  const int nlen = 10;
  const double lenmin = 0;
  const double lenmax = 250;
  hBeamLen = new TH2D("p9BeamLenNOHPRFPUR"+tag,"", nlen, lenmin, lenmax, nPass, Passmin, Passmax); lout->Add(hBeamLen);
  hBeamPCLen = new TH2D("p9BeamPCLenNOHPRFPUR"+tag,"", nlen, lenmin, lenmax, nPass, Passmin, Passmax); lout->Add(hBeamPCLen);
  */

  const int nres = 20;
  const double resmin = -1;
  const double resmax = 1;

  //beam
  const int nbmTheta = 80;
  const double bmThetamin = 0;
  const double bmThetamax = 60;//in deg
  //daughter
  const int ndTheta = 20;
  const double dThetamin = 0;
  const double dThetamax = 180;//in deg

  const int nmomentum = 30;
  const double momentummin = 0;
  const double momentummax = 1.2;

  const int ndedx = 60;
  const double dedxmin = 0;
  const double dedxmax = 30;

  const int nscore = 50;
  const double scoremin = 0;
  const double scoremax = 1;

  const int ncounter = 10;
  const double countermin = -0.5;
  const double countermax = 9.5;

  //NOH, STK are tags for style::Process2DHist

  hBeamThetaRes      = new TH2D("b000BeamThetaResNOH"+tag,"",      nbmTheta, bmThetamin, bmThetamax, 25, -20, 30); lout->Add(hBeamThetaRes);
  hProtonThetaRes    = new TH2D("b001ProtonThetaResNOH"+tag,"",    ndTheta, dThetamin, dThetamax, 25, -20, 30); lout->Add(hProtonThetaRes);
  hProtonMomentumRes = new TH2D("b002ProtonMomentumResNOH"+tag,"", nmomentum, momentummin, momentummax, nres, -0.2, 0.2); lout->Add(hProtonMomentumRes);
  hPiplusThetaRes    = new TH2D("b003PiplusThetaResNOH"+tag,"",        ndTheta, dThetamin, dThetamax, 25, -20, 30); lout->Add(hPiplusThetaRes);
  hPiplusMomentumRes = new TH2D("b004PiplusMomentumResNOH"+tag,"",     nmomentum, momentummin, momentummax, nres, resmin, resmax); lout->Add(hPiplusMomentumRes);

  hRecBeamTheta       = new TH2D("b010bRecBeamThetaSTK"+tag,"",       nbmTheta, bmThetamin, bmThetamax, nevtType, evtTypemin, evtTypemax); lout->Add(hRecBeamTheta);
  hRecProtonTheta     = new TH2D("b011bRecProtonThetaSTK"+tag,"",     ndTheta, dThetamin, dThetamax, nevtType, evtTypemin, evtTypemax); lout->Add(hRecProtonTheta);
  hRecProtonMomentum  = new TH2D("b012bRecProtonMomentumSTK"+tag,"",  nmomentum, momentummin, momentummax, nevtType, evtTypemin, evtTypemax); lout->Add(hRecProtonMomentum);
  hRecPiplusTheta     = new TH2D("b013bRecPiplusThetaSTK"+tag,"",     ndTheta, dThetamin, dThetamax, nevtType, evtTypemin, evtTypemax); lout->Add(hRecPiplusTheta);
  hRecPiplusMomentum  = new TH2D("b014bRecPiplusMomentumSTK"+tag,"",  nmomentum, momentummin, momentummax, nevtType, evtTypemin, evtTypemax); lout->Add(hRecPiplusMomentum);

  hRecProtonLastE2    = new TH2D("b020RecProtonLastE2STK"+tag,"",     ndedx, dedxmin, dedxmax, nevtType, evtTypemin, evtTypemax); lout->Add(hRecProtonLastE2);
  hRecProtonLastE3    = new TH2D("b021RecProtonLastE3STK"+tag,"",     ndedx, dedxmin, dedxmax, nevtType, evtTypemin, evtTypemax); lout->Add(hRecProtonLastE3);
  hRecPiplusLastE2    = new TH2D("b022RecPiplusLastE2STK"+tag,"",     ndedx, dedxmin, dedxmax, nevtType, evtTypemin, evtTypemax); lout->Add(hRecPiplusLastE2);
  hRecPiplusLastE3    = new TH2D("b023RecPiplusLastE3STK"+tag,"",     ndedx, dedxmin, dedxmax, nevtType, evtTypemin, evtTypemax); lout->Add(hRecPiplusLastE3);
  hRecProtonStartE2   = new TH2D("b024RecProtonStartE2STK"+tag,"",     ndedx, dedxmin, dedxmax, nevtType, evtTypemin, evtTypemax); lout->Add(hRecProtonStartE2);
  hRecProtonStartE3   = new TH2D("b025RecProtonStartE3STK"+tag,"",     ndedx, dedxmin, dedxmax, nevtType, evtTypemin, evtTypemax); lout->Add(hRecProtonStartE3);
  hRecPiplusStartE2   = new TH2D("b026RecPiplusStartE2STK"+tag,"",     ndedx, dedxmin, dedxmax, nevtType, evtTypemin, evtTypemax); lout->Add(hRecPiplusStartE2);
  hRecPiplusStartE3   = new TH2D("b027RecPiplusStartE3STK"+tag,"",     ndedx, dedxmin, dedxmax, nevtType, evtTypemin, evtTypemax); lout->Add(hRecPiplusStartE3);
  hRecProtonStartTME  = new TH2D("b028RecProtonStartTMESTK"+tag,"",     ndedx, dedxmin, dedxmax, nevtType, evtTypemin, evtTypemax); lout->Add(hRecProtonStartTME);
  hRecProtonLastTME   = new TH2D("b028RecProtonLastTMESTK"+tag,"",     ndedx, dedxmin, dedxmax, nevtType, evtTypemin, evtTypemax); lout->Add(hRecProtonLastTME);
  hRecPiplusStartTME  = new TH2D("b029RecPiplusStartTMESTK"+tag,"",     ndedx, dedxmin, dedxmax, nevtType, evtTypemin, evtTypemax); lout->Add(hRecPiplusStartTME);
  hRecPiplusLastTME   = new TH2D("b029RecPiplusLastTMESTK"+tag,"",     ndedx, dedxmin, dedxmax, nevtType, evtTypemin, evtTypemax); lout->Add(hRecPiplusLastTME);

  hRecPi0Nshower      = new TH2D("b030bRecPi0NshowerSTK"+tag,"",      ncounter, countermin, countermax, nevtType, evtTypemin, evtTypemax); lout->Add(hRecPi0Nshower);

  hCutNdEdx          = new TH2D("b100CutNdEdxSTK"+tag,"",          50, 0, 50, nevtType, evtTypemin, evtTypemax); lout->Add(hCutNdEdx);
  hCuttrackScore     = new TH2D("b101CuttrackScoreSTK"+tag,"",     nscore, scoremin, scoremax, nevtType, evtTypemin, evtTypemax); lout->Add(hCuttrackScore);
  hCutemScore        = new TH2D("b102CutemScoreSTK"+tag,"",        nscore, scoremin, scoremax, nevtType, evtTypemin, evtTypemax); lout->Add(hCutemScore);
  hCutmichelScore    = new TH2D("b103CutmichelScoreSTK"+tag,"",    nscore, scoremin, scoremax, nevtType, evtTypemin, evtTypemax); lout->Add(hCutmichelScore);
  hCutnHits          = new TH2D("b104CutnHitsSTK"+tag,"",          50, 0, 500, nevtType, evtTypemin, evtTypemax); lout->Add(hCutnHits);
  hCutChi2NDF        = new TH2D("b105CutChi2NDFSTK"+tag,"",        30, 0, 500, nevtType, evtTypemin, evtTypemax); lout->Add(hCutChi2NDF);
  hCutstartE2        = new TH2D("b106CutstartE2STK"+tag,"",        ndedx, dedxmin, dedxmax, nevtType, evtTypemin, evtTypemax); lout->Add(hCutstartE2);
  hCutstartE3        = new TH2D("b107CutstartE3STK"+tag,"",        ndedx, dedxmin, dedxmax, nevtType, evtTypemin, evtTypemax); lout->Add(hCutstartE3);

  hCutbeamID       = new TH1I("c000CutbeamID"+tag,"", 2, -0.5, 1.5); lout->Add(hCutbeamID);
  hCutBeamType     = new TH2D("c001BeamTypeSTK"+tag,  "", 30, -4.5, 25.5, nevtType, evtTypemin, evtTypemax); lout->Add(hCutBeamType);
  hCutBeamPosPass  = new TH2D("c002BeamPosPassSTK"+tag, "", 4, -0.5, 3.5, nevtType, evtTypemin, evtTypemax); lout->Add(hCutBeamPosPass);
  hCutBeamEndZ     = new TH2D("c003BeamEndZSTK"+tag,"",50, 0, 500, nevtType, evtTypemin, evtTypemax); lout->Add(hCutBeamEndZ);
  hCutBeamEndZPass = new TH2D("c004BeamEndZPassSTK"+tag,"",4, -0.5, 3.5, nevtType, evtTypemin, evtTypemax); lout->Add(hCutBeamEndZPass);

  hCutnshower      = new TH2D("c101CutnshowerSTK"+tag,"", ncounter, countermin, countermax, nevtType, evtTypemin, evtTypemax); lout->Add(hCutnshower);
  hCutnmichel      = new TH2D("c102CutnmichelSTK"+tag,"", ncounter, countermin, countermax, nevtType, evtTypemin, evtTypemax); lout->Add(hCutnmichel);
  hCutntrack       = new TH2D("c103CutntrackSTK"+tag,"",  ncounter, countermin, countermax, nevtType, evtTypemin, evtTypemax); lout->Add(hCutntrack);
  hCutnproton      = new TH2D("c104CutnprotonSTK"+tag,"", ncounter, countermin, countermax, nevtType, evtTypemin, evtTypemax); lout->Add(hCutnproton);
  hCutMpi0         = new TH2D("c105CutMpi0STK"+tag,"", 15, 0, 0.3, nevtType, evtTypemin, evtTypemax); lout->Add(hCutMpi0);
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
