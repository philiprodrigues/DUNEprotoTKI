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

void drawTracking(TList *lout, const TString pretag, const TString addCut)
{
  //_____________________________________________________ basic settings _____________________________________________________ 

  const bool kPiZero = pretag.Contains("MPiZero");
  const bool kTrackingProton = !pretag.Contains("TrackingPiPlus");

  cout<<"\n\n                       Running kPiZero "<<kPiZero<<" TrackingProton "<<kTrackingProton<<endl<<endl;

  TTree * tree = (TTree*)gDirectory->Get("mc/tree");
  if(!tree){
    cout<<"no tree!"<<endl;
    exit(1);
  }

  //_____________________________________________________ define cuts _____________________________________________________ 

  const TString varPrec = kTrackingProton?"recProtonmomentum":"recPiPlusmomentum";

  //only look into signal sample, not the whole beam sample
  const TString basCut = "kSignal";

  //reconstructed sample
  const TString recCut = basCut +" && " + varPrec+"!=-999 ";

  //further cut
  const TString selCut = recCut + addCut;

  //_____________________________________________________ Drawing with cuts _____________________________________________________

  //====================== Efficiency ======================

  const int nPtrue = 12;
  const double Ptruemin = 0;
  const double Ptruemax = 1.2;
  TH1D * hPtruebas = new TH1D(pretag+"h0all",basCut,nPtrue, Ptruemin, Ptruemax); lout->Add(hPtruebas);
  TH1D * hPtruerec = new TH1D(pretag+"h1rec",recCut,nPtrue, Ptruemin, Ptruemax); lout->Add(hPtruerec);
  TH1D * hPtruesel = 0x0;
  if(addCut!=""){
    hPtruesel = new TH1D(pretag+"h2sel",selCut,nPtrue, Ptruemin, Ptruemax); lout->Add(hPtruesel);
  }

  const TString trackingparticle = kTrackingProton?"p":"#pi";
  const TString varPtrue = kTrackingProton?"finProtonmomentum":"finPimomentum";
  const TString titN = "N";
  const TString titPtrue = Form("#it{p}_{%s}^{true} (GeV/#it{c})", trackingparticle.Data());
  style::GetHist(varPtrue, titPtrue, titN, tree, hPtruebas);
  style::GetHist(varPtrue, titPtrue, titN, tree, hPtruerec);
  style::GetHist(varPtrue, titPtrue, titN, tree, hPtruesel);

  style::GetEff(hPtruerec, hPtruebas, 1, lout); 
  style::GetEff(hPtruesel, hPtruerec, 1, lout); 

  //====================== 1D Resolution ======================
  
  const int nRes = 15;
  const double Resmin = -1;
  const double Resmax = 1;
  const TString titRes = Form("#it{p}_{%s}^{rec}/#it{p}_{%s}^{true}-1", trackingparticle.Data(), trackingparticle.Data());

  TH1D * hResRec1D = new TH1D(pretag+"hResRec1D",selCut, nRes*3, Resmin, Resmax); lout->Add(hResRec1D);
  const TString varRes = Form("(%s/%s)-1", varPrec.Data(), varPtrue.Data());
  style::GetHist(varRes, titRes, titN, tree, hResRec1D);

  //_____________________________________________________ only draw 2D for beam sample, NOH, PRF are tags for style::Process2DHist  _____________________________________________________
  if(addCut==""){
    
    //new cut! only beam sample!
    const TString cut2d = varPrec+"!=-999 ";

    //====================== 2D Resolution vs p_true ======================
    
    TH2D * hResPtrueRec = new TH2D(pretag+"hResPtrueRecNOHPRF",cut2d, nPtrue, Ptruemin, Ptruemax, nRes, Resmin, Resmax); lout->Add(hResPtrueRec);
    const TString varResPtrue = Form("(%s/%s)-1 : %s", varPrec.Data(), varPtrue.Data(), varPtrue.Data());
    style::GetHist(varResPtrue, titPtrue, titRes, tree, hResPtrueRec);

    //====================== 2D Resolution vs chi2 ======================

    const int nChi2 = 25;
    const double Chi2min = 0;
    const double Chi2max = 500;
    
    const TString varChi2 = "chi2/ndof";
    const TString titChi2 = "#chi^{2}/NDF";
    
    TH2D * hResChi2 = new TH2D(pretag+"hResChi2NOHPRF", cut2d, nChi2, Chi2min, Chi2max, nRes, Resmin, Resmax); lout->Add(hResChi2);
    const TString varResChi2 = Form("(%s/%s)-1 : chi2/ndof", varPrec.Data(), varPtrue.Data());
    style::GetHist(varResChi2, titChi2, titRes, tree, hResChi2);

    //====================== 2D Resolution vs all kinds of dEdx ======================
    const double Emin = 0;
    
    const int nstartE      = kTrackingProton? 30 : 40;
    const double startEmax = kTrackingProton? 30 : 10;
    
    const int nlastE       = kTrackingProton? 40 : 40;
    const double lastEmax  = kTrackingProton? 10 : 10;
    
    //only save startE2 for resolution
    //PRF tag to get profileX to see momentum bais
    TH2D * hResstartE0Rec = new TH2D(pretag+"h5ResstartE0RecNOHPRF", cut2d, nstartE, Emin, startEmax, nRes, Resmin, Resmax); //lout->Add(hResstartE0Rec);
    TH2D * hResstartE1Rec = new TH2D(pretag+"h5ResstartE1RecNOHPRF", cut2d, nstartE, Emin, startEmax, nRes, Resmin, Resmax); //lout->Add(hResstartE1Rec);
    TH2D * hResstartE2Rec = new TH2D(pretag+"h5ResstartE2RecNOHPRF", cut2d, nstartE, Emin, startEmax, nRes, Resmin, Resmax); lout->Add(hResstartE2Rec);
    TH2D * hResstartE3Rec = new TH2D(pretag+"h5ResstartE3RecNOHPRF", cut2d, nstartE, Emin, startEmax, nRes, Resmin, Resmax); //lout->Add(hResstartE3Rec);
    TH2D * hResstartE4Rec = new TH2D(pretag+"h5ResstartE4RecNOHPRF", cut2d, nstartE, Emin, startEmax, nRes, Resmin, Resmax); //lout->Add(hResstartE4Rec);
    TH2D * hResstartE5Rec = new TH2D(pretag+"h5ResstartE5RecNOHPRF", cut2d, nstartE, Emin, startEmax, nRes, Resmin, Resmax); //lout->Add(hResstartE5Rec);
    TH2D * hResstartEsRec[]={hResstartE0Rec, hResstartE1Rec, hResstartE2Rec, hResstartE3Rec, hResstartE4Rec, hResstartE5Rec};
    
    TH2D * hReslastE0Rec = new TH2D(pretag+"h6ReslastE0RecNOHPRF", cut2d, nlastE, Emin, lastEmax, nRes, Resmin, Resmax); //lout->Add(hReslastE0Rec);
    TH2D * hReslastE1Rec = new TH2D(pretag+"h6ReslastE1RecNOHPRF", cut2d, nlastE, Emin, lastEmax, nRes, Resmin, Resmax); //lout->Add(hReslastE1Rec);
    TH2D * hReslastE2Rec = new TH2D(pretag+"h6ReslastE2RecNOHPRF", cut2d, nlastE, Emin, lastEmax, nRes, Resmin, Resmax); //lout->Add(hReslastE2Rec);
    TH2D * hReslastE3Rec = new TH2D(pretag+"h6ReslastE3RecNOHPRF", cut2d, nlastE, Emin, lastEmax, nRes, Resmin, Resmax); //lout->Add(hReslastE3Rec);
    TH2D * hReslastE4Rec = new TH2D(pretag+"h6ReslastE4RecNOHPRF", cut2d, nlastE, Emin, lastEmax, nRes, Resmin, Resmax); //lout->Add(hReslastE4Rec);
    TH2D * hReslastE5Rec = new TH2D(pretag+"h6ReslastE5RecNOHPRF", cut2d, nlastE, Emin, lastEmax, nRes, Resmin, Resmax); //lout->Add(hReslastE5Rec);
    TH2D * hReslastEsRec[]={hReslastE0Rec, hReslastE1Rec, hReslastE2Rec, hReslastE3Rec, hReslastE4Rec, hReslastE5Rec};
    
    //only saving TME10 for resolution
    TH2D * hResstartTruncatedMeanE10Rec = new TH2D(pretag+"h7ResstartTruncatedMeanE10RecNOHPRF", cut2d, nstartE, Emin, startEmax, nRes, Resmin, Resmax); lout->Add(hResstartTruncatedMeanE10Rec);
    TH2D * hResstartTruncatedMeanE20Rec = new TH2D(pretag+"h7ResstartTruncatedMeanE20RecNOHPRF", cut2d, nstartE, Emin, startEmax, nRes, Resmin, Resmax); //lout->Add(hResstartTruncatedMeanE20Rec);
    TH2D * hResstartTruncatedMeanE30Rec = new TH2D(pretag+"h7ResstartTruncatedMeanE30RecNOHPRF", cut2d, nstartE, Emin, startEmax, nRes, Resmin, Resmax); //lout->Add(hResstartTruncatedMeanE30Rec);
    TH2D * hResstartTruncatedMeanE40Rec = new TH2D(pretag+"h7ResstartTruncatedMeanE40RecNOHPRF", cut2d, nstartE, Emin, startEmax, nRes, Resmin, Resmax); //lout->Add(hResstartTruncatedMeanE40Rec);
    TH2D * hResstartTruncatedMeanE50Rec = new TH2D(pretag+"h7ResstartTruncatedMeanE50RecNOHPRF", cut2d, nstartE, Emin, startEmax, nRes, Resmin, Resmax); //lout->Add(hResstartTruncatedMeanE50Rec);
    TH2D * hResstartTruncatedMeanEsRec[]={0x0, hResstartTruncatedMeanE10Rec, hResstartTruncatedMeanE20Rec, hResstartTruncatedMeanE30Rec, hResstartTruncatedMeanE40Rec, hResstartTruncatedMeanE50Rec};

    //====================== dEdx vs. p_true ======================
    //h8 h9 enough to show Bragg peak and motivate ESC cut for proton
    TH2D * hstartE0Ptrue = new TH2D(pretag+"h8PtruestartE0PtrueNOH", cut2d, nPtrue, Ptruemin, Ptruemax, nstartE/2, Emin, startEmax); lout->Add(hstartE0Ptrue);
    TH2D * hstartE1Ptrue = new TH2D(pretag+"h8PtruestartE1PtrueNOH", cut2d, nPtrue, Ptruemin, Ptruemax, nstartE/2, Emin, startEmax); lout->Add(hstartE1Ptrue);
    TH2D * hstartE2Ptrue = new TH2D(pretag+"h8PtruestartE2PtrueNOH", cut2d, nPtrue, Ptruemin, Ptruemax, nstartE/2, Emin, startEmax); lout->Add(hstartE2Ptrue);
    TH2D * hstartE3Ptrue = new TH2D(pretag+"h8PtruestartE3PtrueNOH", cut2d, nPtrue, Ptruemin, Ptruemax, nstartE/2, Emin, startEmax); lout->Add(hstartE3Ptrue);
    TH2D * hstartE4Ptrue = new TH2D(pretag+"h8PtruestartE4PtrueNOH", cut2d, nPtrue, Ptruemin, Ptruemax, nstartE/2, Emin, startEmax); lout->Add(hstartE4Ptrue);
    TH2D * hstartE5Ptrue = new TH2D(pretag+"h8PtruestartE5PtrueNOH", cut2d, nPtrue, Ptruemin, Ptruemax, nstartE/2, Emin, startEmax); lout->Add(hstartE5Ptrue);
    TH2D * hstartEsPtrue[]={hstartE0Ptrue, hstartE1Ptrue, hstartE2Ptrue, hstartE3Ptrue, hstartE4Ptrue, hstartE5Ptrue};
    
    TH2D * hlastE0Ptrue = new TH2D(pretag+"h9PtruelastE0PtrueNOH", cut2d, nPtrue, Ptruemin, Ptruemax, nlastE/2, Emin, lastEmax); lout->Add(hlastE0Ptrue);
    TH2D * hlastE1Ptrue = new TH2D(pretag+"h9PtruelastE1PtrueNOH", cut2d, nPtrue, Ptruemin, Ptruemax, nlastE/2, Emin, lastEmax); lout->Add(hlastE1Ptrue);
    TH2D * hlastE2Ptrue = new TH2D(pretag+"h9PtruelastE2PtrueNOH", cut2d, nPtrue, Ptruemin, Ptruemax, nlastE/2, Emin, lastEmax); lout->Add(hlastE2Ptrue);
    TH2D * hlastE3Ptrue = new TH2D(pretag+"h9PtruelastE3PtrueNOH", cut2d, nPtrue, Ptruemin, Ptruemax, nlastE/2, Emin, lastEmax); lout->Add(hlastE3Ptrue);
    TH2D * hlastE4Ptrue = new TH2D(pretag+"h9PtruelastE4PtrueNOH", cut2d, nPtrue, Ptruemin, Ptruemax, nlastE/2, Emin, lastEmax); lout->Add(hlastE4Ptrue);
    TH2D * hlastE5Ptrue = new TH2D(pretag+"h9PtruelastE5PtrueNOH", cut2d, nPtrue, Ptruemin, Ptruemax, nlastE/2, Emin, lastEmax); lout->Add(hlastE5Ptrue);
    TH2D * hlastEsPtrue[]={hlastE0Ptrue, hlastE1Ptrue, hlastE2Ptrue, hlastE3Ptrue, hlastE4Ptrue, hlastE5Ptrue};
    
    //---------- draw all dEdx stuff together ----------

    const TString varlastEbase="lastE";
    const TString varstartEbase="startE";
    const TString varstartTruncatedMeanEbase="startTruncatedMeanE";
    const TString titlastEbase="end #it{E}";
    const TString titstartEbase="start #it{E}";
    const TString titstartTruncatedMeanEbase="Truncated Mean (40-95%) of start #it{E} with Nodes 2 to 2 +";
    for(unsigned int ii=0; ii<6; ii++){
      const TString varlastE = varlastEbase+Form("%d",ii);
      const TString varstartE = varstartEbase+Form("%d",ii);
      const TString titlastE = titlastEbase+Form("_{%d}",ii)+" (MeV/cm)";
      const TString titstartE = titstartEbase+Form("_{%d}",ii)+" (MeV/cm)";
      
      const TString varResstartE = Form("(%s/%s)-1 : ", varPrec.Data(), varPtrue.Data())+varstartE;
      style::GetHist(varResstartE, titstartE, titRes, tree, hResstartEsRec[ii]);
      
      const TString varReslastE = Form("(%s/%s)-1 : ", varPrec.Data(), varPtrue.Data())+varlastE;
      style::GetHist(varReslastE, titlastE, titRes, tree, hReslastEsRec[ii]);
      
      if(ii){
        const TString varstartTruncatedMeanE = varstartTruncatedMeanEbase+Form("%d",ii*10);
        const TString titstartTruncatedMeanE = titstartTruncatedMeanEbase+Form(" %d",ii*10)+" (MeV/cm)";
        
        const TString varResstartTruncatedMeanE = Form("(%s/%s)-1 : ", varPrec.Data(), varPtrue.Data())+varstartTruncatedMeanE;
        style::GetHist(varResstartTruncatedMeanE, titstartTruncatedMeanE, titRes, tree, hResstartTruncatedMeanEsRec[ii]);
      }
      
      const TString varstartEPtrue =varstartE+":"+varPtrue;
      style::GetHist(varstartEPtrue, titPtrue, titstartE, tree, hstartEsPtrue[ii]);
      
      const TString varlastEPtrue =varlastE+":"+varPtrue;
      style::GetHist(varlastEPtrue, titPtrue, titlastE, tree, hlastEsPtrue[ii]);
    }
    
  }//====================== end of only draw 2D for beam sample
}

int main(int argc, char* argv[])
{
  //_____________________________________________________ basic settings _____________________________________________________ 

  if(argc!=3){
    cout<<"argc!=2 !"<<argc<<endl;
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
  TList * lout= new TList;

  //_____________________________________________________ draw MC truth-matched reconstruction performance for different cuts _____________________________________________________
  const TString inputfn = Form("../anaData/output/outanaData_%s_anaTruth.root", tag.Data());
  TFile *fin = new TFile(inputfn);
  if(!fin->IsOpen()){
    cout<<"fin not open!"<<endl;
    exit(1);
  }
  cout<<"Reading "<<inputfn<<endl;

  drawTracking(lout, tag+"0000raw", "");

  if(kProton){
    //historical cuts
    /*
      seeana010.log:check cut kpi0 0 ndedxcut 16
      seeana010.log:check cut kpi0 0 proton tag Chi2NDF 50.00 nhits 260
      
      seeana110.log:check cut kpi0 1 ndedxcut 6
      seeana110.log:check cut kpi0 1 proton tag startE2 10.00 nhits 260 startE3 9.00
    */
    const TString esccut = "&& (ndEdxCls>=6) && (startE2>10) && (startE3>9)";
    const TString chicut = "&& (ndEdxCls>=6) && (chi2/ndof<31)";//490; 30, 484; 31, 489; 31.5, 493; 35, 504; 
    const TString tmecut10 = "&& (ndEdxCls>=6) && (startTruncatedMeanE10>5.8)";//5, 505; 5.5, 499; 5.8, 491; 5.9, 487
    const TString tmecut20 = "&& (ndEdxCls>=6) && (startTruncatedMeanE20>2)";//3, 401; 2, 412
    drawTracking(lout, tag+"0001ESC",   esccut);
    drawTracking(lout, tag+"0002CHI",   chicut);
    //drawTracking(lout, tag+"0003TME10", tmecut10);
    //drawTracking(lout, tag+"0003TME20", tmecut20);
  }

  //_____________________________________________________ Automatically handling histograms and saving _____________________________________________________

  style::Process2DHist(lout);
  const bool ktext = false;
  const bool kfast = false;
  //no MC scaling because of no overlay from data
  style::DrawHist(lout, 1, 0x0, "output", tag, ktext, kfast);

  TFile * fout= new TFile(Form("output/outdrawTracking_%s.root", tag.Data()),"recreate");
  lout->Write();
  fout->Save();
  fout->Close();

  return 0;
}
