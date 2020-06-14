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

void getHist(const TString var, const TString xtit, const TString ytit, TTree *tree, TH1 *hist)
{
  const TString cut = hist->GetTitle();
  const TString hname = hist->GetName();
  printf("Drawing %s %s %s\n", var.Data(), hname.Data(), cut.Data());
  tree->Draw(var+">>"+hname, cut);

  hist->SetXTitle(xtit);
  hist->SetYTitle(ytit);
}

void drawTracking(TList *lout, const TString pretag, const TString addCut="")
{
  const bool kPiZero = pretag.Contains("MPiZero");
  const bool kTrackingProton = !pretag.Contains("TrackingPiPlus");

  cout<<"\n\n                       Running kPiZero "<<kPiZero<<" TrackingProton "<<kTrackingProton<<endl<<endl;

  TTree * tree = (TTree*)gDirectory->Get("tree");
  if(!tree){
    cout<<"no tree!"<<endl;
    exit(1);
  }

  const TString varPrec = kTrackingProton?"recProtonmomentum":"recPiPlusmomentum";

  const TString selCut = varPrec+"!=-999 " + addCut;

  //========================================================================
  const int nPtrue = 20;
  const double Ptruemin = 0;
  const double Ptruemax = 1.2;
  TH1D * hPtrueall = new TH1D(pretag+"h0all","1",nPtrue, Ptruemin, Ptruemax); lout->Add(hPtrueall);
  TH1D * hPtruerec = new TH1D(pretag+"h1rec",varPrec+"!=-999",nPtrue, Ptruemin, Ptruemax); lout->Add(hPtruerec);
  TH1D * hPtruencl = new TH1D(pretag+"h2ncl",selCut,nPtrue, Ptruemin, Ptruemax); lout->Add(hPtruencl);
  TH1D * hPtruesig = new TH1D(pretag+"h3sig",selCut+"&& kSignal",nPtrue, Ptruemin, Ptruemax); lout->Add(hPtruesig);

  const TString trackingparticle = kTrackingProton?"p":"#pi";
  const TString varPtrue = kTrackingProton?"finProtonmomentum":"finPimomentum";
  const TString ytitN = "N";
  const TString xtitPtrue = Form("#it{p}_{%s}^{true} (GeV/#it{c})", trackingparticle.Data());
  getHist(varPtrue, xtitPtrue, ytitN, tree, hPtrueall);
  getHist(varPtrue, xtitPtrue, ytitN, tree, hPtruerec);
  getHist(varPtrue, xtitPtrue, ytitN, tree, hPtruencl);
  getHist(varPtrue, xtitPtrue, ytitN, tree, hPtruesig);

  TH1D * effPtruerec = style::GetEff(hPtruerec, hPtrueall, 0.7); lout->Add(effPtruerec);
  TH1D * effPtruencl = style::GetEff(hPtruencl, hPtrueall, 0.7); lout->Add(effPtruencl);
  TH1D * effPtruesig = style::GetEff(hPtruesig, hPtrueall, 0.2); lout->Add(effPtruesig);

  //========================================================================
  const int nRes = 15;
  const double Resmin = -1;
  const double Resmax = 1;

  TH2D * hResPtrueRec = new TH2D(pretag+"hResPtrueRecNOHPRF",selCut, nPtrue, Ptruemin, Ptruemax, nRes, Resmin, Resmax); lout->Add(hResPtrueRec);

  const TString varResPtrue = Form("(%s/%s)-1 : %s", varPrec.Data(), varPtrue.Data(), varPtrue.Data());
  const TString ytitRes = Form("#it{p}_{%s}^{rec}/#it{p}_{%s}^{true}-1", trackingparticle.Data(), trackingparticle.Data());
  getHist(varResPtrue, xtitPtrue, ytitRes, tree, hResPtrueRec);

  TH1D * hResRecSignal = new TH1D(pretag+"hResRecSignal",selCut+"&& kSignal", nRes*3, Resmin, Resmax); lout->Add(hResRecSignal);
  TH1D * hResRecBeam = new TH1D(pretag+"hResRecBeam",selCut, nRes*3, Resmin, Resmax); lout->Add(hResRecBeam);
  const TString varRes = Form("(%s/%s)-1", varPrec.Data(), varPtrue.Data());
  const TString xtitRes = ytitRes;
  getHist(varRes, xtitRes, ytitN, tree, hResRecSignal);
  getHist(varRes, xtitRes, ytitN, tree, hResRecBeam);

  if(addCut==""){
  //========================================================================
  const int nChi2 = 20;
  const double Chi2min = 0;
  const double Chi2max = 100;

  const TString varChi2 = "chi2/ndof";
  const TString xtitChi2 = "#chi^{2}/NDF";

  TH2D * hResChi2 = new TH2D(pretag+"hResChi2NOHPRF", selCut, nChi2, Chi2min, Chi2max, nRes, Resmin, Resmax); lout->Add(hResChi2);
  const TString varResChi2 = Form("(%s/%s)-1 : chi2/ndof", varPrec.Data(), varPtrue.Data());
  getHist(varResChi2, xtitChi2, ytitRes, tree, hResChi2);

  //========================================================================
  const int nE = 40;
  const double Emin = 0;
  const double Emax = 40;

  TH2D * hResstartE0Rec = new TH2D(pretag+"h5ResstartE0RecNOHPRF", selCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartE0Rec);
  TH2D * hResstartE1Rec = new TH2D(pretag+"h5ResstartE1RecNOHPRF", selCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartE1Rec);
  TH2D * hResstartE2Rec = new TH2D(pretag+"h5ResstartE2RecNOHPRF", selCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartE2Rec);
  TH2D * hResstartE3Rec = new TH2D(pretag+"h5ResstartE3RecNOHPRF", selCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartE3Rec);
  TH2D * hResstartE4Rec = new TH2D(pretag+"h5ResstartE4RecNOHPRF", selCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartE4Rec);
  TH2D * hResstartE5Rec = new TH2D(pretag+"h5ResstartE5RecNOHPRF", selCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartE5Rec);
  TH2D * hResstartEsRec[]={hResstartE0Rec, hResstartE1Rec, hResstartE2Rec, hResstartE3Rec, hResstartE4Rec, hResstartE5Rec};

  TH2D * hReslastE0Rec = new TH2D(pretag+"h5ReslastE0RecNOHPRF", selCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hReslastE0Rec);
  TH2D * hReslastE1Rec = new TH2D(pretag+"h5ReslastE1RecNOHPRF", selCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hReslastE1Rec);
  TH2D * hReslastE2Rec = new TH2D(pretag+"h5ReslastE2RecNOHPRF", selCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hReslastE2Rec);
  TH2D * hReslastE3Rec = new TH2D(pretag+"h5ReslastE3RecNOHPRF", selCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hReslastE3Rec);
  TH2D * hReslastE4Rec = new TH2D(pretag+"h5ReslastE4RecNOHPRF", selCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hReslastE4Rec);
  TH2D * hReslastE5Rec = new TH2D(pretag+"h5ReslastE5RecNOHPRF", selCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hReslastE5Rec);
  TH2D * hReslastEsRec[]={hReslastE0Rec, hReslastE1Rec, hReslastE2Rec, hReslastE3Rec, hReslastE4Rec, hReslastE5Rec};

  TH2D * hResstartTruncatedMeanE10Rec = new TH2D(pretag+"h5ResstartTruncatedMeanE10RecNOHPRF", selCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartTruncatedMeanE10Rec);
  TH2D * hResstartTruncatedMeanE20Rec = new TH2D(pretag+"h5ResstartTruncatedMeanE20RecNOHPRF", selCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartTruncatedMeanE20Rec);
  TH2D * hResstartTruncatedMeanE30Rec = new TH2D(pretag+"h5ResstartTruncatedMeanE30RecNOHPRF", selCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartTruncatedMeanE30Rec);
  TH2D * hResstartTruncatedMeanE40Rec = new TH2D(pretag+"h5ResstartTruncatedMeanE40RecNOHPRF", selCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartTruncatedMeanE40Rec);
  TH2D * hResstartTruncatedMeanE50Rec = new TH2D(pretag+"h5ResstartTruncatedMeanE50RecNOHPRF", selCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartTruncatedMeanE50Rec);
  TH2D * hResstartTruncatedMeanEsRec[]={0x0, hResstartTruncatedMeanE10Rec, hResstartTruncatedMeanE20Rec, hResstartTruncatedMeanE30Rec, hResstartTruncatedMeanE40Rec, hResstartTruncatedMeanE50Rec};

  const TString varlastEbase="lastE";
  const TString varstartEbase="startE";
  const TString varstartTruncatedMeanEbase="startTruncatedMeanE";
  const TString xtitlastEbase="end #it{E}";
  const TString xtitstartEbase="start #it{E}";
  const TString xtitstartTruncatedMeanEbase="Truncated Mean (40-95%) of start #it{E} with Nodes 2 to 2 +";
  for(unsigned int ii=0; ii<6; ii++){
    const TString varlastE = varlastEbase+Form("%d",ii);
    const TString varstartE = varstartEbase+Form("%d",ii);
    const TString xtitlastE = xtitlastEbase+Form("_{%d}",ii);
    const TString xtitstartE = xtitstartEbase+Form("_{%d}",ii);

    const TString varResstartE = Form("(%s/%s)-1 : ", varPrec.Data(), varPtrue.Data())+varstartE;
    getHist(varResstartE, xtitstartE+" (MeV/cm)", ytitRes, tree, hResstartEsRec[ii]);

    const TString varReslastE = Form("(%s/%s)-1 : ", varPrec.Data(), varPtrue.Data())+varlastE;
    getHist(varReslastE, xtitlastE+" (MeV/cm)", ytitRes, tree, hReslastEsRec[ii]);

    if(ii){
      const TString varstartTruncatedMeanE = varstartTruncatedMeanEbase+Form("%d",ii*10);
      const TString xtitstartTruncatedMeanE = xtitstartTruncatedMeanEbase+Form(" %d",ii*10);

      const TString varResstartTruncatedMeanE = Form("(%s/%s)-1 : ", varPrec.Data(), varPtrue.Data())+varstartTruncatedMeanE;
      getHist(varResstartTruncatedMeanE, xtitstartTruncatedMeanE+" (MeV/cm)", ytitRes, tree, hResstartTruncatedMeanEsRec[ii]);
    }
  }
  }
}

int main(int argc, char* argv[])
{
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

  const TString inputfn = Form("../anaData/output/outanaData_%s_anaTruth.root", tag.Data());
  TFile *fin = new TFile(inputfn);
  if(!fin->IsOpen()){
    cout<<"fin not open!"<<endl;
    exit(1);
  }
  cout<<"Reading "<<inputfn<<endl;

  drawTracking(lout, tag+"0000raw");

  if(kProton){
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
  drawTracking(lout, tag+"0003TME10", tmecut10);
  drawTracking(lout, tag+"0003TME20", tmecut20);
  }

  style::Process2DHist(lout);
  style::DrawHist(lout, "output", tag, false, true);

  TFile * fout= new TFile(Form("output/outdrawTracking_%s.root", tag.Data()),"recreate");
  lout->Write();
  fout->Save();
  fout->Close();

  return 0;
}
