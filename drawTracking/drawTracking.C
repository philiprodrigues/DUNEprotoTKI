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

TH1D * getEff(const TH1D * ha, const TH1D *hb)
{
  const TString hname = ha->GetName();
  TH1D * heff = (TH1D*) ha->Clone(Form("%seff", hname.Data()));
  heff->Sumw2();
  heff->Divide(ha, hb, 1, 1, "B");
  heff->SetMinimum(0);
  heff->SetMaximum(hname.Contains("sig")?0.2:0.7);//1.1);
  heff->SetTitle(Form("%s/%s", ha->GetTitle(), hb->GetTitle()));
  heff->SetYTitle("Eff.");

  return heff;
}

void getHist(const TString var, const TString xtit, const TString ytit, TTree *tree, TH1 *hist)
{
  const TString cut = hist->GetTitle();
  const TString hname = hist->GetName();
  printf("Drawing %s %s %s\n", var.Data(), hname.Data(), cut.Data());
  tree->Draw(var+">>"+hname, cut);

  hist->SetXTitle(xtit);
  hist->SetYTitle(ytit);
}

void drawHist(TList* lout)
{
  TCanvas * c1=new TCanvas("c1","", 600, 400);
  style::PadSetup(c1);
  gPad->SetTopMargin(0.07);
  gPad->SetRightMargin(0.1);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(1);
  gStyle->SetStatY(0.9);

  for(int ii=0; ii<lout->GetSize(); ii++){
    TH1 * hh=dynamic_cast<TH1*> (lout->At(ii));
    if(!hh){
      printf("not a TH1 at %d!\n", ii); exit(1);
    }
    style::ResetStyle(hh);

    TH2D * htmp = dynamic_cast<TH2D *>(hh);
    const bool k2d = htmp;

    const TString hname = hh->GetName();
    if(hname.Contains("eff")||k2d){
      gStyle->SetOptStat(0);
    }
    else{
      gStyle->SetOptStat("enuomr");
      gStyle->SetStatColor(0);
      gStyle->SetStatStyle(0);
      gStyle->SetStatY(0.9);
    }

    hh->UseCurrentStyle();
    gPad->SetGrid(k2d, k2d);

    hh->SetLineColor(kBlack);
    hh->SetLineWidth(2);
    TString dopt("hist");
    if(k2d){
      dopt="colz";
    }
    if(hname.Contains("eff")){
      dopt="E hist";
      hh->SetMarkerStyle(24);
      hh->SetMarkerSize(2);
      hh->SetMarkerColor(kRed);
    }
    /*
    TString hn = hh->GetName();
    hn.ReplaceAll("_TrackingProton","");
    hh->SetName(hn);
    */
    hh->Draw(dopt);

    const TString ptag = hh->GetName();
    c1->Print(Form("output/%s.png", ptag.Data()));
    c1->Print(Form("output/%s.pdf", ptag.Data()));
    c1->Print(Form("output/%s.eps", ptag.Data()));
  }
}

void drawTracking(TList *lout, const TString pretag, const TString ESCcut="")
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

  const TString nclCut = varPrec+"!=-999 && ndEdxCls>=6" + ESCcut;

  //========================================================================
  const int nPtrue = 20;
  const double Ptruemin = 0;
  const double Ptruemax = 1.2;
  TH1D * hPtrueall = new TH1D(pretag+"h0all","1",nPtrue, Ptruemin, Ptruemax); lout->Add(hPtrueall);
  TH1D * hPtruerec = new TH1D(pretag+"h1rec",varPrec+"!=-999",nPtrue, Ptruemin, Ptruemax); lout->Add(hPtruerec);
  TH1D * hPtruencl = new TH1D(pretag+"h2ncl",nclCut,nPtrue, Ptruemin, Ptruemax); lout->Add(hPtruencl);
  TH1D * hPtruesig = new TH1D(pretag+"h3sig",nclCut+"&& kSignal",nPtrue, Ptruemin, Ptruemax); lout->Add(hPtruesig);

  const TString trackingparticle = kTrackingProton?"p":"#pi";
  const TString varPtrue = kTrackingProton?"finProtonmomentum":"finPimomentum";
  const TString ytitN = "N";
  const TString xtitPtrue = Form("#it{p}_{%s}^{true} (GeV/#it{c})", trackingparticle.Data());
  getHist(varPtrue, xtitPtrue, ytitN, tree, hPtrueall);
  getHist(varPtrue, xtitPtrue, ytitN, tree, hPtruerec);
  getHist(varPtrue, xtitPtrue, ytitN, tree, hPtruencl);
  getHist(varPtrue, xtitPtrue, ytitN, tree, hPtruesig);

  TH1D * effPtruerec = getEff(hPtruerec, hPtrueall); lout->Add(effPtruerec);
  TH1D * effPtruencl = getEff(hPtruencl, hPtrueall); lout->Add(effPtruencl);
  TH1D * effPtruesig = getEff(hPtruesig, hPtrueall); lout->Add(effPtruesig);

  //========================================================================
  const int nRes = 15;
  const double Resmin = -1;
  const double Resmax = 1;

  TH2D * hResPtrueRec = new TH2D(pretag+"hResPtrueRec",nclCut, nPtrue, Ptruemin, Ptruemax, nRes, Resmin, Resmax); lout->Add(hResPtrueRec);

  const TString varResPtrue = Form("(%s/%s)-1 : %s", varPrec.Data(), varPtrue.Data(), varPtrue.Data());
  const TString ytitRes = Form("#it{p}_{%s}^{rec}/#it{p}_{%s}^{true}-1", trackingparticle.Data(), trackingparticle.Data());
  getHist(varResPtrue, xtitPtrue, ytitRes, tree, hResPtrueRec);

  TH1D * hResPtrueRecpdf = 0x0; 
  TH1D * hResPtrueReccdf = 0x0; 
  TH2D * hResPtrueRecNor = style::NormalHist(hResPtrueRec, hResPtrueRecpdf, hResPtrueReccdf, 0, true); lout->Add(hResPtrueRecNor); lout->Add(hResPtrueRecpdf); lout->Add(hResPtrueReccdf);

  TH1D * hResRecSignal = new TH1D(pretag+"hResRecSignal",nclCut+"&& kSignal", nRes*3, Resmin, Resmax); lout->Add(hResRecSignal);
  TH1D * hResRecBeam = new TH1D(pretag+"hResRecBeam",nclCut, nRes*3, Resmin, Resmax); lout->Add(hResRecBeam);
  const TString varRes = Form("(%s/%s)-1", varPrec.Data(), varPtrue.Data());
  const TString xtitRes = ytitRes;
  getHist(varRes, xtitRes, ytitN, tree, hResRecSignal);
  getHist(varRes, xtitRes, ytitN, tree, hResRecBeam);

  //========================================================================
  const int nChi2 = 20;
  const double Chi2min = 0;
  const double Chi2max = 100;

  TH1D * hChi2 = new TH1D(pretag+"hChi2", nclCut, nChi2, Chi2min, Chi2max); lout->Add(hChi2);
  const TString varChi2 = "chi2/ndof";
  const TString xtitChi2 = "#chi^{2}/NDF";
  getHist(varChi2, xtitChi2, ytitN, tree, hChi2);

  TH2D * hResChi2 = new TH2D(pretag+"hResChi2", nclCut, nChi2, Chi2min, Chi2max, nRes, Resmin, Resmax); lout->Add(hResChi2);
  const TString varResChi2 = Form("(%s/%s)-1 : chi2/ndof", varPrec.Data(), varPtrue.Data());
  getHist(varResChi2, xtitChi2, ytitRes, tree, hResChi2);
  TH1D * hResChi2pdf = 0x0;
  TH1D * hResChi2cdf = 0x0;
  TH2D * hResChi2Nor = style::NormalHist(hResChi2, hResChi2pdf, hResChi2cdf, 0, true); lout->Add(hResChi2Nor); lout->Add(hResChi2pdf); lout->Add(hResChi2cdf);
  //========================================================================
  const int nE = 40;
  const double Emin = 0;
  const double Emax = 40;
  TH1D * hlastE0 = new TH1D(pretag+"h4lastE0", nclCut, nE, Emin, Emax); lout->Add(hlastE0);
  TH1D * hlastE1 = new TH1D(pretag+"h4lastE1", nclCut, nE, Emin, Emax); lout->Add(hlastE1);
  TH1D * hlastE2 = new TH1D(pretag+"h4lastE2", nclCut, nE, Emin, Emax); lout->Add(hlastE2);
  TH1D * hlastE3 = new TH1D(pretag+"h4lastE3", nclCut, nE, Emin, Emax); lout->Add(hlastE3);
  TH1D * hlastE4 = new TH1D(pretag+"h4lastE4", nclCut, nE, Emin, Emax); lout->Add(hlastE4);
  TH1D * hlastE5 = new TH1D(pretag+"h4lastE5", nclCut, nE, Emin, Emax); lout->Add(hlastE5);
  TH1D * hlastEs[]={hlastE0, hlastE1, hlastE2, hlastE3, hlastE4, hlastE5};

  TH1D * hstartE0 = new TH1D(pretag+"h4startE0", nclCut, nE, Emin, Emax); lout->Add(hstartE0);
  TH1D * hstartE1 = new TH1D(pretag+"h4startE1", nclCut, nE, Emin, Emax); lout->Add(hstartE1);
  TH1D * hstartE2 = new TH1D(pretag+"h4startE2", nclCut, nE, Emin, Emax); lout->Add(hstartE2);
  TH1D * hstartE3 = new TH1D(pretag+"h4startE3", nclCut, nE, Emin, Emax); lout->Add(hstartE3);
  TH1D * hstartE4 = new TH1D(pretag+"h4startE4", nclCut, nE, Emin, Emax); lout->Add(hstartE4);
  TH1D * hstartE5 = new TH1D(pretag+"h4startE5", nclCut, nE, Emin, Emax); lout->Add(hstartE5);
  TH1D * hstartEs[]={hstartE0, hstartE1, hstartE2, hstartE3, hstartE4, hstartE5};

  TH1D * hstartTruncatedMeanE10 = new TH1D(pretag+"h4startTruncatedMeanE10", nclCut, nE, Emin, Emax); lout->Add(hstartTruncatedMeanE10);
  TH1D * hstartTruncatedMeanE20 = new TH1D(pretag+"h4startTruncatedMeanE20", nclCut, nE, Emin, Emax); lout->Add(hstartTruncatedMeanE20);
  TH1D * hstartTruncatedMeanE30 = new TH1D(pretag+"h4startTruncatedMeanE30", nclCut, nE, Emin, Emax); lout->Add(hstartTruncatedMeanE30); 
  TH1D * hstartTruncatedMeanE40 = new TH1D(pretag+"h4startTruncatedMeanE40", nclCut, nE, Emin, Emax); lout->Add(hstartTruncatedMeanE40);
  TH1D * hstartTruncatedMeanE50 = new TH1D(pretag+"h4startTruncatedMeanE50", nclCut, nE, Emin, Emax); lout->Add(hstartTruncatedMeanE50);
  TH1D * hstartTruncatedMeanEs[]={0x0, hstartTruncatedMeanE10, hstartTruncatedMeanE20, hstartTruncatedMeanE30, hstartTruncatedMeanE40, hstartTruncatedMeanE50};

  TH2D * hResstartE0Rec = new TH2D(pretag+"h5ResstartE0Rec", nclCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartE0Rec);
  TH2D * hResstartE1Rec = new TH2D(pretag+"h5ResstartE1Rec", nclCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartE1Rec);
  TH2D * hResstartE2Rec = new TH2D(pretag+"h5ResstartE2Rec", nclCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartE2Rec);
  TH2D * hResstartE3Rec = new TH2D(pretag+"h5ResstartE3Rec", nclCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartE3Rec);
  TH2D * hResstartE4Rec = new TH2D(pretag+"h5ResstartE4Rec", nclCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartE4Rec);
  TH2D * hResstartE5Rec = new TH2D(pretag+"h5ResstartE5Rec", nclCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartE5Rec);
  TH2D * hResstartEsRec[]={hResstartE0Rec, hResstartE1Rec, hResstartE2Rec, hResstartE3Rec, hResstartE4Rec, hResstartE5Rec};

  TH2D * hResstartTruncatedMeanE10Rec = new TH2D(pretag+"h5ResstartTruncatedMeanE10Rec", nclCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartTruncatedMeanE10Rec);
  TH2D * hResstartTruncatedMeanE20Rec = new TH2D(pretag+"h5ResstartTruncatedMeanE20Rec", nclCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartTruncatedMeanE20Rec);
  TH2D * hResstartTruncatedMeanE30Rec = new TH2D(pretag+"h5ResstartTruncatedMeanE30Rec", nclCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartTruncatedMeanE30Rec);
  TH2D * hResstartTruncatedMeanE40Rec = new TH2D(pretag+"h5ResstartTruncatedMeanE40Rec", nclCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartTruncatedMeanE40Rec);
  TH2D * hResstartTruncatedMeanE50Rec = new TH2D(pretag+"h5ResstartTruncatedMeanE50Rec", nclCut, nE, Emin, Emax, nRes, Resmin, Resmax); lout->Add(hResstartTruncatedMeanE50Rec);
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

    getHist(varlastE, xtitlastE+" (MeV/cm)", ytitN, tree, hlastEs[ii]);
    getHist(varstartE, xtitstartE+" (MeV/cm)", ytitN, tree, hstartEs[ii]);

    const TString varResstartE = Form("(%s/%s)-1 : ", varPrec.Data(), varPtrue.Data())+varstartE;
    getHist(varResstartE, xtitstartE+" (MeV/cm)", ytitRes, tree, hResstartEsRec[ii]);

    TH1D * hResstartERecpdf = 0x0;
    TH1D * hResstartEReccdf = 0x0;
    TH2D * hResstartERecNor = style::NormalHist(hResstartEsRec[ii], hResstartERecpdf, hResstartEReccdf, 0, true); lout->Add(hResstartERecNor); lout->Add(hResstartERecpdf); lout->Add(hResstartEReccdf);

    if(ii){
      const TString varstartTruncatedMeanE = varstartTruncatedMeanEbase+Form("%d",ii*10);
      const TString xtitstartTruncatedMeanE = xtitstartTruncatedMeanEbase+Form(" %d",ii*10);

      getHist(varstartTruncatedMeanE, xtitstartTruncatedMeanE+" (MeV/cm)", ytitN, tree, hstartTruncatedMeanEs[ii]);

      const TString varResstartTruncatedMeanE = Form("(%s/%s)-1 : ", varPrec.Data(), varPtrue.Data())+varstartTruncatedMeanE;
      getHist(varResstartTruncatedMeanE, xtitstartTruncatedMeanE+" (MeV/cm)", ytitRes, tree, hResstartTruncatedMeanEsRec[ii]);

      TH1D * hResstartTruncatedMeanERecpdf = 0x0;
      TH1D * hResstartTruncatedMeanEReccdf = 0x0;
      TH2D * hResstartTruncatedMeanERecNor = style::NormalHist(hResstartTruncatedMeanEsRec[ii], hResstartTruncatedMeanERecpdf, hResstartTruncatedMeanEReccdf, 0, true); lout->Add(hResstartTruncatedMeanERecNor); lout->Add(hResstartTruncatedMeanERecpdf); lout->Add(hResstartTruncatedMeanEReccdf);
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

  const TString chi2cut = "&& (chi2/ndof<15)";
  const TString esc2cut = "&& (startE2>10)";
  const TString esc3cut = "&& (startE3>7)";
  const TString tme10cut = "&& (startTruncatedMeanE10>6)";
  drawTracking(lout, tag+"0000raw");
  drawTracking(lout, tag+"0001CHIcut",   chi2cut);
  drawTracking(lout, tag+"0002ESC2",     esc2cut);
  drawTracking(lout, tag+"0004TME10",    tme10cut);
  drawTracking(lout, tag+"0012CHIESC2",  chi2cut+esc2cut);
  drawTracking(lout, tag+"0123CHIESC23", chi2cut+esc2cut+esc3cut);

  drawHist(lout);

  TFile * fout= new TFile(Form("output/outdrawTracking_%s.root", tag.Data()),"recreate");
  lout->Write();
  fout->Save();
  fout->Close();

  return 0;
}
