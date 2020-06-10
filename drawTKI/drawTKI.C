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

TString setX(const TString var, const bool kPiZero, int &nbin, double &xmin, double &xmax, double &lxoff)
{
  TString vn;
  lxoff = 0;
  const double xr = 0.5;
  const TString sfinpi = Form("#pi^{%s}", kPiZero?"0":"+");
  const TString prefix = "p"+sfinpi+"-";

  if(var=="dalphat"){
    nbin = 18;
    xmin = 0;
    xmax = 180;
    vn="#delta#alpha_{T} (degree)";
  }
  else if(var=="dphit"){
    nbin = 18;
    xmin = 0; 
    xmax = 180;
    vn="#delta#phi_{T} (degree)";
    lxoff = xr;
  }
  else if(var=="dpt"){
    nbin = 30;
    xmin = 0;
    xmax = 1.2;
    vn="#delta#it{p}_{T} (GeV/#it{c})";
    lxoff = xr;
  }
  else if(var=="pn"){
    nbin = 30;
    xmin = 0; 
    xmax = 1.2;
    vn="#it{p}_{N} (GeV/#it{c})";
    lxoff = xr;
  }
  else if(var=="iniPimomentum"){
    nbin = 30;
    xmin = 0;
    xmax = 1.3;
    vn="#it{p}_{#pi^{+}}^{in} (GeV/#it{c})";
  }
  else if(var=="finPimomentum"){
    nbin = 20;
    xmin = 0;
    xmax = 1.1;
    vn="#it{p}_{"+sfinpi+"} (GeV/#it{c})";
    lxoff = xr;
  }
  else if(var=="finProtonmomentum"){
    nbin = 30;
    xmin = 0;
    xmax = 1.5;
    vn="#it{p}_{p} (GeV/#it{c})";
    lxoff = xr;
  }
  else if(var=="iniPitheta"){
    nbin = 30;
    xmin = 0;
    xmax = 30;
    vn="#theta_{#pi^{+}}^{in, lab} (degree)";
  }
  else if(var=="finPitheta"){
    nbin = 20;
    xmin = 0;
    xmax = 180;
    vn="#theta_{"+sfinpi+"} (degree)";
    lxoff = xr;
  }
  else if(var=="finProtontheta"){
    nbin = 20;
    xmin = 0;
    xmax = 180;
    vn="#theta_{p} (degree)";
    lxoff = xr;
  }
  else if(var=="fin2Pmom"){
    nbin = 30;
    xmin = 0;
    xmax = 1.5;
    vn="#it{p}_{p}^{subleading} (GeV/#it{c})";
    lxoff = xr;
  }
  else{
    cout<<"setX unknown var! "<<var<<endl;
    exit(1);
  }

  return prefix+vn;
}

TString getUnit(const TString tit)
{
  const TString unit = tit(tit.First("(")+1,tit.Length());
  return unit(0, unit.First(")"));
}

void drawTKI(const TString var, TList *lout, const TString pretag, const bool kPScut, const bool kSLcut, const double ppthres)
{
  const bool kPiZero = pretag.Contains("MPiZero");
  cout<<"\n\n                       Running kPiZero "<<kPiZero<<endl<<endl;
 
  TTree * tree = (TTree*)gDirectory->Get("tree");
  if(!tree){
    cout<<"no tree!"<<endl;
    exit(1);
  }

  int nbin=-999;
  double xmin=-999, xmax=-999, lxoff=-999;
  const TString vn = setX(var, kPiZero, nbin, xmin, xmax, lxoff);
  const TString xunit = getUnit(vn);

  printf("var %s vn %s nbin %d xmin %f xmax %f\n", var.Data(), vn.Data(), nbin, xmin, xmax);

  vector<TString> cuts;//need to match cns
  cuts.push_back("1");
  cuts.push_back("nproton==1 && nneutron==0");
  cuts.push_back("nproton!=1 && nneutron==0");
  cuts.push_back("nproton==1 && nneutron!=0");
  cuts.push_back("nproton!=1 && nneutron!=0");

  TString PScut = Form("&& (finProtonmomentum > %f)", ppthres);
  if(!kPiZero){
    PScut += "&& (finPimomentum > 0.15)";
  }
  const TString SLcut = Form("&& (fin2Pmom < %f)", ppthres);
  vector<TString> cns;
  cns.push_back("all");
  cns.push_back("1p0n");
  cns.push_back("Np0n");
  cns.push_back("1pMn");
  cns.push_back("NpMn");

  THStack * stk = new THStack("s"+var,var); lout->Add(stk);
  TLegend * lg = new TLegend(lxoff+0.15, 0.65, lxoff+0.5, 0.9);
  style::ResetStyle(lg);

  //need system color
  const int col[]={kGray, kRed, kBlue, kOrange, kGreen+3};
  const TString tag = Form("%s_%s_ppthres%.0fPScut%dSLcut%d", var.Data(), pretag.Data(), ppthres*1E3, kPScut, kSLcut);

  double ntotall = -999;
  for(unsigned int ii=0; ii<cuts.size(); ii++){
    const TString ntmp="h"+tag+cns[ii];
    TH1D * htmp = new TH1D(ntmp, ntmp, nbin, xmin, xmax); lout->Add(htmp);
    style::ResetStyle(htmp);
    htmp->SetXTitle(vn);
    htmp->SetYTitle("Events/"+xunit);
    htmp->SetTitle("");

    const TString darg = var+">>"+ntmp;
    TString dcut = "kSignal && "+cuts[ii];
    if(kPScut){
      dcut += PScut;
    }
    if(kSLcut){
      dcut += SLcut;
    }

    const int nev = tree->Draw(darg, dcut);
    printf("ntmp %s ii %d darg %s cut %s nev %d underflow %f overflow %f integral %f\n\n", ntmp.Data(), ii, darg.Data(), dcut.Data(), nev, htmp->GetBinContent(0), htmp->GetBinContent(nbin+1), htmp->Integral(0, 1000));

    htmp->Scale(1,"width");
    if(cns[ii]!="all"){
      htmp->SetFillColor(col[ii]);
      stk->Add(htmp);
      lg->AddEntry(htmp, cns[ii], "f");
    }
    else{
      ntotall = htmp->Integral(0,1000,"width");
    }
  }

  TCanvas * c1=new TCanvas("c1","", 600, 400);
  style::PadSetup(c1);
  gPad->SetTopMargin(0.07);
  gPad->SetRightMargin(0.03);
  gPad->SetBottomMargin(0.15);
  gStyle->SetOptTitle(1);

  style::ResetStyle(stk);
  stk->Draw();
  stk->GetXaxis()->SetTitle(vn);
  stk->GetYaxis()->SetTitle("Event/"+xunit);
  stk->SetTitle(Form("Total %.0f events", ntotall));
  stk->SetMinimum(0);
  stk->Draw();

  TString lheader("full phase space");
  if(kPScut){
    lheader = Form("#it{p}_{p}>%.2f", ppthres);
    if(!kPiZero){
      lheader += ", #it{p}_{#pi^{+}}>0.15";
    }
  }
  if(kSLcut){
    lheader += Form(", #it{p}_{p}^{s.l.}<%.2f", ppthres);
  }
  lg->SetHeader(lheader);
  lg->Draw();

  TString ptag = tag;
  ptag.ReplaceAll("_TrackingProton","");
  c1->Print("output/"+ptag+".png");
  c1->Print("output/"+ptag+".eps");
  c1->Print("output/"+ptag+".pdf");
}

int main(int argc, char* argv[])
{
  if(argc!=2){
    cout<<"argc!=2 !"<<argc<<endl;
    return 0;
  }

  gSystem->Load("libTree");

  style::SetGlobalStyle();

  const bool kPiZero = atoi(argv[1]);
  const TString tag = kPiZero?"MPiZero":"1PiPlus";

  TList * lout= new TList;

  TFile *fin = new TFile(Form("../anaData/output/outanaData_%s_TrackingProton.root", tag.Data()));
  if(!fin->IsOpen()){
    cout<<"fin not open!"<<endl;
    exit(1);
  }

  const TString vars[]={"dalphat","dphit","dpt","pn", "iniPimomentum", "finPimomentum", "finProtonmomentum", "iniPitheta", "finPitheta", "finProtontheta", "fin2Pmom"};
  for(unsigned int ii=0; ii<sizeof(vars)/sizeof(TString); ii++){
    //0.45 -> 100 MeV K.E.
    drawTKI(vars[ii], lout, tag,  0, 0, 0.45);
    drawTKI(vars[ii], lout, tag,  1, 0, 0.45);
    drawTKI(vars[ii], lout, tag,  1, 1, 0.45);
    //0.25 -> 30 MeV K.E.
    drawTKI(vars[ii], lout, tag,  1, 0, 0.25);
    drawTKI(vars[ii], lout, tag,  1, 1, 0.25);
  }

  TFile * fout= new TFile(Form("output/outdrawTKI_%s.root", tag.Data()),"recreate");
  lout->Write();
  fout->Save();
  fout->Close();

  return 0;
}
