/******************************************************
                                                      *
                                                      *
                  Author: Xianguo Lu                  *
                                                      *
                  xianguo.lu@physics.ox.ac.uk         *
                  lu.xianguo@gmail.com                *
                                                      *
                                                      *
                                                      *
*******************************************************/

#include "style.h"

Double_t style::fgkTextSize = 0.05;
Double_t style::fgkTitleSize = 0.05;
Double_t style::fgkMarkerSize = 1;
Double_t style::fgkLineWidth = 2;
Int_t style::fgkTextFont = 42;
Double_t style::fgkLabelOffset = 0.01;
Double_t style::fgkXTitleOffset = 1.25;//1.1;//1.25;
Double_t style::fgkYTitleOffset = 1.1;//1.2;
Double_t style::fgkTickLength = 0.02;

//ClassImp(style);

#ifndef EPSILON
#define EPSILON 1E-12
#endif

const TString gTagPRF="PRF";
const TString gTagPUR="PUR";
const TString gTagEFF="EFF";
const TString gTagNOH="NOH";
const TString gTagSTK="STK";

void style::Process2DHist(TList *lout)
{
  const int nhist = lout->GetSize();//the size will increase
  for(int ii=0; ii<nhist; ii++){
    TH2D * htmp = dynamic_cast<TH2D *>( lout->At(ii) );
    if(htmp){
      const TString tag = htmp->GetName();

      if(tag.Contains(gTagPRF)){
        TH1D * hpro = htmp->ProfileX(Form("%s_profileX", htmp->GetName())); lout->Add(hpro);
        if(tag.Contains(gTagPUR)){
          hpro->SetYTitle("Purity");
        }
        else if(tag.Contains(gTagEFF)){
          hpro->SetYTitle("Efficiency");
        }

        hpro->SetTitle(tag);
      }

      if(tag.Contains(gTagNOH)){
        TH1D * hpdf = 0x0;
        TH1D * hcdf = 0x0;
        const double thres = 5;
        TH2D * hnor = NormalHist(htmp, hpdf, hcdf, thres, true);
        hpdf->SetYTitle("p.d.f.");
        hcdf->SetYTitle("c.d.f.");
        hpdf->SetTitle(tag);
        hcdf->SetTitle(tag);
        hnor->SetTitle(tag);
        lout->Add(hpdf);
        lout->Add(hcdf);
        lout->Add(hnor); 
      }
    }
  }
}

void style::DrawHist(TList *lout, const TString outdir, const TString tag, const bool kfast)
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
    gStyle->SetOptTitle(1);
    if(k2d){
      //gStyle->SetOptStat(0);
      gStyle->SetOptStat("eou");
    }
    else{
      gStyle->SetOptStat("eoum");
      gStyle->SetStatColor(0);
      gStyle->SetStatStyle(0);
      gStyle->SetStatY(0.9);
    }

    const TString tag = hh->GetName();
    cout<<"Printing "<<tag<<endl;

    style::ResetStyle(hh);
    hh->UseCurrentStyle();//can't go after setmarkersize
    hh->SetMaximum(hh->GetBinContent(hh->GetMaximumBin())*1.3);
    hh->SetMinimum(0);
    hh->UseCurrentStyle();
    hh->SetMarkerSize(3);
    TString dopt="text hist";
    if(k2d){
      if(tag.Contains(gTagPRF)){
        dopt = "box";
      }
      else if(tag.Contains(gTagNOH)){
        dopt = "colz";
      }
    }
    else{
      if(tag.Contains(gTagPRF)){
        //hh->SetMaximum(0.1);
        //hh->SetMaximum(0.3);
        dopt = "hist";
      }
      else if(tag.Contains(gTagNOH)){
        dopt = "hist";
      }
    }
    hh->Draw(dopt);

    c1->Print(outdir+"/"+tag+".png");

    if(!kfast){
      c1->Print(outdir+"/"+tag+".pdf");
      c1->Print(outdir+"/"+tag+".eps");
    }
  }

}

TH1D * style::GetCDF(const TH2D *hraw, const TString hname)
{
  TH1D * hcdf = hraw->ProjectionX(hname);
  TH1D * hold = (TH1D*) hcdf->Clone("hold");

  hcdf->Scale(0);

  const Int_t n2 = hold->GetNbinsX()+1;
  const Double_t ntot = hold->Integral(0, n2);

  for(Int_t ii=0; ii<= n2; ii++){
    const Double_t num = hold->Integral(0, ii);
    if(num<EPSILON){
      continue;
    }

    const Double_t nee = TMath::Sqrt(num);

    const Double_t fra = num/(ntot+EPSILON);
    const Double_t fee = nee/(ntot+EPSILON);

    hcdf->SetBinContent(ii, fra);
    hcdf->SetBinError(ii, fee);
  }

  delete hold;

  return hcdf;
}

TH1D* style::ToPDF(const TH1 *hraw, const TString hn)
{
  const Int_t x0 = 0;
  const Int_t x1 = hraw->GetNbinsX()+1;
  const Double_t tmpnt = hraw->Integral(x0, x1);

  TH1D * hist = (TH1D*) hraw->Clone((hn+hraw->GetName())+"pdf");
  hist->Scale(0);

  for(Int_t ib=x0; ib<=x1; ib++){
    const Double_t bw = hraw->GetBinWidth(ib);
    const Double_t cont = hraw->GetBinContent(ib);
    if(cont<EPSILON)
      continue;

    //in case of finit number of bins (i.e. eff not always small), Binomial error is more accurate than Poisson error

    const Double_t eff = cont/tmpnt;
    const Double_t pdf = eff/bw;

    const Double_t dpdf = sqrt(eff*(1-eff)/tmpnt) / bw;
    hist->SetBinContent(ib, pdf);
    hist->SetBinError(ib, dpdf);
  }

  hist->SetEntries(tmpnt);

  return hist;
}


TH2D* style::NormalHist(const TH2D *hraw, TH1D * &hpdf, TH1D * &hcdf,  const Double_t thres, const Bool_t kmax)
{
  hpdf = hraw->ProjectionX(Form("%spdf",hraw->GetName()));

  hcdf = GetCDF(hraw, Form("%scdf",hraw->GetName()));

  TH2D *hh=(TH2D*)hraw->Clone(Form("%snor",hraw->GetName()));
  hh->Scale(0);

  const Int_t x0 = hh->GetXaxis()->GetFirst();
  const Int_t x1 = hh->GetXaxis()->GetLast();
  const Int_t y0 = hh->GetYaxis()->GetFirst();
  const Int_t y1 = hh->GetYaxis()->GetLast();

  Double_t hmax = -1e10;
  Double_t hmin = 1e10;
  Double_t nent = 0;
  for(Int_t ix=x0; ix<=x1; ix++){

    //if option "e" is specified, the errors are computed. if option "o" original axis range of the taget axes will be kept, but only bins inside the selected range will be filled.

    TH1D * sliceh = hraw->ProjectionY(Form("tmpnormalhist%sx%d", hh->GetName(), ix), ix, ix, "oe");
    const Double_t tot = sliceh->GetEntries();

    TH1D * pdfh=0x0;


    if(tot>EPSILON){
      nent += tot;

      Double_t imax = -999;

      if(!kmax){
        pdfh = ToPDF(sliceh,"tmp");
      }
      else{
        imax = sliceh->GetBinContent(sliceh->GetMaximumBin());
      }

      for(Int_t iy=y0; iy<=y1; iy++){
        const Double_t cont = kmax ? sliceh->GetBinContent(iy)/imax : pdfh->GetBinContent(iy);
        const Double_t ierr = kmax ? sliceh->GetBinError(iy)/imax   : pdfh->GetBinError(iy);
        if(tot>thres && cont>0){
          hh->SetBinContent(ix, iy, cont);
          hh->SetBinError(ix,iy, ierr);
          if(cont>hmax) hmax = cont;
          if(cont<hmin) hmin = cont;
        }
      }
    }


    delete pdfh;
    delete sliceh;
  }

  hh->SetEntries(nent);
  hh->SetMinimum(0.99*hmin);
  hh->SetMaximum(1.1*hmax);

  TString xtit(hraw->GetXaxis()->GetTitle());
  if(xtit.Contains("(")){
    xtit=xtit(0, xtit.First('('));
  }

  TString ytit(hraw->GetYaxis()->GetTitle());
  if(ytit.Contains("(")){
    ytit=ytit(0, ytit.First('('));
  }

  hh->SetTitle(Form("f(%s|%s) %s", ytit.Data(), xtit.Data(), hraw->GetTitle()));

  return hh;
}



void style::PrintHist(const TH1D *hh, const TString tag)
{
  printf("\n");
  for(int ii=0; ii<=hh->GetNbinsX()+1; ii++){
    printf("%s bin %d x %f %f content %e error %e\n", tag.Data(), ii, hh->GetXaxis()->GetBinLowEdge(ii), hh->GetXaxis()->GetBinUpEdge(ii), hh->GetBinContent(ii), hh->GetBinError(ii));
  }
  printf("\n");
}

void style::GetHistB0L(const TH1D * h0, const double thres,  int & b0,  int & ll)
{
  b0=0;
  int b1=h0->GetNbinsX()+1;//do use overflow bin

  for(int ii=b0; ii<=b1; ii++){
    const double yy = h0->GetBinContent(ii);
    if(yy>thres){
      b0=ii;
      break;
    }
  }

  for(int ii=b1; ii>=b0; ii--){
    const double yy = h0->GetBinContent(ii);
    if(yy>thres){
      b1=ii;
      break;
    }
  }

  ll=b1-b0+1;
}

TMatrixDSym *style::ExtractShapeOnlyCovar(const TMatrixDSym *full_covar, const TH1D *data_hist, TH1D * &outhist) 
{
  //https://github.com/NUISANCEMC/nuisance/blob/master/src/Statistical/StatUtils.cxx#L1077

  printf("\n\n========================= style::ExtractShapeOnlyCovar ======================\n\n");
  PrintHist(data_hist, data_hist->GetName());
  full_covar->Print();

  if(outhist!=0x0){
    printf("outhist not null!!\n"); exit(1);
  }

  outhist = (TH1D*) data_hist->Clone(Form("shapeonly%s", data_hist->GetName()));

  const double data_scale = 1;

  int nbins = full_covar->GetNrows();
  TMatrixDSym *shape_covar = new TMatrixDSym(nbins);

  // Check nobody is being silly
  if (data_hist->GetNbinsX() != nbins) {
    cout<<"wrong nbins "<<nbins<<" "<<data_hist->GetNbinsX()<<endl; exit(1);
  }

  double total_data = 0;
  double total_covar = 0;

  // Initial loop to calculate some constants
  for (int i = 0; i < nbins; ++i) {
    total_data += data_hist->GetBinContent(i + 1) * data_scale;
    for (int j = 0; j < nbins; ++j) {
      total_covar += (*full_covar)(i, j);
    }
  }

  if (total_data == 0 || total_covar == 0) {
    cout<<"wrong data or covar"<<endl; exit(1);
  }

  cout<<"Norm error = " << (sqrt(total_covar) / total_data)<<endl;;

  // Now loop over and calculate the shape-only matrix
  for (int i = 0; i < nbins; ++i) {
    double data_i = data_hist->GetBinContent(i + 1) * data_scale;

    for (int j = 0; j < nbins; ++j) {
      double data_j = data_hist->GetBinContent(j + 1) * data_scale;

      double norm_term =
          data_i * data_j * total_covar / total_data / total_data;
      double mix_sum1 = 0;
      double mix_sum2 = 0;

      for (int k = 0; k < nbins; ++k) {
        mix_sum1 += (*full_covar)(k, j);
        mix_sum2 += (*full_covar)(i, k);
      }

      double mix_term1 =
          data_i * (mix_sum1 / total_data -
                    total_covar * data_j / total_data / total_data);
      double mix_term2 =
          data_j * (mix_sum2 / total_data -
                    total_covar * data_i / total_data / total_data);

      (*shape_covar)(i, j) =
          (*full_covar)(i, j) - mix_term1 - mix_term2 - norm_term;
    }
  }

  for(int ii=0; ii<nbins; ii++){
    const double errNew = sqrt((*shape_covar)(ii,ii));
    outhist->SetBinError(ii+1, errNew);
  }

  PrintHist(outhist,outhist->GetName());

  return shape_covar;
}

void style::AddStackToLegend(THStack *stk, TLegend *lg, const TString tag, const double drawthres)
{
  TList * hgstack = stk->GetHists();

  double totcount = 0;
  for(int ii=0; ii<hgstack->GetEntries(); ii++){
    TH1D * hh=(TH1D*) hgstack->At(ii);
    totcount += hh->Integral(0,1000, "width");
  }

  for(int ii=0; ii<hgstack->GetEntries(); ii++){
    TH1D * hh=(TH1D*) hgstack->At(ii);
    const Color_t col=hh->GetFillColor();
    hh->SetLineColor(col);
    const double modCount =  hh->Integral(0,1000, "width");
    const double modRatio = modCount/totcount;
    if( modRatio > drawthres ){
      lg->AddEntry(hh, hh->GetTitle()+tag, "f");
    }
    printf("style::AddStackToLegend ii %d %s %e %e %e\n", ii, hh->GetName(), modCount, totcount, modRatio);
  }
}

TMatrixD * style::Hist2Matrix(const TH2D * hist, const double unit)
{
  //hist symmetric
  const int ndim = hist->GetNbinsX();
  TMatrixD *cov=new TMatrixD(ndim+2, ndim+2);//include overflow and underflow
  for(int ix=1; ix<=ndim; ix++){
    for(int jy=1; jy<=ndim; jy++){
      //same index in matrix and hist
      (*cov)[ix][jy] = hist->GetBinContent(ix, jy)/unit;
    }
  }

  return cov;
}

void style::ScaleXaxis(THStack * stk, const double scale)
{
  //THStack * sout = new THStack;
  const TList * ll = stk->GetHists();
  for(Int_t ii=0; ii<ll->GetEntries(); ii++){
    TH1D *hh = (TH1D*)ll->At(ii);
    ScaleXaxis(hh, scale);
    //sout->Add((TH1D*)hh->Clone(hh->GetName()));
  }

  //delete stk;
  //stk = sout;
}

void style::ScaleXaxis(TH1D * h0, const double scale)
{
  PrintHist(h0, "style::ScaleXaxis before: ");

  const int nbin = h0->GetNbinsX();
  double xxs[100]={0};
  double yys[100]={0};
  double ers[100]={0};
  for(int ii=1; ii<=nbin+1; ii++){
    const double iy = h0->GetBinContent(ii)/scale;
    const double ie = h0->GetBinError(ii)/scale;
    const double ix = h0->GetXaxis()->GetBinLowEdge(ii)*scale;
    xxs[ii-1]=ix;
    yys[ii]=iy;
    ers[ii]=ie;
    printf("style::ScaleXaxis test1 %s %d %f %e %e\n", h0->GetName(), ii, ix, iy, ie);
  }

  /*
  const TString hname = h0->GetName();
  const double xmin =  h0->GetXaxis()->GetBinLowEdge(1);
  const double xmax =  h0->GetXaxis()->GetBinLowEdge(nbin+1);
  
  TH1D *hh= new TH1D(hname,"", nbin, xmin, xmax);
  */
  h0->SetBins(nbin, xxs);

  for(int ii=1; ii<=nbin+1; ii++){
    h0->SetBinContent(ii, yys[ii]);
    h0->SetBinError(ii, ers[ii]);
    printf("style::ScaleXaxis test2 %d %s %f %e %e\n", ii, h0->GetName(), h0->GetXaxis()->GetBinLowEdge(ii), yys[ii], ers[ii]);
  }

  PrintHist(h0, "style::ScaleXaxis after: ");
}

TH1D * style::GetMCDataNsigma(TH1D * hmc, TH1D * hdata, const TString tag)
{
  TH1D * hnsig = (TH1D*) hdata->Clone(Form("%sration", hdata->GetName()));
  hnsig->Clear();
  hnsig->SetTitle(tag);

  for(int ii=1; ii<=hnsig->GetNbinsX(); ii++){
    const double errordata = hdata->GetBinError(ii);
    if(errordata<1E-39*EPSILON){
      continue;
    }

    const double dataxx = hdata->GetBinCenter(ii);

    const int mcbin = hmc->GetXaxis()->FindBin(dataxx);
    const double mcxx = hmc->GetBinCenter(mcbin);
    if(fabs(dataxx/mcxx-1)>0.05){
      printf("style::GetMCDataNsigma %s ii %d dataxx %f mcxx %f %f\n", hdata->GetName(), ii, dataxx, mcxx, fabs(dataxx/mcxx-1));
      continue;
    }

    const double countdata = hdata->GetBinContent(ii);

    const double countmc = hmc->GetBinContent(mcbin);

    const double nsigma = (countmc-countdata)/errordata;

    hnsig->SetBinContent(ii, nsigma);
    printf("style::GetMCDataNsigma %s ii %d dataxx %f data count %e error %e mc count %e nsigma %f\n", hdata->GetName(), ii, dataxx, countdata, errordata, countmc, nsigma);
  }

  return hnsig;
}

void style::SetStackColorCB(THStack * hh, const int *cbcol, const int lw)
{
  TList * hgstack = hh->GetHists();
  for(int ii=0; ii<hgstack->GetEntries(); ii++){
    TH1D * tmph=(TH1D*) hgstack->At(ii);
    const int tmpcol = cbcol[ii];
    printf("style::SetStackColorCB ii %d col %d\n", ii, tmpcol);
    tmph->SetFillColor(tmpcol);
    tmph->SetLineColor(tmpcol);
    tmph->SetLineWidth(lw);
  }
}

void style::IniColorCB()
{
  //http://www.somersault1824.com/tips-for-designing-scientific-figures-for-color-blind-readers/
  const Int_t ibase=1001;
  Int_t id=ibase;
  new TColor(id++, 0./255., 0./255., 0./255., "CB1_Black",1.0);
  new TColor(id++, 0./255., 73./255., 73./255., "CB2_Forest",1.0);
  new TColor(id++, 0./255., 146./255., 146./255., "CB3_Teal",1.0);
  new TColor(id++, 255./255., 109./255., 182./255., "CB4_HotPink",1.0);
  new TColor(id++, 255./255., 182./255., 119./255., "CB5_BabyPink",1.0);
  new TColor(id++, 73./255., 0./255., 146./255., "CB6_Purple",1.0);
  new TColor(id++, 0./255., 109./255., 219./255., "CB7_RoyalBlue",1.0);
  new TColor(id++, 182./255., 109./255., 255./255., "CB8_Lilac",1.0);
  new TColor(id++, 109./255., 182./255., 255./255., "CB9_BlueGrey",1.0);
  new TColor(id++, 182./255., 219./255., 255./255., "CB10_SpaceWolves",1.0);
  new TColor(id++, 146./255., 0./255., 0./255., "CB11_Maroon",1.0);
  new TColor(id++, 146./255., 73./255., 0./255., "CB12_Tan",1.0);
  new TColor(id++, 219./255., 209./255., 0./255., "CB13_Orange",1.0);
  new TColor(id++, 36./255., 255./255., 36./255., "CB14_DayGleen",1.0);
  new TColor(id++, 255./255., 255./255., 109./255., "CB15_SunFlower",1.0);

  /*
  for(Int_t ii=0; ii<20; ii++){
    TColor * tmpcol = gROOT->GetColor(ibase+ii);
    if(tmpcol){
      if(ii==0) cout<<"\nstyle::IniColorCB()"<<endl;
      tmpcol->Print();
      if(ii==19) cout<<endl;
    }
  }
  */
}

TGraphAsymmErrors* style::ScaleGraph(const TGraphAsymmErrors * gin, const  TGraphAsymmErrors * gref)
{
  TGraphAsymmErrors *gout = new TGraphAsymmErrors;
  //for each ref, search in, save only when both exist
  for(Int_t iref=0; iref<gref->GetN(); iref++){
    const Double_t xx = gref->GetX()[iref];

    Int_t targetid = -999;
    for(Int_t iin = 0; iin<gin->GetN(); iin++){
      if( fabs(gin->GetX()[iin]-xx)<EPSILON){
        targetid = iin;
        break;
      }
    }
    if(targetid<0){
      continue;
    }

    const Double_t yref = gref->GetY()[iref];

    //only for yref > 0
    if(yref<EPSILON){
      continue;
    }

    const Double_t xin = gin->GetX()[targetid];
    const Double_t exl = gin->GetEXlow()[targetid];
    const Double_t exh = gin->GetEXhigh()[targetid];

    const Double_t yin = gin->GetY()[targetid];
    const Double_t eyl = gin->GetEYlow()[targetid];
    const Double_t eyh = gin->GetEYhigh()[targetid];

    const Int_t iout = gout->GetN();

    gout->SetPoint(iout, xin, yin/yref);
    gout->SetPointError(iout, exl, exh, eyl/yref, eyh/yref);
  }

  gout->SetName(Form("%sDivBy%s", gin->GetName(), gref->GetName()));
  return gout;
}

void style::GraphMinMaxXY(const TGraphAsymmErrors * gr, Double_t & xmin, Double_t & xmax, Double_t &ymin, Double_t &ymax )
{
  for(Int_t ii=0; ii<gr->GetN(); ii++){
    const Double_t yy = gr->GetY()[ii];
    const Double_t xx = gr->GetX()[ii];
    const Double_t ey = gr->GetEYlow()[ii]*2;
    const Double_t ex = gr->GetEXlow()[ii]*2;

    if(yy-ey<ymin)
      ymin = yy-ey;

    if(yy+ey>ymax)
      ymax = yy+ey;

    if(xx-ex<xmin)
      xmin = xx-ex;

    if(xx+ex>xmax)
      xmax = xx+ex;
  }
}

void style::PadSetup(TPad *currentPad, const Double_t currentLeft, const Double_t currentTop, const Double_t currentRight, const Double_t currentBottom)
{
  currentPad->SetTicks(1,1);
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);

  currentPad->SetFillColor(0);//this is the desired one!!!
}

void style::PadSetup(TVirtualPad *currentPad, const Double_t currentLeft, const Double_t currentTop, const Double_t currentRight, const Double_t currentBottom)
{
  currentPad->SetTicks(1,1);
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);

  currentPad->SetFillColor(0);//this is the desired one!!!
}

void style::ToNaturalScale(TGraph* gr)
{
  for(Int_t ip=0; ip<gr->GetN(); ip++){
    Double_t xx=-999, yy=-999;
    gr->GetPoint(ip,xx,yy);
    gr->SetPoint(ip,TMath::Power(10,xx), yy);
  }
}

void style::ToNaturalScale(TGraphErrors* gr)
{
  for(Int_t ip=0; ip<gr->GetN(); ip++){
    Double_t xx=-999, yy=-999;
    gr->GetPoint(ip,xx,yy);
    //printf("%d %e %e\n",ip,xx,yy);
    gr->SetPoint(ip,TMath::Power(10,xx), yy);
  }
}

void style::ToNaturalScale(TGraphAsymmErrors* gr)
{
  for(Int_t ip=0; ip<gr->GetN(); ip++){
    const Double_t xx = gr->GetX()[ip];
    const Double_t exl = gr->GetEXlow()[ip];
    const Double_t exh = gr->GetEXhigh()[ip];

    const Double_t yy = gr->GetY()[ip];
    const Double_t eyl = gr->GetEYlow()[ip];
    const Double_t eyh = gr->GetEYhigh()[ip];

    //printf("%d %e %e\n",ip,xx,yy);
    gr->SetPoint(ip,TMath::Power(10,xx), yy);
    gr->SetPointError(ip, TMath::Power(10, xx) - TMath::Power(10, xx-exl), TMath::Power(10, xx+exh)-TMath::Power(10, xx), eyl, eyh);
  }
}

void style::ToNaturalScale(TAxis *ax)
{
  TAxis* oldx = (TAxis*)ax->Clone("oldx");
  ax->SetLimits(TMath::Power(10,oldx->GetXmin()), TMath::Power(10,oldx->GetXmax()));
  const Int_t nb = oldx->GetNbins();
  Double_t *bins = new Double_t[nb+1];
  bins[0]=TMath::Power(10,oldx->GetXmin());
  for(Int_t ii=1; ii<=nb; ii++){
    bins[ii]=TMath::Power(10,oldx->GetBinUpEdge(ii));
  }
  ax->Set(nb, bins);

  delete oldx;
  delete bins;
}

void style::ToNaturalScale(TH1 *hh)
{
  ToNaturalScale(hh->GetXaxis());
}


void style::AxisStyle(TAxis *ax, Bool_t kcen)
{
  ax->SetTickLength(fgkTickLength);

  ax->SetLabelFont(fgkTextFont);
  ax->SetLabelSize(fgkTextSize);
  ax->SetLabelOffset(fgkLabelOffset);

  ax->SetTitleFont(fgkTextFont);
  //printf("title %f text %f \n", fgkTitleSize, fgkTextSize);
  ax->SetTitleSize(fgkTitleSize);

  kcen = 1;
  ax->CenterTitle(kcen);

  ax->SetNdivisions(505);
}

void style::AxisStyle(TGaxis *ax, Bool_t kcen)
{
  ax->SetTickLength(fgkTickLength);

  ax->SetLabelFont(fgkTextFont);
  ax->SetLabelSize(fgkTextSize);
  ax->SetLabelOffset(fgkLabelOffset);

  ax->SetTitleFont(fgkTextFont);
  //printf("title %f text %f \n", fgkTitleSize, fgkTextSize);
  ax->SetTitleSize(fgkTitleSize);

  kcen = 1;
  ax->CenterTitle(kcen);

  ax->SetNdivisions(505);

  ax->SetTitleOffset(style::fgkXTitleOffset);
}

void style::ResetStyle(TLegend *obj, Double_t mar, const Double_t ff)
{
  if(mar>0){
    obj->SetMargin(mar);
  }
  obj->SetFillStyle(-1);
  obj->SetBorderSize(-1);
  obj->SetTextFont(fgkTextFont);
  obj->SetTextSize(fgkTextSize*ff);
}

void style::ResetStyle(TLatex *obj, const Double_t ff, const bool kndc)
{
  if(kndc) obj->SetNDC();
  obj->SetTextFont(fgkTextFont);
  obj->SetTextSize(fgkTextSize*ff);
}

void style::ResetStyle(TPaveText *obj, const Double_t ff)
{
  obj->SetFillStyle(-1);
  obj->SetBorderSize(-1);
  obj->SetTextFont(fgkTextFont);
  obj->SetTextSize(fgkTextSize*ff);
  obj->SetTextAlign(11);
}

void style::ResetStyle(TGraph *obj)
{
  obj->SetMarkerSize(fgkMarkerSize);
  obj->SetLineWidth(fgkLineWidth);

  AxisStyle(obj->GetXaxis());
  AxisStyle(obj->GetYaxis());

  obj->GetXaxis()->SetTitleOffset(fgkXTitleOffset);
  obj->GetYaxis()->SetTitleOffset(fgkYTitleOffset);
}

void style::ResetStyle(TF1 * obj, Bool_t kcen)
{
  obj->SetMarkerSize(fgkMarkerSize);

  AxisStyle(obj->GetXaxis(), kcen);
  AxisStyle(obj->GetYaxis(), kcen);

  obj->GetXaxis()->SetTitleOffset(fgkXTitleOffset);
  obj->GetYaxis()->SetTitleOffset(fgkYTitleOffset);
}


void style::ResetStyle(TH1 * obj, TVirtualPad* cpad, Bool_t kcen)
{
  if(!obj){
    printf("style::ResetStyle obj null!\n"); exit(1);
  }

  obj->SetMarkerSize(fgkMarkerSize);

  AxisStyle(obj->GetXaxis(), kcen);
  AxisStyle(obj->GetYaxis(), kcen);

  obj->GetXaxis()->SetTitleOffset(fgkXTitleOffset);
  obj->GetYaxis()->SetTitleOffset(fgkYTitleOffset);

  if(cpad){
    TPaletteAxis *palette = (TPaletteAxis*)obj->GetListOfFunctions()->FindObject("palette");
    if(!palette){
      printf("ResetStyle no palette!!\n"); 
      obj->GetListOfFunctions()->Print();
    }
    else{
      palette->SetX1NDC(1-cpad->GetRightMargin()+0.005);
      palette->SetX2NDC(1-cpad->GetRightMargin()/3*2);
      palette->SetY1NDC(cpad->GetBottomMargin());
      palette->SetY2NDC(1-cpad->GetTopMargin());
      palette->SetLabelFont(fgkTextFont);
      palette->SetLabelSize(fgkTextSize);
      palette->SetLabelOffset(fgkLabelOffset);
    }
  }

}

TMatrixD style::GetDiagonal(const TMatrixD err, const bool kinvert)
{
  const Int_t ndim = err.GetNrows();
  TMatrixD diaMatrix(ndim, ndim);
  for(int ii=0; ii<ndim; ii++){
    const double nor = sqrt( err[ii][ii]) ;

    if(kinvert){

      if(fabs(nor)<1E-100){
        diaMatrix[ii][ii] = 0;
      }
      else{
        diaMatrix[ii][ii] = 1/nor;
      }

    }
    else{
      diaMatrix[ii][ii] = nor;
    }
  }
  return diaMatrix;
}


void style::Corr2Err(const TMatrixD corrMatrix, TMatrixD & err)
{
  //err already has diagonal entries
  const TMatrixD diaMatrix = GetDiagonal(err);

  printf("\n\nstyle::Corr2Err\n");
  printf("Original Error Matrix (should expect only diagonal entries:\n");
  err.Print();
  err = diaMatrix * corrMatrix * diaMatrix;
  printf("Error Matrix:\n");
  err.Print();
  printf("Correlation Matrix:\n");
  corrMatrix.Print();
  printf("<-style::Corr2Err\n\n");
}

void style::Err2Corr(const TMatrixD err, TMatrixD & corrMatrix)
{
  if(err.GetNrows()!=corrMatrix.GetNrows()){
    printf("style::Err2Corr %d %d\n", err.GetNrows(), corrMatrix.GetNrows());
    exit(1);
  }

  const TMatrixD diaInv = GetDiagonal(err, true);

  corrMatrix = diaInv * err * diaInv;

  printf("\n\nstyle::Err2Corr\n");
  printf("Error Matrix:\n");
  err.Print();
  printf("Correlation Matrix:\n");
  corrMatrix.Print();
  printf("<-style::Err2Corr\n\n");
}

double style::GetChi2(TH1D * rawdata, const TMatrixD rawcov, const int noff, const int ndim, TH1D * rawtheory, const double unit)
{
  //test->
  //printf("\n\n\n super test Setting theory to 0 !!\n\n");
  //rawtheory->Scale(0);
  //<--

  printf("\n\nstyle::GetChi2\n\n");

  TH2D *hdummy = 0x0;
  TMatrixD covMatrix = CleanErrorMatrix(rawdata, rawcov, noff, ndim, hdummy, unit);
  delete hdummy;

  TMatrixD diffMatrix(ndim,1);
  for(int ii=0; ii<ndim; ii++){
    const int preidii = ii+noff;

    const double dataxx = rawdata->GetBinCenter(preidii); 
    const double datayy = rawdata->GetBinContent(preidii)/unit;

    //update error according to the covariance which might change depending on the compoents considered
    rawdata->SetBinError(preidii, sqrt(covMatrix[ii][ii])*unit);

    const int theorybin = rawtheory->FindBin(dataxx);
    const double theoryyy =  rawtheory->GetBinContent(theorybin)/unit;

    diffMatrix[ii]=(datayy-theoryyy);
  }

  //_______________________________________________________________________________________________________
  //_______________________________________________________________________________________________________

  TMatrixD diffTranspose=diffMatrix;
  diffTranspose.Transpose(diffTranspose);

  printf("Data Theory Difference:\n");
  diffMatrix.Print();
  diffTranspose.Print();
  printf("<<--Data Theory Difference\n");
  //_______________________________________________________________________________________________________
  //_______________________________________________________________________________________________________
  
  TMatrixD covInvert = covMatrix.Invert();
  //covInvert.Print();

  TMatrixD covdiff = covInvert*diffMatrix;
  covdiff.Print();

  TMatrixD chi2Matrix = diffTranspose*covInvert*diffMatrix;
  chi2Matrix.Print();

  const double chsq = chi2Matrix[0][0];

  return chsq;
}

TMatrixD style::CleanErrorMatrix(TH1D * rawhist, const TMatrixD rawcov, const Int_t noff, const Int_t ndim, TH2D * & covHist, const double unit)
{
  printf("\n\nstyle::CleanErrorMatrix\n\n");

  if(covHist){
    printf("covHist already exist!\n"); exit(1);
  }  

  double xb[100];
  GetBins(rawhist, xb, noff, ndim);
  covHist = new TH2D(Form("covHist%s", rawhist->GetName()),"",ndim, xb, ndim, xb);
  covHist->SetXTitle(rawhist->GetXaxis()->GetTitle());
  covHist->SetYTitle(rawhist->GetXaxis()->GetTitle());

  TMatrixD covMatrix(ndim, ndim);
  for(int ii=0; ii<ndim; ii++){
    const int preidii = ii+noff;

    for(int jj=0; jj<ndim; jj++){
      const int preidjj = jj+noff;

      covMatrix[ii][jj] = rawcov[preidii][preidjj]/unit/unit;

      covHist->SetBinContent(ii+1, jj+1,  covMatrix[ii][jj]);
    }
  }

  printf("style::CleanErrorMatrix cleaned Covariance:\n");
  covMatrix.Print();
  printf("<<-cleaned Covariance\n");

  return covMatrix;
}

void style::GetBins(TH1D * hh, double *xb, const int noff, const int ndim)
{
  printf("\n\nstyle::GetBins\n\n");
  for(int ii=0; ii<=ndim; ii++){
    xb[ii]=hh->GetXaxis()->GetBinLowEdge(ii+noff);
    printf("%f ", xb[ii]);
  }
  printf("\n");
}

TH1D * style::GetStackedSum(THStack *stk, const char * name, const int col, const int lsty, const int lwid, const int fsty)
{
  const TList * ll = stk->GetHists();
  TH1D * hout = (TH1D*)ll->At(0)->Clone(Form("%ssum",name?name:stk->GetName()));
  hout->SetDirectory(0);
  for(Int_t ii=1; ii<ll->GetEntries(); ii++){
    hout->Add((TH1D*)ll->At(ii));
  }

  hout->SetLineColor(col);
  hout->SetLineStyle(lsty);
  hout->SetLineWidth(lwid);
  hout->SetFillStyle(fsty);

  return hout;
}

void style::ResetStyle(THStack * obj, Bool_t kcen)
{
  //cout<<"\n\n********************* Will Crash! no ResetStyle for THStack! ******************\n\n"<<endl;

  return;

  AxisStyle(obj->GetXaxis(), kcen);
  AxisStyle(obj->GetYaxis(), kcen);

  obj->GetXaxis()->SetTitleOffset(fgkXTitleOffset);
  obj->GetYaxis()->SetTitleOffset(fgkYTitleOffset);
}


void style::SetGlobalStyle(const Int_t lStat, Bool_t kcolor)
{
  // Set gStyle
  // From plain

  IniColorCB();

  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(-1);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(-1);
  gStyle->SetLegendBorderSize(-1);
  //
  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(fgkTextFont);
  gStyle->SetStatFont(fgkTextFont);
  gStyle->SetStatFontSize(fgkTextSize);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  gStyle->SetTickLength(fgkTickLength,"xy");
  gStyle->SetEndErrorSize(3);
  gStyle->SetLabelSize(fgkTextSize,"xyz");
  gStyle->SetLabelFont(fgkTextFont,"xyz"); 
  gStyle->SetLabelOffset(fgkLabelOffset,"xyz");
  gStyle->SetTitleFont(fgkTextFont,"xyz");  
  gStyle->SetTitleFont(fgkTextFont,"");  
  gStyle->SetTitleFontSize(fgkTitleSize);
  gStyle->SetTitleOffset(fgkXTitleOffset,"x");  
  gStyle->SetTitleOffset(fgkYTitleOffset,"y");  
  gStyle->SetTitleOffset(1.0,"z");  
  gStyle->SetTitleSize(fgkTitleSize,"xyz");  
  gStyle->SetTitleSize(fgkTitleSize,"");  
  gStyle->SetMarkerSize(fgkMarkerSize); 
  gStyle->SetPalette(1,0); 
  if (lStat){
    gStyle->SetOptTitle(1);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    }
  else {
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
  }

  TGaxis::SetMaxDigits(3);
  gStyle->SetTitleBorderSize(-1);

  if(kcolor){
    SetColor();
  }

  gROOT->ForceStyle();
}

  
void style::SetColor()
{
  gStyle->SetHistFillColor(0);
  //gStyle->SetFillColor(0);//it conflicts with color palette
  gStyle->SetFrameFillColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetTitleFillColor(0);

  //---

  const Int_t nRGBs = 5;
  const Int_t nCont = 100;
  Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[nRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[nRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  const Int_t cgc = TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);

  const Int_t nc = nCont;
  Int_t colors[nc];
  gStyle->SetNumberContours(nc); 
  for(Int_t ii=0; ii<nc; ii++){
    colors[ii]=cgc+ii;
  }
  gStyle->SetPalette(nc, colors);
}

void style::BinLogX(TAxis* axis)
{
  //void XGLUtils::BinLogX(TAxis *axis)
  //
  // Method for the correct logarithmic binning of histograms
  // copied and modified from AliTPCcalibBase

  const Int_t bins = axis->GetNbins();

  const Double_t from = axis->GetXmin();
  const Double_t to = axis->GetXmax();
  if (from<EPSILON) return;
  Double_t *new_bins = new Double_t[bins + 1];

  new_bins[0] = from;
  const Double_t factor = TMath::Power(to/from, 1./bins);

  for (int i = 1; i <= bins; i++) {
   new_bins[i] = factor * new_bins[i-1];
  }
  axis->Set(bins, new_bins);
  delete [] new_bins;

}

void style::UpdateLogX(TH1 *obj)
{
  TAxis * axis = obj->GetXaxis();
  double xmin = axis->GetXmin();
  double xmax = axis->GetXmax();
  axis->SetLimits(TMath::Power(10,xmin), TMath::Power(10,xmax));
  BinLogX(axis);
}

void style::UpdateLogX(TGraphAsymmErrors *obj)
{
  const Int_t np = obj->GetN();
  for(Int_t ii=0; ii<np; ii++){
    Double_t xx =-999, yy=-999;
    obj->GetPoint(ii,xx,yy);
    Double_t ex0 = obj->GetErrorXlow(ii);
    Double_t ex1 = obj->GetErrorXhigh(ii);
    const Double_t ey0 = obj->GetErrorYlow(ii);
    const Double_t ey1 = obj->GetErrorYhigh(ii);

    ex0 = TMath::Power(10, xx) - TMath::Power(10, xx-ex0);
    ex1 = TMath::Power(10, xx+ex1) - TMath::Power(10, xx);
    xx = TMath::Power(10, xx);

    obj->SetPoint(ii, xx, yy);
    obj->SetPointError(ii, ex0, ex1, ey0, ey1);
  }
}


void style::UpdateLogX(TGraphErrors *obj)
{
  const Int_t np = obj->GetN();
  for(Int_t ii=0; ii<np; ii++){
    Double_t xx =-999, yy=-999;
    obj->GetPoint(ii,xx,yy);
    Double_t ex = obj->GetErrorX(ii);
    const Double_t ey = obj->GetErrorY(ii);

    ex = TMath::Power(10, xx) - TMath::Power(10, xx-ex);
    xx = TMath::Power(10, xx);

    obj->SetPoint(ii, xx, yy);
    obj->SetPointError(ii, ex, ey);
  }
}


void style::DividegPad(int nx, int ny, double l, double r, double t, double b)
{
  gStyle->SetPadTopMargin(0.);
  gStyle->SetPadBottomMargin(0.);
  gStyle->SetPadLeftMargin(0.);
  gStyle->SetPadRightMargin(0.); 
 
   int ix, iy, n=0;
   double x1, x2, y1, y2;
   double dx = ((1-r)*(1-l))/((1-r)*(1-l)*(nx-2)-r+2-l);
   double dl = dx/(1-l);
   double dy = ((1-t)*(1-b))/((1-t)*(1-b)*(ny-2)-t+2-b);
   double db = dy/(1-b);
   char *name  = new char [strlen(gPad->GetName())+6];

   y1 = 0;
   y2 = db;
   for (iy=0; iy<ny; iy++) {
      x1 = 0;
      x2 = dl;
      for (ix=0;ix<nx;ix++) {
         if (x1 > x2) continue;
         n++;
         sprintf(name,"%s_%d",gPad->GetName(),n);
         TPad *pad = new TPad(name,name,x1,y1,x2,y2,0);
         if (ix==0)    pad->SetLeftMargin(l);
         if (ix==nx-1) pad->SetRightMargin(r);
         if (iy==ny-1) pad->SetTopMargin(t);
         if (iy==0)    pad->SetBottomMargin(b);
         x1 = x2;
         if (ix==nx-2) x2 = 1;
         else          x2 = x1+dx;
         pad->SetNumber(n);
         pad->Draw();
         pad->SetTicks();
      }
      y1 = y2;
      if (iy==ny-2) y2 = 1;
      else          y2 = y1+dy;
   }
}


void style::DividePad(TPad *pin, int nx, int ny, double l, double r, double t, double b, TList * lout, Int_t &npad, const Int_t lastny)
{
  pin->SetTopMargin(0.);
  pin->SetBottomMargin(0.);
  pin->SetLeftMargin(0.);
  pin->SetRightMargin(0.); 
 
   int ix, iy, n=0;
   double x1, x2, y1, y2;
   double dx = nx>=2? ( ((1-r)*(1-l))/((1-r)*(1-l)*(nx-2)-r+2-l) ) : 0;
   double dl = nx>=2? dx/(1-l) : 1;
   double dy = ((1-t)*(1-b))/((1-t)*(1-b)*(ny-2)-t+2-b);
   double db = dy/(1-b);
   char *name  = new char [strlen(pin->GetName())+6];

   //printf("test dx %f dl %f dy %f db %f\n", dx, dl, dy, db);

   y1 = 0;
   y2 = db;
   for (iy=0; iy<lastny+1; iy++) {
     //printf("test1 iy %d\n", iy);
      x1 = 0;
      x2 = dl;
      for (ix=0;ix<nx;ix++) {
         if (x1 > x2) continue;
         //printf("test2 ix %d\n", ix);

         n++;
         sprintf(name,"%s_%d",pin->GetName(),n);
         TPad *pad = new TPad(name,name,x1,y1,x2,y2,0);lout->Add(pad);
         PadSetup(pad);

  pad->SetTopMargin(0.);
  pad->SetBottomMargin(0.);
  pad->SetLeftMargin(0.);
  pad->SetRightMargin(0.); 
 
         if (ix==0)    pad->SetLeftMargin(l);
         if (ix==nx-1) pad->SetRightMargin(r);
         if (iy==lastny) pad->SetTopMargin(t);
         if (iy==0)    pad->SetBottomMargin(b);
         x1 = x2;
         if (ix==nx-2) x2 = 1;
         else          x2 = x1+dx;
         pad->SetNumber(npad);npad++;
         pad->Draw();
         pad->SetTicks();
      }
      y1 = y2;
      if (iy==lastny-1) y2 = 1;
      else          y2 = y1+dy;
   }
}

TGraphAsymmErrors * style::GetInverse(const TGraphAsymmErrors * graw)
{
  TGraphAsymmErrors * ginv = (TGraphAsymmErrors*) graw->Clone(Form("%sInv", graw->GetName()));

  const Int_t np = graw->GetN();

  Int_t ipoint = 0;
  for(Int_t ii=0; ii<np; ii++){
    const Double_t xx = graw->GetX()[ii];
    const Double_t exl = graw->GetEXlow()[ii];
    const Double_t exh = graw->GetEXhigh()[ii];

    const Double_t tmpyy = graw->GetY()[ii];
    const Double_t tmpeyl = graw->GetEYlow()[ii];
    const Double_t tmpeyh = graw->GetEYhigh()[ii];

    if(tmpyy<EPSILON)
      continue;

    const Double_t yy = 1/tmpyy;
    const Double_t eyl = yy * tmpeyl/tmpyy;
    const Double_t eyh = yy * tmpeyh/tmpyy;

    ginv->SetPoint(ipoint, xx, yy);
    ginv->SetPointError(ipoint, exl, exh, eyl, eyh);
    ipoint++;
  }

  ginv->Set(ipoint);

  //have to reset title, at least, after Set                                                                                                                                                                                                                 
  ginv->GetYaxis()->SetTitle(Form("inverse of %s",graw->GetYaxis()->GetTitle()));
  ginv->GetXaxis()->SetTitle(graw->GetXaxis()->GetTitle());

  return ginv;
}


TH1D * style::Func2Hist(TF1 * ff, const Bool_t klogx)
{
  TH1D * hh=0x0;
  TString hname(Form("hist%s%d", ff->GetName(), gRandom->Integer(1000)));
  if(gDirectory->Get(hname)){
    hname = Form("hist%s%d", ff->GetName(), gRandom->Integer(1000));
  }

  hh=new TH1D(hname,"",500, ff->GetXmin(), ff->GetXmax());
  if(klogx){
    BinLogX(hh->GetXaxis());
  }

  hh->SetLineColor(ff->GetLineColor());
  hh->SetLineStyle(ff->GetLineStyle());

  for(Int_t ii=1; ii<=hh->GetNbinsX(); ii++){
    const Double_t xi = hh->GetBinCenter(ii);
    const Double_t yi = ff->Eval(xi);
    hh->SetBinContent(ii, yi);
  }

  hh->GetXaxis()->SetTitle(ff->GetXaxis()->GetTitle());
  hh->GetYaxis()->SetTitle(ff->GetYaxis()->GetTitle());
  return hh;
}



