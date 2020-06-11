#ifndef _ANAUTILS_H_
#define _ANAUTILS_H_


using namespace std;

namespace AnaUtils
{

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
  gkNucleus,

  //15
  gkNotCaredType
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
  else if(pdg==11){
    type = gkElectron;
  }
  else if(pdg==-13){
    type = gkMuPlus;
  }
  else if(pdg==13){
    type = gkMuMinus;
  }
  else if(pdg==321||pdg==-321||pdg==310||pdg==130){
    type = gkKaon;
  }
  else if(pdg==2112){
    type = gkNeutron;
  }
  else if(pdg==22){
    type = gkGamma;
  }
  else if(pdg==14 || pdg==-14 || pdg==12 || pdg==-12){
    type = gkNeutrino;
  }
  else if(pdg==3122||pdg==3212||pdg==3222){
    type = gkHyperon;
  }
  else if(pdg>9999){
    type = gkNucleus;
  }
  else{
    cout<<"getParticleType unknown pdg "<<pdg<<endl; exit(1);
    //type = gkNotCaredType;
  }

  return type;
}

void getProfileX(TList *lout)
{
  const int nhist = lout->GetSize();//the size will increase
  for(int ii=0; ii<nhist; ii++){
    TH2D * htmp = dynamic_cast<TH2D *>( lout->At(ii) );
    if(htmp){
      TH1D * hpro = htmp->ProfileX(Form("%s_profileX", htmp->GetName())); lout->Add(hpro);

      TH1D * hpdf = 0x0;
      TH1D * hcdf = 0x0;
      const double thres = 5;
      TH2D * hnor = style::NormalHist(htmp, hpdf, hcdf, thres, true);
      hpdf->SetYTitle("p.d.f.");
      hcdf->SetYTitle("c.d.f.");
      hpro->SetYTitle("purity");
      hpdf->SetTitle(htmp->GetName());
      hcdf->SetTitle(htmp->GetName());
      hpro->SetTitle(htmp->GetName());
      lout->Add(hpdf);
      lout->Add(hcdf);
      lout->Add(hnor); 
    }
  }
}

void drawHist(TList *lout, const TString outdir, const TString tag, const bool kfast=false)
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
      gStyle->SetOptStat(0);
      //gStyle->SetOptStat("enou");
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
      if(tag.Contains("profileX")){
        dopt = "box";
      }
      else if(tag.Contains("nor")){
        dopt = "colz";
      }
    }
    else{
      if(tag.Contains("profileX")){
        hh->SetMaximum(0.25);
        dopt = "hist";
      }
      else if(tag.Contains("pdf") || tag.Contains("cdf")){
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

}


#endif

