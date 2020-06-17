#ifndef _ANAUTILS_H_
#define _ANAUTILS_H_

#include "AnaIO.h"

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


TLorentzVector * GetPiZero(const int ksig, const vector<TLorentzVector> & shws,  const bool kprint, const bool kfill)
{
  const int shsize = shws.size();
  if(kfill){
    AnaIO::hRecPi0Nshower->Fill(shsize, !ksig);
  }

  vector<TLorentzVector> piarr;
  int maxEidx = -999;
  double pi0Emax = -999;
  if(shsize>=2){
    for(int ii = 0; ii<shsize; ii++){
      for(int kk = ii+1; kk<shsize; kk++){
        const TLorentzVector tmpgg = shws[ii]+shws[kk];
        if(kprint){
          printf("\n============ %d %d %d \n", shsize, ii, kk);
          tmpgg.Print();        
        }

        if(tmpgg.E() > pi0Emax){
          pi0Emax = tmpgg.E();
          maxEidx = piarr.size();
        }

        piarr.push_back(tmpgg);
      }
    }
  }  

  TLorentzVector * outcopy = 0x0;
  if(maxEidx>=0){
    outcopy = new TLorentzVector(piarr[maxEidx]);
  }

  return outcopy;
}

int GetParticleType(const int pdg)
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



double GetTruncatedMean(vector<double> array, const int nsample)
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

double GetTruncatedMean(const vector<double> &tmparr, const unsigned int nsample0, const unsigned int nsample1, const double lowerFrac, const double upperFrac)
{
  //for proton Bragg peak use 0.4-0.95. Seen by CDF of startE using signal proton samples in drawTracking

  if(nsample1>=tmparr.size()){
    return -999;
  }

  vector<double> array;
  for(unsigned int ii=nsample0; ii<=nsample1; ii++){
    array.push_back(tmparr[ii]);
  }

  std::sort(array.begin(), array.end());

  double sum =0.0;

  const int iter0 = array.size()*lowerFrac;
  const int iter1 = array.size()*upperFrac;

  for(int ii=iter0; ii< iter1; ii++){
    sum += array[ii];
  }
  return sum / ( (iter1-iter0)+1E-10);
}


void GetdEdx(const vector<double> &arraydEdx, vector<double> &startE, vector<double> &endE, const unsigned int padding=0)
{
  const unsigned int ncls = arraydEdx.size();

  if(ncls<=padding){
    return;
  }

  //start from [2] because [0] and [1] in both start and last are weird
  for(unsigned int kk=padding; kk<ncls; kk++){
    startE.push_back(arraydEdx[kk]);

    const double endpe = arraydEdx[ncls-1-kk];
    endE.push_back(endpe);
  }
}


int GetFillEventType()
{
  int filleventtype = -999;
  if(AnaIO::kSignal){
    filleventtype = 0;
  }
  else if(AnaIO::true_beam_PDG==211){
    filleventtype = 1;
  }
  else{
    filleventtype = 2;
  }
  
  return filleventtype;
}


void PrintLegend()
{
  TCanvas * c1 = 0x0;

  vector<TString> parType;
  parType.push_back("p");
  parType.push_back("#pi^{+}");
  parType.push_back("EM shower");
  parType.push_back("others");
  c1 = style::DrawLegend(parType, "f");
  c1->Print("output/legend_parType.eps");
  c1->Print("output/legend_parType.pdf");
  c1->Print("output/legend_parType.png");

  vector<TString> evtType;
  evtType.push_back("signal");
  evtType.push_back("background");
  evtType.push_back("#mu^{+} beam");
  c1 = style::DrawLegend(evtType, "f");
  c1->Print("output/legend_evtType.eps");
  c1->Print("output/legend_evtType.pdf");
  c1->Print("output/legend_evtType.png");
}

//=== end
}


#endif

