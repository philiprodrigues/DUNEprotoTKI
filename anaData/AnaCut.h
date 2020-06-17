#ifndef _ANACUT_H_
#define _ANACUT_H_

#include "AnaIO.h"
#include "AnaUtils.h"

using namespace std;

namespace AnaCut
{

int GetTruthFromRec(const int recidx, int & pdg, double & momentum)
{
  pdg = -999;
  momentum = -999;
  int trueidx = -999;

  if(AnaIO::reco_daughter_PFP_true_byHits_ID && AnaIO::reco_daughter_PFP_true_byHits_ID->size() ){//in data this is size = 0 but not null

    const int truthID = (*AnaIO::reco_daughter_PFP_true_byHits_ID)[recidx];
    int counter = 0;
    for(unsigned int ii = 0; ii<AnaIO::true_beam_daughter_ID->size(); ii++){
      if(  (*AnaIO::true_beam_daughter_ID)[ii] == truthID ){
        if(counter){
          printf("GetTruthFromRec truthID found again %d %d\n", recidx, truthID); exit(1);
        }
        trueidx = ii;
        counter++;
      }
    }
    
    if(trueidx>=0){
      pdg = (*AnaIO::true_beam_daughter_PDG)[trueidx];
      
      const TVector3 vec((*AnaIO::true_beam_daughter_startPx)[trueidx], (*AnaIO::true_beam_daughter_startPy)[trueidx], (*AnaIO::true_beam_daughter_startPz)[trueidx]);
      momentum = vec.Mag();
    }
    
  }

  return trueidx;
}


int GetNTrack(const bool kpi0, const bool ksig, int & nproton, int & nshower, int & nmichel, TLorentzVector *  & leadingPi0, const bool kprint=false, const bool kfill=false, const int kdebug = 0)
{
  /*
root [11] beamana->Scan("reco_daughter_PFP_ID:reco_daughter_PFP_true_byHits_ID:reco_daughter_allTrack_ID:reco_daughter_allShower_ID:true_beam_daughter_ID:true_beam_daughter_reco_byHits_PFP_ID:true_beam_daughter_reco_byHits_allTrack_ID:true_beam_daughter_reco_byHits_allShower_ID")
***********************************************************************************************************************
*    Row   * Instance * reco_daug * reco_daug * reco_daug * reco_daug * true_beam * true_beam * true_beam * true_beam *
***********************************************************************************************************************
*       13 *        0 *        86 *     40841 *        86 *        86 *     40800 *           *           *           *
*       13 *        1 *       226 *     40855 *       225 *       225 *     40841 *           *           *           *

So
(reco_daughter_PFP_ID = reco_daughter_allTrack_ID = reco_daughter_allShower_ID)
                             !=
 (reco_daughter_PFP_true_byHits_ID = true_beam_daughter_ID)

All 
true_beam_daughter_reco_byHits_PFP_ID:
true_beam_daughter_reco_byHits_allTrack_ID:
true_beam_daughter_reco_byHits_allShower_ID
not set
   */

  const int recsize = AnaIO::reco_daughter_PFP_ID->size();
  const int mbr     = AnaIO::reco_daughter_allTrack_momByRange_proton->size();
  if(recsize!=mbr){
    printf("recsize != mbr %d %d\n", recsize, mbr); exit(1); 
  }

  /*
    est rec 1/2 trueID 15 pdg 22 truemomentum 0.002512 recp -999.000000 resolution -397725.603302 nhits 7 trackScore 0.016739 emScore 0.983240 michelScore 0.289565 sum 1.289545
  */
  int ntracks = 0;
  nshower = 0;
  nmichel = 0;
  nproton = 0;
 
  vector<TLorentzVector> showerArray;

  static bool kShowCut = true;

  for(int ii=0; ii<recsize; ii++){

    int pdg = -999;
    double truemomentum = -999;
    const int trueidx = GetTruthFromRec(ii, pdg, truemomentum);

    int fillstktype = -999;
    if(pdg==2212){//proton
      fillstktype = 0;
    }
    else if(pdg==211){//pi+ put it beyond others
      fillstktype = 1;
    }
    else if(pdg==22){//gamma
      fillstktype = 4;
    }
    else if(pdg==-999){//shower no true
      fillstktype = 2;
    }
    else{//all others
      fillstktype = 3;
    }

    //---> need to be done before any counting!!!
    const vector<double>  dedxarray = (*AnaIO::reco_daughter_allTrack_calibrated_dEdX_SCE)[ii];
    const int NdEdx = dedxarray.size();
    if(kfill){
      AnaIO::hCutNdEdx->Fill(NdEdx, fillstktype);
    }
    const int ndedxcut = kpi0 ? 6 : 16; //need E[3], at least 4 cls
    if(kShowCut){
      printf("check cut kpi0 %d ndedxcut %d\n", kpi0, ndedxcut);
    }
    if(NdEdx<ndedxcut){
      //do not count it as anything
      continue;
    }
    //<--- need to be done before any counting!!!

    //track and em scores at 0.5 look reasonable as baseline (after 0.5 the proton fraction look flat)
    //michel is not found (all below 0.5, those above 0.5 are shower and other_type)
    const double trackScore  = (*AnaIO::reco_daughter_PFP_trackScore)[ii];
    if(trackScore>0.5){
      if( (*AnaIO::reco_daughter_allTrack_ID)[ii]==-1 ){
        printf("bad track ID! %d\n", ii); exit(1);
      }
      ntracks++;
    }
    const double emScore     = (*AnaIO::reco_daughter_PFP_emScore)[ii];
    if(emScore>0.5){
      if( (*AnaIO::reco_daughter_allShower_ID)[ii]==-1 ){
        printf("bad shower ID! %d\n", ii); exit(1);
      }
      nshower++;
      
      const TVector3 showerDir( (*AnaIO::reco_daughter_allShower_dirX)[ii], (*AnaIO::reco_daughter_allShower_dirY)[ii], (*AnaIO::reco_daughter_allShower_dirZ)[ii] );
      const TVector3 showerMomentum = showerDir.Unit()*(*AnaIO::reco_daughter_allShower_energy)[ii] * 1E-3; //MeV to GeV
      const TLorentzVector showerLv( showerMomentum, showerMomentum.Mag() );
      showerArray.push_back(showerLv);
    }
    const double michelScore = (*AnaIO::reco_daughter_PFP_michelScore)[ii];
    if(michelScore>0.5){
      nmichel++;
    }

    //nhit 260 is just clean-up
    const int nhits          = (*AnaIO::reco_daughter_PFP_nHits)[ii];
    /*//they are indeed different: 236 80
    if(nhits != NdEdx){
      printf("nhits!=dedx size %d %d\n", nhits, NdEdx); exit(1);
    }
    */

    //========== proton tagging now!
  
    const double startE2 = dedxarray[2];
    const double startE3 = dedxarray[3];
    const double lastE2 = dedxarray[NdEdx-1-2];
    const double lastE3 = dedxarray[NdEdx-1-3];

    const double chi2 = (*AnaIO::reco_daughter_allTrack_Chi2_proton)[ii];
    const double ndof = (*AnaIO::reco_daughter_allTrack_Chi2_ndof)[ii];
    const double Chi2NDF = chi2/(ndof+1E-10);

    bool isSelProton = false;
    const double cutNH = 260;
    if(kpi0){
      //all signal protons have nhit below 260
      //all signal protons have startE3 > 9
      //with Chi2NDF cut  44/191 = 23% purity; without 46/197 = 23% purity, slightly higher efficiency
      //test if(startE2>10 && nhits<260 && startE3>9 && Chi2NDF<50){
      
      const double cutSE2 = 10;
      const double cutSE3 = 9;
      isSelProton = (startE2>cutSE2 && nhits<cutNH && startE3>cutSE3);
      if(kShowCut){
        printf("check cut kpi0 %d proton tag startE2 %.2f nhits %.0f startE3 %.2f\n", kpi0, cutSE2, cutNH, cutSE3); 
      }
    }
    else{
      const double cutCHI = 50;
      isSelProton = (Chi2NDF<cutCHI && nhits<cutNH);
      if(kShowCut){
        printf("check cut kpi0 %d proton tag Chi2NDF %.2f nhits %.0f\n", kpi0, cutCHI, cutNH);
      }
    }

    double recp        = -999;
    if(isSelProton){
      nproton++;
      recp = (*AnaIO::reco_daughter_allTrack_momByRange_proton)[ii];
    }
    else{
      //temporary pi+ momentum
      recp = (*AnaIO::reco_daughter_allTrack_momByRange_muon)[ii];
    }
    //========== proton tagging done!
   

    if(kprint){
      printf("test sig %d rii %d/%d tii %4d pdg %4d sE2 %6.1f sE3 %6.1f lE2 %6.1f lE3 %6.1f nhi %4d tkS %6.1f emS %6.1f miS %6.1f sum %6.1f chf %6.1f\n", ksig, ii, recsize, trueidx, pdg, startE2, startE3, lastE2, lastE3, nhits, trackScore, emScore, michelScore, trackScore+emScore+michelScore, Chi2NDF);
    }

    if(kfill){
      if(pdg==2212 && isSelProton){
        AnaIO::hProtonMomentumRes->Fill(truemomentum, recp/truemomentum-1);
      }
      else if(pdg==211 && !isSelProton){
        AnaIO::hPiMomentumRes->Fill(truemomentum, recp/truemomentum-1);
      }

      if(!kdebug || 
         (kdebug==1 && isSelProton) ||
         (kdebug==2 && !isSelProton)
         ){
      AnaIO::hCutstartE2->Fill(startE2, fillstktype);
      AnaIO::hCutstartE3->Fill(startE3, fillstktype);
      AnaIO::hCutlastE2->Fill(lastE2, fillstktype);
      AnaIO::hCutlastE3->Fill(lastE3, fillstktype);

      AnaIO::hCutnHits->Fill(nhits, fillstktype);
      AnaIO::hCutChi2NDF->Fill(Chi2NDF, fillstktype);
      AnaIO::hCuttrackScore->Fill(trackScore, fillstktype);
      AnaIO::hCutemScore->Fill(emScore, fillstktype);
      AnaIO::hCutmichelScore->Fill(michelScore, fillstktype);
      }
    }

    kShowCut = false;
  }

  //only fill nshowers when doing fill
  leadingPi0 = AnaUtils::GetPiZero(ksig, showerArray, false, kfill);

  if(kprint){
    printf("\n\n");
  }

  return ntracks;
}


bool CutTopology(const bool kpi0)
{
    /*
      PiPlus: nshower = 0, nproton=1 -> 75/218=34% in signal sample; 63/467 = 13% purity, 63/218 = 29% efficiency, p*e = 63*63/218/467 = 3.9%
      PiPlus: nshower = 0, nproton=1,  noly nhit<260 and startE3>9 and Chi2NDF<50 is proton, no pion-analysis pre-cuts! -> 61/380 = 16% 
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 and Chi2NDF<50 is proton, no pion-analysis pre-cuts! -> 44/191 = 23% purity, e*p = 44.*44./191./218. = 4.6%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9                is proton, no pion-analysis pre-cuts! -> 46/197 = 23% purity, slightly higher efficiency
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=4 instead of 3) is proton, no pion-analysis pre-cuts! -> 49/199 = 25%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=5 instead of 3) is proton, no pion-analysis pre-cuts! -> 48/201 = 24%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=6 instead of 3) is proton, no pion-analysis pre-cuts! -> 50/202 = 25%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=7 instead of 3) is proton, no pion-analysis pre-cuts! -> 50/200 = 25%, e*p = 5.7%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=8 instead of 3) is proton, no pion-analysis pre-cuts! -> 50/194 = 26%, e*p = 5.9%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=9 instead of 3) is proton, no pion-analysis pre-cuts! -> 47/176 = 27%, e*p = 5.8%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=10 instead of 3)is proton, no pion-analysis pre-cuts! -> 46/168 = 27%, e*p = 5.8%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=13 instead of 3)is proton, no pion-analysis pre-cuts! -> 43/159 = 27%, e*p = 5.3%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=14 instead of 3)is proton, no pion-analysis pre-cuts! -> 42/149 = 28%, e*p = 5.4%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=15 instead of 3)is proton, no pion-analysis pre-cuts! -> 40/135 = 30%, e*p = 5.4%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=17 instead of 3)is proton, no pion-analysis pre-cuts! -> 39/126 = 31%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=20 instead of 3)is proton, no pion-analysis pre-cuts! -> 34/114 = 30%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=25 instead of 3)is proton, no pion-analysis pre-cuts! -> 25/82 = 30%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=30 instead of 3)is proton, no pion-analysis pre-cuts! -> 22/71 = 31%

      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and chi2<20 (NdEdx>=16 instead of 3, no startE2, E3 cuts at all)is proton, no pion-analysis pre-cuts! -> 39/122 = 32%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and chi2<23 (NdEdx>=16 instead of 3, no startE2, E3 cuts at all)is proton, no pion-analysis pre-cuts! -> 41/123 = 33%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and chi2<25 (NdEdx>=16 instead of 3, no startE2, E3 cuts at all)is proton, no pion-analysis pre-cuts! -> 42/125 = 34%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and chi2<30 (NdEdx>=16 instead of 3, no startE2, E3 cuts at all)is proton, no pion-analysis pre-cuts! -> 42/127 = 33%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and chi2<35 (NdEdx>=16 instead of 3, no startE2, E3 cuts at all)is proton, no pion-analysis pre-cuts! -> 43/128 = 34%, e*p = 6.6%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and chi2<40 (NdEdx>=16 instead of 3, no startE2, E3 cuts at all)is proton, no pion-analysis pre-cuts! -> 43/127 = 34%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and chi2<60 (NdEdx>=16 instead of 3, no startE2, E3 cuts at all)is proton, no pion-analysis pre-cuts! -> 44/138 = 32%

      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and startE3>9 (NdEdx>=16 instead of 3)is proton, no pion-analysis pre-cuts! -> 40/127 = 32%, e*p = 5.8%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and chi2<50 (NdEdx>=16 instead of 3, no startE2, E3 cuts at all)is proton, HAS pion-analysis pre-cuts! -> 39/111 = 35%, e*p = 6.3%
      PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and chi2<50 (NdEdx>=16 instead of 3, no startE2, E3 cuts at all)is proton, no pion-analysis pre-cuts! -> 44/132 = 33%, e*p = 6.7%

      (*) PiPlus New Def:  2) p_lading 0.45-1, p_sl<0.45. PiPlus: nshower = 0, nproton=1, ntrack = 2,  noly nhit<260 and chi2<50 (NdEdx>=16 instead of 3, no startE2, E3 cuts at all)is proton, no pion-analysis pre-cuts! -> 43/132 = 33%, e*p = 43.*43./132./220. = 6.4%

      <-- //p_leading 0.45-1, p_sublieading < 0.45, no cut on pi+
      <-- //p-pi+: 220 events

      PiZero: nshower = 2, nmichel=0, nproton=1 -> 32/260=12% in signal sample; 27/68 = 40% purity, 27/260 = 10% efficiency, p*e = 27*27/68/260 = 4.1%
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1 -> ; 19/40 = 48% purity
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 is proton -> 18/37=49% purity
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 is proton -> 18/32=56% purity, p*e = 18.*18./32./260. = 3.9%
      PiZero: nshower = 2, nmichel=0, nproton=1,           noly nhit<260 and startE3>9 is proton -> 21/47=45% purity -> so it is important to have ntrack=1

      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 and Chi2NDF<50 is proton, no pion-analysis pre-cuts! -> 18/34 = 53% purity
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 is proton, no pion-analysis pre-cuts! -> 20/36=56% purity, eff*purity = 20.*20./36./260. = 4.3%
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=4 instead of 3 ) is proton, no pion-analysis pre-cuts! -> 19/33=58%
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=5 instead of 3 ) is proton, no pion-analysis pre-cuts! -> 21/34=62%
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=7 instead of 3 ) is proton, no pion-analysis pre-cuts! -> 20/33=61%
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=8 instead of 3 ) is proton, no pion-analysis pre-cuts! -> 19/31=61%, e*p = 4.5%
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=10 instead of 3 )is proton, no pion-analysis pre-cuts! -> 13/20=65%, e*p = 3.3%
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=20 instead of 3 )is proton, no pion-analysis pre-cuts! -> 9/10=90%,  e*p = 3.1%

      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and chi2<25 ( NdEdx>=6 instead of 3 but no startE2,E3 at all) is proton, no pion-analysis pre-cuts! -> 20/33=61%
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and chi2<30 ( NdEdx>=6 instead of 3 but no startE2,E3 at all) is proton, no pion-analysis pre-cuts! -> 20/33=61%
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and chi2<35 ( NdEdx>=6 instead of 3 but no startE2,E3 at all) is proton, no pion-analysis pre-cuts! -> 20/34=59%

      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=6 instead of 3 ) is proton, HAS pion-analysis pre-cuts! -> 21/32=66%, e*p= 5.3%
      PiZero: nshower = 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=6 instead of 3 ) is proton, no pion-analysis pre-cuts! -> 22/35=63%, e*p = 5.3%

      (*) PiZero: nshower >= 2, nmichel=0, nproton=1, ntrack=1, noly nhit<260 and startE3>9 ( NdEdx>=6 instead of 3 ) is proton, no pion-analysis pre-cuts! -> 30/47=64%, e*p = 7.4%, eff = 30/260 = 12%, this is the SAME for both 1) p_leading>0.45, p_sl<0.45, p+>0.15 and 2) p_lading 0.45-1, p_sl<0.45. Eff*pur for new def 2) is 30.*30./47./242 = 7.9%
      
      <-- //p_leading 0.45-1, p_sublieading < 0.45, no cut on pi+
      <-- //p-pi0: 242 events
     */

    int cutnproton = 0;
    int cutnshower = 0;
    int cutnmichel = 0;
    //vector<TLorentzVector> piZeroArray;
    TLorentzVector * leadingPi0 = 0x0;
 
    //debug=1 proton candidate; debug=2 non-proton candidate
    const int kdebug = 0;
    if(kdebug==0){
      //only filling
      AnaIO::nTrack = GetNTrack(kpi0, AnaIO::kSignal, cutnproton, cutnshower, cutnmichel, leadingPi0, false, true, 0);
    }
    else{//just no filling, still no debug
      AnaIO::nTrack = GetNTrack(kpi0, AnaIO::kSignal, cutnproton, cutnshower, cutnmichel, leadingPi0, false, false, 0);
    }

    const int filleventtype = AnaUtils::GetFillEventType();
    //----------- do cuts

    AnaIO::hCutnshower->Fill(cutnshower, filleventtype);
    if(kpi0){
      if(cutnshower<2){
        return false;
      }
    }
    else{
      if(cutnshower>0){
        return false;
      }
    }

    AnaIO::hCutnmichel->Fill(cutnmichel, filleventtype);
    if(kpi0){
      if(cutnmichel>0){
        return false;
      }
    }
    else{
      //not found in signal, need Michel to tell 2-proton events from p-pi; ignore michel count for the momentum, things will change after implementing it.
        /*
          if(cutnmichel!=1){
          continue;
          }
        */
    }

    AnaIO::hCutntrack->Fill(AnaIO::nTrack, filleventtype);
    if(kpi0){
      //with this 56% purity, without 45%
      if(AnaIO::nTrack!=1){
        return false;
      }
    }
    else{
      if(AnaIO::nTrack!=2){
        return false;
      }
    }

    AnaIO::hCutnproton->Fill(cutnproton, filleventtype);
    if(cutnproton!=1){
      return false;
    }
    
    //===================================================
    //only checking after-selection distributions
    if(kpi0){
      if(!leadingPi0){
        printf("leadingpi0 null!!\n"); exit(1);
      }
      AnaIO::hCutMpi0->Fill(leadingPi0->M(), filleventtype);
    }

    TLorentzVector *dummypi0 = 0x0;
    //no getting of pizeroarray
    if(kdebug==0){
      //just for printing, no filling
      GetNTrack(kpi0, AnaIO::kSignal, cutnproton, cutnshower, cutnmichel, dummypi0, true, false, 0);
    }
    else{
      //debug mode: major background is -999 shower mocking protons, no efficient cut variables found
      //print, fill, and debug
      GetNTrack(kpi0, AnaIO::kSignal, cutnproton, cutnshower, cutnmichel, dummypi0, true, true, kdebug);
    }

    return true;
}

bool CutBeamID(const std::vector<int> & pidCandidates)
{
  //https://github.com/calcuttj/PionStudies/blob/master/rDataFrame/eventSelection.h
  //line 149

  //https://github.com/calcuttj/PionStudies/blob/master/rDataFrame/eventSelection.C
  //line 128
  // auto data_all = data_frame
  //    .Define("beamPID", data_beam_PID, {"data_BI_PDG_candidates"})
  //line 157
  // .Filter("beamPID == true"); 

  return ( (std::find(pidCandidates.begin(), pidCandidates.end(), 211)) != pidCandidates.end());
};

bool Manual_beamPos_mc(const double beam_startX, const double beam_startY, const double beam_startZ, const double beam_dirX, const double beam_dirY,   const double beam_dirZ, const double true_dirX,   const double true_dirY, const double true_dirZ,   const double true_startX, const double true_startY, const double true_startZ) 
{
  //https://github.com/calcuttj/PionStudies/blob/master/rDataFrame/eventSelection.h

  //For MC from Owen Goodwins studies
  const double xlow = -3.,  xhigh = 7.,  ylow = -8.,  yhigh = 7.;
  const double zlow = 27.5,  zhigh = 32.5,  coslow = 0.93;

  const double projectX = (true_startX + -1*true_startZ*(true_dirX/true_dirZ) );
  const double projectY = (true_startY + -1*true_startZ*(true_dirY/true_dirZ) );
  const double cos = true_dirX*beam_dirX + true_dirY*beam_dirY + true_dirZ*beam_dirZ;

  if ( (beam_startX - projectX) < xlow )
    return false;
  
  if ( (beam_startX - projectX) > xhigh )
    return false;

  if ( (beam_startY - projectY) < ylow )
    return false;

  if ( (beam_startY - projectY) > yhigh )
    return false;
  
  if (beam_startZ < zlow || zhigh < beam_startZ)
    return false;
  
  if ( cos < coslow)
    return false;

  return true;
}


bool CutMCBeamPos()
{
  //https://github.com/calcuttj/PionStudies/blob/master/rDataFrame/eventSelection.C
  /*
.Define("passBeamCut", manual_beamPos_mc, 
            {"reco_beam_startX", "reco_beam_startY", "reco_beam_startZ",
             "reco_beam_trackDirX", "reco_beam_trackDirY", "reco_beam_trackDirZ",
             "true_beam_startDirX", "true_beam_startDirY", "true_beam_startDirZ",
             "true_beam_startX", "true_beam_startY", "true_beam_startZ"})
   */

  return Manual_beamPos_mc(AnaIO::reco_beam_startX, AnaIO::reco_beam_startY, AnaIO::reco_beam_startZ,
                           AnaIO::reco_beam_trackDirX, AnaIO::reco_beam_trackDirY, AnaIO::reco_beam_trackDirZ,
                           AnaIO::true_beam_startDirX, AnaIO::true_beam_startDirY, AnaIO::true_beam_startDirZ,
                           AnaIO::true_beam_startX, AnaIO::true_beam_startY, AnaIO::true_beam_startZ);
}


bool Manual_beamPos_data(const int event,            const double data_startX,
                         const double data_startY,   const double data_startZ,
                         const double data_dirX,     const double data_dirY,
                         const double data_dirZ,     const double data_BI_X,
                         const double data_BI_Y,     const double data_BI_dirX,
                         const double data_BI_dirY,  const double data_BI_dirZ,
                         const int data_BI_nMomenta, const int data_BI_nTracks) 
{

  //For Data from Owen Goodwin
  const double data_xlow = 0., data_xhigh = 10., data_ylow= -5.;
  const double data_yhigh= 10., data_zlow=30., data_zhigh=35., data_coslow=.93;

  const double deltaX = data_startX - data_BI_X;
  const double deltaY = data_startY - data_BI_Y;
  const double cos = data_BI_dirX*data_dirX + data_BI_dirY*data_dirY +
               data_BI_dirZ*data_dirZ;

  if(data_BI_nMomenta != 1 || data_BI_nTracks != 1)
    return false;

  if( (deltaX < data_xlow) || (deltaX > data_xhigh) )
    return false;

  if ( (deltaY < data_ylow) || (deltaY > data_yhigh) )
    return false;

  if ( (data_startZ < data_zlow) || (data_startZ > data_zhigh) )
    return false;

  if (cos < data_coslow)
    return false;

  return true;

};

bool CutDataBeamPos()
{
  //https://github.com/calcuttj/PionStudies/blob/master/rDataFrame/eventSelection.C
  /*
.Define("passBeamCut", manual_beamPos_data,
            {"event","reco_beam_startX", "reco_beam_startY", "reco_beam_startZ",
            "reco_beam_trackDirX", "reco_beam_trackDirY", "reco_beam_trackDirZ",
            "data_BI_X", "data_BI_Y", "data_BI_dirX", "data_BI_dirY",
            "data_BI_dirZ", "data_BI_nMomenta", "data_BI_nTracks"})
   */

  return Manual_beamPos_data(-999, AnaIO::reco_beam_startX, AnaIO::reco_beam_startY, AnaIO::reco_beam_startZ,
                             AnaIO::reco_beam_trackDirX, AnaIO::reco_beam_trackDirY, AnaIO::reco_beam_trackDirZ,
                             AnaIO::data_BI_X, AnaIO::data_BI_Y, AnaIO::data_BI_dirX, AnaIO::data_BI_dirY,
                             AnaIO::data_BI_dirZ, AnaIO::data_BI_nMomenta, AnaIO::data_BI_nTracks);
  
}


bool CutBeamAllInOne(const bool kmc)
{
  //0. true beam particle //"In data the beam Instrumentation is able to filter for these events but is not inside the MC" so read data file has this
  if(kmc){
    if(AnaIO::true_beam_PDG != 211 && AnaIO::true_beam_PDG != -13){
      return false;
    }
  }
  else{
    const bool data_beamID = CutBeamID((*AnaIO::data_BI_PDG_candidates));
    AnaIO::hCutbeamID->Fill(data_beamID);
    if(!data_beamID){
      return false;
    }
  }

  const int filleventtype = AnaUtils::GetFillEventType();

  //for pi+ it is worse to include this pion-ana-cut (e*p: 6.7->6.3%), for pi0 it is the same e*p, better purity, worse efficiency.
  if(1){//====================================================== switch off pion analysis cuts
        //x. primary beam type 
    AnaIO::hCutBeamType->Fill(AnaIO::reco_beam_type, filleventtype);
    //no effect, shadowed by TMeanStart cut
    if(AnaIO::reco_beam_type!=13){//13: Pandora "track like"
      return false;
    }
    
    //1. beam position MC cut, need MC truth, how is it possible in analysis?
    const bool kBeamPosPass = kmc ? CutMCBeamPos() : CutDataBeamPos();
    AnaIO::hCutBeamPosPass->Fill(kBeamPosPass, filleventtype);
    if(!kBeamPosPass){
      return false;
    }
    //-> now signal purity 138/3537 = 3.9%, 2283 pi+ bea, 801 e+ beam
    
    //2. APA3 
    AnaIO::hCutBeamEndZ->Fill(AnaIO::reco_beam_endZ, filleventtype);
    AnaIO::hCutBeamEndZPass->Fill(!(AnaIO::reco_beam_endZ>=226), filleventtype);
    if(AnaIO::reco_beam_endZ>=226){
      return false;
    }
    //-> now signal purity 135/3143 = 4.3%, 2102 pi+ beam, 801 e+ beam
  }//====================================================== switch off pion analysis cuts

  return true;
}



bool CutBeamdEdx(const double varSignal)
{
  vector<double> startE, lastE;
  AnaUtils::GetdEdx( *AnaIO::reco_beam_calibrated_dEdX, startE, lastE, 2);
  AnaIO::nBeamdEdxCls = startE.size();
  AnaIO::hnBeamdEdxCls->Fill(AnaIO::nBeamdEdxCls);
  if(AnaIO::nBeamdEdxCls<6){
    return false;
  }

  //no bragg peak
  AnaIO::beamStartE0 = startE[0];
  AnaIO::beamStartE1 = startE[1];
  AnaIO::beamStartE2 = startE[2];
  AnaIO::beamStartE3 = startE[3];
  AnaIO::beamStartE4 = startE[4];
  AnaIO::beamStartE5 = startE[5];
  
  AnaIO::beamTMeanStart = AnaUtils::GetTruncatedMean(startE, AnaIO::nBeamdEdxCls-6);
  
  AnaIO::hSignalVsStartE0->Fill(AnaIO::beamStartE0, varSignal);
  AnaIO::hSignalVsStartE1->Fill(AnaIO::beamStartE1, varSignal);
  AnaIO::hSignalVsStartE2->Fill(AnaIO::beamStartE2, varSignal);
  AnaIO::hSignalVsStartE3->Fill(AnaIO::beamStartE3, varSignal);
  AnaIO::hSignalVsStartE4->Fill(AnaIO::beamStartE4, varSignal);
  AnaIO::hSignalVsStartE5->Fill(AnaIO::beamStartE5, varSignal);
  
  AnaIO::hSignalVsTMeanStart->Fill(AnaIO::beamTMeanStart, varSignal);
  
  //has Bragg Peak
  AnaIO::beamLastE0 = lastE[0];
  AnaIO::beamLastE1 = lastE[1];
  AnaIO::beamLastE2 = lastE[2];
  AnaIO::beamLastE3 = lastE[3];
  AnaIO::beamLastE4 = lastE[4];
  AnaIO::beamLastE5 = lastE[5];
  
  AnaIO::beamTMeanLast = AnaUtils::GetTruncatedMean(lastE, 6);
  
  AnaIO::hSignalVsLastE0->Fill(AnaIO::beamLastE0, varSignal);
  AnaIO::hSignalVsLastE1->Fill(AnaIO::beamLastE1, varSignal);
  AnaIO::hSignalVsLastE2->Fill(AnaIO::beamLastE2, varSignal);
  AnaIO::hSignalVsLastE3->Fill(AnaIO::beamLastE3, varSignal);
  AnaIO::hSignalVsLastE4->Fill(AnaIO::beamLastE4, varSignal);
  AnaIO::hSignalVsLastE5->Fill(AnaIO::beamLastE5, varSignal);
  
  //both need to tune for different energy
  //it is so clean that no need to cut on last since there is no Bragg peak form proton any more
  if(AnaIO::beamTMeanStart>2.8){
    return false;
  }
  
  AnaIO::hSignalVsTMeanLast->Fill(AnaIO::beamTMeanLast, varSignal);
  
  AnaIO::hSigAfterVsStartE0->Fill(AnaIO::beamStartE0, varSignal);
  AnaIO::hSigAfterVsStartE1->Fill(AnaIO::beamStartE1, varSignal);
  AnaIO::hSigAfterVsStartE2->Fill(AnaIO::beamStartE2, varSignal);
  
  AnaIO::hSigAfterVsLastE0->Fill(AnaIO::beamLastE0, varSignal);
  AnaIO::hSigAfterVsLastE1->Fill(AnaIO::beamLastE1, varSignal);
  AnaIO::hSigAfterVsLastE2->Fill(AnaIO::beamLastE2, varSignal);
  
  AnaIO::hSigAfterVsStartE3->Fill(AnaIO::beamStartE3, varSignal);
  AnaIO::hSigAfterVsLastE3->Fill(AnaIO::beamLastE3, varSignal);
  AnaIO::hSigAfterVsStartE4->Fill(AnaIO::beamStartE4, varSignal);
  AnaIO::hSigAfterVsLastE4->Fill(AnaIO::beamLastE4, varSignal);
  AnaIO::hSigAfterVsStartE5->Fill(AnaIO::beamStartE5, varSignal);
  AnaIO::hSigAfterVsLastE5->Fill(AnaIO::beamLastE5, varSignal);
  
  return true;
}


//end of namespace
}

#endif
