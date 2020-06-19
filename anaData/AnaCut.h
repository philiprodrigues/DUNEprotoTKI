#ifndef _ANACUT_H_
#define _ANACUT_H_

using namespace std;

namespace AnaCut
{

int GetTruthFromRec(const int recidx, int & pdg, TLorentzVector *& momRefBeam)
{
  pdg = -999;
  momRefBeam = 0x0;
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
      
      momRefBeam = new TLorentzVector(AnaUtils::GetMomentumRefBeam(true, trueidx, pdg==2212));
    }
    
  }

  return trueidx;
}


//int GetNTrack(const bool kpi0, const bool ksig, int & nproton, int & nshower, int & nmichel, TLorentzVector * & leadingProton, TLorentzVector *& leadingPiplus, TLorentzVector *  & leadingPi0, const bool kprint=false, const bool kfill=false)
int GetNTrack(const bool kpi0, const bool ksig, int & nproton, int & nshower, int & nmichel, TLorentzVector *  & leadingPi0, const bool kprint, const bool kfill)
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

  //proton and pion array only keeps the leading one, ie, only [0] will be returned
  //vector<TLorentzVector> protonArray; 
  //vector<TLorentzVector> pionArray; //that is non-proton actually
  vector<TLorentzVector> showerArray;

  static bool kPrintCutInfo = true;

  for(int ii=0; ii<recsize; ii++){

    //__________________________________________ Get Truth information __________________________________________
    int pdg = -999;
    TLorentzVector * truthMomRefBeam = 0x0;
    const int trueidx = GetTruthFromRec(ii, pdg, truthMomRefBeam);

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

    //__________________________________________ Cut PFP before counting  __________________________________________

    //---> need to be done before any counting!!!
    const vector<double>  dedxarray = (*AnaIO::reco_daughter_allTrack_calibrated_dEdX_SCE)[ii];
    const int NdEdx = dedxarray.size();
    if(kfill){
      style::FillInRange(AnaIO::hCutNdEdx, NdEdx, fillstktype);
    }
    const int ndedxcut = kpi0 ? 6 : 16; //need E[3], at least 4 cls
    if(kPrintCutInfo){
      printf("check cut kpi0 %d ndedxcut %d\n", kpi0, ndedxcut);
    }

    if(NdEdx<ndedxcut){
      //do not count it as anything
      continue;
    }
    //<--- need to be done before any counting!!!

    //__________________________________________ Count without continue  __________________________________________

    //track and em scores at 0.5 look reasonable as baseline (after 0.5 the proton fraction look flat)
    //michel is not found (all below 0.5, those above 0.5 are shower and other_type)
    const double trackScore  = (*AnaIO::reco_daughter_PFP_trackScore_collection)[ii];
    if(trackScore>0.5 && (*AnaIO::reco_daughter_allTrack_ID)[ii]!=-1 ){
      if( (*AnaIO::reco_daughter_allTrack_ID)[ii]==-1 ){
        printf("bad track ID! %d\n", ii); exit(1);
      }
      ntracks++;
    }
    const double emScore     = (*AnaIO::reco_daughter_PFP_emScore_collection)[ii];
    //already requiring good shower
    if(emScore>0.5 && (*AnaIO::reco_daughter_allShower_ID)[ii]!=-1){
      if( (*AnaIO::reco_daughter_allShower_ID)[ii]==-1 ){
        printf("bad shower ID! %d\n", ii); exit(1);
      }
      nshower++;
      
      if(AnaIO::reco_daughter_allShower_energy){
        const TVector3 showerDir( (*AnaIO::reco_daughter_allShower_dirX)[ii], (*AnaIO::reco_daughter_allShower_dirY)[ii], (*AnaIO::reco_daughter_allShower_dirZ)[ii] );
        const TVector3 showerMomentum = showerDir.Unit()*(*AnaIO::reco_daughter_allShower_energy)[ii] * 1E-3; //MeV to GeV
        const TLorentzVector showerLv( showerMomentum, showerMomentum.Mag() );
        showerArray.push_back(showerLv);
      }
      else{
        printf("shower energy null!!\n"); exit(1);
      }
    }
    const double michelScore = (*AnaIO::reco_daughter_PFP_michelScore_collection)[ii];
    if(michelScore>0.5){
      nmichel++;
    }


    //--------------------------------------------> need to move inside trackScore !! ->
    //nhit 260 is just clean-up
    const int nhits          = (*AnaIO::reco_daughter_PFP_nHits)[ii];
    /*//they are indeed different: 236 80
    if(nhits != NdEdx){
      printf("nhits!=dedx size %d %d\n", nhits, NdEdx); exit(1);
    }
    */

    //========== proton tagging now!
  
    const double startE2 = NdEdx<3?-999:dedxarray[2];
    const double startE3 = NdEdx<4?-999:dedxarray[3];
    const double lastE2  = NdEdx<3?-999:dedxarray[NdEdx-1-2];
    const double lastE3  = NdEdx<4?-999:dedxarray[NdEdx-1-3];

    const double chi2 = (*AnaIO::reco_daughter_allTrack_Chi2_proton)[ii];
    const double ndof = (*AnaIO::reco_daughter_allTrack_Chi2_ndof)[ii];
    const double Chi2NDF = chi2/(ndof+1E-10);

    int recParticleType = -999;
    const double cutNH = 260;
    if(kpi0){
      //all signal protons have nhit below 260
      //all signal protons have startE3 > 9
      //with Chi2NDF cut  44/191 = 23% purity; without 46/197 = 23% purity, slightly higher efficiency
      //test if(startE2>10 && nhits<260 && startE3>9 && Chi2NDF<50){
      const double cutSE2 = 10;
      const double cutSE3 = 9;
      if(startE2>cutSE2 && nhits<cutNH && startE3>cutSE3){
        recParticleType = AnaUtils::gkProton;
      }

      if(kPrintCutInfo){
        printf("check cut kpi0 %d proton tag startE2 %.2f nhits %.0f startE3 %.2f\n", kpi0, cutSE2, cutNH, cutSE3); 
      }
    }
    else{
      const double cutCHI = 50;
      if(Chi2NDF<cutCHI && nhits<cutNH){
        recParticleType = AnaUtils::gkProton;
      }
      else{
        recParticleType = AnaUtils::gkPiPlus;
      }

      if(kPrintCutInfo){
        printf("check cut kpi0 %d proton tag Chi2NDF %.2f nhits %.0f\n", kpi0, cutCHI, cutNH);
      }
    }

    if(recParticleType==AnaUtils::gkProton){
      nproton++;
    }
    //========== proton tagging done!
    //--------------------------------------------> need to move inside trackScore !! <---
   
    //__________________________________________ Get Reco kinematics  __________________________________________
    const TLorentzVector recMomRefBeam = AnaUtils::GetMomentumRefBeam(false, ii, recParticleType==AnaUtils::gkProton);

    //__________________________________________ Print and Fill  __________________________________________
    if(kprint){
      printf("test sig %d rii %d/%d tii %4d pdg %4d sE2 %6.1f sE3 %6.1f lE2 %6.1f lE3 %6.1f nhi %4d tkS %6.1f emS %6.1f miS %6.1f sum %6.1f chf %6.1f\n", ksig, ii, recsize, trueidx, pdg, startE2, startE3, lastE2, lastE3, nhits, trackScore, emScore, michelScore, trackScore+emScore+michelScore, Chi2NDF);
    }

    if(kfill){
      const double momentumRes = truthMomRefBeam? recMomRefBeam.P()/truthMomRefBeam->P()-1 : -999;
      const double thetaRes    = truthMomRefBeam? (recMomRefBeam.Theta()-truthMomRefBeam->Theta())*TMath::RadToDeg() : -999;

      if(recParticleType==AnaUtils::gkProton){
        if(pdg==2212){
          style::FillInRange(AnaIO::hProtonMomentumRes, truthMomRefBeam->P(), momentumRes);
          style::FillInRange(AnaIO::hProtonThetaRes, truthMomRefBeam->Theta()*TMath::RadToDeg(), thetaRes);
        }
        style::FillInRange(AnaIO::hRecProtonMomentum, recMomRefBeam.P(), fillstktype);
        style::FillInRange(AnaIO::hRecProtonTheta, recMomRefBeam.Theta()*TMath::RadToDeg(), fillstktype);
        style::FillInRange(AnaIO::hRecProtonLastE2, lastE2, fillstktype);
        style::FillInRange(AnaIO::hRecProtonLastE3, lastE3, fillstktype);
      }
      else if(recParticleType==AnaUtils::gkPiPlus){
        if(pdg==211){
          style::FillInRange(AnaIO::hPiMomentumRes, truthMomRefBeam->P(), momentumRes);
          style::FillInRange(AnaIO::hPiThetaRes, truthMomRefBeam->Theta()*TMath::RadToDeg(), thetaRes);
        }
        style::FillInRange(AnaIO::hRecPiplusMomentum, recMomRefBeam.P(), fillstktype);
        style::FillInRange(AnaIO::hRecPiplusTheta, recMomRefBeam.Theta()*TMath::RadToDeg(), fillstktype);
        style::FillInRange(AnaIO::hRecPiplusLastE2, lastE2, fillstktype);
        style::FillInRange(AnaIO::hRecPiplusLastE3, lastE3, fillstktype);
      }
           
      style::FillInRange(AnaIO::hCutnHits, nhits, fillstktype);
      style::FillInRange(AnaIO::hCutChi2NDF, Chi2NDF, fillstktype);
      style::FillInRange(AnaIO::hCuttrackScore, trackScore, fillstktype);
      style::FillInRange(AnaIO::hCutemScore, emScore, fillstktype);
      style::FillInRange(AnaIO::hCutmichelScore, michelScore, fillstktype);

      style::FillInRange(AnaIO::hCutstartE2, startE2, fillstktype);
      style::FillInRange(AnaIO::hCutstartE3, startE3, fillstktype);
    }

    kPrintCutInfo = false;
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
  int cutnproton = 0;
  int cutnshower = 0;
  int cutnmichel = 0;
  TLorentzVector * leadingPi0 = 0x0;  

  const bool kFillBefore = true;
  if(kFillBefore){
    //no print, fill
    AnaIO::nTrack = GetNTrack(kpi0, AnaIO::kSignal, cutnproton, cutnshower, cutnmichel, leadingPi0, false, true);
  }
  else{
    //no print, no fill
    AnaIO::nTrack = GetNTrack(kpi0, AnaIO::kSignal, cutnproton, cutnshower, cutnmichel, leadingPi0, false, false);
  }
  
  const int filleventtype = AnaUtils::GetFillEventType();
  //----------- do cuts
  
  style::FillInRange(AnaIO::hCutnshower, cutnshower, filleventtype);
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
  
  style::FillInRange(AnaIO::hCutnmichel, cutnmichel, filleventtype);
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
  
  style::FillInRange(AnaIO::hCutntrack, AnaIO::nTrack, filleventtype);
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
  
  style::FillInRange(AnaIO::hCutnproton, cutnproton, filleventtype);
  if(cutnproton!=1){
    return false;
  }
  
  //===================================================
  //only checking after-selection distributions
  if(kpi0){
    if(leadingPi0){
      style::FillInRange(AnaIO::hCutMpi0, leadingPi0->M(), filleventtype);
    }
    else{
      printf("leadingpi0 null!!\n"); exit(1);
    }
  }
  
  TLorentzVector *dummypi0 = 0x0;
  //no getting of pizero
  if(kFillBefore){
    //print, no fill
    GetNTrack(kpi0, AnaIO::kSignal, cutnproton, cutnshower, cutnmichel, dummypi0, true, false);
  }
  else{
    //debug mode: major background is -999 shower mocking protons, no efficient cut variables found
    //print, fill
    GetNTrack(kpi0, AnaIO::kSignal, cutnproton, cutnshower, cutnmichel, dummypi0, true, true);
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
    const bool mc_beampass =(AnaIO::true_beam_PDG == 211 || AnaIO::true_beam_PDG == -13);
    style::FillInRange(AnaIO::hCutbeamID, mc_beampass);
    if(!mc_beampass){
      return false;
    }
  }
  else{
    const bool data_beamID = CutBeamID((*AnaIO::data_BI_PDG_candidates));
    style::FillInRange(AnaIO::hCutbeamID, data_beamID);
    if(!data_beamID){
      return false;
    }
  }

  const int filleventtype = AnaUtils::GetFillEventType();

  //for pi+ it is worse to include this pion-ana-cut (e*p: 6.7->6.3%), for pi0 it is the same e*p, better purity, worse efficiency.
  //x. primary beam type 
  style::FillInRange(AnaIO::hCutBeamType, AnaIO::reco_beam_type, filleventtype);
  //no effect, shadowed by TMeanStart cut
  if(AnaIO::reco_beam_type!=13){//13: Pandora "track like"
    return false;
  }
  
  //1. beam position MC cut, need MC truth, how is it possible in analysis?
  const bool kBeamPosPass = kmc ? CutMCBeamPos() : CutDataBeamPos();
  style::FillInRange(AnaIO::hCutBeamPosPass, kBeamPosPass, filleventtype);
  if(!kBeamPosPass){
    return false;
  }
  //-> now signal purity 138/3537 = 3.9%, 2283 pi+ bea, 801 e+ beam
  
  //2. APA3 
  style::FillInRange(AnaIO::hCutBeamEndZ, AnaIO::reco_beam_endZ, filleventtype);
  style::FillInRange(AnaIO::hCutBeamEndZPass, !(AnaIO::reco_beam_endZ>=226), filleventtype);
  if(AnaIO::reco_beam_endZ>=226){
    return false;
  }
  //-> now signal purity 135/3143 = 4.3%, 2102 pi+ beam, 801 e+ beam

  //___________________________ fill beam kinematics _________________________

  const TVector3 recBeam = AnaUtils::GetRecBeamDir();//currently only use dir

  if(kmc){
    const TVector3 truthBeam = AnaUtils::GetTruthBeamFull();
    const double beamthetaRes = (recBeam.Theta()-truthBeam.Theta())*TMath::RadToDeg();//use absolute difference 
    style::FillInRange(AnaIO::hBeamThetaRes, truthBeam.Theta()*TMath::RadToDeg(), beamthetaRes);
  }

  style::FillInRange(AnaIO::hRecBeamTheta, recBeam.Theta()*TMath::RadToDeg(), filleventtype);

  return true;
}



bool CutBeamdEdx(const double varSignal)
{
  vector<double> startE, lastE;
  AnaUtils::GetdEdx( *AnaIO::reco_beam_calibrated_dEdX, startE, lastE, 2);
  AnaIO::nBeamdEdxCls = startE.size();
  style::FillInRange(AnaIO::hnBeamdEdxCls, AnaIO::nBeamdEdxCls);
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
  
  style::FillInRange(AnaIO::hSignalVsStartE0, AnaIO::beamStartE0, varSignal);
  style::FillInRange(AnaIO::hSignalVsStartE1, AnaIO::beamStartE1, varSignal);
  style::FillInRange(AnaIO::hSignalVsStartE2, AnaIO::beamStartE2, varSignal);
  style::FillInRange(AnaIO::hSignalVsStartE3, AnaIO::beamStartE3, varSignal);
  style::FillInRange(AnaIO::hSignalVsStartE4, AnaIO::beamStartE4, varSignal);
  style::FillInRange(AnaIO::hSignalVsStartE5, AnaIO::beamStartE5, varSignal);
  
  style::FillInRange(AnaIO::hSignalVsTMeanStart, AnaIO::beamTMeanStart, varSignal);
  
  //has Bragg Peak
  AnaIO::beamLastE0 = lastE[0];
  AnaIO::beamLastE1 = lastE[1];
  AnaIO::beamLastE2 = lastE[2];
  AnaIO::beamLastE3 = lastE[3];
  AnaIO::beamLastE4 = lastE[4];
  AnaIO::beamLastE5 = lastE[5];
  
  AnaIO::beamTMeanLast = AnaUtils::GetTruncatedMean(lastE, 6);
  
  style::FillInRange(AnaIO::hSignalVsLastE0, AnaIO::beamLastE0, varSignal);
  style::FillInRange(AnaIO::hSignalVsLastE1, AnaIO::beamLastE1, varSignal);
  style::FillInRange(AnaIO::hSignalVsLastE2, AnaIO::beamLastE2, varSignal);
  style::FillInRange(AnaIO::hSignalVsLastE3, AnaIO::beamLastE3, varSignal);
  style::FillInRange(AnaIO::hSignalVsLastE4, AnaIO::beamLastE4, varSignal);
  style::FillInRange(AnaIO::hSignalVsLastE5, AnaIO::beamLastE5, varSignal);
  
  //both need to tune for different energy
  //it is so clean that no need to cut on last since there is no Bragg peak form proton any more
  if(AnaIO::beamTMeanStart>2.8){
    return false;
  }
  
  style::FillInRange(AnaIO::hSignalVsTMeanLast, AnaIO::beamTMeanLast, varSignal);
  
  style::FillInRange(AnaIO::hSigAfterVsStartE0, AnaIO::beamStartE0, varSignal);
  style::FillInRange(AnaIO::hSigAfterVsStartE1, AnaIO::beamStartE1, varSignal);
  style::FillInRange(AnaIO::hSigAfterVsStartE2, AnaIO::beamStartE2, varSignal);
  
  style::FillInRange(AnaIO::hSigAfterVsLastE0, AnaIO::beamLastE0, varSignal);
  style::FillInRange(AnaIO::hSigAfterVsLastE1, AnaIO::beamLastE1, varSignal);
  style::FillInRange(AnaIO::hSigAfterVsLastE2, AnaIO::beamLastE2, varSignal);
  
  style::FillInRange(AnaIO::hSigAfterVsStartE3, AnaIO::beamStartE3, varSignal);
  style::FillInRange(AnaIO::hSigAfterVsLastE3, AnaIO::beamLastE3, varSignal);
  style::FillInRange(AnaIO::hSigAfterVsStartE4, AnaIO::beamStartE4, varSignal);
  style::FillInRange(AnaIO::hSigAfterVsLastE4, AnaIO::beamLastE4, varSignal);
  style::FillInRange(AnaIO::hSigAfterVsStartE5, AnaIO::beamStartE5, varSignal);
  style::FillInRange(AnaIO::hSigAfterVsLastE5, AnaIO::beamLastE5, varSignal);
  
  return true;
}


//end of namespace
}

#endif
