
////////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb  2 10:51:23 2022 by ROOT version 6.22/07
// from TTree Events/Events
// found on file: WZ_UL_2016_MC.root (an MC file for UL2016)
////////////////////////////////////////////////////////////

#ifndef VLLAna_h
#define VLLAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

//Headers needed by this particular selector
#include <vector>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include "TLorentzVector.h"
#include <fstream>
#include <iostream>
#include "TString.h"
#include <bitset>
#include <time.h>

#include "setup/nlohmann/json.hpp"
using json=nlohmann::json;

class VLLAna : public TSelector {
public:
  bool _run3;
private:
  //for Run3
  using iterator     = Int_t;
  using int_or_char  = UChar_t;
  using int_or_short = Short_t;
  using int_or_ushort = UShort_t;

  //for Run2
  //using iterator     = UInt_t;
  //using int_or_char  = Int_t;
  //using int_or_short = Int_t;
  //using int_or_ushort = Int_t;
  
public :
  TTreeReader     fReader_common;      //reads the common branches
  TTreeReader     fReader_commonMC;    //reads the MC branches
  TTreeReader     fReader_Run2;        //reads the Run3 branches
  TTreeReader     fReader_Run2MC;      //reads the Run3 branches
  TTreeReader     fReader_Run3;        //reads the Run3 branches
  TTreeReader     fReader_Run3MC;      //reads the Run3 branches
  TTreeReader     fReader_2016;        //reads the Special 2016 branches
  TTreeReader     fReader_2017;        //reads the Special 2017 branches
  TTreeReader     fReader_2018;        //reads the Special 2018 branches
  TTreeReader     fReader_2022;        //reads the Special 2018 branches  
  TTreeReader     fReader_2017_2018;        //2017&2018 trigger specific branches
  TTreeReader     fReader_2018_2022;        //2018&2022 specific
  TTreeReader     fReader_commonMCSpecial; //reads special branches available in few samples

  TTree          *fChain = 0;    //!pointer to the analyzed TTree or TChain

  
  // Readers to access the data (delete the ones you do not need).
  // Unnecessary branches can be commented out.

  // READ ONLY NECESSARY BRANCHES TO SPEED UP
  
  //#######################################################################
  //
  //                          Rules
  // 1. Run2 and Run3 same branch name and same type  : fReader_common
  // 2. Run2 and Run3 same branch name and diff type  : fReader_Run2 & fReader_Run3
  // 3. Exclusive Run2 branches only                  : fReader_Run2
  // 4. Exclusive Run3 branches only                  : fReader_Run3
  // 5. Yearwise dependent branches                   : fReader_<Year>, eg, fReader_2016
  // 6. Special branches for specific MC samples      : fReader_commonMCspecial
  //
  // 7. General Rules                     MC branches : <readerName>MC
  //
  //########################################################################

  TTreeReaderValue<UInt_t> run = {fReader_common, "run"};
  TTreeReaderValue<UInt_t> luminosityBlock = {fReader_common, "luminosityBlock"};
  TTreeReaderValue<ULong64_t> event = {fReader_common, "event"};

  //Event Flags
  //Bad Charged Candidate flags
  TTreeReaderValue<Bool_t> Flag_BadChargedCandidateFilter = {fReader_common, "Flag_BadChargedCandidateFilter"};
  TTreeReaderValue<Bool_t> Flag_BadChargedCandidateSummer16Filter = {fReader_common, "Flag_BadChargedCandidateSummer16Filter"};

  //BadMuons flags
  TTreeReaderValue<Bool_t> Flag_BadPFMuonDzFilter         = {fReader_common, "Flag_BadPFMuonDzFilter"};
  TTreeReaderValue<Bool_t> Flag_BadPFMuonFilter           = {fReader_common, "Flag_BadPFMuonFilter"};
  TTreeReaderValue<Bool_t> Flag_BadPFMuonSummer16Filter   = {fReader_common, "Flag_BadPFMuonSummer16Filter"};

  //CSC-chamber-related flags
  TTreeReaderValue<Bool_t> Flag_CSCTightHalo2015Filter    = {fReader_common, "Flag_CSCTightHalo2015Filter"};
  TTreeReaderValue<Bool_t> Flag_CSCTightHaloFilter        = {fReader_common, "Flag_CSCTightHaloFilter"};
  TTreeReaderValue<Bool_t> Flag_CSCTightHaloTrkMuUnvetoFilter = {fReader_common, "Flag_CSCTightHaloTrkMuUnvetoFilter"};

  //EcalDeadCell flags
  TTreeReaderValue<Bool_t> Flag_EcalDeadCellBoundaryEnergyFilter = {fReader_common, "Flag_EcalDeadCellBoundaryEnergyFilter"};
  TTreeReaderValue<Bool_t> Flag_EcalDeadCellTriggerPrimitiveFilter = {fReader_common, "Flag_EcalDeadCellTriggerPrimitiveFilter"};

  //Hcal flags
  TTreeReaderValue<Bool_t> Flag_HBHENoiseFilter           = {fReader_common, "Flag_HBHENoiseFilter"};
  TTreeReaderValue<Bool_t> Flag_HBHENoiseIsoFilter        = {fReader_common, "Flag_HBHENoiseIsoFilter"};
  TTreeReaderValue<Bool_t> Flag_HcalStripHaloFilter       = {fReader_common, "Flag_HcalStripHaloFilter"};
  //MET-filter
  TTreeReaderValue<Bool_t> Flag_METFilters                = {fReader_common, "Flag_METFilters"};
  TTreeReaderValue<Bool_t> Flag_chargedHadronTrackResolutionFilter= {fReader_common, "Flag_chargedHadronTrackResolutionFilter"};
  
  //EcalBadCalibration flags
  TTreeReaderValue<Bool_t> Flag_ecalBadCalibFilter        = {fReader_common, "Flag_ecalBadCalibFilter"};
  TTreeReaderValue<Bool_t> Flag_ecalLaserCorrFilter       = {fReader_common, "Flag_ecalLaserCorrFilter"};
  TTreeReaderValue<Bool_t> Flag_eeBadScFilter             = {fReader_common, "Flag_eeBadScFilter"};
  //Beam-halo flags
  TTreeReaderValue<Bool_t> Flag_globalSuperTightHalo2016Filter= {fReader_common, "Flag_globalSuperTightHalo2016Filter"};
  TTreeReaderValue<Bool_t> Flag_globalTightHalo2016Filter = {fReader_common, "Flag_globalTightHalo2016Filter"};

  //Good quality vertices flags
  TTreeReaderValue<Bool_t> Flag_goodVertices              = {fReader_common, "Flag_goodVertices"};

  //Hcal noise filters
  TTreeReaderValue<Bool_t> Flag_hcalLaserEventFilter      = {fReader_common, "Flag_hcalLaserEventFilter"};
  TTreeReaderValue<Bool_t> Flag_hfNoisyHitsFilter         = {fReader_common, "Flag_hfNoisyHitsFilter"};

  //Track related filters
  TTreeReaderValue<Bool_t> Flag_muonBadTrackFilter        = {fReader_common, "Flag_muonBadTrackFilter"};
  TTreeReaderValue<Bool_t> Flag_trkPOGFilters             = {fReader_common, "Flag_trkPOGFilters"};
  TTreeReaderValue<Bool_t> Flag_trkPOG_logErrorTooManyClusters= {fReader_common, "Flag_trkPOG_logErrorTooManyClusters"};
  TTreeReaderValue<Bool_t> Flag_trkPOG_manystripclus53X = {fReader_common, "Flag_trkPOG_manystripclus53X"};
  TTreeReaderValue<Bool_t> Flag_trkPOG_toomanystripclus53X = {fReader_common, "Flag_trkPOG_toomanystripclus53X"};

  
  //Electron
  TTreeReaderValue<iterator> nElectron          = {fReader_common, "nElectron"};
  TTreeReaderArray<int_or_char> Electron_cutBased      = {fReader_common, "Electron_cutBased"};
  TTreeReaderArray<int_or_short> Electron_genPartIdx   = {fReader_commonMC, "Electron_genPartIdx"};
  TTreeReaderArray<int_or_short> Electron_jetIdx       = {fReader_common, "Electron_jetIdx"};
  TTreeReaderArray<int_or_char> Electron_tightCharge   = {fReader_common, "Electron_tightCharge"};
  TTreeReaderArray<Int_t> Electron_charge = {fReader_common, "Electron_charge"};
  TTreeReaderArray<Float_t> Electron_deltaEtaSC = {fReader_common, "Electron_deltaEtaSC"};
  TTreeReaderArray<Float_t> Electron_dxy = {fReader_common, "Electron_dxy"};
  TTreeReaderArray<Float_t> Electron_dxyErr = {fReader_common, "Electron_dxyErr"};
  TTreeReaderArray<Float_t> Electron_dz = {fReader_common, "Electron_dz"};
  TTreeReaderArray<Float_t> Electron_dzErr = {fReader_common, "Electron_dzErr"};
  TTreeReaderArray<Float_t> Electron_eta = {fReader_common, "Electron_eta"};
  TTreeReaderArray<UChar_t> Electron_genPartFlav = {fReader_commonMC, "Electron_genPartFlav"};
  TTreeReaderArray<Float_t> Electron_ip3d = {fReader_common, "Electron_ip3d"};
  TTreeReaderArray<Float_t> Electron_jetPtRelv2 = {fReader_common, "Electron_jetPtRelv2"};
  TTreeReaderArray<Float_t> Electron_mass = {fReader_common, "Electron_mass"};
  TTreeReaderArray<Float_t> Electron_mvaTTH = {fReader_common, "Electron_mvaTTH"};
  TTreeReaderArray<Int_t> Electron_pdgId = {fReader_common, "Electron_pdgId"};
  TTreeReaderArray<Float_t> Electron_pfRelIso03_all = {fReader_common, "Electron_pfRelIso03_all"};
  TTreeReaderArray<Float_t> Electron_pfRelIso03_chg = {fReader_common, "Electron_pfRelIso03_chg"};
  TTreeReaderArray<Float_t> Electron_phi = {fReader_common, "Electron_phi"};
  TTreeReaderArray<Float_t> Electron_pt = {fReader_common, "Electron_pt"};
  TTreeReaderArray<Float_t> Electron_r9 = {fReader_common, "Electron_r9"};
  TTreeReaderArray<Float_t> Electron_scEtOverPt = {fReader_common, "Electron_scEtOverPt"};
  TTreeReaderArray<Float_t> Electron_sieie = {fReader_common, "Electron_sieie"};
  TTreeReaderArray<Float_t> Electron_sip3d = {fReader_common, "Electron_sip3d"};
  TTreeReaderArray<Int_t> Electron_vidNestedWPBitmap = {fReader_common, "Electron_vidNestedWPBitmap"};

  //Unused electron branches
  //TTreeReaderArray<int_or_short> Electron_photonIdx    = {fReader_common, "Electron_photonIdx"};
  //TTreeReaderArray<UChar_t> Electron_cleanmask = {fReader_Run2, "Electron_cleanmask"};
  //TTreeReaderArray<Bool_t> Electron_convVeto = {fReader_common, "Electron_convVeto"};
  //TTreeReaderArray<Bool_t> Electron_cutBased_HEEP = {fReader_common, "Electron_cutBased_HEEP"};
  //TTreeReaderArray<Float_t> Electron_dEscaleDown = {fReader_Run2, "Electron_dEscaleDown"};
  //TTreeReaderArray<Float_t> Electron_dEscaleUp = {fReader_Run2, "Electron_dEscaleUp"};
  //TTreeReaderArray<Float_t> Electron_dEsigmaDown = {fReader_Run2, "Electron_dEsigmaDown"};
  //TTreeReaderArray<Float_t> Electron_dEsigmaUp = {fReader_Run2, "Electron_dEsigmaUp"};
  //TTreeReaderArray<Float_t> Electron_dr03EcalRecHitSumEt = {fReader_common, "Electron_dr03EcalRecHitSumEt"};
  //TTreeReaderArray<Float_t> Electron_dr03HcalDepth1TowerSumEt = {fReader_common, "Electron_dr03HcalDepth1TowerSumEt"};
  //TTreeReaderArray<Float_t> Electron_dr03TkSumPt = {fReader_common, "Electron_dr03TkSumPt"};
  //TTreeReaderArray<Float_t> Electron_dr03TkSumPtHEEP = {fReader_common, "Electron_dr03TkSumPtHEEP"};
  //TTreeReaderArray<Float_t> Electron_eCorr = {fReader_Run2, "Electron_eCorr"};
  //TTreeReaderArray<Float_t> Electron_eInvMinusPInv = {fReader_common, "Electron_eInvMinusPInv"};
  //TTreeReaderArray<Float_t> Electron_energyErr = {fReader_common, "Electron_energyErr"};
  //TTreeReaderArray<Float_t> Electron_hoe = {fReader_common, "Electron_hoe"};
  //TTreeReaderArray<Bool_t> Electron_isPFcand = {fReader_common, "Electron_isPFcand"};
  //TTreeReaderArray<UChar_t> Electron_jetNDauCharged = {fReader_common, "Electron_jetNDauCharged"};
  //TTreeReaderArray<Float_t> Electron_jetRelIso = {fReader_common, "Electron_jetRelIso"};
  //TTreeReaderArray<UChar_t> Electron_lostHits = {fReader_common, "Electron_lostHits"};
  //TTreeReaderArray<Float_t> Electron_miniPFRelIso_all = {fReader_common, "Electron_miniPFRelIso_all"};
  //TTreeReaderArray<Float_t> Electron_miniPFRelIso_chg = {fReader_common, "Electron_miniPFRelIso_chg"};
  //TTreeReaderArray<Float_t> Electron_mvaFall17V2Iso = {fReader_Run2, "Electron_mvaFall17V2Iso"};
  //TTreeReaderArray<Bool_t> Electron_mvaFall17V2Iso_WP80 = {fReader_Run2, "Electron_mvaFall17V2Iso_WP80"};
  //TTreeReaderArray<Bool_t> Electron_mvaFall17V2Iso_WP90 = {fReader_Run2, "Electron_mvaFall17V2Iso_WP90"};
  //TTreeReaderArray<Bool_t> Electron_mvaFall17V2Iso_WPL = {fReader_Run2, "Electron_mvaFall17V2Iso_WPL"};
  //TTreeReaderArray<Float_t> Electron_mvaFall17V2noIso = {fReader_Run2, "Electron_mvaFall17V2noIso"};
  //TTreeReaderArray<Bool_t> Electron_mvaFall17V2noIso_WP80 = {fReader_Run2, "Electron_mvaFall17V2noIso_WP80"};
  //TTreeReaderArray<Bool_t> Electron_mvaFall17V2noIso_WP90 = {fReader_Run2, "Electron_mvaFall17V2noIso_WP90"};
  //TTreeReaderArray<Bool_t> Electron_mvaFall17V2noIso_WPL = {fReader_Run2, "Electron_mvaFall17V2noIso_WPL"};
  //TTreeReaderArray<UChar_t> Electron_seedGain = {fReader_common, "Electron_seedGain"};
  //TTreeReaderArray<Int_t> Electron_vidNestedWPBitmapHEEP = {fReader_common, "Electron_vidNestedWPBitmapHEEP"};
  
  //Extra Run3 Electron branches
  //TTreeReaderArray<Short_t> Electron_fsrPhotonIdx = {fReader_Run3, "Electron_fsrPhotonIdx"};
  //TTreeReaderArray<Float_t> Electron_mvaHZZIso = {fReader_Run3, "Electron_mvaHZZIso"};
  //TTreeReaderArray<Float_t> Electron_mvaIso = {fReader_Run3, "Electron_mvaIso"};
  //TTreeReaderArray<Bool_t> Electron_mvaIso_WP80 = {fReader_Run3, "Electron_mvaIso_WP80"};
  //TTreeReaderArray<Bool_t> Electron_mvaIso_WP90 = {fReader_Run3, "Electron_mvaIso_WP90"};
  //TTreeReaderArray<Float_t> Electron_mvaNoIso = {fReader_Run3, "Electron_mvaNoIso"};
  //TTreeReaderArray<Bool_t> Electron_mvaNoIso_WP80 = {fReader_Run3, "Electron_mvaNoIso_WP80"};
  //TTreeReaderArray<Bool_t> Electron_mvaNoIso_WP90 = {fReader_Run3, "Electron_mvaNoIso_WP90"};
  //TTreeReaderArray<Char_t> Electron_seediEtaOriX = {fReader_Run3, "Electron_seediEtaOriX"};
  //TTreeReaderArray<Int_t> Electron_seediPhiOriY = {fReader_Run3, "Electron_seediPhiOriY"};
  //TTreeReaderArray<Short_t> Electron_svIdx = {fReader_Run3, "Electron_svIdx"};

  //Muon
  TTreeReaderValue<iterator> nMuon              = {fReader_common, "nMuon"};
  //TTreeReaderArray<int_or_short> Muon_fsrPhotonIdx        = {fReader_common, "Muon_fsrPhotonIdx"};
  TTreeReaderArray<int_or_short> Muon_genPartIdx          = {fReader_commonMC, "Muon_genPartIdx"};
  TTreeReaderArray<int_or_short> Muon_jetIdx              = {fReader_common, "Muon_jetIdx"};
  TTreeReaderArray<int_or_char> Muon_nStations            = {fReader_common, "Muon_nStations"};
  TTreeReaderArray<int_or_char> Muon_nTrackerLayers       = {fReader_common, "Muon_nTrackerLayers"};
  TTreeReaderArray<int_or_char> Muon_tightCharge          = {fReader_common, "Muon_tightCharge"};
  TTreeReaderArray<Int_t> Muon_charge = {fReader_common, "Muon_charge"};
  //TTreeReaderArray<UChar_t> Muon_cleanmask = {fReader_Run2, "Muon_cleanmask"};
  TTreeReaderArray<Float_t> Muon_dxy = {fReader_common, "Muon_dxy"};
  TTreeReaderArray<Float_t> Muon_dxyErr = {fReader_common, "Muon_dxyErr"};
  TTreeReaderArray<Float_t> Muon_dxybs = {fReader_common, "Muon_dxybs"};
  TTreeReaderArray<Float_t> Muon_dz = {fReader_common, "Muon_dz"};
  TTreeReaderArray<Float_t> Muon_dzErr = {fReader_common, "Muon_dzErr"};
  TTreeReaderArray<Float_t> Muon_eta = {fReader_common, "Muon_eta"};
  TTreeReaderArray<UChar_t> Muon_genPartFlav = {fReader_commonMC, "Muon_genPartFlav"};
  //TTreeReaderArray<UChar_t> Muon_highPtId = {fReader_common, "Muon_highPtId"};
  //TTreeReaderArray<Bool_t> Muon_highPurity = {fReader_common, "Muon_highPurity"};
  //TTreeReaderArray<Bool_t> Muon_inTimeMuon = {fReader_common, "Muon_inTimeMuon"};
  TTreeReaderArray<Float_t> Muon_ip3d = {fReader_common, "Muon_ip3d"};
  //TTreeReaderArray<Bool_t> Muon_isGlobal = {fReader_common, "Muon_isGlobal"};
  //TTreeReaderArray<Bool_t> Muon_isPFcand = {fReader_common, "Muon_isPFcand"};
  //TTreeReaderArray<Bool_t> Muon_isStandalone = {fReader_common, "Muon_isStandalone"};
  //TTreeReaderArray<Bool_t> Muon_isTracker = {fReader_common, "Muon_isTracker"};
  //TTreeReaderArray<UChar_t> Muon_jetNDauCharged = {fReader_common, "Muon_jetNDauCharged"};
  TTreeReaderArray<Float_t> Muon_jetPtRelv2 = {fReader_common, "Muon_jetPtRelv2"};
  TTreeReaderArray<Float_t> Muon_jetRelIso = {fReader_common, "Muon_jetRelIso"};
  TTreeReaderArray<Bool_t> Muon_looseId = {fReader_common, "Muon_looseId"};
  TTreeReaderArray<Float_t> Muon_mass = {fReader_common, "Muon_mass"};
  TTreeReaderArray<Bool_t> Muon_mediumId = {fReader_common, "Muon_mediumId"};
  //TTreeReaderArray<Bool_t> Muon_mediumPromptId = {fReader_common, "Muon_mediumPromptId"};
  //TTreeReaderArray<UChar_t> Muon_miniIsoId = {fReader_common, "Muon_miniIsoId"};
  //TTreeReaderArray<Float_t> Muon_miniPFRelIso_all = {fReader_common, "Muon_miniPFRelIso_all"};
  //TTreeReaderArray<Float_t> Muon_miniPFRelIso_chg = {fReader_common, "Muon_miniPFRelIso_chg"};
  //TTreeReaderArray<UChar_t> Muon_multiIsoId = {fReader_common, "Muon_multiIsoId"};
  TTreeReaderArray<UChar_t> Muon_mvaId = {fReader_Run2, "Muon_mvaId"};
  //TTreeReaderArray<Float_t> Muon_mvaLowPt = {fReader_common, "Muon_mvaLowPt"};
  //TTreeReaderArray<UChar_t> Muon_mvaLowPtId = {fReader_Run2, "Muon_mvaLowPtId"};
  TTreeReaderArray<Float_t> Muon_mvaTTH = {fReader_common, "Muon_mvaTTH"};
  TTreeReaderArray<Int_t> Muon_pdgId = {fReader_common, "Muon_pdgId"};
  // TTreeReaderArray<UChar_t> Muon_pfIsoId = {fReader_common, "Muon_pfIsoId"};
  TTreeReaderArray<Float_t> Muon_pfRelIso03_all = {fReader_common, "Muon_pfRelIso03_all"};
  TTreeReaderArray<Float_t> Muon_pfRelIso03_chg = {fReader_common, "Muon_pfRelIso03_chg"};
  TTreeReaderArray<Float_t> Muon_pfRelIso04_all = {fReader_common, "Muon_pfRelIso04_all"};
  TTreeReaderArray<Float_t> Muon_phi = {fReader_common, "Muon_phi"};
  TTreeReaderArray<Float_t> Muon_pt = {fReader_common, "Muon_pt"};
  //TTreeReaderArray<Float_t> Muon_ptErr = {fReader_common, "Muon_ptErr"};
  //TTreeReaderArray<UChar_t> Muon_puppiIsoId = {fReader_common, "Muon_puppiIsoId"};
  //TTreeReaderArray<Float_t> Muon_segmentComp = {fReader_common, "Muon_segmentComp"};
  TTreeReaderArray<Float_t> Muon_sip3d = {fReader_common, "Muon_sip3d"};
  //TTreeReaderArray<Bool_t> Muon_softId = {fReader_common, "Muon_softId"};
  //TTreeReaderArray<Float_t> Muon_softMva = {fReader_common, "Muon_softMva"};
  //TTreeReaderArray<Bool_t> Muon_softMvaId = {fReader_common, "Muon_softMvaId"};
  TTreeReaderArray<Bool_t> Muon_tightId = {fReader_common, "Muon_tightId"};
  TTreeReaderArray<UChar_t> Muon_tkIsoId = {fReader_common, "Muon_tkIsoId"};
  TTreeReaderArray<Float_t> Muon_tkRelIso = {fReader_common, "Muon_tkRelIso"};
  //TTreeReaderArray<Bool_t> Muon_triggerIdLoose = {fReader_common, "Muon_triggerIdLoose"};
  //TTreeReaderArray<Float_t> Muon_tunepRelPt = {fReader_common, "Muon_tunepRelPt"};
  //Extra Run3 branches
  //TTreeReaderArray<Float_t> Muon_bsConstrainedChi2 = {fReader_Run3, "Muon_bsConstrainedChi2"};
  //TTreeReaderArray<Float_t> Muon_bsConstrainedPt = {fReader_Run3, "Muon_bsConstrainedPt"};
  //TTreeReaderArray<Float_t> Muon_bsConstrainedPtErr = {fReader_Run3, "Muon_bsConstrainedPtErr"};
  //TTreeReaderArray<Float_t> Muon_mvaMuID = {fReader_Run3, "Muon_mvaMuID"};
  //TTreeReaderArray<UChar_t> Muon_mvaMuID_WP = {fReader_Run3, "Muon_mvaMuID_WP"};
  //TTreeReaderArray<Short_t> Muon_svIdx = {fReader_Run3, "Muon_svIdx"};
  
  //Tau
  TTreeReaderValue<iterator> nTau               = {fReader_common, "nTau"};
  TTreeReaderArray<int_or_short> Tau_charge               = {fReader_common, "Tau_charge"};
  TTreeReaderArray<int_or_char> Tau_decayMode             = {fReader_common, "Tau_decayMode"};
  TTreeReaderArray<int_or_short> Tau_genPartIdx           = {fReader_commonMC, "Tau_genPartIdx"};
  TTreeReaderArray<int_or_short> Tau_jetIdx               = {fReader_common, "Tau_jetIdx"};
  //TTreeReaderArray<Float_t> Tau_chargedIso          = {fReader_common, "Tau_chargedIso"};
  //TTreeReaderArray<UChar_t> Tau_cleanmask           = {fReader_Run2, "Tau_cleanmask"};
  TTreeReaderArray<Float_t> Tau_dxy                 = {fReader_common, "Tau_dxy"};
  TTreeReaderArray<Float_t> Tau_dz                  = {fReader_common, "Tau_dz"};
  TTreeReaderArray<Float_t> Tau_eta                 = {fReader_common, "Tau_eta"};
  TTreeReaderArray<UChar_t> Tau_genPartFlav         = {fReader_commonMC, "Tau_genPartFlav"};
  //TTreeReaderArray<Bool_t> Tau_idAntiEleDeadECal    = {fReader_common, "Tau_idAntiEleDeadECal"};
  //TTreeReaderArray<UChar_t> Tau_idAntiMu            = {fReader_common, "Tau_idAntiMu"};
  //TTreeReaderArray<Bool_t> Tau_idDecayModeOldDMs       = {fReader_common, "Tau_idDecayModeOldDMs"};
  TTreeReaderArray<UChar_t> Tau_idDeepTau2017v2p1VSe   = {fReader_common, "Tau_idDeepTau2017v2p1VSe"};
  TTreeReaderArray<UChar_t> Tau_idDeepTau2017v2p1VSjet = {fReader_common, "Tau_idDeepTau2017v2p1VSjet"};
  TTreeReaderArray<UChar_t> Tau_idDeepTau2017v2p1VSmu  = {fReader_common, "Tau_idDeepTau2017v2p1VSmu"};
  //TTreeReaderArray<Float_t> Tau_leadTkDeltaEta         = {fReader_common, "Tau_leadTkDeltaEta"};
  //TTreeReaderArray<Float_t> Tau_leadTkDeltaPhi         = {fReader_common, "Tau_leadTkDeltaPhi"};
  //TTreeReaderArray<Float_t> Tau_leadTkPtOverTauPt      = {fReader_common, "Tau_leadTkPtOverTauPt"};
  TTreeReaderArray<Float_t> Tau_mass                   = {fReader_common, "Tau_mass"};
  TTreeReaderArray<Float_t> Tau_neutralIso             = {fReader_common, "Tau_neutralIso"};
  TTreeReaderArray<Float_t> Tau_phi                    = {fReader_common, "Tau_phi"};
  //TTreeReaderArray<Float_t> Tau_photonsOutsideSignalCone = {fReader_common, "Tau_photonsOutsideSignalCone"};
  TTreeReaderArray<Float_t> Tau_pt                       = {fReader_common, "Tau_pt"};
  TTreeReaderArray<Float_t> Tau_puCorr                   = {fReader_common, "Tau_puCorr"};
  //TTreeReaderArray<Float_t> Tau_rawDeepTau2017v2p1VSe    = {fReader_common, "Tau_rawDeepTau2017v2p1VSe"};
  //TTreeReaderArray<Float_t> Tau_rawDeepTau2017v2p1VSjet  = {fReader_common, "Tau_rawDeepTau2017v2p1VSjet"};
  //TTreeReaderArray<Float_t> Tau_rawDeepTau2017v2p1VSmu   = {fReader_common, "Tau_rawDeepTau2017v2p1VSmu"};
  //TTreeReaderArray<Float_t> Tau_rawIso                   = {fReader_common, "Tau_rawIso"};
  //TTreeReaderArray<Float_t> Tau_rawIsodR03               = {fReader_common, "Tau_rawIsodR03"};
  //Extra Run3 branches
  TTreeReaderArray<Short_t> Tau_decayModePNet            = {fReader_Run3, "Tau_decayModePNet"};
  //TTreeReaderArray<Short_t> Tau_eleIdx                   = {fReader_Run3, "Tau_eleIdx"};
  TTreeReaderArray<Bool_t> Tau_idDecayModeNewDMs         = {fReader_Run3, "Tau_idDecayModeNewDMs"};
  TTreeReaderArray<UChar_t> Tau_idDeepTau2018v2p5VSe     = {fReader_Run3, "Tau_idDeepTau2018v2p5VSe"};
  TTreeReaderArray<UChar_t> Tau_idDeepTau2018v2p5VSjet   = {fReader_Run3, "Tau_idDeepTau2018v2p5VSjet"};
  TTreeReaderArray<UChar_t> Tau_idDeepTau2018v2p5VSmu    = {fReader_Run3, "Tau_idDeepTau2018v2p5VSmu"};
  TTreeReaderArray<Short_t> Tau_muIdx                    = {fReader_Run3, "Tau_muIdx"};
  TTreeReaderArray<UChar_t> Tau_nSVs                     = {fReader_Run3, "Tau_nSVs"};
  TTreeReaderArray<Float_t> Tau_probDM0PNet              = {fReader_Run3, "Tau_probDM0PNet"};
  TTreeReaderArray<Float_t> Tau_probDM10PNet             = {fReader_Run3, "Tau_probDM10PNet"};
  TTreeReaderArray<Float_t> Tau_probDM11PNet             = {fReader_Run3, "Tau_probDM11PNet"};
  TTreeReaderArray<Float_t> Tau_probDM1PNet              = {fReader_Run3, "Tau_probDM1PNet"};
  TTreeReaderArray<Float_t> Tau_probDM2PNet              = {fReader_Run3, "Tau_probDM2PNet"};
  TTreeReaderArray<Float_t> Tau_ptCorrPNet               = {fReader_Run3, "Tau_ptCorrPNet"};
  TTreeReaderArray<Float_t> Tau_qConfPNet                = {fReader_Run3, "Tau_qConfPNet"};
  //TTreeReaderArray<Float_t> Tau_rawDeepTau2018v2p5VSe    = {fReader_Run3, "Tau_rawDeepTau2018v2p5VSe"};
  //TTreeReaderArray<Float_t> Tau_rawDeepTau2018v2p5VSjet  = {fReader_Run3, "Tau_rawDeepTau2018v2p5VSjet"};
  //TTreeReaderArray<Float_t> Tau_rawDeepTau2018v2p5VSmu   = {fReader_Run3, "Tau_rawDeepTau2018v2p5VSmu"};
  //TTreeReaderArray<Float_t> Tau_rawPNetVSe               = {fReader_Run3, "Tau_rawPNetVSe"};
  //TTreeReaderArray<Float_t> Tau_rawPNetVSjet             = {fReader_Run3, "Tau_rawPNetVSjet"};
  //TTreeReaderArray<Float_t> Tau_rawPNetVSmu              = {fReader_Run3, "Tau_rawPNetVSmu"};
  TTreeReaderArray<Short_t> Tau_svIdx1                   = {fReader_Run3, "Tau_svIdx1"};
  TTreeReaderArray<Short_t> Tau_svIdx2                   = {fReader_Run3, "Tau_svIdx2"};

  //Jet
  TTreeReaderValue<iterator> nJet               = {fReader_common, "nJet"};
  //TTreeReaderArray<int_or_short> Jet_electronIdx1         = {fReader_common, "Jet_electronIdx1"};
  //TTreeReaderArray<int_or_short> Jet_electronIdx2         = {fReader_common, "Jet_electronIdx2"};
  TTreeReaderArray<int_or_short> Jet_genJetIdx            = {fReader_commonMC, "Jet_genJetIdx"};
  TTreeReaderArray<int_or_char> Jet_hadronFlavour         = {fReader_common, "Jet_hadronFlavour"};
  TTreeReaderArray<int_or_char> Jet_jetId                 = {fReader_common, "Jet_jetId"};
  //TTreeReaderArray<int_or_short> Jet_muonIdx1             = {fReader_common, "Jet_muonIdx1"};
  //TTreeReaderArray<int_or_short> Jet_muonIdx2             = {fReader_common, "Jet_muonIdx2"};
  TTreeReaderArray<int_or_char> Jet_nElectrons            = {fReader_common, "Jet_nElectrons"};
  TTreeReaderArray<int_or_char> Jet_nMuons                = {fReader_common, "Jet_nMuons"};
  TTreeReaderArray<int_or_short> Jet_partonFlavour        = {fReader_common, "Jet_partonFlavour"};
  TTreeReaderArray<Float_t> Jet_area            = {fReader_common, "Jet_area"};
  //TTreeReaderArray<Float_t> Jet_bRegCorr        = {fReader_Run2, "Jet_bRegCorr"};
  //TTreeReaderArray<Float_t> Jet_bRegRes         = {fReader_Run2, "Jet_bRegRes"};
  //TTreeReaderArray<Float_t> Jet_btagCSVV2       = {fReader_Run2, "Jet_btagCSVV2"};
  //TTreeReaderArray<Float_t> Jet_btagDeepB       = {fReader_Run2, "Jet_btagDeepB"};
  //TTreeReaderArray<Float_t> Jet_btagDeepCvB     = {fReader_Run2, "Jet_btagDeepCvB"};
  //TTreeReaderArray<Float_t> Jet_btagDeepCvL     = {fReader_Run2, "Jet_btagDeepCvL"};
  TTreeReaderArray<Float_t> Jet_btagDeepFlavB   = {fReader_common, "Jet_btagDeepFlavB"};
  TTreeReaderArray<Float_t> Jet_btagDeepFlavCvB = {fReader_common, "Jet_btagDeepFlavCvB"};
  TTreeReaderArray<Float_t> Jet_btagDeepFlavCvL = {fReader_common, "Jet_btagDeepFlavCvL"};
  TTreeReaderArray<Float_t> Jet_btagDeepFlavQG  = {fReader_common, "Jet_btagDeepFlavQG"};
  //TTreeReaderArray<Float_t> Jet_cRegCorr        = {fReader_Run2, "Jet_cRegCorr"};
  //TTreeReaderArray<Float_t> Jet_cRegRes         = {fReader_Run2, "Jet_cRegRes"};
  TTreeReaderArray<Float_t> Jet_chEmEF          = {fReader_common, "Jet_chEmEF"};
  //TTreeReaderArray<Float_t> Jet_chFPV0EF        = {fReader_Run2, "Jet_chFPV0EF"};
  TTreeReaderArray<Float_t> Jet_chHEF           = {fReader_common, "Jet_chHEF"};
  //TTreeReaderArray<UChar_t> Jet_cleanmask       = {fReader_Run2, "Jet_cleanmask"};
  TTreeReaderArray<Float_t> Jet_eta             = {fReader_common, "Jet_eta"};
  // TTreeReaderArray<Int_t> Jet_hfadjacentEtaStripsSize = {fReader_common, "Jet_hfadjacentEtaStripsSize"};
  //TTreeReaderArray<Int_t> Jet_hfcentralEtaStripSize   = {fReader_common, "Jet_hfcentralEtaStripSize"};
  //TTreeReaderArray<Float_t> Jet_hfsigmaEtaEta   = {fReader_common, "Jet_hfsigmaEtaEta"};
  //TTreeReaderArray<Float_t> Jet_hfsigmaPhiPhi   = {fReader_common, "Jet_hfsigmaPhiPhi"};
  TTreeReaderArray<Float_t> Jet_mass            = {fReader_common, "Jet_mass"};
  TTreeReaderArray<Float_t> Jet_muEF            = {fReader_common, "Jet_muEF"};
  //TTreeReaderArray<Float_t> Jet_muonSubtrFactor = {fReader_common, "Jet_muonSubtrFactor"};
  TTreeReaderArray<UChar_t> Jet_nConstituents   = {fReader_common, "Jet_nConstituents"};
  TTreeReaderArray<Float_t> Jet_neEmEF          = {fReader_common, "Jet_neEmEF"};
  TTreeReaderArray<Float_t> Jet_neHEF           = {fReader_common, "Jet_neHEF"};
  TTreeReaderArray<Float_t> Jet_phi             = {fReader_common, "Jet_phi"};
  TTreeReaderArray<Float_t> Jet_pt              = {fReader_common, "Jet_pt"};
  TTreeReaderArray<Int_t> Jet_puId              = {fReader_Run2, "Jet_puId"};
  TTreeReaderArray<Float_t> Jet_puIdDisc        = {fReader_Run2, "Jet_puIdDisc"};
  TTreeReaderArray<Float_t> Jet_qgl             = {fReader_Run2, "Jet_qgl"};
  TTreeReaderArray<Float_t> Jet_rawFactor       = {fReader_common, "Jet_rawFactor"};
  //Extra Run3 branches
  TTreeReaderArray<Float_t> Jet_PNetRegPtRawCorr         = {fReader_Run3, "Jet_PNetRegPtRawCorr"};
  TTreeReaderArray<Float_t> Jet_PNetRegPtRawCorrNeutrino = {fReader_Run3, "Jet_PNetRegPtRawCorrNeutrino"};
  TTreeReaderArray<Float_t> Jet_PNetRegPtRawRes      = {fReader_Run3, "Jet_PNetRegPtRawRes"};
  TTreeReaderArray<Float_t> Jet_btagPNetB            = {fReader_Run3, "Jet_btagPNetB"};
  TTreeReaderArray<Float_t> Jet_btagPNetCvB          = {fReader_Run3, "Jet_btagPNetCvB"};
  TTreeReaderArray<Float_t> Jet_btagPNetCvL          = {fReader_Run3, "Jet_btagPNetCvL"};
  TTreeReaderArray<Float_t> Jet_btagPNetQvG          = {fReader_Run3, "Jet_btagPNetQvG"};
  TTreeReaderArray<Float_t> Jet_btagPNetTauVJet      = {fReader_Run3, "Jet_btagPNetTauVJet"};
  TTreeReaderArray<Float_t> Jet_btagRobustParTAK4B   = {fReader_Run3, "Jet_btagRobustParTAK4B"};
  TTreeReaderArray<Float_t> Jet_btagRobustParTAK4CvB = {fReader_Run3, "Jet_btagRobustParTAK4CvB"};
  TTreeReaderArray<Float_t> Jet_btagRobustParTAK4CvL = {fReader_Run3, "Jet_btagRobustParTAK4CvL"};
  TTreeReaderArray<Float_t> Jet_btagRobustParTAK4QG  = {fReader_Run3, "Jet_btagRobustParTAK4QG"};
  TTreeReaderArray<UChar_t> Jet_nSVs                 = {fReader_Run3, "Jet_nSVs"};
  //TTreeReaderArray<Short_t> Jet_svIdx1               = {fReader_Run3, "Jet_svIdx1"};
  //TTreeReaderArray<Short_t> Jet_svIdx2               = {fReader_Run3, "Jet_svIdx2"};

  
  //FatJet
  TTreeReaderValue<iterator> nFatJet            = {fReader_common, "nFatJet"};
  TTreeReaderArray<int_or_short> FatJet_electronIdx3SJ = {fReader_common, "FatJet_electronIdx3SJ"};
  TTreeReaderArray<int_or_short> FatJet_genJetAK8Idx   = {fReader_commonMC, "FatJet_genJetAK8Idx"};
  TTreeReaderArray<int_or_char> FatJet_hadronFlavour   = {fReader_common, "FatJet_hadronFlavour"};
  TTreeReaderArray<int_or_char> FatJet_jetId           = {fReader_common, "FatJet_jetId"};
  TTreeReaderArray<int_or_short> FatJet_muonIdx3SJ     = {fReader_common, "FatJet_muonIdx3SJ"};
  TTreeReaderArray<int_or_short> FatJet_subJetIdx1     = {fReader_common, "FatJet_subJetIdx1"};
  TTreeReaderArray<int_or_short> FatJet_subJetIdx2     = {fReader_common, "FatJet_subJetIdx2"};
  TTreeReaderArray<Float_t> FatJet_area = {fReader_common, "FatJet_area"};
  TTreeReaderArray<Float_t> FatJet_btagCSVV2 = {fReader_Run2, "FatJet_btagCSVV2"};
  TTreeReaderArray<Float_t> FatJet_btagDDBvLV2 = {fReader_common, "FatJet_btagDDBvLV2"};
  TTreeReaderArray<Float_t> FatJet_btagDDCvBV2 = {fReader_common, "FatJet_btagDDCvBV2"};
  TTreeReaderArray<Float_t> FatJet_btagDDCvLV2 = {fReader_common, "FatJet_btagDDCvLV2"};
  TTreeReaderArray<Float_t> FatJet_btagDeepB = {fReader_common, "FatJet_btagDeepB"};
  TTreeReaderArray<Float_t> FatJet_btagHbb = {fReader_common, "FatJet_btagHbb"};
  TTreeReaderArray<Float_t> FatJet_deepTagMD_H4qvsQCD = {fReader_Run2, "FatJet_deepTagMD_H4qvsQCD"};
  TTreeReaderArray<Float_t> FatJet_deepTagMD_HbbvsQCD = {fReader_Run2, "FatJet_deepTagMD_HbbvsQCD"};
  TTreeReaderArray<Float_t> FatJet_deepTagMD_TvsQCD = {fReader_Run2, "FatJet_deepTagMD_TvsQCD"};
  TTreeReaderArray<Float_t> FatJet_deepTagMD_WvsQCD = {fReader_Run2, "FatJet_deepTagMD_WvsQCD"};
  TTreeReaderArray<Float_t> FatJet_deepTagMD_ZHbbvsQCD = {fReader_Run2, "FatJet_deepTagMD_ZHbbvsQCD"};
  TTreeReaderArray<Float_t> FatJet_deepTagMD_ZHccvsQCD = {fReader_Run2, "FatJet_deepTagMD_ZHccvsQCD"};
  TTreeReaderArray<Float_t> FatJet_deepTagMD_ZbbvsQCD = {fReader_Run2, "FatJet_deepTagMD_ZbbvsQCD"};
  TTreeReaderArray<Float_t> FatJet_deepTagMD_ZvsQCD = {fReader_Run2, "FatJet_deepTagMD_ZvsQCD"};
  TTreeReaderArray<Float_t> FatJet_deepTagMD_bbvsLight = {fReader_Run2, "FatJet_deepTagMD_bbvsLight"};
  TTreeReaderArray<Float_t> FatJet_deepTagMD_ccvsLight = {fReader_Run2, "FatJet_deepTagMD_ccvsLight"};
  TTreeReaderArray<Float_t> FatJet_deepTag_H = {fReader_Run2, "FatJet_deepTag_H"};
  TTreeReaderArray<Float_t> FatJet_deepTag_QCD = {fReader_Run2, "FatJet_deepTag_QCD"};
  TTreeReaderArray<Float_t> FatJet_deepTag_QCDothers = {fReader_Run2, "FatJet_deepTag_QCDothers"};
  TTreeReaderArray<Float_t> FatJet_deepTag_TvsQCD = {fReader_Run2, "FatJet_deepTag_TvsQCD"};
  TTreeReaderArray<Float_t> FatJet_deepTag_WvsQCD = {fReader_Run2, "FatJet_deepTag_WvsQCD"};
  TTreeReaderArray<Float_t> FatJet_deepTag_ZvsQCD = {fReader_Run2, "FatJet_deepTag_ZvsQCD"};
  TTreeReaderArray<Float_t> FatJet_eta = {fReader_common, "FatJet_eta"};
  TTreeReaderArray<Float_t> FatJet_lsf3 = {fReader_common, "FatJet_lsf3"};
  TTreeReaderArray<Float_t> FatJet_mass = {fReader_common, "FatJet_mass"};
  TTreeReaderArray<Float_t> FatJet_msoftdrop = {fReader_common, "FatJet_msoftdrop"};
  TTreeReaderArray<Float_t> FatJet_n2b1 = {fReader_common, "FatJet_n2b1"};
  TTreeReaderArray<Float_t> FatJet_n3b1 = {fReader_common, "FatJet_n3b1"};
  TTreeReaderArray<UChar_t> FatJet_nBHadrons = {fReader_common, "FatJet_nBHadrons"};
  TTreeReaderArray<UChar_t> FatJet_nCHadrons = {fReader_common, "FatJet_nCHadrons"};
  TTreeReaderArray<UChar_t> FatJet_nConstituents = {fReader_common, "FatJet_nConstituents"};
  TTreeReaderArray<Float_t> FatJet_particleNetMD_QCD = {fReader_Run2, "FatJet_particleNetMD_QCD"};
  TTreeReaderArray<Float_t> FatJet_particleNetMD_Xbb = {fReader_Run2, "FatJet_particleNetMD_Xbb"};
  TTreeReaderArray<Float_t> FatJet_particleNetMD_Xcc = {fReader_Run2, "FatJet_particleNetMD_Xcc"};
  TTreeReaderArray<Float_t> FatJet_particleNetMD_Xqq = {fReader_Run2, "FatJet_particleNetMD_Xqq"};
  TTreeReaderArray<Float_t> FatJet_particleNet_H4qvsQCD = {fReader_Run2, "FatJet_particleNet_H4qvsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_HbbvsQCD = {fReader_Run2, "FatJet_particleNet_HbbvsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_HccvsQCD = {fReader_Run2, "FatJet_particleNet_HccvsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_QCD = {fReader_common, "FatJet_particleNet_QCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_TvsQCD = {fReader_Run2, "FatJet_particleNet_TvsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_WvsQCD = {fReader_Run2, "FatJet_particleNet_WvsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_ZvsQCD = {fReader_Run2, "FatJet_particleNet_ZvsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_mass = {fReader_Run2, "FatJet_particleNet_mass"};
  TTreeReaderArray<Float_t> FatJet_phi = {fReader_common, "FatJet_phi"};
  TTreeReaderArray<Float_t> FatJet_pt = {fReader_common, "FatJet_pt"};
  TTreeReaderArray<Float_t> FatJet_rawFactor = {fReader_common, "FatJet_rawFactor"};
  TTreeReaderArray<Float_t> FatJet_tau1 = {fReader_common, "FatJet_tau1"};
  TTreeReaderArray<Float_t> FatJet_tau2 = {fReader_common, "FatJet_tau2"};
  TTreeReaderArray<Float_t> FatJet_tau3 = {fReader_common, "FatJet_tau3"};
  TTreeReaderArray<Float_t> FatJet_tau4 = {fReader_common, "FatJet_tau4"};
  //Extra Run3 branches
  TTreeReaderArray<Float_t> FatJet_particleNetWithMass_H4qvsQCD = {fReader_Run3, "FatJet_particleNetWithMass_H4qvsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNetWithMass_HbbvsQCD = {fReader_Run3, "FatJet_particleNetWithMass_HbbvsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNetWithMass_HccvsQCD = {fReader_Run3, "FatJet_particleNetWithMass_HccvsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNetWithMass_QCD    = {fReader_Run3, "FatJet_particleNetWithMass_QCD"};
  TTreeReaderArray<Float_t> FatJet_particleNetWithMass_TvsQCD = {fReader_Run3, "FatJet_particleNetWithMass_TvsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNetWithMass_WvsQCD = {fReader_Run3, "FatJet_particleNetWithMass_WvsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNetWithMass_ZvsQCD = {fReader_Run3, "FatJet_particleNetWithMass_ZvsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_QCD0HF = {fReader_Run3, "FatJet_particleNet_QCD0HF"};
  TTreeReaderArray<Float_t> FatJet_particleNet_QCD1HF = {fReader_Run3, "FatJet_particleNet_QCD1HF"};
  TTreeReaderArray<Float_t> FatJet_particleNet_QCD2HF = {fReader_Run3, "FatJet_particleNet_QCD2HF"};
  TTreeReaderArray<Float_t> FatJet_particleNet_XbbVsQCD = {fReader_Run3, "FatJet_particleNet_XbbVsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_XccVsQCD = {fReader_Run3, "FatJet_particleNet_XccVsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_XggVsQCD = {fReader_Run3, "FatJet_particleNet_XggVsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_XqqVsQCD = {fReader_Run3, "FatJet_particleNet_XqqVsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_XteVsQCD = {fReader_Run3, "FatJet_particleNet_XteVsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_XtmVsQCD = {fReader_Run3, "FatJet_particleNet_XtmVsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_XttVsQCD = {fReader_Run3, "FatJet_particleNet_XttVsQCD"};
  TTreeReaderArray<Float_t> FatJet_particleNet_massCorr = {fReader_Run3, "FatJet_particleNet_massCorr"};

  //MET
  TTreeReaderValue<Float_t> MET_MetUnclustEnUpDeltaX = {fReader_common, "MET_MetUnclustEnUpDeltaX"};
  TTreeReaderValue<Float_t> MET_MetUnclustEnUpDeltaY = {fReader_common, "MET_MetUnclustEnUpDeltaY"};
  TTreeReaderValue<Float_t> MET_covXX                = {fReader_common, "MET_covXX"};
  TTreeReaderValue<Float_t> MET_covXY                = {fReader_common, "MET_covXY"};
  TTreeReaderValue<Float_t> MET_covYY                = {fReader_common, "MET_covYY"};
  TTreeReaderValue<Float_t> MET_fiducialGenPhi       = {fReader_commonMC, "MET_fiducialGenPhi"};
  TTreeReaderValue<Float_t> MET_fiducialGenPt        = {fReader_commonMC, "MET_fiducialGenPt"};
  TTreeReaderValue<Float_t> MET_phi                  = {fReader_common, "MET_phi"};
  TTreeReaderValue<Float_t> MET_pt                   = {fReader_common, "MET_pt"};
  TTreeReaderValue<Float_t> MET_significance         = {fReader_common, "MET_significance"};
  TTreeReaderValue<Float_t> MET_sumEt                = {fReader_common, "MET_sumEt"};
  TTreeReaderValue<Float_t> MET_sumPtUnclustered     = {fReader_common, "MET_sumPtUnclustered"};

  //PuppiMET
  TTreeReaderValue<Float_t> PuppiMET_phi         = {fReader_common, "PuppiMET_phi"};
  TTreeReaderValue<Float_t> PuppiMET_phiJERDown  = {fReader_common, "PuppiMET_phiJERDown"};
  TTreeReaderValue<Float_t> PuppiMET_phiJERUp    = {fReader_common, "PuppiMET_phiJERUp"};
  TTreeReaderValue<Float_t> PuppiMET_phiJESDown  = {fReader_common, "PuppiMET_phiJESDown"};
  TTreeReaderValue<Float_t> PuppiMET_phiJESUp    = {fReader_common, "PuppiMET_phiJESUp"};
  TTreeReaderValue<Float_t> PuppiMET_phiUnclusteredDown = {fReader_common, "PuppiMET_phiUnclusteredDown"};
  TTreeReaderValue<Float_t> PuppiMET_phiUnclusteredUp   = {fReader_common, "PuppiMET_phiUnclusteredUp"};
  TTreeReaderValue<Float_t> PuppiMET_pt          = {fReader_common, "PuppiMET_pt"};
  TTreeReaderValue<Float_t> PuppiMET_ptJERDown   = {fReader_common, "PuppiMET_ptJERDown"};
  TTreeReaderValue<Float_t> PuppiMET_ptJERUp     = {fReader_common, "PuppiMET_ptJERUp"};
  TTreeReaderValue<Float_t> PuppiMET_ptJESDown   = {fReader_common, "PuppiMET_ptJESDown"};
  TTreeReaderValue<Float_t> PuppiMET_ptJESUp     = {fReader_common, "PuppiMET_ptJESUp"};
  TTreeReaderValue<Float_t> PuppiMET_ptUnclusteredDown = {fReader_common, "PuppiMET_ptUnclusteredDown"};
  TTreeReaderValue<Float_t> PuppiMET_ptUnclusteredUp   = {fReader_common, "PuppiMET_ptUnclusteredUp"};
  TTreeReaderValue<Float_t> PuppiMET_sumEt       = {fReader_common, "PuppiMET_sumEt"};

  //PVs
  TTreeReaderValue<iterator> nOtherPV           = {fReader_common, "nOtherPV"};
  TTreeReaderValue<int_or_char> PV_npvs                   = {fReader_common, "PV_npvs"};
  TTreeReaderValue<int_or_char> PV_npvsGood               = {fReader_common, "PV_npvsGood"};
  TTreeReaderArray<Float_t> OtherPV_z    = {fReader_common, "OtherPV_z"};
  TTreeReaderArray<Float_t> PSWeight     = {fReader_common, "PSWeight"};
  TTreeReaderValue<Float_t> PV_chi2      = {fReader_common, "PV_chi2"};
  TTreeReaderValue<Float_t> PV_ndof      = {fReader_common, "PV_ndof"};
  TTreeReaderValue<Float_t> PV_score     = {fReader_common, "PV_score"};
  TTreeReaderValue<Float_t> PV_x         = {fReader_common, "PV_x"};
  TTreeReaderValue<Float_t> PV_y         = {fReader_common, "PV_y"};
  TTreeReaderValue<Float_t> PV_z         = {fReader_common, "PV_z"};
  //Extra Run3 branches
  TTreeReaderArray<Float_t> OtherPV_score= {fReader_Run3, "OtherPV_score"};

  //TrigObj
  TTreeReaderValue<iterator> nTrigObj             = {fReader_common, "nTrigObj"};
  TTreeReaderArray<int_or_ushort> TrigObj_id      = {fReader_common, "TrigObj_id"};
  TTreeReaderArray<int_or_short> TrigObj_l1charge = {fReader_common, "TrigObj_l1charge"};
  TTreeReaderArray<Float_t> TrigObj_eta      = {fReader_common, "TrigObj_eta"};
  TTreeReaderArray<Int_t> TrigObj_filterBits = {fReader_common, "TrigObj_filterBits"};
  TTreeReaderArray<Int_t> TrigObj_l1iso      = {fReader_common, "TrigObj_l1iso"};
  TTreeReaderArray<Float_t> TrigObj_l1pt     = {fReader_common, "TrigObj_l1pt"};
  TTreeReaderArray<Float_t> TrigObj_l1pt_2   = {fReader_common, "TrigObj_l1pt_2"};
  TTreeReaderArray<Float_t> TrigObj_l2pt     = {fReader_common, "TrigObj_l2pt"};
  TTreeReaderArray<Float_t> TrigObj_phi      = {fReader_common, "TrigObj_phi"};
  TTreeReaderArray<Float_t> TrigObj_pt       = {fReader_common, "TrigObj_pt"};


  //special-branches[Not available in PYTHIA-only samples]




  

  //FSR Photons
  //TTreeReaderValue<iterator> nFsrPhoton         = {fReader_common, "nFsrPhoton"};
  //TTreeReaderArray<int_or_short> FsrPhoton_muonIdx        = {fReader_common, "FsrPhoton_muonIdx"};
  //TTreeReaderArray<Float_t> FsrPhoton_dROverEt2 = {fReader_common, "FsrPhoton_dROverEt2"};
  //TTreeReaderArray<Float_t> FsrPhoton_eta = {fReader_common, "FsrPhoton_eta"};
  //TTreeReaderArray<Float_t> FsrPhoton_phi = {fReader_common, "FsrPhoton_phi"};
  //TTreeReaderArray<Float_t> FsrPhoton_pt = {fReader_common, "FsrPhoton_pt"};
  //TTreeReaderArray<Float_t> FsrPhoton_relIso03 = {fReader_common, "FsrPhoton_relIso03"};
  //Extra Run3 branches
  //TTreeReaderArray<Short_t> FsrPhoton_electronIdx = {fReader_Run3, "FsrPhoton_electronIdx"};

  
  //HTXS_Higgs
  //TTreeReaderValue<Float_t> HTXS_Higgs_pt = {fReader_common, "HTXS_Higgs_pt"};
  //TTreeReaderValue<Float_t> HTXS_Higgs_y = {fReader_common, "HTXS_Higgs_y"};
  //TTreeReaderValue<UChar_t> HTXS_njets25 = {fReader_common, "HTXS_njets25"};
  //TTreeReaderValue<UChar_t> HTXS_njets30 = {fReader_common, "HTXS_njets30"};
  //TTreeReaderValue<Int_t> HTXS_stage1_1_cat_pTjet25GeV = {fReader_common, "HTXS_stage1_1_cat_pTjet25GeV"};
  //TTreeReaderValue<Int_t> HTXS_stage1_1_cat_pTjet30GeV = {fReader_common, "HTXS_stage1_1_cat_pTjet30GeV"};
  //TTreeReaderValue<Int_t> HTXS_stage1_1_fine_cat_pTjet25GeV = {fReader_common, "HTXS_stage1_1_fine_cat_pTjet25GeV"};
  // TTreeReaderValue<Int_t> HTXS_stage1_1_fine_cat_pTjet30GeV = {fReader_common, "HTXS_stage1_1_fine_cat_pTjet30GeV"};
  // TTreeReaderValue<Int_t> HTXS_stage1_2_cat_pTjet25GeV = {fReader_common, "HTXS_stage1_2_cat_pTjet25GeV"};
  //TTreeReaderValue<Int_t> HTXS_stage1_2_cat_pTjet30GeV = {fReader_common, "HTXS_stage1_2_cat_pTjet30GeV"};
  //TTreeReaderValue<Int_t> HTXS_stage1_2_fine_cat_pTjet25GeV = {fReader_common, "HTXS_stage1_2_fine_cat_pTjet25GeV"};
  //TTreeReaderValue<Int_t> HTXS_stage1_2_fine_cat_pTjet30GeV = {fReader_common, "HTXS_stage1_2_fine_cat_pTjet30GeV"};
  // TTreeReaderValue<Int_t> HTXS_stage_0 = {fReader_common, "HTXS_stage_0"};
  //TTreeReaderValue<Int_t> HTXS_stage_1_pTjet25 = {fReader_common, "HTXS_stage_1_pTjet25"};
  //TTreeReaderValue<Int_t> HTXS_stage_1_pTjet30 = {fReader_common, "HTXS_stage_1_pTjet30"};

  //IsoTrack
  //TTreeReaderValue<iterator> nIsoTrack          = {fReader_common, "nIsoTrack"};
  //TTreeReaderArray<int_or_short> IsoTrack_charge          = {fReader_common, "IsoTrack_charge"};
  //TTreeReaderArray<int_or_short> IsoTrack_fromPV          = {fReader_common, "IsoTrack_fromPV"};
  //TTreeReaderArray<Float_t> IsoTrack_dxy              = {fReader_common, "IsoTrack_dxy"};
  //TTreeReaderArray<Float_t> IsoTrack_dz               = {fReader_common, "IsoTrack_dz"};
  //TTreeReaderArray<Float_t> IsoTrack_eta              = {fReader_common, "IsoTrack_eta"};
  //TTreeReaderArray<Bool_t> IsoTrack_isFromLostTrack   = {fReader_common, "IsoTrack_isFromLostTrack"};
  //TTreeReaderArray<Bool_t> IsoTrack_isHighPurityTrack = {fReader_common, "IsoTrack_isHighPurityTrack"};
  //TTreeReaderArray<Bool_t> IsoTrack_isPFcand          = {fReader_common, "IsoTrack_isPFcand"};
  //TTreeReaderArray<Float_t> IsoTrack_miniPFRelIso_all = {fReader_common, "IsoTrack_miniPFRelIso_all"};
  //TTreeReaderArray<Float_t> IsoTrack_miniPFRelIso_chg = {fReader_common, "IsoTrack_miniPFRelIso_chg"};
  //TTreeReaderArray<Int_t> IsoTrack_pdgId              = {fReader_common, "IsoTrack_pdgId"};
  //TTreeReaderArray<Float_t> IsoTrack_pfRelIso03_all   = {fReader_common, "IsoTrack_pfRelIso03_all"};
  //TTreeReaderArray<Float_t> IsoTrack_pfRelIso03_chg   = {fReader_common, "IsoTrack_pfRelIso03_chg"};
  //TTreeReaderArray<Float_t> IsoTrack_phi              = {fReader_common, "IsoTrack_phi"};
  //TTreeReaderArray<Float_t> IsoTrack_pt               = {fReader_common, "IsoTrack_pt"};
  
  //LowPtElectron
  //TTreeReaderValue<iterator> nLowPtElectron     = {fReader_common, "nLowPtElectron"};
  //TTreeReaderArray<int_or_char> LowPtElectron_convWP      = {fReader_common, "LowPtElectron_convWP"};
  //TTreeReaderArray<int_or_short> LowPtElectron_genPartIdx = {fReader_commonMC, "LowPtElectron_genPartIdx"};
  //TTreeReaderArray<Float_t> LowPtElectron_ID            = {fReader_common, "LowPtElectron_ID"};
  //TTreeReaderArray<Int_t> LowPtElectron_charge          = {fReader_common, "LowPtElectron_charge"};
  //TTreeReaderArray<Bool_t> LowPtElectron_convVeto       = {fReader_common, "LowPtElectron_convVeto"};
  //TTreeReaderArray<Float_t> LowPtElectron_convVtxRadius = {fReader_common, "LowPtElectron_convVtxRadius"};
  //TTreeReaderArray<Float_t> LowPtElectron_deltaEtaSC    = {fReader_common, "LowPtElectron_deltaEtaSC"};
  //TTreeReaderArray<Float_t> LowPtElectron_dxy           = {fReader_common, "LowPtElectron_dxy"};
  //TTreeReaderArray<Float_t> LowPtElectron_dxyErr        = {fReader_common, "LowPtElectron_dxyErr"};
  //TTreeReaderArray<Float_t> LowPtElectron_dz            = {fReader_common, "LowPtElectron_dz"};
  //TTreeReaderArray<Float_t> LowPtElectron_dzErr         = {fReader_common, "LowPtElectron_dzErr"};
  //TTreeReaderArray<Float_t> LowPtElectron_eInvMinusPInv = {fReader_common, "LowPtElectron_eInvMinusPInv"};
  //TTreeReaderArray<Float_t> LowPtElectron_embeddedID    = {fReader_Run2, "LowPtElectron_embeddedID"};
  //TTreeReaderArray<Float_t> LowPtElectron_energyErr     = {fReader_common, "LowPtElectron_energyErr"};
  //TTreeReaderArray<Float_t> LowPtElectron_eta           = {fReader_common, "LowPtElectron_eta"};
  //TTreeReaderArray<UChar_t> LowPtElectron_genPartFlav   = {fReader_commonMC, "LowPtElectron_genPartFlav"};
  //TTreeReaderArray<Float_t> LowPtElectron_hoe           = {fReader_common, "LowPtElectron_hoe"};
  //TTreeReaderArray<UChar_t> LowPtElectron_lostHits      = {fReader_common, "LowPtElectron_lostHits"};
  //TTreeReaderArray<Float_t> LowPtElectron_mass             = {fReader_common, "LowPtElectron_mass"};
  //TTreeReaderArray<Float_t> LowPtElectron_miniPFRelIso_all = {fReader_common, "LowPtElectron_miniPFRelIso_all"};
  //TTreeReaderArray<Float_t> LowPtElectron_miniPFRelIso_chg = {fReader_common, "LowPtElectron_miniPFRelIso_chg"};
  //TTreeReaderArray<Int_t> LowPtElectron_pdgId           = {fReader_common, "LowPtElectron_pdgId"};
  //TTreeReaderArray<Float_t> LowPtElectron_phi           = {fReader_common, "LowPtElectron_phi"};
  //TTreeReaderArray<Float_t> LowPtElectron_pt            = {fReader_common, "LowPtElectron_pt"};
  //TTreeReaderArray<Float_t> LowPtElectron_ptbiased      = {fReader_common, "LowPtElectron_ptbiased"};
  //TTreeReaderArray<Float_t> LowPtElectron_r9            = {fReader_common, "LowPtElectron_r9"};
  //TTreeReaderArray<Float_t> LowPtElectron_scEtOverPt    = {fReader_common, "LowPtElectron_scEtOverPt"};
  //TTreeReaderArray<Float_t> LowPtElectron_sieie         = {fReader_common, "LowPtElectron_sieie"};
  //TTreeReaderArray<Float_t> LowPtElectron_unbiased      = {fReader_common, "LowPtElectron_unbiased"};
  //Extra Run3 branches
  //TTreeReaderArray<Short_t> LowPtElectron_electronIdx   = {fReader_Run3, "LowPtElectron_electronIdx"};
  //TTreeReaderArray<Short_t> LowPtElectron_photonIdx     = {fReader_Run3, "LowPtElectron_photonIdx"};

  //Photon
  //TTreeReaderValue<iterator> nPhoton                      = {fReader_common, "nPhoton"};
  //TTreeReaderArray<int_or_char> Photon_cutBased           = {fReader_common, "Photon_cutBased"};
  //TTreeReaderArray<int_or_short> Photon_electronIdx       = {fReader_common, "Photon_electronIdx"};
  //TTreeReaderArray<int_or_short> Photon_genPartIdx        = {fReader_commonMC, "Photon_genPartIdx"};
  //TTreeReaderArray<int_or_short> Photon_jetIdx            = {fReader_common, "Photon_jetIdx"};
  //TTreeReaderArray<Int_t> Photon_charge                   = {fReader_Run2, "Photon_charge"};
  //TTreeReaderArray<UChar_t> Photon_cleanmask              = {fReader_Run2, "Photon_cleanmask"};
  //TTreeReaderArray<Int_t> Photon_cutBased_Fall17V1Bitmap  = {fReader_Run2, "Photon_cutBased_Fall17V1Bitmap"};
  //TTreeReaderArray<Float_t> Photon_dEscaleDown            = {fReader_Run2, "Photon_dEscaleDown"};
  //TTreeReaderArray<Float_t> Photon_dEscaleUp              = {fReader_Run2, "Photon_dEscaleUp"};
  //TTreeReaderArray<Float_t> Photon_dEsigmaDown            = {fReader_Run2, "Photon_dEsigmaDown"};
  //TTreeReaderArray<Float_t> Photon_dEsigmaUp              = {fReader_Run2, "Photon_dEsigmaUp"};
  //TTreeReaderArray<Float_t> Photon_eCorr                  = {fReader_Run2, "Photon_eCorr"};
  //TTreeReaderArray<Bool_t> Photon_electronVeto            = {fReader_common, "Photon_electronVeto"};
  //TTreeReaderArray<Float_t> Photon_energyErr              = {fReader_common, "Photon_energyErr"};
  //TTreeReaderArray<Float_t> Photon_eta                    = {fReader_common, "Photon_eta"};
  //TTreeReaderArray<UChar_t> Photon_genPartFlav            = {fReader_commonMC, "Photon_genPartFlav"};
  //TTreeReaderArray<Float_t> Photon_hoe                    = {fReader_common, "Photon_hoe"};
  //TTreeReaderArray<Bool_t> Photon_isScEtaEB               = {fReader_common, "Photon_isScEtaEB"};
  //TTreeReaderArray<Bool_t> Photon_isScEtaEE               = {fReader_common, "Photon_isScEtaEE"};
  //TTreeReaderArray<Float_t> Photon_mass                   = {fReader_Run2, "Photon_mass"};
  //TTreeReaderArray<Float_t> Photon_mvaID                  = {fReader_common, "Photon_mvaID"};
  //TTreeReaderArray<Float_t> Photon_mvaID_Fall17V1p1       = {fReader_Run2, "Photon_mvaID_Fall17V1p1"};
  //TTreeReaderArray<Bool_t> Photon_mvaID_WP80              = {fReader_common, "Photon_mvaID_WP80"};
  //TTreeReaderArray<Bool_t> Photon_mvaID_WP90              = {fReader_common, "Photon_mvaID_WP90"};
  //TTreeReaderArray<Int_t> Photon_pdgId                    = {fReader_Run2, "Photon_pdgId"};
  //TTreeReaderArray<Float_t> Photon_pfRelIso03_all         = {fReader_Run2, "Photon_pfRelIso03_all"};
  //TTreeReaderArray<Float_t> Photon_pfRelIso03_chg         = {fReader_Run2, "Photon_pfRelIso03_chg"};
  //TTreeReaderArray<Float_t> Photon_phi                    = {fReader_common, "Photon_phi"};
  //TTreeReaderArray<Bool_t> Photon_pixelSeed               = {fReader_common, "Photon_pixelSeed"};
  //TTreeReaderArray<Float_t> Photon_pt                     = {fReader_common, "Photon_pt"};
  //TTreeReaderArray<Float_t> Photon_r9                     = {fReader_common, "Photon_r9"};
  //TTreeReaderArray<UChar_t> Photon_seedGain               = {fReader_common, "Photon_seedGain"};
  //TTreeReaderArray<Float_t> Photon_sieie                  = {fReader_common, "Photon_sieie"};
  //TTreeReaderArray<Int_t> Photon_vidNestedWPBitmap        = {fReader_common, "Photon_vidNestedWPBitmap"};
  //Extra Run3 branches
  //TTreeReaderArray<Float_t> Photon_ecalPFClusterIso       = {fReader_Run3, "Photon_ecalPFClusterIso"};
  //TTreeReaderArray<Float_t> Photon_energyRaw              = {fReader_Run3, "Photon_energyRaw"};
  //TTreeReaderArray<Float_t> Photon_esEffSigmaRR           = {fReader_Run3, "Photon_esEffSigmaRR"};
  //TTreeReaderArray<Float_t> Photon_esEnergyOverRawE       = {fReader_Run3, "Photon_esEnergyOverRawE"};
  //TTreeReaderArray<Float_t> Photon_etaWidth               = {fReader_Run3, "Photon_etaWidth"};
  //TTreeReaderArray<Float_t> Photon_haloTaggerMVAVal       = {fReader_Run3, "Photon_haloTaggerMVAVal"};
  //TTreeReaderArray<Bool_t> Photon_hasConversionTracks     = {fReader_Run3, "Photon_hasConversionTracks"};
  //TTreeReaderArray<Float_t> Photon_hcalPFClusterIso       = {fReader_Run3, "Photon_hcalPFClusterIso"};
  //TTreeReaderArray<Float_t> Photon_hoe_PUcorr             = {fReader_Run3, "Photon_hoe_PUcorr"};
  //TTreeReaderArray<Float_t> Photon_pfChargedIso           = {fReader_Run3, "Photon_pfChargedIso"};
  //TTreeReaderArray<Float_t> Photon_pfChargedIsoPFPV       = {fReader_Run3, "Photon_pfChargedIsoPFPV"};
  //TTreeReaderArray<Float_t> Photon_pfChargedIsoWorstVtx   = {fReader_Run3, "Photon_pfChargedIsoWorstVtx"};
  //TTreeReaderArray<Float_t> Photon_pfPhoIso03             = {fReader_Run3, "Photon_pfPhoIso03"};
  //TTreeReaderArray<Float_t> Photon_pfRelIso03_all_quadratic = {fReader_Run3, "Photon_pfRelIso03_all_quadratic"};
  //TTreeReaderArray<Float_t> Photon_pfRelIso03_chg_quadratic = {fReader_Run3, "Photon_pfRelIso03_chg_quadratic"};
  //TTreeReaderArray<Float_t> Photon_phiWidth               = {fReader_Run3, "Photon_phiWidth"};
  //TTreeReaderArray<Float_t> Photon_s4                     = {fReader_Run3, "Photon_s4"};
  //TTreeReaderArray<Char_t> Photon_seediEtaOriX            = {fReader_Run3, "Photon_seediEtaOriX"};
  //TTreeReaderArray<Int_t> Photon_seediPhiOriY             = {fReader_Run3, "Photon_seediPhiOriY"};
  //TTreeReaderArray<Float_t> Photon_sieip                  = {fReader_Run3, "Photon_sieip"};
  //TTreeReaderArray<Float_t> Photon_sipip                  = {fReader_Run3, "Photon_sipip"};
  //TTreeReaderArray<Float_t> Photon_trkSumPtHollowConeDR03 = {fReader_Run3, "Photon_trkSumPtHollowConeDR03"};
  //TTreeReaderArray<Float_t> Photon_trkSumPtSolidConeDR04  = {fReader_Run3, "Photon_trkSumPtSolidConeDR04"};
  //TTreeReaderArray<Float_t> Photon_x_calo                 = {fReader_Run3, "Photon_x_calo"};
  //TTreeReaderArray<Float_t> Photon_y_calo                 = {fReader_Run3, "Photon_y_calo"};
  //TTreeReaderArray<Float_t> Photon_z_calo                 = {fReader_Run3, "Photon_z_calo"};

  //RawMET
  //TTreeReaderValue<Float_t> RawMET_phi        = {fReader_common, "RawMET_phi"};
  //TTreeReaderValue<Float_t> RawMET_pt         = {fReader_common, "RawMET_pt"};
  //TTreeReaderValue<Float_t> RawMET_sumEt      = {fReader_common, "RawMET_sumEt"};
  //TTreeReaderValue<Float_t> RawPuppiMET_phi   = {fReader_common, "RawPuppiMET_phi"};
  //TTreeReaderValue<Float_t> RawPuppiMET_pt    = {fReader_common, "RawPuppiMET_pt"};
  //TTreeReaderValue<Float_t> RawPuppiMET_sumEt = {fReader_common, "RawPuppiMET_sumEt"};

  //SV
  //TTreeReaderValue<iterator> nSV             = {fReader_common, "nSV"};
  //TTreeReaderArray<int_or_short> SV_charge   = {fReader_common, "SV_charge"};  
  //TTreeReaderArray<Float_t> SV_chi2    = {fReader_common, "SV_chi2"};
  //TTreeReaderArray<Float_t> SV_dlen    = {fReader_common, "SV_dlen"};
  //TTreeReaderArray<Float_t> SV_dlenSig = {fReader_common, "SV_dlenSig"};
  //TTreeReaderArray<Float_t> SV_dxy     = {fReader_common, "SV_dxy"};
  //TTreeReaderArray<Float_t> SV_dxySig  = {fReader_common, "SV_dxySig"};
  //TTreeReaderArray<Float_t> SV_eta     = {fReader_common, "SV_eta"};
  //TTreeReaderArray<Float_t> SV_mass    = {fReader_common, "SV_mass"};
  //TTreeReaderArray<Float_t> SV_ndof    = {fReader_common, "SV_ndof"};
  //TTreeReaderArray<UChar_t> SV_ntracks = {fReader_common, "SV_ntracks"};
  //TTreeReaderArray<Float_t> SV_pAngle  = {fReader_common, "SV_pAngle"};
  //TTreeReaderArray<Float_t> SV_phi     = {fReader_common, "SV_phi"};
  //TTreeReaderArray<Float_t> SV_pt      = {fReader_common, "SV_pt"};
  //TTreeReaderArray<Float_t> SV_x       = {fReader_common, "SV_x"};
  //TTreeReaderArray<Float_t> SV_y       = {fReader_common, "SV_y"};
  //TTreeReaderArray<Float_t> SV_z       = {fReader_common, "SV_z"};

  //SoftActivityJet
  //TTreeReaderValue<iterator> nSoftActivityJet   = {fReader_common, "nSoftActivityJet"};
  //TTreeReaderValue<Float_t> SoftActivityJetHT    = {fReader_common, "SoftActivityJetHT"};
  //TTreeReaderValue<Float_t> SoftActivityJetHT10  = {fReader_common, "SoftActivityJetHT10"};
  //TTreeReaderValue<Float_t> SoftActivityJetHT2   = {fReader_common, "SoftActivityJetHT2"};
  //TTreeReaderValue<Float_t> SoftActivityJetHT5   = {fReader_common, "SoftActivityJetHT5"};
  //TTreeReaderValue<Int_t> SoftActivityJetNjets10 = {fReader_common, "SoftActivityJetNjets10"};
  //TTreeReaderValue<Int_t> SoftActivityJetNjets2  = {fReader_common, "SoftActivityJetNjets2"};
  //TTreeReaderValue<Int_t> SoftActivityJetNjets5  = {fReader_common, "SoftActivityJetNjets5"};
  //TTreeReaderArray<Float_t> SoftActivityJet_eta  = {fReader_common, "SoftActivityJet_eta"};
  //TTreeReaderArray<Float_t> SoftActivityJet_phi  = {fReader_common, "SoftActivityJet_phi"};
  //TTreeReaderArray<Float_t> SoftActivityJet_pt   = {fReader_common, "SoftActivityJet_pt"};

  //SubJet
  //TTreeReaderValue<iterator> nSubJet                  = {fReader_common, "nSubJet"};
  //TTreeReaderArray<int_or_char> SubJet_hadronFlavour  = {fReader_common, "SubJet_hadronFlavour"};
  //TTreeReaderArray<Float_t> SubJet_btagCSVV2 = {fReader_Run2, "SubJet_btagCSVV2"};
  //TTreeReaderArray<Float_t> SubJet_btagDeepB = {fReader_common, "SubJet_btagDeepB"};
  //TTreeReaderArray<Float_t> SubJet_eta       = {fReader_common, "SubJet_eta"};
  //TTreeReaderArray<Float_t> SubJet_mass      = {fReader_common, "SubJet_mass"};
  //TTreeReaderArray<Float_t> SubJet_n2b1      = {fReader_common, "SubJet_n2b1"};
  //TTreeReaderArray<Float_t> SubJet_n3b1      = {fReader_common, "SubJet_n3b1"};
  //TTreeReaderArray<UChar_t> SubJet_nBHadrons = {fReader_common, "SubJet_nBHadrons"};
  //TTreeReaderArray<UChar_t> SubJet_nCHadrons = {fReader_common, "SubJet_nCHadrons"};
  //TTreeReaderArray<Float_t> SubJet_phi       = {fReader_common, "SubJet_phi"};
  //TTreeReaderArray<Float_t> SubJet_pt        = {fReader_common, "SubJet_pt"};
  //TTreeReaderArray<Float_t> SubJet_rawFactor = {fReader_common, "SubJet_rawFactor"};
  //TTreeReaderArray<Float_t> SubJet_tau1      = {fReader_common, "SubJet_tau1"};
  //TTreeReaderArray<Float_t> SubJet_tau2      = {fReader_common, "SubJet_tau2"};
  //TTreeReaderArray<Float_t> SubJet_tau3      = {fReader_common, "SubJet_tau3"};
  //TTreeReaderArray<Float_t> SubJet_tau4      = {fReader_common, "SubJet_tau4"};

  //boostedTau
  //TTreeReaderValue<iterator> nboostedTau        = {fReader_common, "nboostedTau"};
  //TTreeReaderArray<int_or_short> boostedTau_genPartIdx    = {fReader_commonMC, "boostedTau_genPartIdx"};
  //TTreeReaderArray<int_or_short> boostedTau_jetIdx        = {fReader_common, "boostedTau_jetIdx"};
  //TTreeReaderArray<int_or_short> boostedTau_rawAntiEleCat2018 = {fReader_common, "boostedTau_rawAntiEleCat2018"};
  //TTreeReaderArray<Int_t> boostedTau_charge              = {fReader_common, "boostedTau_charge"};
  //TTreeReaderArray<Float_t> boostedTau_chargedIso        = {fReader_common, "boostedTau_chargedIso"};
  //TTreeReaderArray<Int_t> boostedTau_decayMode           = {fReader_common, "boostedTau_decayMode"};
  //TTreeReaderArray<Float_t> boostedTau_eta               = {fReader_common, "boostedTau_eta"};
  //TTreeReaderArray<UChar_t> boostedTau_genPartFlav       = {fReader_commonMC, "boostedTau_genPartFlav"};
  //TTreeReaderArray<UChar_t> boostedTau_idAntiEle2018     = {fReader_common, "boostedTau_idAntiEle2018"};
  //TTreeReaderArray<UChar_t> boostedTau_idAntiMu          = {fReader_common, "boostedTau_idAntiMu"};
  //TTreeReaderArray<UChar_t> boostedTau_idMVAnewDM2017v2  = {fReader_common, "boostedTau_idMVAnewDM2017v2"};
  //TTreeReaderArray<UChar_t> boostedTau_idMVAoldDM2017v2  = {fReader_common, "boostedTau_idMVAoldDM2017v2"};
  //TTreeReaderArray<UChar_t> boostedTau_idMVAoldDMdR032017v2 = {fReader_Run2, "boostedTau_idMVAoldDMdR032017v2"};
  //TTreeReaderArray<Float_t> boostedTau_leadTkDeltaEta    = {fReader_common, "boostedTau_leadTkDeltaEta"};
  //TTreeReaderArray<Float_t> boostedTau_leadTkDeltaPhi    = {fReader_common, "boostedTau_leadTkDeltaPhi"};
  //TTreeReaderArray<Float_t> boostedTau_leadTkPtOverTauPt = {fReader_common, "boostedTau_leadTkPtOverTauPt"};
  //TTreeReaderArray<Float_t> boostedTau_mass              = {fReader_common, "boostedTau_mass"};
  //TTreeReaderArray<Float_t> boostedTau_neutralIso        = {fReader_common, "boostedTau_neutralIso"};
  //TTreeReaderArray<Float_t> boostedTau_phi               = {fReader_common, "boostedTau_phi"};
  //TTreeReaderArray<Float_t> boostedTau_photonsOutsideSignalCone = {fReader_common, "boostedTau_photonsOutsideSignalCone"};
  //TTreeReaderArray<Float_t> boostedTau_pt                = {fReader_common, "boostedTau_pt"};
  //TTreeReaderArray<Float_t> boostedTau_puCorr            = {fReader_common, "boostedTau_puCorr"};
  //TTreeReaderArray<Float_t> boostedTau_rawAntiEle2018    = {fReader_common, "boostedTau_rawAntiEle2018"};
  //TTreeReaderArray<Float_t> boostedTau_rawIso            = {fReader_common, "boostedTau_rawIso"};
  //TTreeReaderArray<Float_t> boostedTau_rawIsodR03        = {fReader_common, "boostedTau_rawIsodR03"};
  //TTreeReaderArray<Float_t> boostedTau_rawMVAnewDM2017v2 = {fReader_common, "boostedTau_rawMVAnewDM2017v2"};
  //TTreeReaderArray<Float_t> boostedTau_rawMVAoldDM2017v2 = {fReader_common, "boostedTau_rawMVAoldDM2017v2"};
  //TTreeReaderArray<Float_t> boostedTau_rawMVAoldDMdR032017v2 = {fReader_Run2, "boostedTau_rawMVAoldDMdR032017v2"};
  
  //Beamspot
  //TTreeReaderValue<Float_t> BeamSpot_sigmaZ      = {fReader_Run3, "BeamSpot_sigmaZ"};
  //TTreeReaderValue<Float_t> BeamSpot_sigmaZError = {fReader_Run3, "BeamSpot_sigmaZError"};
  //TTreeReaderValue<Char_t> BeamSpot_type         = {fReader_Run3, "BeamSpot_type"};
  //TTreeReaderValue<Float_t> BeamSpot_z           = {fReader_Run3, "BeamSpot_z"};
  //TTreeReaderValue<Float_t> BeamSpot_zError      = {fReader_Run3, "BeamSpot_zError"}; 

  //CaloMET
  //TTreeReaderValue<Float_t> CaloMET_phi    = {fReader_common, "CaloMET_phi"};
  //TTreeReaderValue<Float_t> CaloMET_pt     = {fReader_common, "CaloMET_pt"};
  //TTreeReaderValue<Float_t> CaloMET_sumEt  = {fReader_common, "CaloMET_sumEt"};

  //TrackMET
  //TTreeReaderValue<Float_t> TkMET_phi   = {fReader_common, "TkMET_phi"};
  //TTreeReaderValue<Float_t> TkMET_pt    = {fReader_common, "TkMET_pt"};
  //TTreeReaderValue<Float_t> TkMET_sumEt = {fReader_common, "TkMET_sumEt"};

  //ChsMET
  //TTreeReaderValue<Float_t> ChsMET_phi    = {fReader_common, "ChsMET_phi"};
  //TTreeReaderValue<Float_t> ChsMET_pt     = {fReader_common, "ChsMET_pt"};
  //TTreeReaderValue<Float_t> ChsMET_sumEt  = {fReader_common, "ChsMET_sumEt"};

  //CorrT1METJet
  //TTreeReaderValue<iterator> nCorrT1METJet      = {fReader_common, "nCorrT1METJet"};
  //TTreeReaderArray<Float_t> CorrT1METJet_area = {fReader_common, "CorrT1METJet_area"};
  //TTreeReaderArray<Float_t> CorrT1METJet_eta  = {fReader_common, "CorrT1METJet_eta"};
  //TTreeReaderArray<Float_t> CorrT1METJet_muonSubtrFactor = {fReader_common, "CorrT1METJet_muonSubtrFactor"};
  //TTreeReaderArray<Float_t> CorrT1METJet_phi  = {fReader_common, "CorrT1METJet_phi"};
  //TTreeReaderArray<Float_t> CorrT1METJet_rawPt= {fReader_common, "CorrT1METJet_rawPt"};

  //DeepMETResolution
  //TTreeReaderValue<Float_t> DeepMETResolutionTune_phi = {fReader_common, "DeepMETResolutionTune_phi"};
  //TTreeReaderValue<Float_t> DeepMETResolutionTune_pt  = {fReader_common, "DeepMETResolutionTune_pt"};
  //TTreeReaderValue<Float_t> DeepMETResponseTune_phi   = {fReader_common, "DeepMETResponseTune_phi"};
  //TTreeReaderValue<Float_t> DeepMETResponseTune_pt    = {fReader_common, "DeepMETResponseTune_pt"};
  

  
  //-------------------------------//
  //    Extra Run3 branches        //
  //-------------------------------//
  
  //TTreeReaderValue<UInt_t> bunchCrossing = {fReader_Run3, "bunchCrossing"};

  //-------------------------------//
  //     Trigger branches          //
  //-------------------------------//

  //used in analysis
  TTreeReaderValue<Bool_t> HLT_IsoMu24           = {fReader_common, "HLT_IsoMu24"}; //used in analysis
  TTreeReaderValue<Bool_t> HLT_IsoMu27           = {fReader_common, "HLT_IsoMu27"}; //used in analysis
  TTreeReaderValue<Bool_t> HLT_IsoTkMu24         = {fReader_2016  , "HLT_IsoTkMu24"};//Only 2016, used in analysis  
  TTreeReaderValue<Bool_t> HLT_Ele30_WPTight_Gsf = {fReader_common, "HLT_Ele30_WPTight_Gsf"};//used in analysis
  TTreeReaderValue<Bool_t> HLT_Ele27_WPTight_Gsf = {fReader_common, "HLT_Ele27_WPTight_Gsf"};
  TTreeReaderValue<Bool_t> HLT_Ele32_WPTight_Gsf = {fReader_2018_2022, "HLT_Ele32_WPTight_Gsf"};//common for 2018/2022
  
  //Common trigger branch for year 2017 and 2018
  TTreeReaderValue<Bool_t> HLT_Ele32_WPTight_Gsf_L1DoubleEG = {fReader_2017_2018, "HLT_Ele32_WPTight_Gsf_L1DoubleEG"};

  //Others
  //HLT-Single-Muon
  TTreeReaderValue<Bool_t> HLT_IsoMu20           = {fReader_common, "HLT_IsoMu20"};
  TTreeReaderValue<Bool_t> HLT_IsoMu30           = {fReader_common, "HLT_IsoMu30"};
  
  //non-isolated muon trigger
  TTreeReaderValue<Bool_t> HLT_Mu12              = {fReader_common, "HLT_Mu12"};
  TTreeReaderValue<Bool_t> HLT_Mu15              = {fReader_common, "HLT_Mu15"};
  TTreeReaderValue<Bool_t> HLT_Mu17              = {fReader_common, "HLT_Mu17"};
  TTreeReaderValue<Bool_t> HLT_Mu18_Mu9          = {fReader_common, "HLT_Mu18_Mu9"};
  TTreeReaderValue<Bool_t> HLT_Mu19              = {fReader_common, "HLT_Mu19"};
  TTreeReaderValue<Bool_t> HLT_Mu20              = {fReader_common, "HLT_Mu20"};
  TTreeReaderValue<Bool_t> HLT_Mu23_Mu12         = {fReader_common, "HLT_Mu23_Mu12"};
  TTreeReaderValue<Bool_t> HLT_Mu27              = {fReader_common, "HLT_Mu27"};
  TTreeReaderValue<Bool_t> HLT_Mu37_TkMu27       = {fReader_common, "HLT_Mu37_TkMu27"};
  TTreeReaderValue<Bool_t> HLT_Mu55              = {fReader_common, "HLT_Mu55"};
  TTreeReaderValue<Bool_t> HLT_Mu8               = {fReader_common, "HLT_Mu8"};
  TTreeReaderValue<Bool_t> HLT_Mu50              = {fReader_common, "HLT_Mu50"};
  
  //HLT-Single-Electron
  TTreeReaderValue<Bool_t> HLT_Ele20_WPTight_Gsf = {fReader_common, "HLT_Ele20_WPTight_Gsf"};
  TTreeReaderValue<Bool_t> HLT_Ele28_WPTight_Gsf = {fReader_common, "HLT_Ele28_WPTight_Gsf"};
  TTreeReaderValue<Bool_t> HLT_Ele35_WPTight_Gsf = {fReader_common, "HLT_Ele35_WPTight_Gsf"};
  TTreeReaderValue<Bool_t> HLT_Ele38_WPTight_Gsf = {fReader_common, "HLT_Ele38_WPTight_Gsf"};
  TTreeReaderValue<Bool_t> HLT_Ele40_WPTight_Gsf = {fReader_common, "HLT_Ele40_WPTight_Gsf"};
  
  //LooseWP
  TTreeReaderValue<Bool_t> HLT_Ele15_WPLoose_Gsf = {fReader_common, "HLT_Ele15_WPLoose_Gsf"};
  TTreeReaderValue<Bool_t> HLT_Ele17_WPLoose_Gsf = {fReader_common, "HLT_Ele17_WPLoose_Gsf"};
  TTreeReaderValue<Bool_t> HLT_Ele20_WPLoose_Gsf = {fReader_common, "HLT_Ele20_WPLoose_Gsf"};

  //Double-electron-trigger
  TTreeReaderValue<Bool_t> HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL = {fReader_common, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"};
  //cross-trigger
  TTreeReaderValue<Bool_t> HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL = {fReader_common, "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"};
  TTreeReaderValue<Bool_t> HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ = {fReader_common, "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"};

 
  //#########################################################################
  //Branches which are only present in the MC files are read using fReader_MC
  //#########################################################################

  //GenPart
  TTreeReaderValue<iterator> nGenPart           = {fReader_commonMC, "nGenPart"};
  TTreeReaderArray<int_or_short> GenPart_genPartIdxMother = {fReader_commonMC, "GenPart_genPartIdxMother"};
  TTreeReaderArray<int_or_ushort> GenPart_statusFlags     = {fReader_commonMC, "GenPart_statusFlags"};
  TTreeReaderArray<Float_t> GenPart_eta  = {fReader_commonMC, "GenPart_eta"};
  TTreeReaderArray<Float_t> GenPart_mass = {fReader_commonMC, "GenPart_mass"};
  TTreeReaderArray<Int_t> GenPart_pdgId  = {fReader_commonMC, "GenPart_pdgId"};
  TTreeReaderArray<Float_t> GenPart_phi  = {fReader_commonMC, "GenPart_phi"};
  TTreeReaderArray<Float_t> GenPart_pt   = {fReader_commonMC, "GenPart_pt"};
  TTreeReaderArray<Int_t> GenPart_status = {fReader_commonMC, "GenPart_status"};

  

  //GenDressedLepton
  //TTreeReaderValue<iterator> nGenDressedLepton  = {fReader_commonMC, "nGenDressedLepton"};
  //TTreeReaderArray<Float_t> GenDressedLepton_eta      = {fReader_commonMC, "GenDressedLepton_eta"};
  //TTreeReaderArray<Bool_t> GenDressedLepton_hasTauAnc = {fReader_commonMC, "GenDressedLepton_hasTauAnc"};
  //TTreeReaderArray<Float_t> GenDressedLepton_mass     = {fReader_commonMC, "GenDressedLepton_mass"};
  //TTreeReaderArray<Int_t> GenDressedLepton_pdgId      = {fReader_commonMC, "GenDressedLepton_pdgId"};
  //TTreeReaderArray<Float_t> GenDressedLepton_phi      = {fReader_commonMC, "GenDressedLepton_phi"};
  //TTreeReaderArray<Float_t> GenDressedLepton_pt       = {fReader_commonMC, "GenDressedLepton_pt"};

  //GenIsolatedPhoton
  //TTreeReaderValue<iterator> nGenIsolatedPhoton = {fReader_commonMC, "nGenIsolatedPhoton"};
  //TTreeReaderArray<Float_t> GenIsolatedPhoton_eta    = {fReader_commonMC, "GenIsolatedPhoton_eta"};
  //TTreeReaderArray<Float_t> GenIsolatedPhoton_mass   = {fReader_commonMC, "GenIsolatedPhoton_mass"};
  //TTreeReaderArray<Float_t> GenIsolatedPhoton_phi    = {fReader_commonMC, "GenIsolatedPhoton_phi"};
  //TTreeReaderArray<Float_t> GenIsolatedPhoton_pt     = {fReader_commonMC, "GenIsolatedPhoton_pt"};

  //GenJetAK8
  TTreeReaderValue<iterator> nGenJetAK8         = {fReader_commonMC, "nGenJetAK8"};
  TTreeReaderArray<int_or_short> GenJetAK8_partonFlavour  = {fReader_commonMC, "GenJetAK8_partonFlavour"};
  TTreeReaderArray<Float_t> GenJetAK8_eta            = {fReader_commonMC, "GenJetAK8_eta"};
  TTreeReaderArray<UChar_t> GenJetAK8_hadronFlavour  = {fReader_commonMC, "GenJetAK8_hadronFlavour"};
  TTreeReaderArray<Float_t> GenJetAK8_mass           = {fReader_commonMC, "GenJetAK8_mass"};
  TTreeReaderArray<Float_t> GenJetAK8_phi            = {fReader_commonMC, "GenJetAK8_phi"};
  TTreeReaderArray<Float_t> GenJetAK8_pt             = {fReader_commonMC, "GenJetAK8_pt"};

  //GenJet
  TTreeReaderValue<iterator> nGenJet            = {fReader_commonMC, "nGenJet"};
  TTreeReaderArray<int_or_short> GenJet_partonFlavour     = {fReader_commonMC, "GenJet_partonFlavour"};
  TTreeReaderArray<Float_t> GenJet_eta            = {fReader_commonMC, "GenJet_eta"};
  TTreeReaderArray<UChar_t> GenJet_hadronFlavour  = {fReader_commonMC, "GenJet_hadronFlavour"};
  TTreeReaderArray<Float_t> GenJet_mass           = {fReader_commonMC, "GenJet_mass"};
  TTreeReaderArray<Float_t> GenJet_phi            = {fReader_commonMC, "GenJet_phi"};
  TTreeReaderArray<Float_t> GenJet_pt             = {fReader_commonMC, "GenJet_pt"};

  //SubGenJetAK8
  //TTreeReaderArray<Float_t> SubGenJetAK8_eta      = {fReader_commonMC, "SubGenJetAK8_eta"};
  //TTreeReaderArray<Float_t> SubGenJetAK8_mass     = {fReader_commonMC, "SubGenJetAK8_mass"};
  //TTreeReaderArray<Float_t> SubGenJetAK8_phi      = {fReader_commonMC, "SubGenJetAK8_phi"};
  //TTreeReaderArray<Float_t> SubGenJetAK8_pt       = {fReader_commonMC, "SubGenJetAK8_pt"};
  
  //GenMET
  //TTreeReaderValue<Float_t> GenMET_phi = {fReader_commonMC, "GenMET_phi"};
  //TTreeReaderValue<Float_t> GenMET_pt  = {fReader_commonMC, "GenMET_pt"};


  //GenVisTau
  //TTreeReaderValue<iterator> nGenVisTau         = {fReader_commonMC, "nGenVisTau"};
  //TTreeReaderArray<int_or_short> GenVisTau_charge         = {fReader_commonMC, "GenVisTau_charge"};
  //TTreeReaderArray<int_or_short> GenVisTau_genPartIdxMother = {fReader_commonMC, "GenVisTau_genPartIdxMother"};
  //TTreeReaderArray<int_or_char> GenVisTau_status          = {fReader_commonMC, "GenVisTau_status"};
  //TTreeReaderArray<Float_t> GenVisTau_eta = {fReader_commonMC, "GenVisTau_eta"};
  //TTreeReaderArray<Float_t> GenVisTau_mass = {fReader_commonMC, "GenVisTau_mass"};
  //TTreeReaderArray<Float_t> GenVisTau_phi = {fReader_commonMC, "GenVisTau_phi"};
  //TTreeReaderArray<Float_t> GenVisTau_pt = {fReader_commonMC, "GenVisTau_pt"};
  
  //GenVertex
  //TTreeReaderValue<Float_t> GenVtx_t0 = {fReader_commonMC, "GenVtx_t0"};
  //TTreeReaderValue<Float_t> GenVtx_x = {fReader_commonMC, "GenVtx_x"};
  //TTreeReaderValue<Float_t> GenVtx_y = {fReader_commonMC, "GenVtx_y"};
  //TTreeReaderValue<Float_t> GenVtx_z = {fReader_commonMC, "GenVtx_z"};

  //Generator
  TTreeReaderValue<Float_t> Generator_binvar   = {fReader_commonMC, "Generator_binvar"};
  TTreeReaderValue<Int_t> Generator_id1        = {fReader_commonMC, "Generator_id1"};
  TTreeReaderValue<Int_t> Generator_id2        = {fReader_commonMC, "Generator_id2"};
  TTreeReaderValue<Float_t> Generator_scalePDF = {fReader_commonMC, "Generator_scalePDF"};
  TTreeReaderValue<Float_t> Generator_weight   = {fReader_commonMC, "Generator_weight"};
  TTreeReaderValue<Float_t> Generator_x1       = {fReader_commonMC, "Generator_x1"};
  TTreeReaderValue<Float_t> Generator_x2       = {fReader_commonMC, "Generator_x2"};
  TTreeReaderValue<Float_t> Generator_xpdf1    = {fReader_commonMC, "Generator_xpdf1"};
  TTreeReaderValue<Float_t> Generator_xpdf2    = {fReader_commonMC, "Generator_xpdf2"};

  //L1PreFiringWeight
  TTreeReaderValue<Float_t> L1PreFiringWeight_Dn          = {fReader_Run2MC, "L1PreFiringWeight_Dn"};
  TTreeReaderValue<Float_t> L1PreFiringWeight_ECAL_Dn     = {fReader_Run2MC, "L1PreFiringWeight_ECAL_Dn"};
  TTreeReaderValue<Float_t> L1PreFiringWeight_ECAL_Nom    = {fReader_Run2MC, "L1PreFiringWeight_ECAL_Nom"};
  TTreeReaderValue<Float_t> L1PreFiringWeight_ECAL_Up     = {fReader_Run2MC, "L1PreFiringWeight_ECAL_Up"};
  TTreeReaderValue<Float_t> L1PreFiringWeight_Muon_Nom    = {fReader_Run2MC, "L1PreFiringWeight_Muon_Nom"};
  TTreeReaderValue<Float_t> L1PreFiringWeight_Muon_StatDn = {fReader_Run2MC, "L1PreFiringWeight_Muon_StatDn"};
  TTreeReaderValue<Float_t> L1PreFiringWeight_Muon_StatUp = {fReader_Run2MC, "L1PreFiringWeight_Muon_StatUp"};
  TTreeReaderValue<Float_t> L1PreFiringWeight_Muon_SystDn = {fReader_Run2MC, "L1PreFiringWeight_Muon_SystDn"};
  TTreeReaderValue<Float_t> L1PreFiringWeight_Muon_SystUp = {fReader_Run2MC, "L1PreFiringWeight_Muon_SystUp"};
  TTreeReaderValue<Float_t> L1PreFiringWeight_Nom         = {fReader_Run2MC, "L1PreFiringWeight_Nom"};
  TTreeReaderValue<Float_t> L1PreFiringWeight_Up          = {fReader_Run2MC, "L1PreFiringWeight_Up"};

  //LHE-branches are not available any pythia-only generated samples
  //These branches are tracked by fReader: <reader_name>MCSpecial 
  //LHEPart
  TTreeReaderValue<iterator> nLHEPart              = {fReader_commonMCSpecial, "nLHEPart"};
  TTreeReaderArray<Float_t> LHEPart_eta         = {fReader_commonMCSpecial, "LHEPart_eta"};
  TTreeReaderArray<Float_t> LHEPart_incomingpz  = {fReader_commonMCSpecial, "LHEPart_incomingpz"};
  TTreeReaderArray<Float_t> LHEPart_mass        = {fReader_commonMCSpecial, "LHEPart_mass"};
  TTreeReaderArray<Int_t> LHEPart_pdgId         = {fReader_commonMCSpecial, "LHEPart_pdgId"};
  TTreeReaderArray<Float_t> LHEPart_phi         = {fReader_commonMCSpecial, "LHEPart_phi"};
  TTreeReaderArray<Float_t> LHEPart_pt          = {fReader_commonMCSpecial, "LHEPart_pt"};
  TTreeReaderArray<Int_t> LHEPart_spin          = {fReader_commonMCSpecial, "LHEPart_spin"};
  TTreeReaderArray<Int_t> LHEPart_status        = {fReader_commonMCSpecial, "LHEPart_status"};

  //LHEPdfWeight
  TTreeReaderValue<iterator> nLHEPdfWeight         = {fReader_commonMCSpecial, "nLHEPdfWeight"};
  TTreeReaderArray<Float_t> LHEPdfWeight        = {fReader_commonMCSpecial, "LHEPdfWeight"};

  //LHEReweightingWeight
  //TTreeReaderValue<iterator> nLHEReweightingWeight = {fReader_commonMCSpecial, "nLHEReweightingWeight"};
  //TTreeReaderArray<Float_t> LHEReweightingWeight = {fReader_commonMCSpecial, "LHEReweightingWeight"};
  
  //LHEScaleWeight
  TTreeReaderValue<iterator> nLHEScaleWeight       = {fReader_commonMCSpecial, "nLHEScaleWeight"};
  TTreeReaderArray<Float_t> LHEScaleWeight       = {fReader_commonMCSpecial, "LHEScaleWeight"};
  TTreeReaderValue<Float_t> LHE_AlphaS           = {fReader_commonMCSpecial, "LHE_AlphaS"};
  TTreeReaderValue<Float_t> LHE_HT               = {fReader_commonMCSpecial, "LHE_HT"};
  TTreeReaderValue<Float_t> LHE_HTIncoming       = {fReader_commonMCSpecial, "LHE_HTIncoming"};
  TTreeReaderValue<UChar_t> LHE_Nb               = {fReader_commonMCSpecial, "LHE_Nb"};
  TTreeReaderValue<UChar_t> LHE_Nc               = {fReader_commonMCSpecial, "LHE_Nc"};
  TTreeReaderValue<UChar_t> LHE_Nglu             = {fReader_commonMCSpecial, "LHE_Nglu"};
  TTreeReaderValue<UChar_t> LHE_Njets            = {fReader_commonMCSpecial, "LHE_Njets"};
  TTreeReaderValue<UChar_t> LHE_NpLO             = {fReader_commonMCSpecial, "LHE_NpLO"};
  TTreeReaderValue<UChar_t> LHE_NpNLO            = {fReader_commonMCSpecial, "LHE_NpNLO"};
  TTreeReaderValue<UChar_t> LHE_Nuds             = {fReader_commonMCSpecial, "LHE_Nuds"};
  TTreeReaderValue<Float_t> LHE_Vpt              = {fReader_commonMCSpecial, "LHE_Vpt"};
  TTreeReaderValue<Float_t> LHEWeight_originalXWGTUP = {fReader_commonMCSpecial, "LHEWeight_originalXWGTUP"};

  //PileUp
  //TTreeReaderValue<Float_t> Pileup_gpudensity    = {fReader_commonMC, "Pileup_gpudensity"};
  TTreeReaderValue<Int_t> Pileup_nPU             = {fReader_commonMC, "Pileup_nPU"};
  TTreeReaderValue<Float_t> Pileup_nTrueInt      = {fReader_commonMC, "Pileup_nTrueInt"};
  //TTreeReaderValue<Float_t> Pileup_pudensity     = {fReader_commonMC, "Pileup_pudensity"};
  //TTreeReaderValue<Int_t> Pileup_sumEOOT         = {fReader_commonMC, "Pileup_sumEOOT"};
  //TTreeReaderValue<Int_t> Pileup_sumLOOT         = {fReader_commonMC, "Pileup_sumLOOT"};

  //PS & genWeight
  //TTreeReaderValue<Float_t> genWeight            = {fReader_commonMC, "genWeight"};
  //TTreeReaderValue<iterator> nPSWeight          = {fReader_commonMC, "nPSWeight"};
  
  //FixedGridRhoFastjet [Need for JES uncertainties]
  TTreeReaderValue<Float_t> fixedGridRhoFastjetAll = {fReader_Run2MC, "fixedGridRhoFastjetAll"};
  //TTreeReaderValue<Float_t> fixedGridRhoFastjetCentral = {fReader_Run2MC, "fixedGridRhoFastjetCentral"};
  //TTreeReaderValue<Float_t> fixedGridRhoFastjetCentralCalo = {fReader_Run2MC, "fixedGridRhoFastjetCentralCalo"};
  //TTreeReaderValue<Float_t> fixedGridRhoFastjetCentralChargedPileUp = {fReader_Run2MC, "fixedGridRhoFastjetCentralChargedPileUp"};
  //TTreeReaderValue<Float_t> fixedGridRhoFastjetCentralNeutral = {fReader_Run2MC, "fixedGridRhoFastjetCentralNeutral"};
  
  //Run3 equivalent branch name [Naming convension changes]
  //TTreeReaderValue<Float_t> Rho_fixedGridRhoAll = {fReader_Run3MC, "Rho_fixedGridRhoAll"};
  TTreeReaderValue<Float_t> Rho_fixedGridRhoFastjetAll = {fReader_Run3MC, "Rho_fixedGridRhoFastjetAll"};
  //TTreeReaderValue<Float_t> Rho_fixedGridRhoFastjetCentral = {fReader_Run3MC, "Rho_fixedGridRhoFastjetCentral"};
  //TTreeReaderValue<Float_t> Rho_fixedGridRhoFastjetCentralCalo = {fReader_Run3MC, "Rho_fixedGridRhoFastjetCentralCalo"};
  //TTreeReaderValue<Float_t> Rho_fixedGridRhoFastjetCentralChargedPileUp = {fReader_Run3MC, "Rho_fixedGridRhoFastjetCentralChargedPileUp"};
  //TTreeReaderValue<Float_t> Rho_fixedGridRhoFastjetCentralNeutral = {fReader_Run3MC, "Rho_fixedGridRhoFastjetCentralNeutral"};


  //##################################################################################################
  
  VLLAna(TTree * /*tree*/ =0) { }
  virtual ~VLLAna() { }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();

  //User defined functions are declared here
  void SetHstFileName(const char *HstFileName){ _HstFileName  = HstFileName;}
  void SetSkimFileName(const char *SkimFileName){ _SkimFileName = SkimFileName;}
  void SetSumFileName(const char *SumFileName){ _SumFileName  = SumFileName;}
  void SetSample(int sample){_sample=sample;}
  void SetLep(int lep){_lep=lep;}
  void SetVerbose(int verbose){ _verbosity = verbose; }
  void SetData(int data){_data=data;}
  void SetYear(int year){_year = year;}
  void SetEra(TString era){_era=era;}
  void SetMCwt(int mcwt){_mcwt=mcwt;}
  void SetGenEventSumW(float value){_geneventsumw = value;}
  void BookHistograms();
  
  
  // My functions
  bool PassEventFlags(bool isData, int year);
  bool PassLeptonTrigger(int year,int LeptonType);
  int getBit ( int bitvar, int nthbit );
  float delta_phi(float phi1, float phi2);
  int  MotherID(int partindex, int momindex);
  int electronCustomID(Int_t bitmap,int quality, int skipCut);
  bool LeptonCustomIDCut(int year, float sip3d, float deepjet);
  float getInvMass(TLorentzVector a,TLorentzVector b);
  float getMOSSF();
  int getOSSFN();
  bool isOSSF(int i, int j);
  float OSSFMass(int i, int j);
  int NonBestOnZLeptonIndex();
  bool TriggerObjectMatching(int id,TLorentzVector t);
  bool MuonCleaning(TLorentzVector t, int opt);
  bool TaulepCleaning(TLorentzVector t);
  bool TaujetCleaning(TLorentzVector t);
  bool LepjetCleaning(TLorentzVector t);
  void TruthLevelInfo();
  
  
  //CorrectionFunctions
  float muonIDSF(float pt, float eta, string era, string mode);
  float muonIsoSF(float pt, float eta, string era, string mode);
  float egammaIDSF(float pt, float eta, string era, string mode);
  float egammaIDSF_Run3(float pt, float eta, float phi,string era, string mode);
  float LeptonIDSF(int id,float pt,float eta, float phi, string era, string mode);
  float LeptonISOSF(int id,float pt,float eta,string era, string mode);
  float jetSF(float pt, float eta, string era, string mode);
  float jetRF(float pt, float eta, float phi, int matchedjetidx, float genpt, float geneta, float genphi, float rho, string era, string mode);
  float btagWPSFfromPOG(float pt, float eta, int flav, string era, string mode);
  float btagMCEff(int MCsample, float pt, float eta, int flav, string era);
  float pileupWeight(float nTrueInt,string era, string mode);
  json load_json_data(int year);
  int checkJson(bool isData, int runno, int lsno);
  void CreatePassFailOutputJsonFiles(bool validRunLSFlag, int runno, int lsno);
  float muon_RochesterSF(int charge,float pt, float phi,float eta, bool isMC, bool genmatched, float genpt, int nl, string era,string mode);
  pair<float,float> METXYCorr_Met_MetPhi(float uncormet, float uncormet_phi, int runnb, string year, bool isData, int npv);

  //AnalysisTree
  void Analysis_LJJ(TTree* outputTree,string era, string jetsystvar);
  //void InitializeAnalysisTreeBranch(const std::string &outputFileName, TFile*& outputFile, TTree*& outputTree);
  void InitializeAnalysisTreeBranch(const std::string &outputFileName,TTree*& outputTree);
  
  
public:
  struct Hists {
    //Histograms are declared here.
    //BookHistogram()
    TH1F *nevt,*trigobj[8], *nlep[8], *ptlep[9],*etmiss[4];
    TH2F *trigobj2d[3];

    TH1F *channel1L2J[1];
    
    TH1F *genstudy[2];
    
    //TH1F *channelLJJ[16][1][2][15];
    
    //TH1F *mucut, *elcut, *taucut,*phocut;
    //TH1F *genpltmu[5],*genpltele[5],*genplttau[5],*genpltjet[9],*genPart[3];
    //TH1F *evtwt[4];
    //TH1F *tauacc[8],*vllclass,*region[2];
    
    //1L2J VLL Analysis Histograms
    //Book1L2JAnalysisHistogram()
    //TH1F *vllplotl2j[37][25],*testplotl2j[11],*l2jevtwt[4];
    //TH2F *test2Dplotl2j[3];
    
    //2L VLL Analysis Histograms
    //Book2LAnalysisHistogram()
    //TH1F *vllplot2l[30][12],*landscape2L,*testrac[5];

    
    // JANUARY 27 2023
    // PLOTTING ALL VARIABLES IN 1L2J FINAL STATE
    //TH1F *vllvarsl2j[63][3];
    //TH1F *morel2jvar[5];
  };
  
  struct Lepton {//The struct 'Lepton' can store the following variables:
    TLorentzVector v;
    int id;
    int ind;
    float wt;
    int flavor;
    int charge;
    bool lepcleaning;
    bool taucleaning;
    bool muoncleaning;
    int momid;
    int grandmomid;
    int genmatch;
    int jetmatch;
    int status;
    int pdgid;
    int partonflav;
    UChar_t hadronflav;
    int hadronflavor;
    bool bflav;
    float sip3d;
    float deepjet;
    float iso;
  };
  
  //PUBLIC FUNCTIONS
  //FUNCTIONS THAT DEPENDS ON struct<Lepton>
  void Sort(vector<Lepton>lepvec);
  float btagIDSF(int MCSample, vector<Lepton>Jet,string era, string mode);
  vector<float> btagIDSFv2(int MCSample, vector<Lepton>Jet,string era, string flav, string mode);
  vector<int>BestJetPairFinder(vector<Lepton>jetarray);

  
  //functions defined in different header file for dedicated study
  void GenStudy();
  float get_WMass();
  

protected:
  Hists h;

private:
  //Global variables go here. Make them global only if necessary.
  TFile *_HstFile,*_mvaFile,*_Analysis1L2J_HstFile,*_skimFile;
  TTree *tree,*skimTree,*tree1;
  
  TTree *outputTree1,*outputTree2,*outputTree3,*outputTree4,*outputTree5;
  //TFile *outputFile1,*outputFile2;
  vector<TFile*> outputFileVector;
  
  const char *_HstFileName,*_SkimFileName;
  const char *_SumFileName;
  bool produce_histos,produce_trees,channel_2L,channel_1L2J;
  
  int _verbosity,_exclude,_sample;
  float _geneventsumw;
  
  int nEvtTotal,nEvtRan,nEvtValidRunLS,nEvtTrigger,nEvtTree;
  float genEventsumw;
  json jsondata,passing_json,failing_json;
  
  int n_4L,n_3L1T,n_3L,n_2L2T,n_2L1T,n_1L3T,n_1L2T;
  int n_2L,n_1L2J,n_muJJ_passglobalsel,n_muJJ_passtreesel;
  
  //global variables to control different analysis settings
  TString AnalysisWP;
  bool apply_LHEHTcut,apply_LHEVptcut;
  
  
  int _data, _lep, _year,_mcwt;
  bool GoodEvt, GoodEvt2016, GoodEvt2017, GoodEvt2018,triggerRes,trigger2016,trigger2017,trigger2018,passeleTrigger,passmuTrigger,HEMVeto;
  float metpt, metphi,evwt,prob,evtwt,prob1,puppimetpt,puppimetphi,lheht;
  TString _era;
  
  //Lepton array
  vector<Lepton> llep, taus, Muon, loosemuon, Electron,jets,LooseLep,BTaggedJet;
  vector<Lepton> genPart,genMuon,genVisTau,genTau,truthtaus,genHadronicTau,genllep,genElectron,genjet;  
  
  
  time_t start, end;
  
  //MVA VARIABLES
  float lep0_pt,lep0_mt,jet0_pt,jet1_pt,jet0_mt,jet1_mt,dijet_mass,dijet_pt,dijet_mt;
  float deltaR_jet01,deltaPhi_metjet0,deltaPhi_metjet1,deltaPhi_metlep0,deltaPhi_jet0lep0;
  float deltaPhi_jet1lep0,deltaPhi_dijetlep0,deltaPhi_metdijet;
  float event_MET,event_METPhi,event_HT,event_ST;
  int nEvt,lep0_flavor,n_Jet,n_bJet;
  float genEventSumw;
  
  //EXTRA FOR ANALYSIS TREE
  float lep0_iso,lep0_phi,lep0_eta,lep0_sip3d,lep0_deepjet,jet0_phi,jet0_eta,jet1_phi,jet1_eta;
  bool lep0_tight;
  //SF FOR ANALYSIS TREE
  float event_btagsf,event_lepIDSF,event_lepIDSF_up,event_lepIDSF_down;
  float event_lepISOSF,event_lepISOSF_up,event_lepISOSF_down;
  float event_btagsfUncorrelated_up,event_btagsfUncorrelated_down,event_btagsfbcUncorrelated_up,event_btagsfbcUncorrelated_down;
  float event_btagsfCorrelated_up,event_btagsfCorrelated_down,event_btagsfbcCorrelated_up,event_btagsfbcCorrelated_down;
  
  //Added EXTRA Variables for wjets separation
  float deltaPhi_ljjsysmet,deltaPhi_ljjsyslep0,deltaPhi_ljjsysjet0,deltaPhi_ljjsysjet1,ljjsys_mass;
  float deepjetQG_jet0,deepjetQG_jet1,dijet_deltaeta,ljjsys_PT;
  int avgnpart,npart_jet0,npart_jet1;
  float event_avgCvsBscore,event_avgCvsLscore,event_avgQGscore,event_Rpt,event_zstar;

  float GenWeight,PU_N,PU_TrueN,L1PreFireWt_nom,L1PreFireWt_up,L1PreFireWt_down;
  float PileUpWt_nom,PileUpWt_up,PileUpWt_down,PDF_up,PDF_down,qcd_scale_up,qcd_scale_down;
  bool event_passmuTrigger, event_passeleTrigger,event_HEMVeto;
  
  //Extra Variables for 2L Analysis
  //float lep0_iso03,lep1_iso03,lep1_iso,lep1_phi,lep1_eta,lep1_sip3d,lep1_deepjet,lep1_pt,lep1_mt,lep1_flavor;
  //float lep0_charge,lep1_charge,n_lep;
  //float dilep_mass,dilep_pt,dilep_mt;
  //float deltaR_lep01,deltaPhi_metlep1,deltaPhi_lep01,deltaPhi_metdilep;
  //float event_LT,M_OS,M_SS,M_OSSF,M_OSOF,M_SSSF,M_SSOF;
  //bool passEleTrigger,passMuTrigger;

  //check
  //int n_2LOF,n_2LSF,n_WW,n_ZZ,n_HH,n_WZ,n_WH,n_ZH,n_1t,n_2t;
  
  ClassDef(VLLAna,0);
  
};

#endif

#ifdef VLLAna_cxx
void VLLAna::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the reader is initialized.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  cout<<"INIT function..."<<endl;
  cout<<"Data/Year/Sample: "<<_data<<"/"<<_year<<"/"<<_sample<<endl;

  fReader_common.SetTree(tree); //same name and same type branches Run2 and Run3
  if(_year<2020)fReader_Run2.SetTree(tree);//common branches of diff type + extra Run2 branches
  else if(_year>2020)fReader_Run3.SetTree(tree);//common branches of diff type + extra Run3 branches
  

  //MC branches
  if(_data == 0){                  
    fReader_commonMC.SetTree(tree);
    if(_year>2020)fReader_Run3MC.SetTree(tree);
    if(_year<2020)fReader_Run2MC.SetTree(tree);
    //special
    if((_sample<9) || (_sample>19))
      fReader_commonMCSpecial.SetTree(tree);
  }
  
  if(_year == 2016)                   //OnlyFor2016
    fReader_2016   .SetTree(tree);
  else if(_year == 2017)              //OnlyFor2017
    fReader_2017   .SetTree(tree);
  else if(_year == 2018)              //OnlyFor2018
    fReader_2018   .SetTree(tree);
  else if(_year == 2022)              //OnlyFor2022
    fReader_2022   .SetTree(tree);
  
  //Common
  if(_year == 2017 || _year==2018)fReader_2017_2018.SetTree(tree);
  if(_year == 2018 || _year==2022)fReader_2018_2022.SetTree(tree);
}

Bool_t VLLAna::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  
  return kTRUE;
}


#endif // #ifdef nano9Ana_cxx
