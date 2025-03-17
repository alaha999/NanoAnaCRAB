#include "correction.h"

using namespace std;

auto muonjson2022 = correction::CorrectionSet::from_file("setup/POG/MUO/2022_Summer22/muon_Z.json.gz");
auto muoniso2022  = muonjson2022->at("NUM_TightPFIso_DEN_MediumID");
auto muonid2022   = muonjson2022->at("NUM_MediumID_DEN_TrackerMuons");  
  
auto muonjson2022EE = correction::CorrectionSet::from_file("setup/POG/MUO/2022_Summer22EE/muon_Z.json.gz");
auto muoniso2022EE  = muonjson2022EE->at("NUM_TightPFIso_DEN_MediumID");
auto muonid2022EE   = muonjson2022EE->at("NUM_MediumID_DEN_TrackerMuons");
  
auto muonjson2023 = correction::CorrectionSet::from_file("setup/POG/MUO/2023_Summer23/muon_Z.json.gz");
auto muoniso2023  = muonjson2023->at("NUM_TightPFIso_DEN_MediumID");
auto muonid2023   = muonjson2023->at("NUM_MediumID_DEN_TrackerMuons");

auto muonjson2023BPix = correction::CorrectionSet::from_file("setup/POG/MUO/2023_Summer23BPix/muon_Z.json.gz");
auto muoniso2023BPix  = muonjson2023BPix->at("NUM_TightPFIso_DEN_MediumID");
auto muonid2023BPix   = muonjson2023BPix->at("NUM_MediumID_DEN_TrackerMuons");  
  
auto muonjson2018 = correction::CorrectionSet::from_file("setup/POG/MUO/2018_UL/muon_Z.json.gz");
auto muoniso2018  = muonjson2018->at("NUM_TightRelIso_DEN_MediumID");
auto muonid2018   = muonjson2018->at("NUM_MediumID_DEN_genTracks");
  
auto muonjson2017 = correction::CorrectionSet::from_file("setup/POG/MUO/2017_UL/muon_Z.json.gz");
auto muoniso2017  = muonjson2017->at("NUM_TightRelIso_DEN_MediumID");
auto muonid2017   = muonjson2017->at("NUM_MediumID_DEN_genTracks");
  
auto muonjson2016postVFP = correction::CorrectionSet::from_file("setup/POG/MUO/2016postVFP_UL/muon_Z.json.gz");
auto muoniso2016postVFP  = muonjson2016postVFP->at("NUM_TightRelIso_DEN_MediumID");
auto muonid2016postVFP   = muonjson2016postVFP->at("NUM_MediumID_DEN_genTracks");
  
auto muonjson2016preVFP = correction::CorrectionSet::from_file("setup/POG/MUO/2016preVFP_UL/muon_Z.json.gz");
auto muoniso2016preVFP  = muonjson2016preVFP->at("NUM_TightRelIso_DEN_MediumID");
auto muonid2016preVFP   = muonjson2016preVFP->at("NUM_MediumID_DEN_genTracks");

//electrons
auto egammajson2022      = correction::CorrectionSet::from_file("setup/POG/EGM/2022_Summer22/electron.json.gz");
auto egammaID2022        = egammajson2022->at("Electron-ID-SF");
auto egammajson2022EE    = correction::CorrectionSet::from_file("setup/POG/EGM/2022_Summer22EE/electron.json.gz");
auto egammaID2022EE      = egammajson2022EE->at("Electron-ID-SF");
auto egammajson2023      = correction::CorrectionSet::from_file("setup/POG/EGM/2023_Summer23/electron.json.gz");
auto egammaID2023        = egammajson2023->at("Electron-ID-SF");
auto egammajson2023BPix  = correction::CorrectionSet::from_file("setup/POG/EGM/2023_Summer23BPix/electron.json.gz");
auto egammaID2023BPix    = egammajson2023BPix->at("Electron-ID-SF");



void testSF(){
  cout<<"Test file for SF"<<endl;
  
  float pt =50;
  float eta=1.5;
  string mode="nominal";
  
  std::vector<correction::Variable::Type>values;
  values.emplace_back(fabs(eta));
  values.emplace_back(pt);
  values.emplace_back(mode);

  //muons
  cout<<"---------------------"<<endl;
  cout<<"2022_Summer22: muon ID/ISO SF      : "<<muonid2022->evaluate(values)<<"/"<<muoniso2022->evaluate(values)<<endl;
  cout<<"2022_Summer22EE: muon ID/ISO SF    : "<<muonid2022EE->evaluate(values)<<"/"<<muoniso2022EE->evaluate(values)<<endl;
  cout<<"2023_Summer23: muon ID/ISO SF      : "<<muonid2023->evaluate(values)<<"/"<<muoniso2023->evaluate(values)<<endl;
  cout<<"2023_Summer23BPix: muon ID/ISO SF  : "<<muonid2023BPix->evaluate(values)<<"/"<<muoniso2023BPix->evaluate(values)<<endl;

  //electrons
  cout<<"---------------------"<<endl;
  cout<<"2022_Summer22 ele IDSF    : "<<egammaID2022->evaluate({"2022Re-recoBCD","sf","Medium",eta,pt})<<endl;
  cout<<"2022EE_Summer22 ele IDSF  : "<<egammaID2022EE->evaluate({"2022Re-recoE+PromptFG","sf","Medium",eta,pt})<<endl;
  float phi=1.0;
  cout<<"2023_Summer23 ele IDSF    : "<<egammaID2023->evaluate({"2023PromptC","sf","Medium",eta,pt,phi})<<endl;
  cout<<"2023BPix_Summer23 ele IDSF: "<<egammaID2023BPix->evaluate({"2023PromptD","sf","Medium",eta,pt,phi})<<endl;
  
}
