#include "correction.h"

auto egammajson2022      = correction::CorrectionSet::from_file("setup/POG/EGM/2022_Summer22/electron.json.gz");
auto egammaID2022        = egammajson2022->at("Electron-ID-SF");
auto egammajson2022EE    = correction::CorrectionSet::from_file("setup/POG/EGM/2022_Summer22EE/electron.json.gz");
auto egammaID2022EE      = egammajson2022EE->at("Electron-ID-SF");
auto egammajson2023      = correction::CorrectionSet::from_file("setup/POG/EGM/2023_Summer23/electron.json.gz");
auto egammaID2023        = egammajson2023->at("Electron-ID-SF");
auto egammajson2023BPix  = correction::CorrectionSet::from_file("setup/POG/EGM/2023_Summer23BPix/electron.json.gz");
auto egammaID2023BPix    = egammajson2023BPix->at("Electron-ID-SF");

auto egammajson2018      = correction::CorrectionSet::from_file("setup/POG/EGM/2018_UL/electron.json.gz");
auto egammaID2018        = egammajson2018->at("UL-Electron-ID-SF");
//auto egammaScalejson2018 = correction::CorrectionSet::from_file("setup/POG/EGM/2018_UL/EGM_ScaleUnc.json.gz");
//auto egammaScale2018     = egammaScalejson2018->at("UL-EGM_ScaleUnc");

auto egammajson2017      = correction::CorrectionSet::from_file("setup/POG/EGM/2017_UL/electron.json.gz");
auto egammaID2017        = egammajson2017->at("UL-Electron-ID-SF");
//auto egammaScalejson2017 = correction::CorrectionSet::from_file("setup/POG/EGM/2017_UL/EGM_ScaleUnc.json.gz");
//auto egammaScale2017     = egammaScalejson2017->at("UL-EGM_ScaleUnc");

auto egammajson2016postVFP      = correction::CorrectionSet::from_file("setup/POG/EGM/2016postVFP_UL/electron.json.gz");
auto egammaID2016postVFP        = egammajson2016postVFP->at("UL-Electron-ID-SF");
//auto egammaScalejson2016postVFP = correction::CorrectionSet::from_file("setup/POG/EGM/2016postVFP_UL/EGM_ScaleUnc.json.gz");
//auto egammaScale2016postVFP     = egammaScalejson2016postVFP->at("UL-EGM_ScaleUnc");

auto egammajson2016preVFP      = correction::CorrectionSet::from_file("setup/POG/EGM/2016preVFP_UL/electron.json.gz");
auto egammaID2016preVFP        = egammajson2016preVFP->at("UL-Electron-ID-SF");
//auto egammaScalejson2016preVFP = correction::CorrectionSet::from_file("setup/POG/EGM/2016preVFP_UL/EGM_ScaleUnc.json.gz");
//auto egammaScale2016preVFP     = egammaScalejson2016preVFP->at("UL-EGM_ScaleUnc");

float VLLAna::egammaIDSF(float pt, float eta,string era, string mode){

  //https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/EGM_2016postVFP_UL_electron.html
  
  string mode_;
  if(      mode=="nom"  ) mode_="sf";
  else if( mode=="up"   ) mode_="sfup";
  else if( mode=="down" ) mode_="sfdown";
  
  std::vector<correction::Variable::Type>  values;
  values.emplace_back(era);
  values.emplace_back(mode_);
  values.emplace_back("Medium");
  values.emplace_back(eta);//Do we need to add deltaEtaSC??
  values.emplace_back(pt);
  
  float sf = 1.0;
  if( pt >10.){
    if( era=="2018" ){
      sf = egammaID2018->evaluate(values); 
    }
    if( era=="2017" ){
      sf = egammaID2017->evaluate(values); 
    }
    if( era=="2016postVFP" ){
      sf = egammaID2016postVFP->evaluate(values); 
    }
    if( era=="2016preVFP" ){
      sf = egammaID2016preVFP->evaluate(values); 
    }
  }
  
  return sf;
  
}

float VLLAna::egammaIDSF_Run3(float pt, float eta, float phi,string era, string mode){
  
  //https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/summaries/EGM_2022_Summer22_electron.html
  
  string mode_;
  if(      mode=="nom"  ) mode_="sf";
  else if( mode=="up"   ) mode_="sfup";
  else if( mode=="down" ) mode_="sfdown";

  string era_="";
  if(era=="2022")era_="2022Re-recoBCD";
  if(era=="2022EE")era_="2022Re-recoE+PromptFG";
  if(era=="2023")era_="2023PromptC";
  if(era=="2023BPix")era_="2023PromptD";
  
  std::vector<correction::Variable::Type>  values;
  values.emplace_back(era);
  values.emplace_back(mode_);
  values.emplace_back("Medium");
  values.emplace_back(eta);//Do we need to add deltaEtaSC??
  values.emplace_back(pt);
  //Phi is needed for 2023 campaign
  if(era=="2023")values.emplace_back(phi);
  if(era=="2023BPix")values.emplace_back(phi);
  
  float sf = 1.0;
  if( pt >10.){
    if     (era=="2022")sf=egammaID2022->evaluate(values);
    else if(era=="2022EE")sf=egammaID2022EE->evaluate(values);
    else if(era=="2023")sf=egammaID2023->evaluate(values);
    else if(era=="2023BPix")sf=egammaID2023BPix->evaluate(values);
  }
  
  return sf;
  
}
