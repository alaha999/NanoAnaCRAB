#include "correction.h"
#include "TRandom3.h"

auto json2018 = correction::CorrectionSet::from_file("setup/POG/JME/2018_UL/jet_jerc.json.gz"); 
auto jecsflabel2018 = "Summer19UL18_V5_MC_Total_AK4PFchs";
auto jecsfjson2018  = json2018->at(jecsflabel2018);
auto jersflabel2018 = "Summer19UL18_JRV2_MC_ScaleFactor_AK4PFchs";
auto jersfjson2018  = json2018->at(jersflabel2018);
auto ptreslabel2018 = "Summer19UL18_JRV2_MC_PtResolution_AK4PFchs";
auto ptresjson2018  = json2018->at(ptreslabel2018);

auto json2017 = correction::CorrectionSet::from_file("setup/POG/JME/2017_UL/jet_jerc.json.gz"); 
auto jecsflabel2017 = "Summer19UL17_V5_MC_Total_AK4PFchs";
auto jecsfjson2017  = json2017->at(jecsflabel2017);
auto jersflabel2017 = "Summer19UL17_JRV2_MC_ScaleFactor_AK4PFchs";
auto jersfjson2017  = json2017->at(jersflabel2017);
auto ptreslabel2017 = "Summer19UL17_JRV2_MC_PtResolution_AK4PFchs";
auto ptresjson2017  = json2017->at(ptreslabel2017);

auto json2016postVFP = correction::CorrectionSet::from_file("setup/POG/JME/2016postVFP_UL/jet_jerc.json.gz"); 
auto jecsflabel2016postVFP = "Summer19UL16_V7_MC_Total_AK4PFchs";
auto jecsfjson2016postVFP  = json2016postVFP->at(jecsflabel2016postVFP);
auto jersflabel2016postVFP = "Summer20UL16_JRV3_MC_ScaleFactor_AK4PFchs";
auto jersfjson2016postVFP  = json2016postVFP->at(jersflabel2016postVFP);
auto ptreslabel2016postVFP = "Summer20UL16_JRV3_MC_PtResolution_AK4PFchs";
auto ptresjson2016postVFP  = json2016postVFP->at(ptreslabel2016postVFP);

auto json2016preVFP = correction::CorrectionSet::from_file("setup/POG/JME/2016preVFP_UL/jet_jerc.json.gz"); 
auto jecsflabel2016preVFP = "Summer19UL16APV_V7_MC_Total_AK4PFchs";
auto jecsfjson2016preVFP  = json2016preVFP->at(jecsflabel2016preVFP);
auto jersflabel2016preVFP = "Summer20UL16APV_JRV3_MC_ScaleFactor_AK4PFchs";
auto jersfjson2016preVFP  = json2016preVFP->at(jersflabel2016preVFP);
auto ptreslabel2016preVFP = "Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs";
auto ptresjson2016preVFP  = json2016preVFP->at(ptreslabel2016preVFP);



float VLLAna::jetSF(float pt, float eta, string era, string mode){
  
  std::vector<correction::Variable::Type>  values;
  values.emplace_back(eta);
  values.emplace_back(pt);
  
  float unc = 0.;
  float sf  = 1.0;

  if( era=="2018" ){
    unc=jecsfjson2018->evaluate(values);
  }
  if( era=="2017" ){
    unc=jecsfjson2017->evaluate(values);
  }
  if( era=="2016postVFP" ){
    unc=jecsfjson2016postVFP->evaluate(values);
  }
  if( era=="2016preVFP" ){
    unc=jecsfjson2016preVFP->evaluate(values);
  }
  
  if (      mode =="up"  ) {sf = 1.0+unc;}
  else if ( mode =="down") {sf = 1.0-unc;}
  else                     {sf =1.;}
  
  return sf;
}

float VLLAna::jetRF(float pt, float eta, float phi, int matchedjetidx, float genpt, float geneta, float genphi, float rho, string era, string mode){
  
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
  
  std::vector<correction::Variable::Type>  jersfvalues;
  std::vector<correction::Variable::Type>  ptresvalues;
  jersfvalues.emplace_back(eta);
  jersfvalues.emplace_back(mode);
  //
  ptresvalues.emplace_back(eta);
  ptresvalues.emplace_back(pt);
  ptresvalues.emplace_back(rho); //NanoAODv9 branch: fixedGridRhoFastjetAll
  //
  float sjer=1.0;
  float sigmajer=0.;
  float cjer =1.0;
  //
  if( era=="2018" ){
    sjer  = jersfjson2018->evaluate(jersfvalues);
    sigmajer = ptresjson2018->evaluate(ptresvalues);
  }
  if( era=="2017" ){
    sjer  = jersfjson2017->evaluate(jersfvalues);
    sigmajer = ptresjson2017->evaluate(ptresvalues);
  }
  if( era=="2016postVFP" ){
    sjer  = jersfjson2016postVFP->evaluate(jersfvalues);
    sigmajer = ptresjson2016postVFP->evaluate(ptresvalues);
  }
  if( era=="2016preVFP" ){
    sjer  = jersfjson2016preVFP->evaluate(jersfvalues);
    sigmajer = ptresjson2016preVFP->evaluate(ptresvalues);
  }
  
  //genjet matching
  TLorentzVector recojetv;
  recojetv.SetPtEtaPhiM(pt,eta,phi,0);

  bool isGenMatched = false;
  
  if(matchedjetidx>0){
    TLorentzVector genjetv;
    genjetv.SetPtEtaPhiM(genpt,geneta,genphi,0);
    float dr = recojetv.DeltaR(genjetv);
    
    if(dr<0.2 && (abs(pt - genpt)<(3*pt*sigmajer)))isGenMatched=true;
  }
  
  //Hybrid method(Usual scaling + stochastic scaling)
  if(isGenMatched){
    cjer = 1.+(sjer-1.)*((pt-genpt)/pt);
  }
  else{
    TRandom3* gRandom = new TRandom3();
    gRandom->SetSeed(0);
    float n = gRandom->Gaus( 0., sigmajer );
    cjer = 1.+n*(sqrt(max( sjer*sjer-1., 0. )));
  }
  if(cjer<0.){cjer=0.;} //cjer is truncated at zero as required
  

  return cjer;
}
