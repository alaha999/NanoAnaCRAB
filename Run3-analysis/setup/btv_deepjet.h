#include "btv_deepjet_MCEff.h"
#include "correction.h"

auto btvjson2018      = correction::CorrectionSet::from_file("setup/POG/BTV/2018_UL/btagging.json.gz");
auto bcjet2018        = btvjson2018->at("deepJet_comb");
auto lightjet2018     = btvjson2018->at("deepJet_incl");

auto btvjson2017      = correction::CorrectionSet::from_file("setup/POG/BTV/2017_UL/btagging.json.gz");
auto bcjet2017        = btvjson2017->at("deepJet_comb");
auto lightjet2017     = btvjson2017->at("deepJet_incl");

auto btvjson2016postVFP      = correction::CorrectionSet::from_file("setup/POG/BTV/2016postVFP_UL/btagging.json.gz");
auto bcjet2016postVFP        = btvjson2016postVFP->at("deepJet_comb");
auto lightjet2016postVFP     = btvjson2016postVFP->at("deepJet_incl");

auto btvjson2016preVFP      = correction::CorrectionSet::from_file("setup/POG/BTV/2016preVFP_UL/btagging.json.gz");
auto bcjet2016preVFP        = btvjson2016preVFP->at("deepJet_comb");
auto lightjet2016preVFP     = btvjson2016preVFP->at("deepJet_incl");


float VLLAna::btagWPSFfromPOG(float pt, float eta, int flav, string era, string mode){

  string mode_;
  if(      mode=="nom"              ) mode_="central";
  else if( mode=="upUncorrelated"   ) mode_="up_uncorrelated";
  else if( mode=="upCorrelated"     ) mode_="up_correlated";
  else if( mode=="downUncorrelated" ) mode_="down_uncorrelated";
  else if( mode=="downCorrelated"   ) mode_="down_correlated";
  
  string wp = "M";   //Working Point of b-tagging algorithm
  
  std::vector<correction::Variable::Type>  values;
  values.emplace_back(mode_);
  values.emplace_back(wp);
  values.emplace_back(flav);
  values.emplace_back(abs(eta));
  values.emplace_back(pt);

  float sf =1.0;
  
  if(pt>20.0 && fabs(eta)<2.5){
    if( era=="2018" ){
      if(flav==5 || flav==4){ sf = bcjet2018->evaluate(values);    }
      else                  { sf = lightjet2018->evaluate(values); }
    }
    if( era=="2017" ){
      if(flav==5 || flav==4){ sf = bcjet2017->evaluate(values);    }
      else                  { sf = lightjet2017->evaluate(values); }
    }
    if( era=="2016postVFP" ){
      if(flav==5 || flav==4){ sf = bcjet2016postVFP->evaluate(values);    }
      else                  { sf = lightjet2016postVFP->evaluate(values); }
    }
    if( era=="2016preVFP" ){
      if(flav==5 || flav==4){ sf = bcjet2016preVFP->evaluate(values);    }    
      else                  { sf = lightjet2016preVFP->evaluate(values); }
    }     
  }

  return sf;
}


float VLLAna::btagIDSF(int MCSample, vector<Lepton>Jet,string era, string mode){

  
  float probability_mc = 1.0;
  float probability_data = 1.0;
  
  
  //btag SF is basically an event reweighting procedure
  for (int i = 0; i < (int)Jet.size(); i++){
    float jet_prob_mc = 1.0;
    float jet_prob_data = 1.0;
    float jet_eff = 1.0; //MC efficiency
    
    
    float jetpt  = Jet.at(i).v.Pt();
    float jeteta = Jet.at(i).v.Eta();
    int jetflav  = Jet.at(i).hadronflavor;
    
    //get MC efficiency
    jet_eff = btagMCEff(MCSample,jetpt,jeteta,jetflav,era);
    

    //SFfromPOG
    float SFfromPOG = btagWPSFfromPOG(jetpt,jeteta,jetflav,era,mode);

    //DeepJetWPThreshold (https://btv-wiki.docs.cern.ch/ScaleFactors/#useful-links)
    //  Year          L           M          T
    //  2018         0.0490     0.2783     0.7100
    //  2017         0.0532     0.3040     0.7476
    //  2016preVFP   0.0508     0.2598     0.6502
    //  2016postVFP  0.0480     0.2489     0.6377
    
    float WPth = 0.;
    if(era =="2018")WPth=0.2783;
    else if (era =="2017")WPth=0.3040;
    else if (era =="2016preVFP")WPth=0.2598;
    else if (era =="2016postVFP")WPth=0.2489;
    
    // check if jet is tagged or not
    if ( Jet_btagDeepB[Jet.at(i).ind] > WPth){
      jet_prob_data = jet_eff*SFfromPOG;
      jet_prob_mc   = jet_eff;
    }
    else {    
      jet_prob_data = 1-jet_eff*SFfromPOG;
      jet_prob_mc   = 1-jet_eff;
    }
    
    probability_data *= jet_prob_data;
    probability_mc   *= jet_prob_mc;
  }
  
  float scaleFactor = 1.0;
  if (probability_mc > 0.0) scaleFactor = probability_data/probability_mc;
  
  return scaleFactor;
  
}


// return scale factor but bc and light separately
vector<float> VLLAna::btagIDSFv2(int MCSample, vector<Lepton>Jet,string era, string flav, string mode){
  
  
  float probability_mc_bc = 1.0;
  float probability_data_bc = 1.0;
  float probability_mc_light = 1.0;
  float probability_data_light = 1.0;  
  
  //btag SF is basically an event reweighting procedure
  for (int i = 0; i < (int)Jet.size(); i++){
    float jet_prob_mc_bc = 1.0;
    float jet_prob_data_bc = 1.0;
    float jet_prob_mc_light = 1.0;
    float jet_prob_data_light = 1.0;    
    float jet_eff = 1.0; //MC efficiency
    
    
    float jetpt  = Jet.at(i).v.Pt();
    float jeteta = Jet.at(i).v.Eta();
    int jetflav  = Jet.at(i).hadronflavor;
    
    bool bcflav    = (jetflav==5 || jetflav==4);
    bool lightflav = jetflav==0;

    //get MC efficiency
    jet_eff = btagMCEff(MCSample,jetpt,jeteta,jetflav,era);
    
    

    //SFfromPOG
    float SFfromPOG = btagWPSFfromPOG(jetpt,jeteta,jetflav,era,mode);
    
    //DeepJetWPThreshold (https://btv-wiki.docs.cern.ch/ScaleFactors/#useful-links)
    //  Year          L           M          T
    //  2018         0.0490     0.2783     0.7100
    //  2017         0.0532     0.3040     0.7476
    //  2016preVFP   0.0508     0.2598     0.6502
    //  2016postVFP  0.0480     0.2489     0.6377
    
    float WPth = 0.;
    if(era =="2018")WPth=0.2783;
    else if (era =="2017")WPth=0.3040;
    else if (era =="2016preVFP")WPth=0.2598;
    else if (era =="2016postVFP")WPth=0.2489;
    
    // check if jet is tagged or not
    if ( Jet_btagDeepB[Jet.at(i).ind] > WPth){
      if(bcflav)jet_prob_data_bc = jet_eff*SFfromPOG,jet_prob_mc_bc   = jet_eff;
      else if(lightflav)jet_prob_data_light = jet_eff*SFfromPOG,jet_prob_mc_light   = jet_eff;
    }
    else {      
      if(bcflav)jet_prob_data_bc = 1-jet_eff*SFfromPOG,jet_prob_mc_bc   = 1-jet_eff;
      else if(lightflav)jet_prob_data_light = 1- jet_eff*SFfromPOG,jet_prob_mc_light   = 1-jet_eff;
    }
    
    probability_data_bc *= jet_prob_data_bc;
    probability_mc_bc   *= jet_prob_mc_bc;
    
    probability_data_light *= jet_prob_data_light;
    probability_mc_light   *= jet_prob_mc_light;
  }
  
  //float scaleFactor = 1.0;
  //if(flav=="bc"){
  //  if (probability_mc_bc > 0.0) scaleFactor = probability_data_bc/probability_mc_bc;vecsf.emplace_back(scaleFactor);
  // }
  //else if(flav=="light"){
  //  if (probability_mc_light > 0.0) scaleFactor = probability_data_light/probability_mc_light;vecsf.emplace_back(scaleFactor);
  //}

  float bcSF=1.0;
  float bclight=1.0;
  vector<float>vecsf;
  if (probability_mc_bc > 0.0)    bcSF    = probability_data_bc/probability_mc_bc;
  if (probability_mc_light > 0.0) bclight = probability_data_light/probability_mc_light;

  vecsf.emplace_back(bcSF);
  vecsf.emplace_back(bclight);
  
  //return scaleFactor;
  return vecsf;
}
