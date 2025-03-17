#include "correction.h"

auto pujson2018 = correction::CorrectionSet::from_file("setup/POG/LUM/2018_UL/puWeights.json.gz");
auto pu2018 = pujson2018->at("Collisions18_UltraLegacy_goldenJSON");

auto pujson2017 = correction::CorrectionSet::from_file("setup/POG/LUM/2017_UL/puWeights.json.gz");
auto pu2017 = pujson2017->at("Collisions17_UltraLegacy_goldenJSON");

auto pujson2016preVFP = correction::CorrectionSet::from_file("setup/POG/LUM/2016preVFP_UL/puWeights.json.gz");
auto pu2016preVFP = pujson2016preVFP->at("Collisions16_UltraLegacy_goldenJSON");

auto pujson2016postVFP = correction::CorrectionSet::from_file("setup/POG/LUM/2016postVFP_UL/puWeights.json.gz");
auto pu2016postVFP = pujson2016postVFP->at("Collisions16_UltraLegacy_goldenJSON");


float VLLAna::pileupWeight(float nTrueInt,string era, string mode){
  string mode_;
  if      (mode=="nom"  ) mode_="nominal";
  else if (mode=="up"   ) mode_="up";
  else if (mode=="down" ) mode_="down";
  
  std::vector<correction::Variable::Type> values;
  values.emplace_back(nTrueInt);
  values.emplace_back(mode_);
  
  //cout<<"era/mode/nTrueInt="<<"/"<<era<<"/"<<mode<<"/"<<nTrueInt<<endl;
  
  float pileupWt = 1.0;
  if      (era=="2018"        ) pileupWt = pu2018->evaluate(values);
  else if (era=="2017"        ) pileupWt = pu2017->evaluate(values);
  else if (era=="2016preVFP"  ) pileupWt = pu2016preVFP->evaluate(values);
  else if (era=="2016postVFP" ) pileupWt = pu2016postVFP->evaluate(values);
  
  return pileupWt;
}
