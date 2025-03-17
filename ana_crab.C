#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include<vector>
#include<TString.h>

void ana_crab( TString ifname , TString ofname, TString data, TString year, TString lep, TString era, TString MC)
{
  gROOT->Time();
  
  const char *hstfilename;
  const char *skimfilename;
  
  TChain *chain = new TChain("Events");
  VLLAna m_selec;
  //cout<<"Declared chains"<<endl;

  //Handle multiple input root files if needed.
  //Add all the root files one by one in chain
  std::vector<TString> input;
  std::string s(ifname.Data());
  size_t start = 0;
  size_t end = s.find(",");
  
  while (end != std::string::npos) {
    input.push_back(TString(s.substr(start, end - start).c_str()));
    start = end + 1;
    end = s.find(",", start);
  }
  
  input.push_back(TString(s.substr(start).c_str()));

  TString redirector="root://cms-xrd-global.cern.ch//";
  for (const auto& file : input) {
    std::cout << "Processing file: " << redirector+file << std::endl;
    chain->Add(redirector+file);
  }
  
  //define output files
  skimfilename = ofname;
  hstfilename = "hstFile.root";
  //cout<<"Output files are "<<hstfilename<<endl;
  m_selec.SetHstFileName(hstfilename);
  m_selec.SetSkimFileName(skimfilename);

  cout<<"ana_crab.C| HstFile: "<<hstfilename<<"| SkimFile: "<<skimfilename<<endl;
  

  //define verbose, data/mc,year,dataset, mc name/sample(int)
  m_selec.SetVerbose(1);
  //  m_selec.SetData(boost::lexical_cast<int>(data)); //0 - not running over data, 1 - running over data
  if(data=="0")
    m_selec.SetData(0); //0 - not running over data, 1 - running over data
  if(data=="1")
    m_selec.SetData(1); //0 - not running over data, 1 - running over data
  
  //  m_selec.SetSample(0); //_sample = 0-WZ, 1-ZZ, 2-TTW, 3-TTZ, 4-TTH, 5-WWW, 6-WWZ, 7-DY, 8-TTDilep, 9-TTSlep, 10-TTBarSlep
  if(year=="2016")
    m_selec.SetYear(2016);
  if(year=="2017")
    m_selec.SetYear(2017);
  if(year=="2018")
    m_selec.SetYear(2018);
  m_selec.SetEra(era);
  if(lep=="el")
    m_selec.SetLep(0);
  if(lep=="mu")
    m_selec.SetLep(1);
  if(era=="preVFP")
    m_selec.SetEra("preVFP");
  if(era=="postVFP")
    m_selec.SetEra("postVFP");
  
  //cout<<"Set the options"<<endl;

  float MC2016_normwt[13] = {1,1.02,1,1.14,1,1,1,1,1,1,1,1,1}; //_sample = 0-WZ, 1-ZZ, 2-TTW, 3-TTZ, 4-TTH, 5-WWW, 6-WWZ, 7-WZZ, 8-ZZZ, 9-signal, 10-DY, 11-TTDilep, 12-rare
  float MC2017_normwt[13] = {1.03,1.00,1,1.40,1,1,1,1,1,1,1,1,1}; //_sample = 0-WZ, 1-ZZ, 2-TTW, 3-TTZ, 4-TTH, 5-WWW, 6-WWZ, 7-WZZ, 8-ZZZ, 9-signal, 10-DY, 11-TTDilep, 12-rare
  float MC2018_normwt[13] = {1.,1.,1,1.,1,1,1,1,1,1,1,1,1}; //_sample = 0-WZ, 1-ZZ, 2-TTW, 3-TTZ, 4-TTH, 5-WWW, 6-WWZ, 7-WZZ, 8-ZZZ, 9-signal, 10-DY, 11-TTDilep, 12-rare

  if(data=="0"){
    int _sample = 80;
    
    //SampleNumbering
    // DY: 1-9
    // QCD: 10-20
    // ZZ : 20-29
    // WZ : 30-39
    // WW : 40-49
    // SingleTop: 50-59
    // TTBar: 60-69
    // HTBinnedWJets: 70-79
    // VLL: 80-90

    if(MC=="DYJetsToLL_M10to50")_sample = 1;
    if(MC=="DYJetsToLL_M50")_sample = 2;

    if(MC=="QCD_MuEnriched_20to30")_sample = 10;
    if(MC=="QCD_MuEnriched_30to50")_sample = 11;
    if(MC=="QCD_MuEnriched_50to80")_sample = 12;
    if(MC=="QCD_MuEnriched_80to120")_sample = 13;
    if(MC=="QCD_MuEnriched_120to170")_sample = 14;
    if(MC=="QCD_MuEnriched_170to300")_sample = 15;
    if(MC=="QCD_MuEnriched_300to470")_sample = 16;
    if(MC=="QCD_MuEnriched_470to600")_sample = 17;
    if(MC=="QCD_MuEnriched_600to800")_sample = 18;
    if(MC=="QCD_MuEnriched_800to1000")_sample = 19;

    if(MC=="ZZ_ZZTo2Q2L")_sample = 20;
    if(MC=="ZZ_ZZTo2Q2Nu")_sample = 21;
    if(MC=="ZZ_ZZTo4L")_sample = 22;
    if(MC=="ZZ_ZZTo2L2Nu")_sample = 23;
    

    if(MC=="WZ_WZTo2Q2L")_sample = 30;
    if(MC=="WZ_WZTo3LNu")_sample = 31;
    if(MC=="WZ_WZTo1L1Nu2Q")_sample = 32;
    
    
    if(MC=="WW_WWTo1L1Nu2Q")_sample = 40;
    if(MC=="WW_WWto2L2Nu")_sample = 41;
    if(MC=="WW_WWto4Q")_sample = 42;

    if(MC=="SingleTop_t-channel_Top_InclusiveDecays")_sample = 50;
    if(MC=="SingleTop_t-channel_AntiTop_InclusiveDecays")_sample = 51;
    if(MC=="SingleTop_s-channel_LeptonDecays")_sample = 52;
    if(MC=="SingleTop_tW_AntiTop_InclusiceDecays")_sample = 53;
    if(MC=="SingleTop_tW_Top_InclusiveDecays")_sample = 54;

    if(MC=="TTBar_TTToSemiLeptonic")_sample = 60;
    if(MC=="TTBar_TTTo2L2Nu")_sample = 61;
    
    if(MC=="HTbinnedWJets_Inclusive")_sample = 70;
    if(MC=="HTbinnedWJets_70to100")_sample = 70;
    if(MC=="HTbinnedWJets_100to200")_sample = 70;
    if(MC=="HTbinnedWJets_200to400")_sample = 70;
    if(MC=="HTbinnedWJets_400to600")_sample = 70;
    if(MC=="HTbinnedWJets_600to800")_sample = 70;
    if(MC=="HTbinnedWJets_800to1200")_sample = 70;
    if(MC=="HTbinnedWJets_1200to2500")_sample = 70;
    if(MC=="HTbinnedWJets_2500toInf")_sample = 70;
    
    if(MC=="VLL_M100")_sample = 80;
    if(MC=="VLL_M125")_sample = 81;
    if(MC=="VLL_M150")_sample = 82;
    if(MC=="VLL_M200")_sample = 83;
    if(MC=="VLL_M250")_sample = 84;
    if(MC=="VLL_M300")_sample = 85;
    if(MC=="VLL_M350")_sample = 86;
    if(MC=="VLL_M400")_sample = 87;
    
    m_selec.SetSample(_sample);
    
    //if(year=="2016")
    //  m_selec.SetMCwt(MC2016_normwt[_sample]);
    //if(year=="2017")
    //  m_selec.SetMCwt(MC2017_normwt[_sample]);
    //if(year=="2018")
    //  m_selec.SetMCwt(MC2018_normwt[_sample]);
  }
  
  chain->Process(&m_selec);
  gROOT->Time();
}
