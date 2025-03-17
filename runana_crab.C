#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
void runana_crab(TString ifname="input_file" ,TString ofname="outputfile" ,TString data="data",TString year="year",TString lep="lep",TString era="era",TString MC="mc")
{
  TString anastring =".x ana_crab.C(\""+ifname+"\",\""+ofname+"\",\""+data+"\",\""+year+"\",\""+lep+"\",\""+era+"\",\""+MC+"\")";
  //gSystem->Load("VLLAna_C.so");
  gSystem->CompileMacro("VLLAna.C", "k");
  gROOT->ProcessLine(anastring);
}
