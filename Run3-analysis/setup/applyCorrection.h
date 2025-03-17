#include "muon.h"
#include "egamma.h"
//#include "jetmet.h"
//#include "btv_deepjet.h"
//#include "puweight.h"
//#include "muonRochesterCorr.h"
#include "lumiJsonFilter.h"

float VLLAna::LeptonIDSF(int id,float pt,float eta, float phi, string era, string mode){

  //comments: phi variable is only used for electron ID case and 2023, 2023BPix campaign
  
  bool is_muon=(abs(id)==13);
  bool is_ele =(abs(id)==11);
  
  float sf=1.0;
  
  if(is_muon){
    sf = muonIDSF(pt,eta,era,mode);    
  }
  if(is_ele){
    if(era=="2023" || era=="2023BPix")sf = egammaIDSF_Run3(pt,eta,phi,era,mode);
    else sf = egammaIDSF(pt,eta,era,mode);
  }
  
  return sf;
}
  
float VLLAna::LeptonISOSF(int id,float pt,float eta,string era, string mode){
  
  bool is_muon=(abs(id)==13);
  bool is_ele =(abs(id)==11);
  
  float sf=1.0;
  
  if(is_muon){
    sf = muonIsoSF(pt,eta,era,mode);    
  }
  if(is_ele){
    sf=1.0; //Isolation is inbuilt in electron ID; so no separate isolation SF 
  }
  
  return sf;
}
  
