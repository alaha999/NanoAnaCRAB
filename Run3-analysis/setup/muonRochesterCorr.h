#include "RoccoR/RoccoR.cc"
#include "TRandom3.h"

RoccoR rc2018("setup/RoccoR/RoccoR2018UL.txt");
RoccoR rc2017("setup/RoccoR/RoccoR2017UL.txt");
RoccoR rc2016preVFP("setup/RoccoR/RoccoR2016aUL.txt");
RoccoR rc2016postVFP("setup/RoccoR/RoccoR2016bUL.txt");


float VLLAna:: muon_RochesterSF(int charge,float pt, float phi,float eta, bool isMC, bool genmatched, float genpt, int nl, string era,string mode){
  
  //SF for momentum of each muon

  float SF=1.0;
  float SFerror=1.0;
  int s,m;

  TRandom3 random_gen;
  float u=random_gen.Rndm();
  
  //2018
  if(era=="2018"){
    if(!isMC){
      SF=rc2018.kScaleDT(charge,pt,eta,phi,s=0,m=0);
      SFerror=0;
    }
    else{
      if(genmatched){
	SF=rc2018.kSpreadMC(charge,pt,eta,phi,genpt,s=0,m=0);
	SFerror=rc2018.kSpreadMCerror(charge,pt,eta,phi,genpt);
      }
      else{
	SF=rc2018.kSmearMC(charge,pt,eta,phi,nl,u,s=0,m=0);
	SFerror=rc2018.kSmearMCerror(charge,pt,eta,phi,nl,u);
      }
    }
  }
  
  //2017
  else if(era=="2017"){
    if(!isMC){
      SF=rc2017.kScaleDT(charge,pt,eta,phi,s=0,m=0);
      SFerror=0;
    }
    else{
      if(genmatched){
	SF=rc2017.kSpreadMC(charge,pt,eta,phi,genpt,s=0,m=0);
	SFerror=rc2017.kSpreadMCerror(charge,pt,eta,phi,genpt);
      }
      else{
	SF=rc2017.kSmearMC(charge,pt,eta,phi,nl,u,s=0,m=0);
	SFerror=rc2017.kSmearMCerror(charge,pt,eta,phi,nl,u);
      }
    }
  }
  
  //2016preVFP
  else if(era=="2016preVFP"){
    if(!isMC){
      SF=rc2016preVFP.kScaleDT(charge,pt,eta,phi,s=0,m=0);
      SFerror=0;
    }
    else{
      if(genmatched){
	SF=rc2016preVFP.kSpreadMC(charge,pt,eta,phi,genpt,s=0,m=0);
	SFerror=rc2016preVFP.kSpreadMCerror(charge,pt,eta,phi,genpt);
      }
      else{
	SF=rc2016preVFP.kSmearMC(charge,pt,eta,phi,nl,u,s=0,m=0);
	SFerror=rc2016preVFP.kSmearMCerror(charge,pt,eta,phi,nl,u);
      }
    }
  }
  //
  else if(era=="2016post"){
    if(!isMC){
      SF=rc2016postVFP.kScaleDT(charge,pt,eta,phi,s=0,m=0);
      SFerror=0;
    }
    else{
      if(genmatched){
	SF=rc2016postVFP.kSpreadMC(charge,pt,eta,phi,genpt,s=0,m=0);
	SFerror=rc2016postVFP.kSpreadMCerror(charge,pt,eta,phi,genpt);
      }
      else{
	SF=rc2016postVFP.kSmearMC(charge,pt,eta,phi,nl,u,s=0,m=0);
	SFerror=rc2016postVFP.kSmearMCerror(charge,pt,eta,phi,nl,u);
      }
    }
  }
  
  float sfvalue=1.0;
  if     (mode=="nom" )sfvalue=SF;
  else if(mode=="up"  )sfvalue=SF+SFerror;
  else if(mode=="down")sfvalue=SF-SFerror;

  return sfvalue;
}

