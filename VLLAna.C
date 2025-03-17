#define VLLAna_cxx
#include "VLLAna.h"
#include <TH2.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
using namespace std;


void VLLAna::Begin(TTree * /*tree*/)
{
  //cout<<"Inside Begin()"<<endl;
  TString option = GetOption();

}

void VLLAna::SlaveBegin(TTree *tree /*tree*/)
{
  time(&start);
  
  //cout<<"Inside SlaveBegin()"<<endl;
  TString option = GetOption();
  nEvtTotal = 0;
  nEvtRan = 0;
  nEvtTrigger=0;
  nEvtSkim=0;
  genEventsumw=0;
  //Initialize event counter
  n_4Lskim   = 0;
  n_3Lskim   = 0;
  n_2Lskim   = 0;
  n_1Lskim   = 0;
  
  cout<<"HstFileName: "<<_HstFileName<<endl;
  //Histograms
  //Create the histogram file
  _HstFile = new TFile(_HstFileName,"recreate");
  BookHistograms();
  
  //SkimFile
  //_skimFile = new TFile(_SkimFileName,"recreate");
  //skimTree = tree->CloneTree(0);
  
  //Deactivate Branches
  //tree->SetBranchStatus("*",0);
  //Activate Branches
  //ActivateBranch(tree);   //Look at skimmerHelper.h file
  //ReadBranch();
  
}


void VLLAna::SlaveTerminate()
{
  _HstFile->cd();
  _HstFile->Write();
  _HstFile->Close();

  //create a TH1F histogram and write in skimroot file                                                                                                        
  TH1F *hCount = new TH1F("hCount", "hCount;;",5, 0.5, 5.5);
  hCount->SetBinContent(1,genEventsumw);
  hCount->SetBinContent(2,nEvtTotal);
  hCount->SetBinContent(3,nEvtRan);
  hCount->SetBinContent(4,nEvtTrigger);
  hCount->SetBinContent(5,nEvtSkim);
  hCount->GetXaxis()->SetBinLabel(1,"genEventSumW");
  hCount->GetXaxis()->SetBinLabel(2,"nEvtGen");
  hCount->GetXaxis()->SetBinLabel(3,"nEvtRan");
  hCount->GetXaxis()->SetBinLabel(4,"nEvtTrigger");
  hCount->GetXaxis()->SetBinLabel(5,"nEvtSkim");

  
  _skimFile->cd();
  skimTree->Write();
  hCount->Write();
  _skimFile->Close();
  cout<<"   Done!  "<<endl;
  
  //Output to screen
  cout<<"Total events           = "<<nEvtTotal<<endl;
  cout<<"Total events ran       = "<<nEvtRan<<endl;
  cout<<"Total Triggered events = "<< nEvtTrigger <<endl;
  cout<<"Total events skimmed   = "<<nEvtSkim<<endl;
  
  cout<<"-----------------------------------------------"<<endl;
  cout<<"               Event Selection                 "<<endl;
  cout<<"-----------------------------------------------"<<endl;
  cout<<">>>Skim4L     = "<<n_4Lskim<<endl;
  cout<<">>>Skim3L     = "<<n_3Lskim<<endl;
  cout<<">>>Skim2L     = "<<n_2Lskim<<endl;
  cout<<">>>Skim1L     = "<<n_1Lskim<<endl;
  cout<<"__________________________________________________"<<endl;
  

  //Open the text output file
  ofstream fout(_SumFileName);
  //Put text output in the summary file.
  fout<<"Total events ran      = "<<nEvtRan<<endl;
  fout<<"Total events          = "<<nEvtTotal<<endl;
  fout<<"Total events skimmed  = "<<nEvtSkim<<endl;
  
  time(&end);

  double time_taken = double(end-start);
  cout<<"Time taken by the programe is= "<<fixed<<time_taken<<setprecision(5);
  cout<<"sec"<<endl;

  cout<<"HST FILE  ="<<_HstFileName<<endl;
  cout<<"SKIM FILE ="<<_SkimFileName<<endl;
  
}


void VLLAna::Terminate()
{
  //   cout<<"Inside Terminate()"<<endl;
}

Bool_t VLLAna::Process(Long64_t entry)
{
  fReader.SetLocalEntry(entry);
  if(_data == 0)
    fReader_MC  .SetLocalEntry(entry);
  if(_data == 1)
    fReader_Data.SetLocalEntry(entry);
  if(_year==2016)
    fReader_2016.SetLocalEntry(entry);
  if(_year==2017)
    fReader_2017.SetLocalEntry(entry);
  if(_year==2018)
    fReader_2018.SetLocalEntry(entry);
  if(_year==2018 || _year==2017)
    fReader_1718.SetLocalEntry(entry);
  if( ((_sample < 9) || (_sample>19)) && _data==0)
    fReader_special.SetLocalEntry(entry);
  

  //Verbosity determines the number of processed events after which the root prompt is supposed to display a status update.
  if(_verbosity==0 && nEvtTotal%10000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;
  else if(_verbosity>0 && nEvtTotal%10000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;
  
  nEvtTotal++;
  h.nevt->Fill(0);
  
  if(_data==0)genEventsumw+= *Generator_weight;
  
  
  //METFilter
  //Add TWIKI:
  h.metfilter[0]->Fill(*Flag_goodVertices);
  h.metfilter[1]->Fill(*Flag_globalSuperTightHalo2016Filter);
  h.metfilter[2]->Fill(*Flag_HBHENoiseFilter);
  h.metfilter[3]->Fill(*Flag_HBHENoiseIsoFilter);
  h.metfilter[4]->Fill(*Flag_EcalDeadCellTriggerPrimitiveFilter);
  h.metfilter[5]->Fill(*Flag_BadPFMuonFilter);  
  h.metfilter[6]->Fill(*Flag_eeBadScFilter);
  h.metfilter[7]->Fill(*Flag_ecalBadCalibFilter);
  h.metfilter[8]->Fill(*Flag_METFilters);
  
  GoodEvt2018 = (_year==2018 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && *Flag_BadPFMuonDzFilter && *Flag_ecalBadCalibFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2017 = (_year==2017 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && *Flag_BadPFMuonDzFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2016 = (_year==2016 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  
  h.metfilter[7]->Fill(GoodEvt2016);
  h.metfilter[8]->Fill(GoodEvt2017);
  h.metfilter[9]->Fill(GoodEvt2018);
  
  GoodEvt = GoodEvt2018 && GoodEvt2017 && GoodEvt2016;
  
  h.metfilter[10]->Fill(GoodEvt);
  
  if(GoodEvt){
    
    nEvtRan++; //only good events
    h.nevt->Fill(1);
    
    //Triggers for Data
    if(_data==1){
      trigger2018 = (_year==2018 ? (_lep==1 ? *HLT_IsoMu24==1 : _lep==0 && (*HLT_Ele32_WPTight_Gsf)) : 1);
      //trigger2017 = (_year==2017 ? (_lep==1 ? *HLT_IsoMu27==1 : _lep==0 && (*HLT_Ele32_WPTight_Gsf||*HLT_Ele32_WPTight_Gsf_L1DoubleEG)) : 1);
      trigger2017 = (_year==2017 ? (_lep==1 ? *HLT_IsoMu27==1 : _lep==0 && *HLT_Ele32_WPTight_Gsf_L1DoubleEG) : 1);
      trigger2016 = (_year==2016 ? (_lep==1 ? (*HLT_IsoMu24==1 || *HLT_IsoTkMu24==1) : _lep==0 && *HLT_Ele27_WPTight_Gsf) : 1);
      
      triggerRes = trigger2018 && trigger2017 && trigger2016;
      
    }
    
    //Triggers for MC
    if(_data==0){
      
      trigger2018 = true;trigger2017 = true;trigger2016 = true;
      
      if(_year==2018)trigger2018 = (*HLT_IsoMu24 ==1                     ) || (*HLT_Ele32_WPTight_Gsf == 1            ); //pass mu OR ele trigger
      if(_year==2017)trigger2017 = (*HLT_IsoMu27 ==1                     ) || (*HLT_Ele32_WPTight_Gsf_L1DoubleEG == 1 ); //pass mu OR ele trigger
      if(_year==2016)trigger2016 = (*HLT_IsoMu24 ==1 ||*HLT_IsoTkMu24==1 ) || (*HLT_Ele27_WPTight_Gsf == 1            ); //pass mu OR ele trigger
      
      triggerRes = trigger2018 && trigger2017 && trigger2016;
    }
    
    triggerRes=true; //Set this flag True if you don't apply trigger
    if(triggerRes){
      nEvtTrigger++; //only triggered events
      h.nevt->Fill(2);
      
      
      //-----------------------------------------------------------------------------------------------------------------------
      //                                              GEN OBJECTS BLOCK BEGINS                                                |
      //-----------------------------------------------------------------------------------------------------------------------
      
      // Filled up with code when Gen studies are needed
      
      //-----------------------------------------------------------------------------------------------------------------------
      //                                              RECO OBJECTS BLOCK BEGINS                                               |
      //-----------------------------------------------------------------------------------------------------------------------
      
      //clear array from previous event
      llep.clear();Muon.clear();Electron.clear();
      loosemuon.clear();LooseLep.clear();
      BTaggedJet.clear();jets.clear();taus.clear();
      
      //Per Event Count Object
      int nmu=0, nel=0,ntau=0,njet=0,nbjet=0;
      
      //RecoMuon
      h.nlep[0]->Fill(*nMuon);
      for(unsigned int i=0; i< (*nMuon); i++){
	Lepton temp; temp.v.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],0.105); 
	temp.id = -13*Muon_charge[i]; temp.ind = i;  temp.charge = Muon_charge[i];
	temp.sip3d=Muon_sip3d[i];temp.deepjet=Jet_btagDeepFlavB[Muon_jetIdx[i]];
	h.ptlep[0]->Fill(temp.v.Pt());
	
	//Selections
	bool passcut= temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4 && Muon_mediumId[i]&& Muon_pfRelIso04_all[i]<1.0;
	bool is_promptmuon= fabs(Muon_dxy[i])<0.05 && fabs(Muon_dz[i])<0.1;	
	
	bool analysis_loosemuon= passcut && is_promptmuon && LeptonCustomIDCut(_year,temp.sip3d,temp.deepjet);
	bool analysis_tightmuon= analysis_loosemuon && Muon_pfRelIso04_all[i]<0.15;	
	
	//Define Loose Muon
	if(analysis_loosemuon)LooseLep.push_back(temp),loosemuon.push_back(temp);
	
	//Define Analysis Level Muon
	if(analysis_loosemuon){
	  nmu++;
	  h.ptlep[1]->Fill(temp.v.Pt());		
	  Muon.push_back(temp);
	  llep.push_back(temp);
	}
      }
      
      //RecoElectron
      h.nlep[2]->Fill(*nElectron);
      for(unsigned int i=0; i< (*nElectron); i++){
	Lepton temp; temp.v.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],0.000511); 
	temp.id = -11*Electron_charge[i]; temp.ind = i; temp.charge = Electron_charge[i];
	temp.muoncleaning=MuonCleaning(temp.v,0);//Cleaning against selected muon: 0=loosemuon,1=analysis muon
	temp.sip3d=Electron_sip3d[i];temp.deepjet=Jet_btagDeepFlavB[Electron_jetIdx[i]];
	h.ptlep[2]->Fill(temp.v.Pt());
	
	bool isprompt = false;
	if(fabs(temp.v.Eta())<=1.479){
	  if(fabs(Electron_dxy[i])<0.05 && fabs(Electron_dz[i])<0.1)
	    isprompt = true;
	}
	if(fabs(temp.v.Eta())>1.479){
	  if(fabs(Electron_dxy[i])<0.1 && fabs(Electron_dz[i])<0.2)
	    isprompt = true;
	}

	bool eleMediumIDNoIso  = electronCustomID(Electron_vidNestedWPBitmap[i],3,7);
	bool eleMediumID       = electronCustomID(Electron_vidNestedWPBitmap[i],3,-1);
	bool looseIsolation    = Electron_pfRelIso03_chg[i]<1.0;	
	bool passcut           = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4 && eleMediumIDNoIso && looseIsolation && temp.muoncleaning;
	bool analysis_looseele = passcut && isprompt && LeptonCustomIDCut(_year,temp.sip3d,temp.deepjet); //loose iso: <1.0
	bool analysis_tightele = analysis_looseele && eleMediumID;//Pt,Eta,dxy,dz,MediumIDWithIso
	
	//Define Loose Electron array
	if(analysis_looseele)LooseLep.push_back(temp);
	
	//Define Analysis Level Electrons
	if(analysis_looseele){
	  nel++;
	  h.ptlep[3]->Fill(temp.v.Pt());
	  Electron.push_back(temp);
	  llep.push_back(temp);
	}
      }
      
      //RecoTau
      h.nlep[4]->Fill(*nTau);
      for(unsigned int i=0; i< (*nTau); i++){
	if(Tau_decayMode[i]&&(Tau_decayMode[i]<3||Tau_decayMode[i]>9)){
	  //Tau energy scale correction
	  float tlv_corr = 1.;
	  if(_year==2016){
	    if(Tau_decayMode[i]==0)
	      tlv_corr = 0.994;
	    if(Tau_decayMode[i]==1)
	      tlv_corr = 0.995;
	    if(Tau_decayMode[i]>9)
	      tlv_corr = 1;
	  }
	  if(_year==2017){
	    if(Tau_decayMode[i]==0)
	      tlv_corr = 1.007;
	    if(Tau_decayMode[i]==1)
	      tlv_corr = 0.998;
	    if(Tau_decayMode[i]==10)
	      tlv_corr = 1.001;
	    if(Tau_decayMode[i]==11)
	      tlv_corr = 0.999;
	  }
	  if(_year==2018){
	    if(Tau_decayMode[i]==0)
	      tlv_corr = 0.987;
	    if(Tau_decayMode[i]==1)
	      tlv_corr = 0.995;
	    if(Tau_decayMode[i]==10)
	      tlv_corr = 0.998;
	    if(Tau_decayMode[i]==11)
	      tlv_corr = 1;
	  }
	  Lepton temp; temp.v.SetPtEtaPhiM(Tau_pt[i],Tau_eta[i],Tau_phi[i],1.77);
	  h.ptlep[4]->Fill(temp.v.Pt());
	  temp.v *= tlv_corr; //energy correction
	  temp.id = -15*Tau_charge[i]; temp.ind = i; temp.charge = Tau_charge[i];
	  temp.lepcleaning = TaulepCleaning(temp.v); //haven't found any lepton in the direction of Tau
	  h.ptlep[5]->Fill(temp.v.Pt());
	  
	  bool passcut=temp.v.Pt()>20 && fabs(temp.v.Eta())<2.3;
	  passcut = passcut && temp.lepcleaning && fabs(Tau_dz[i])<0.2;
	  bool DeepTauID= Tau_idDeepTau2017v2p1VSe[i]>15 && Tau_idDeepTau2017v2p1VSmu[i]>3 && Tau_idDeepTau2017v2p1VSjet[i]>127;
	  //VL AntiEle, VL AntiMu, VT AntiJet
	  
	  bool analysisCut= passcut && DeepTauID;
	  
	  if(analysisCut){
	    ntau++;
	    taus.push_back(temp);
	    h.ptlep[6]->Fill(temp.v.Pt());
	  }
	}
      }
      
      //Jets and b-tagged Jets (DeepJet Medium WP)
      for(unsigned int i=0; i< (*nJet); i++){
	Lepton temp; temp.v.SetPtEtaPhiM(Jet_pt[i],Jet_eta[i],Jet_phi[i],Jet_mass[i]);
	temp.ind = i;
	temp.taucleaning=TaujetCleaning(temp.v);//No jet in the direction of tau
	temp.lepcleaning=LepjetCleaning(temp.v);
	if(_data==0)temp.hadronflavor = Jet_hadronFlavour[i];
	else temp.hadronflavor = 999;
	//h.ptjet[0]->Fill(temp.v.Pt());
	bool passcut =temp.v.Pt()>30 && fabs(temp.v.Eta())<2.4 && temp.lepcleaning; //PF Tight Jet ID: bit 2
	if(passcut){
	  if(_year == 2016 ? Jet_jetId[i]>=1 : Jet_jetId[i]>=2){
	    jets.push_back(temp);
	    njet++;
	    h.ptlep[7]->Fill(temp.v.Pt());
	    if(_year ==2018 ? Jet_btagDeepFlavB[i]>=0.2783: (_year == 2017) ? Jet_btagDeepFlavB[i]>=0.3040 :(_year==2016 && _era=="preVFP") ? Jet_btagDeepFlavB[i]>=0.2598:(_year==2016 && _era=="postVFP") ? Jet_btagDeepFlavB[i]>=0.2489:1.0)
	      //B tag discriminator medium working point(float variable)//2016preVFP: 0.2598 and 2016postVFP:2489
	      {
		BTaggedJet.push_back(temp); nbjet++; h.ptlep[8]->Fill(temp.v.Pt());
	      }
	  }
	}
      }
      
      //Fill Histograms of object
      h.nlep[1]->Fill(nmu);
      h.nlep[3]->Fill(nel);
      h.nlep[5]->Fill(ntau);
      h.nlep[6]->Fill(njet);
      h.nlep[7]->Fill(nbjet);
      
      //Sort objects based on pt
      Sort(Muon);
      Sort(Electron);
      Sort(llep);
      Sort(taus);
      Sort(jets);
      Sort(BTaggedJet);
      
      //MET
      metpt = *MET_pt;
      metphi = *MET_phi;
      
      h.etmiss->Fill(metpt);//fill a histogram with the missing Et of the event.
      
      //  cout<<"jet & met stored"<<endl;

      // Object Selection is completed with different set of criterias on object
      //-----------------------------------------------------------------------------------------------------------------------
      //                                              RECO OBJECTS BLOCK ENDS                                                 |
      //-----------------------------------------------------------------------------------------------------------------------
      
      


      //-----------------------------------------------------------------------------------------------------------------------
      //                                        EVENT SELECTION                                                               |
      //-----------------------------------------------------------------------------------------------------------------------
      /*
	>>>> Let's do inclusive skim for 4L, 3L, 2L, and 1L.
	>> 4L: atleast 4 light leptons
	>> 3L: atleast 3 light leptons
	>> 2L: atleast 2 light leptons
	>> 1L: atleast 1 light leptons
	
	>>> Current Criteria
	>> Atleast 1 triggerable light leptons
	>> 
	>> Analysis Level LOOSE Muons, LOOSE Electrons
      */
      
      
      //Lepton pt per year
      float trigmupt =0;float trigelept =0;
      if(_year==2016)trigmupt=26,trigelept=30;
      if(_year==2017)trigmupt=29,trigelept=35;
      if(_year==2018)trigmupt=26,trigelept=35;
      
      bool is_4Lskim   = false; //4L inclusive
      bool is_3Lskim   = false; //3L inclusive
      bool is_2Lskim   = false; //2L inclusive
      bool is_1Lskim   = false; //1L inclusive
      bool is_1L2Jskim = false; //1L2J skim

      //4L inclusive skim
      if((int)llep.size()>3){ //4L
	bool passTrigger = false;
	if(abs(llep.at(0).id)==13)passTrigger = llep.at(0).v.Pt()>trigmupt;
	else if(abs(llep.at(0).id)==11)passTrigger = llep.at(0).v.Pt()>trigelept;
	
	if(passTrigger)is_4Lskim=true,n_4Lskim++;
      }      

      //3L inclusive skim
      if((int)llep.size()>2){ //3L
	bool passTrigger = false;
	if(abs(llep.at(0).id)==13)passTrigger = llep.at(0).v.Pt()>trigmupt;
	else if(abs(llep.at(0).id)==11)passTrigger = llep.at(0).v.Pt()>trigelept;
	
	if(passTrigger)is_3Lskim=true,n_3Lskim++;
      }

      //2L inclusive skim
      if((int)llep.size()>1){ //2L
	bool passTrigger = false;
	if(abs(llep.at(0).id)==13)passTrigger = llep.at(0).v.Pt()>trigmupt;
	else if(abs(llep.at(0).id)==11)passTrigger = llep.at(0).v.Pt()>trigelept;
	
	if(passTrigger)is_2Lskim=true,n_2Lskim++;
      }
      
      //1L inclusive skim
      if((int)llep.size()>0){ //1L
	bool passTrigger = false;
	if(abs(llep.at(0).id)==13)passTrigger = llep.at(0).v.Pt()>trigmupt;
	else if(abs(llep.at(0).id)==11)passTrigger = llep.at(0).v.Pt()>trigelept;
	
	if(passTrigger)is_1Lskim=true,n_1Lskim++;
      }
      
      //1L2J skim
      is_1L2Jskim = is_1Lskim && (int)jets.size()>1;
      
      //-----------------------------------------------------------------------------------------------------------------------
      //                                        SKIM TTREE                                                                    |
      //-----------------------------------------------------------------------------------------------------------------------
      
      bool keepThisEvent=false;
      if(is_1L2Jskim)keepThisEvent=true; //change it to is_3Lskim, is_2Lskim etc [#1L2J skim]
      
      if(keepThisEvent){
	nEvtSkim++;
	
	//This function activates all the branches needed downstream. In TTreeReader
	//if you don't call the branch, in skimTree those branch will contain empty values
	//leading to segmentation fault or unphysical answer.
	ReadBranch();    
	
	//Fill the tree
	skimTree->Fill();
      }
      
      //=====================================================DON'T TOUCH THESE BRACKET ENDS=========================================//
    }//end of trigger block
  }// end of goodEvt block
  return kTRUE;
}// end of processLoop


//============================================================================================================================================//
//   BELOW THIS LINE ARE USEFUL FUNCTIONS THAT WE CALL IN THE PROCESS LOOP
//============================================================================================================================================//

//Sort Lepton Vector in Pt
void VLLAna::Sort(vector<Lepton>lepvec)
{
  
  for(int i=0; i<(int)lepvec.size()-1; i++){
    for(int j=i+1; j<(int)lepvec.size(); j++){
      if(lepvec[i].v.Pt() < lepvec[j].v.Pt() ) swap(lepvec.at(i),lepvec.at(j));
    }
  }
}

//Debugging MC Function
void VLLAna::TruthLevelInfo()
{
  if(nEvtTotal>20000 && nEvtTotal<20010){
    for(unsigned int i=0; i< (*nGenPart); i++){
      Lepton temp; temp.v.SetPtEtaPhiM(GenPart_pt[i],GenPart_eta[i],GenPart_phi[i],GenPart_mass[i]);temp.status = GenPart_status[i]; temp.ind = i;temp.pdgid = GenPart_pdgId[i];temp.momid=MotherID(i,GenPart_genPartIdxMother[i]);
      if(abs(temp.pdgid)<100 && abs(temp.momid)<100 && temp.status==1){
	cout<<i<<"|"<<temp.pdgid<<"|"<<temp.momid<<"|"<<endl;
	cout<<"===========Event="<<nEvtTotal<<"==========="<<endl;
      }
    }
  }
}


//To Access bits from Bit Variables
int VLLAna::getBit ( int bitvar, int nthbit ){
  int bittrue = 0;
  bittrue = ((bitvar&(1<<nthbit))>>nthbit);
  return bittrue;
}

//Electron ID from BitMap
//use to produce Electron Medium ID without Isolation cuts
int VLLAna::electronCustomID(Int_t bitmap,int quality, int skipCut){  
  int qualityPass=1;
  for(int i=0;i<10;i++){
    if(i==skipCut)continue;
    if(((bitmap>>i*3) & 0x7)<quality){qualityPass=0;break;}
  }
  
  return qualityPass;
}

bool VLLAna::LeptonCustomIDCut(int year, float sip3d, float deepjet){

  bool passCut=false;
  
  if      (year==2018 && sip3d<9 && deepjet<0.3)passCut=true;
  else if (year==2017 && sip3d<12 && deepjet<0.4)passCut=true;
  else if (year==2016 && sip3d<10 && deepjet<0.6)passCut=true;
  //else if (year=="2016postVFP" && sip3d<10 && deepjet<0.6)passCut=true;
  
  return passCut;
}

bool VLLAna::TriggerObjectMatching(int id,TLorentzVector t){
  bool result=false;float dr=999.0,drmin=999.0; float trigobjpt=0;
  int n_HLTid13=0;int n_matchedHLT=0;int n_HLTid13_drmatched=0;
  int trigobjindex=-1;
  for(unsigned int i=0; i< (*nTrigObj); i++){        
    if(abs(id)==abs(TrigObj_id[i])){
      float deta= abs(t.Eta()-TrigObj_eta[i]);
      float dphi= delta_phi(t.Phi(),TrigObj_phi[i]);
      dr = sqrt(pow(deta,2)+pow(dphi,2));
      
      //counting how many id==13 trigobj
      n_HLTid13++;
      h.trigobj[2]->Fill(abs(t.Pt()-TrigObj_pt[i]));
      h.trigobj[3]->Fill(dr);
      h.trigobj2d[0]->Fill(t.Pt(),TrigObj_pt[i]);
      //How many trig 13 objects are close by?
      if(dr<0.2){
	n_HLTid13_drmatched++;
	h.trigobj[4]->Fill(abs(t.Pt()-TrigObj_pt[i]));
	h.trigobj[5]->Fill(dr);
	h.trigobj2d[1]->Fill(t.Pt(),TrigObj_pt[i]);
      }            
      //cout<<"dR(muon,trigobj) and Pt="<<dr<<","<<TrigObj_pt[i]<<endl;
      
      //Find the minimum dr
      if(dr<drmin){drmin=dr,trigobjindex=i,trigobjpt=TrigObj_pt[i],n_matchedHLT++;	
      }
    }
  }
  
  h.trigobj[0]->Fill(n_HLTid13);
  h.trigobj[1]->Fill(n_HLTid13_drmatched);

  float pt_thres = 0;int firebit =0;
  if(abs(id)==13 && (_year==2018 || _year==2016))pt_thres=24.0,firebit = getBit(TrigObj_filterBits[trigobjindex],3); //ISOMu24 Trigger(2016,2018)
  else if(abs(id)==11 && (_year==2018 || _year==2017))pt_thres=32.0,firebit = getBit(TrigObj_filterBits[trigobjindex],1);//Ele32Tight Trigger(2017,2018)
  else if(abs(id)==13 && _year==2017) pt_thres=27.0,firebit = getBit(TrigObj_filterBits[trigobjindex],3); //ISOMu27 Trigger(2017)
  else if(abs(id)==11 && _year==2016) pt_thres=27.0,firebit = getBit(TrigObj_filterBits[trigobjindex],1); //Ele27Tight Trigger(2016)
  
  if(drmin<0.2 && (trigobjpt > pt_thres)){
    result =true;
    //cout<<"drMin="<<drmin<<endl;
    //cout<<"event no, trigobjpt and lepton pt="<<nEvtTotal<<","<<trigobjpt<<","<<t.Pt()<<endl;
    h.trigobj[6]->Fill((trigobjpt-t.Pt())/t.Pt());
    h.trigobj[7]->Fill(drmin);
    h.trigobj2d[2]->Fill(t.Pt(),trigobjpt);
  }
  return result;
}


void VLLAna::BookHistograms()
{
  //  cout<<"Inside BookHist()"<<endl;
  h.nevt = new TH1F("nEvents","0-Total events, 1-Total events ran, 2-Total events with trigger applied",5,-1,4);
  h.metfilter[0] = new TH1F("METfilter_goodVertices","METfilter_goodVertices",5,-1,4);
  h.metfilter[1] = new TH1F("METfilter_globalSuperTightHalo2016Filter","METfilter_globalSuperTightHalo2016Filter",5,-1,4);
  h.metfilter[2] = new TH1F("METfilter_HBHENoiseFilter","METfilter_HBHENoiseFilter",5,-1,4);
  h.metfilter[3] = new TH1F("METfilter_HBHENoiseIsoFilter","METfilter_HBHENoiseIsoFilter",5,-1,4);
  h.metfilter[4] = new TH1F("METfilter_EcalDeadCellTriggerPrimitiveFilter","METfilter_EcalDeadCellTriggerPrimitiveFilter",5,-1,4);
  h.metfilter[5] = new TH1F("METfilter_BadPFMuonFilter","METfilter_BadPFMuonFilter",5,-1,4);
  h.metfilter[6] = new TH1F("METfilter_eeBadScFilter","METfilter_eeBadScFilter",5,-1,4);
  h.metfilter[7] = new TH1F("METfilter_GoodEvt2016","METfilter_GoodEvt2016",5,-1,4);
  h.metfilter[8] = new TH1F("METfilter_GoodEvt2017","METfilter_GoodEvt2017",5,-1,4);
  h.metfilter[9] = new TH1F("METfilter_GoodEvt2018","METfilter_GoodEvt2018",5,-1,4);
  h.metfilter[10] = new TH1F("METfilter_GoodEvt","METfilter_GoodEvt",5,-1,4);

  //nTrigObj Histograms
  h.trigobj[0] = new TH1F("ntrigobj_id13","No of trig obj with abs(id)=13",10,0,10);
  h.trigobj[1] = new TH1F("ntrigobj_id13dr","No of trig obj with abs(id)=13 and within dr<0.2",10,0,10);
  h.trigobj[2] = new TH1F("dpt_alltrigobj","dPt between all trig obj(id=13) and lepton",100,0,10);
  h.trigobj[3] = new TH1F("dr_alltrigobj","dR between all trig obj(id=13) and lepton",100,0,10);
  h.trigobj[4] = new TH1F("dpt_allmatchedtrigobj","dPt between all trig obj(id=13) with dR<0.2 and lepton",100,0,10);
  h.trigobj[5] = new TH1F("dr_allmatchedtrigobj","dR between all trig obj(id=13) with dR<0.2 and lepton",100,0,10);
  h.trigobj[6] = new TH1F("resolution_matched","RecoPt-HLTobjPt/RecoPt",4000,-2,2);
  h.trigobj[7] = new TH1F("dr_matchedtrigobj","dR between matched HLT obj(id=13) and lepton",100,0,10);
  for(int i=0; i<8; i++) h.trigobj[i]->Sumw2();
  h.trigobj2d[0] = new TH2F("RecoPt_vs_HLTobjectPt_id13","Reco pT vs HLT object pT",500,0,500,500,0,500);
  h.trigobj2d[1] = new TH2F("RecoPt_vs_HLTobjectPt_allmatched","Reco pT vs HLT object pT all matched",500,0,500,500,0,500);
  h.trigobj2d[2] = new TH2F("RecoPt_vs_HLTobjectPt_matched","Reco pT vs HLT object pT (matched)",500,0,500,500,0,500);
  //Basic plots starts
  //lepton(ele,mu,tau)
  h.nlep[0] = new TH1F("nMu","Number of Muon Candidates",20,0,20);
  h.nlep[1] = new TH1F("ngoodMu","Number of good Muon Candidates",20,0,20);
  h.nlep[2] = new TH1F("nEl","Number of Electron Candidates",20,0,20);
  h.nlep[3] = new TH1F("ngoodEl","Number of good Electron Candidates",20,0,20);
  h.nlep[4] = new TH1F("nTau","Number of Tau Candidates",20,0,20);
  h.nlep[5] = new TH1F("ngoodTau","Number of good Tau Candidates",20,0,20);
  h.nlep[6] = new TH1F("ngoodJets","Number of good jet Candidates",20,0,20);
  h.nlep[7] = new TH1F("nbtaggedjets","Number of b-tagged jet Candidates",20,0,20);
  
  for(int i=0; i<8; i++) h.nlep[i]->Sumw2();
  
  //lepton properties
  h.ptlep[0] = new TH1F("Mu_pt","Muon candidate p_{T}",1000,0,1000);
  h.ptlep[1] = new TH1F("goodMu_pt","Muon p_{T}",1000,0,1000);
  h.ptlep[2] = new TH1F("El_pt","Electron candidate p_{T}",1000,0,1000);
  h.ptlep[3] = new TH1F("goodEl_pt","Electron p_{T}",1000,0,1000);
  h.ptlep[4] = new TH1F("Tau_pt","Tau candidate p_{T}",1000,0,1000);
  h.ptlep[5] = new TH1F("Tau_corrpt","Tau candidate p_{T}",1000,0,1000);
  h.ptlep[6] = new TH1F("goodTau_pt","Tau p_{T}",1000,0,1000);
  h.ptlep[7] = new TH1F("goodJets_pt","good jets p_{T}",1000,0,1000);
  h.ptlep[8] = new TH1F("goodbJets_pt","b-tagged jets p_{T}",1000,0,1000);
  for(int i=0; i<9; i++) h.ptlep[i]->Sumw2();
  
  //missing ET of the event
  h.etmiss= new TH1F("etmiss_beforeMTcut","Missing E_{T}",1000,0,1000); h.etmiss->Sumw2();
  
  //object cutflow
  //h.mucut = new TH1F("mucut","0=pt,1=eta,2=dxydz,3=loose,4=Medium,5=Iso:0.25",10,0,10);h.mucut->Sumw2();
  //h.elcut = new TH1F("elcut","0=pt,1=eta,2=dxydz,3=loose,4=Medium,5=muoncleaning",10,0,10);h.elcut->Sumw2();
  //h.taucut = new TH1F("taucut","0=All,1=NewDM,2=DeepCSV,3=AntiEle,4=AntiMu,5=Prompt,6=LepClean,7=Eta,8=PT",10,0,10);h.taucut->Sumw2();
  //h.phocut = new TH1F("phocut","phocuts",10,0,10);
  
  //evtwt
  //h.evtwt[0] = new TH1F("evtwt"        ,"EventWeight",100,0,1);
  //h.evtwt[1] = new TH1F("evtwt_mu"        ,"EventWeight after mumumu selection",100,0,1);
  //h.evtwt[2] = new TH1F("triggereff"        ,"Trigger Efficiency",100,0,1);
  //h.evtwt[3] = new TH1F("totwt"        ,"Trigger Eff x evtwt",100,0,1);
  //for(int i=0;i<4;i++) h.evtwt[i]->Sumw2();

  //Channel 1L2J
  // h.channel1L2J[0] = new TH1F("channel1L2J_Count","channel 1L2J count: 1=incl,2=mu,3=ele",10,0,10);
  //for(int i=0; i<1; i++) h.channel1L2J[i]->Sumw2();
  
  
}
//end of BOOK HISTOS


//=========================================================================================================================================//
//Custom Functions
//===================
float VLLAna::delta_phi(float phi1, float phi2)
{
  //Calculate the correct deltaPhi=phi1-phi2
  phi1 = TVector2::Phi_0_2pi(phi1);
  phi2 = TVector2::Phi_0_2pi(phi2);
  float dphi = fabs(phi1 - phi2);
  if(dphi>TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
  return dphi;
}

int VLLAna::MotherID(int partindex, int momindex)
{
  int parid =GenPart_pdgId[partindex];
  int momid = GenPart_pdgId[momindex];
  while(parid==momid){
    partindex=momindex;
    momindex=GenPart_genPartIdxMother[momindex];
    parid =GenPart_pdgId[partindex];
    momid = GenPart_pdgId[momindex];
  }
  return momid;
}

//DISTINCT OSSFN PAIR//
int VLLAna::getOSSFN()
{
  if((int)llep.size()==3){
    Lepton lep0=llep.at(0), lep1=llep.at(1),lep2=llep.at(2); 
    {
      if(abs(lep0.id)==abs(lep1.id) && lep0.charge != lep1.charge ) return 1;  
      if(abs(lep0.id)==abs(lep2.id) && lep0.charge != lep2.charge) return 1;
      if(abs(lep1.id)==abs(lep2.id) && lep1.charge != lep2.charge) return 1;
      else return 0;
    }
  }
  else return 0;
}

float VLLAna::getInvMass(TLorentzVector a,TLorentzVector b)
{
  return (a+b).M();
}
//is OSSF?
bool VLLAna::isOSSF(int i, int j)
{
  bool result= false;
  Lepton lep_i=llep.at(i),lep_j=llep.at(j);
  if(abs(lep_i.id)==abs(lep_j.id) && lep_i.charge != lep_j.charge ){
    result=true;
  }
  return result;
}

float VLLAna::OSSFMass(int i, int j)
{
  float mass=-1.0;
  Lepton lep_i=llep.at(i),lep_j=llep.at(j);
  if(abs(lep_i.id)==abs(lep_j.id) && lep_i.charge != lep_j.charge ){
    mass=getInvMass(lep_i.v,lep_j.v);
  }
  return mass;
}

//BEST OSSF PAIR
float VLLAna::getMOSSF()
{
  float bestMOSSF=-1.0, ZMass=91.0, dmMin=999.0;
  for(int i=0; i<(int)llep.size(); i++){
    for(int j=0; j<(int)llep.size();j++){
      if(j!=i && j>i){
	if(abs(llep.at(j).id)==abs(llep.at(i).id) && llep.at(j).charge!=llep.at(i).charge){
	  bool onZpair=getInvMass(llep.at(i).v,llep.at(j).v)>76 && getInvMass(llep.at(i).v,llep.at(j).v)<106;
	  if(onZpair){
	    if(abs(ZMass-getInvMass(llep.at(i).v,llep.at(j).v))<=dmMin){
	      bestMOSSF=getInvMass(llep.at(i).v,llep.at(j).v);
	      dmMin=abs(ZMass-getInvMass(llep.at(i).v,llep.at(j).v));
	    }
	  }
	}
      }
    }
  }
  return bestMOSSF;
}

//Non best on Z lepton index//
int VLLAna::NonBestOnZLeptonIndex()
{
  if((int)llep.size()==3){
    float LightLeptonBestMOSSF=getMOSSF();
    Lepton lep0=llep.at(0),lep1=llep.at(1),lep2=llep.at(2);
    float pair01mass=-1.0,pair02mass=-1.0,pair12mass=-1.0;
    if(isOSSF(0,1))  pair01mass=getInvMass(lep0.v,lep1.v);
    if(isOSSF(0,2))  pair02mass=getInvMass(lep0.v,lep2.v);
    if(isOSSF(1,2))  pair12mass=getInvMass(lep1.v,lep2.v);
    
    //NON Z LEPTON INDEX//
    if(pair01mass==LightLeptonBestMOSSF) return 2;
    if(pair02mass==LightLeptonBestMOSSF) return 1;
    if(pair12mass==LightLeptonBestMOSSF) return 0;
    else return 0;
  }
  else return 0;
}

/////////////////////////////////
//  Object Cleaning Functions  //
/////////////////////////////////

//Muon cleaning: Clean an electron against Muon
bool VLLAna::MuonCleaning(TLorentzVector t, int opt)
{
  bool result1=false;
  vector<Lepton> lepton;
  if(opt==0)lepton=loosemuon;
  if(opt==1)lepton=Muon;
  
  if((int)lepton.size()==0)
    result1=true;
  if((int)lepton.size()>0){
    for(int i=0; i<(int)lepton.size(); i++){
      if(t.DeltaR(lepton[i].v)>0.5)
        result1=true;
      else{
        result1=false;
        break;
      }
    }
  }
  return result1;
}


// Clean a tau against light lepton
bool VLLAna::TaulepCleaning(TLorentzVector t)
{
  bool result1=false;
  if((int)LooseLep.size()==0)
    result1=true;
  if((int)LooseLep.size()>0){
    for(int i=0; i<(int)LooseLep.size(); i++){
      if(t.DeltaR(LooseLep[i].v)>0.5)
        result1=true;
      else{
        result1=false;
        break;
      }
    }
  }
  return result1;
}

// Clean a jet against tau
bool VLLAna::TaujetCleaning(TLorentzVector t)
{
  bool result1=false;
  if((int)taus.size()==0)
    result1=true;
  if((int)taus.size()>0){
    for(int i=0; i<(int)taus.size(); i++){
      if(t.DeltaR(taus[i].v)>0.4)
        result1=true;
      else{
        result1=false;
        break;
      }
    }
  }
  return result1;
}

//Clean a jet against light lepton
bool VLLAna::LepjetCleaning(TLorentzVector t)
{
  bool result1=false;
  if((int)LooseLep.size()==0)
    result1=true;
  if((int)LooseLep.size()>0){
    for(int i=0; i<(int)LooseLep.size(); i++){
      if(t.DeltaR(LooseLep[i].v)>0.4)
        result1=true;
      else{
        result1=false;
        break;
      }
    }
  }
  return result1;
}


void VLLAna::ReadBranch(){
  //ONLY READ
  //Muon,Electron,Jets,Taus,MET,PuppiMET,GenJet,GenPart,GenVisTau,nTrigObj,Flags,HLT_Path
  // Rules:
  // TTreeReaderValue branch: *variable
  // TTreeReaderArray branch: variable[i] where i is looping over no of that object
  
  *run;
  *luminosityBlock;
  *event;
  
  //Branches(TTreeValue) #TTreeArray needs looping
  
  //HLT
  *HLT_IsoMu24;
  *HLT_IsoMu27;
  *HLT_Ele27_WPTight_Gsf;
  if(_year==2017 || _year==2018)*HLT_Ele32_WPTight_Gsf_L1DoubleEG;
  if(_year==2018               )*HLT_Ele32_WPTight_Gsf;
  if(_year==2016               )*HLT_IsoTkMu24;
  
  //
  if(_data==0){
    *L1PreFiringWeight_Dn;
    *L1PreFiringWeight_Up;
    *L1PreFiringWeight_Nom;
    *Generator_weight;
    *Pileup_nPU;
    *Pileup_nTrueInt;
    *fixedGridRhoFastjetAll;
  }
  if(_data==0 && ((_sample <9) || (_sample >19))){//These branches are not available in Non-Pythia samples such as QCD 
    *nLHEPdfWeight;
    *nLHEScaleWeight;    
    *LHE_HT;
    *LHE_Vpt;
    for(unsigned int i=0;i<*nLHEPdfWeight;i++)LHEPdfWeight[i];
    for(unsigned int i=0;i<*nLHEScaleWeight;i++)LHEScaleWeight[i];    
  }
  
  //Flags
  *Flag_HBHENoiseFilter;
  *Flag_HBHENoiseIsoFilter;
  *Flag_CSCTightHaloFilter;
  *Flag_CSCTightHaloTrkMuUnvetoFilter;
  *Flag_CSCTightHalo2015Filter;
  *Flag_globalTightHalo2016Filter;
  *Flag_globalSuperTightHalo2016Filter;
  *Flag_HcalStripHaloFilter;
  *Flag_hcalLaserEventFilter;
  *Flag_EcalDeadCellTriggerPrimitiveFilter;
  *Flag_EcalDeadCellBoundaryEnergyFilter;
  *Flag_ecalBadCalibFilter;
  *Flag_goodVertices;
  *Flag_eeBadScFilter;
  *Flag_ecalLaserCorrFilter;
  *Flag_trkPOGFilters;
  *Flag_chargedHadronTrackResolutionFilter;
  *Flag_muonBadTrackFilter;
  *Flag_BadChargedCandidateFilter;
  *Flag_BadPFMuonFilter;
  *Flag_BadPFMuonDzFilter;
  //*Flag_hfNoisyHitsFilter;
  *Flag_BadChargedCandidateSummer16Filter;
  *Flag_BadPFMuonSummer16Filter;
  *Flag_trkPOG_manystripclus53X;
  *Flag_trkPOG_toomanystripclus53X;
  *Flag_trkPOG_logErrorTooManyClusters;
  *Flag_METFilters;
  
  //PFMET
  *MET_MetUnclustEnUpDeltaX;
  *MET_MetUnclustEnUpDeltaY;
  *MET_covXX;
  *MET_covXY;
  *MET_covYY;
  *MET_phi;
  *MET_pt;
  *MET_significance;
  *MET_sumEt;
  *MET_sumPtUnclustered;
  //if(_data==0){
  //  *MET_fiducialGenPhi;
  //  *MET_fiducialGenPt;
  //}
  //PuppiMET
  *PuppiMET_phi;
  *PuppiMET_phiJERDown;
  *PuppiMET_phiJERUp;
  *PuppiMET_phiJESDown;
  *PuppiMET_phiJESUp;
  *PuppiMET_phiUnclusteredDown;
  *PuppiMET_phiUnclusteredUp;
  *PuppiMET_pt;
  *PuppiMET_ptJERDown;
  *PuppiMET_ptJERUp;
  *PuppiMET_ptJESDown;
  *PuppiMET_ptJESUp;
  *PuppiMET_ptUnclusteredDown;
  *PuppiMET_ptUnclusteredUp;
  *PuppiMET_sumEt;
  
  //GenMET
  //if(_data==0){
  //  *GenMET_phi;
  //  *GenMET_pt;
  // }
  
  //Muon Branches
  *nMuon;
  for(unsigned int i=0;i<(*nMuon);i++){
    Muon_dxy[i];
    Muon_dz[i];
    Muon_eta[i];
    Muon_jetPtRelv2[i];
    Muon_jetRelIso[i];
    Muon_mass[i];
    Muon_pfRelIso04_all[i];
    Muon_phi[i];
    Muon_pt[i];
    Muon_sip3d[i];
    Muon_charge[i];
    Muon_jetIdx[i];
    Muon_pdgId[i];
    Muon_looseId[i];
    Muon_mediumId[i];
    Muon_tightId[i];
    Muon_nTrackerLayers[i];
    if(_data==0)Muon_genPartIdx[i],Muon_genPartFlav[i];
  }
  
  //Electron
  *nElectron;
  for(unsigned int i=0; i<(*nElectron);i++){
    Electron_dEscaleDown[i];
    Electron_dEscaleUp[i];
    Electron_deltaEtaSC[i];
    Electron_dxy[i];
    Electron_dz[i];
    Electron_eta[i];
    Electron_jetPtRelv2[i];
    Electron_jetRelIso[i];
    Electron_mass[i];
    Electron_pfRelIso03_all[i];
    Electron_phi[i];
    Electron_pt[i];
    Electron_sip3d[i];
    Electron_charge[i];
    Electron_cutBased[i];
    Electron_jetIdx[i];
    Electron_pdgId[i];
    Electron_vidNestedWPBitmap[i];
    if(_data==0)Electron_genPartIdx[i],Electron_genPartFlav[i];
  }
  
  //Taus
  *nTau;
  for(unsigned int i=0; i<(*nTau);i++){
    Tau_dxy[i];
    Tau_dz[i];
    Tau_eta[i];
    Tau_mass[i];
    Tau_phi[i];
    Tau_pt[i];
    Tau_charge[i];
    Tau_decayMode[i];
    Tau_jetIdx[i];
    Tau_idDeepTau2017v2p1VSe[i];
    Tau_idDeepTau2017v2p1VSjet[i];
    Tau_idDeepTau2017v2p1VSmu[i];
    //if(_data==0)Tau_genPartIdx[i],Tau_genPartFlav[i];
  }

  //Jets
  *nJet;
  for(unsigned int i=0; i<(*nJet);i++){
    Jet_btagDeepFlavB[i];
    Jet_btagDeepFlavCvB[i];
    Jet_btagDeepFlavCvL[i];
    Jet_btagDeepFlavQG[i];
    Jet_eta[i];
    Jet_mass[i];
    Jet_phi[i];
    Jet_pt[i];
    Jet_puIdDisc[i];
    Jet_qgl[i];
    Jet_jetId[i];
    Jet_puId[i];
    Jet_nConstituents[i];
    if(_data==0){Jet_genJetIdx[i];Jet_hadronFlavour[i];Jet_partonFlavour[i];}
  }

  //PV
  //Taus
  //*nOtherPV;
  //for(unsigned int i=0; i<(*nOtherPV);i++){
  //}
  *PV_npvs;
  *PV_npvsGood;
  
  //BELOW ARE ONLY FOR MC
  if(_data==0){
    //GenJet
    *nGenJet;
    for(unsigned int i=0; i<(*nGenJet);i++){    
      GenJet_eta[i];
      GenJet_mass[i];
      GenJet_phi[i];
      GenJet_pt[i];
      GenJet_partonFlavour[i];
      GenJet_hadronFlavour[i];
    }
    
    //GenPart
    *nGenPart;
    for(unsigned int i=0; i<(*nGenPart);i++){    
      GenPart_eta[i];
      GenPart_mass[i];
      GenPart_phi[i];
      GenPart_pt[i];
      GenPart_genPartIdxMother[i];
      GenPart_pdgId[i];
      GenPart_status[i];
      GenPart_statusFlags[i];
    }
    //GenVisTau
    *nGenVisTau;
    for(unsigned int i=0; i<(*nGenVisTau);i++){    
      GenVisTau_eta[i];
      GenVisTau_mass[i];
      GenVisTau_phi[i];
      GenVisTau_pt[i];
      GenVisTau_charge[i];
      GenVisTau_genPartIdxMother[i];
      GenVisTau_status[i];
    }
  }
  //nTrigObj
  *nTrigObj;
  for(unsigned int i=0; i<(*nTrigObj);i++){
    TrigObj_pt[i];
    TrigObj_eta[i];
    TrigObj_phi[i];
    TrigObj_l1pt[i];
    TrigObj_l1pt_2[i];
    TrigObj_l2pt[i];
    TrigObj_id[i];
    TrigObj_l1iso[i];
    TrigObj_l1charge[i];
    TrigObj_filterBits[i];
  }
  
  //END
}
