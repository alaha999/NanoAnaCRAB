#define VLLAna_cxx
#include "VLLAna.h"
#include <TH2.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;


//Load Custom Header Files
#include "setup/applyCorrection.h"



void VLLAna::Begin(TTree * /*tree*/)
{
  //iiserlogo();
  //cout<<"Inside Begin()"<<endl;
  TString option = GetOption();
}

void VLLAna::SlaveBegin(TTree * /*tree*/)
{
  time(&start);
  
  //cout<<"Inside SlaveBegin()"<<endl;
  TString option = GetOption();
  nEvtTotal = 0;
  nEvtRan = 0;
  nEvtValidRunLS=0;
  nEvtTrigger=0;
  
  //Declare Analysis Working Point (default is Tight selections)
  //OPTIONS: TIGHT, LOOSE
  
  //AnalysisWP="LOOSE";
  AnalysisWP="TIGHT";
    
  //Initialize event counter
  n_4L   = 0;
  n_3L1T = 0;
  n_3L   = 0;
  n_2L2T = 0;
  n_2L1T = 0;
  n_1L3T = 0;
  n_1L2T = 0;
  n_2L   = 0;
  n_1L2J = 0;
  
  n_muJJ_passglobalsel  = 0;
  n_muJJ_passtreesel = 0;
  
  genEventsumw=0;


  string skimfname = string(_SkimFileName);
  size_t pos = skimfname.find(".root");
  if (pos != std::string::npos) {
    skimfname = skimfname.substr(0, pos); // Trim the ".root"
  }
  
  string filename1 = skimfname+"_nominal" + ".root";
  //string filename2 = skimfname+"_jecup"   + ".root";
  //string filename3 = skimfname+"_jecdown" + ".root";
  //string filename4 = skimfname+"_jerup"   + ".root";
  //string filename5 = skimfname+"_jerdown" + ".root";
  
  InitializeAnalysisTreeBranch(filename1,outputTree1);
  //InitializeAnalysisTreeBranch(filename2,outputTree2);
  //InitializeAnalysisTreeBranch(filename3,outputTree3);
  //InitializeAnalysisTreeBranch(filename4,outputTree4);
  //InitializeAnalysisTreeBranch(filename5,outputTree5);
  
  
  //Histograms                                                                                                                                               
  //Create the histogram file
  string hstfilename = skimfname+"_nominal_hst"+".root";
  _HstFile = new TFile(hstfilename.c_str(),"recreate");
  BookHistograms();
  
  
  //Load GoldenJson Data
  //jsondata is global variable of nlohmann::json namespace
  //defined in header file
  jsondata=load_json_data(_year);
  
}


void VLLAna::SlaveTerminate()  
{
  cout<<" "<<endl;
  cout<<"<<< EVENTS PROCESSED! COLLECTING FILES >>>"<<endl;
  
  _HstFile->cd();
  _HstFile->Write();
  _HstFile->Close();
  
  cout<<"genEventsumW from skim file:" << _geneventsumw <<endl;
  //create a TH1F histogram and write in skimroot file
  TH1F *hCount = new TH1F("hCount", "hCount;;",7, 0.5, 7.5);
  hCount->SetBinContent(1,genEventsumw);
  hCount->SetBinContent(2,nEvtTotal);
  hCount->SetBinContent(3,nEvtRan);
  hCount->SetBinContent(4,nEvtValidRunLS);
  hCount->SetBinContent(5,nEvtTrigger);
  hCount->SetBinContent(6,n_muJJ_passtreesel);
  
  //overwrite nEvtTotal for skim samples
  if(_geneventsumw>0)hCount->SetBinContent(2,_geneventsumw);

  //label the axis
  hCount->GetXaxis()->SetBinLabel(1,"genEventSumW");
  hCount->GetXaxis()->SetBinLabel(2,"nEvtGen");
  hCount->GetXaxis()->SetBinLabel(3,"nEvtRan");
  hCount->GetXaxis()->SetBinLabel(4,"nEvtValidRunLS");
  hCount->GetXaxis()->SetBinLabel(5,"nEvtTrigger");
  hCount->GetXaxis()->SetBinLabel(6,"n1L2Jskim");
  hCount->GetXaxis()->SetBinLabel(7,"genEventsumw_skim");
  
  //iterate over files and save
  for (TFile* file : outputFileVector) {
    cout << "Processing file: " << file->GetName() <<endl;
    //cout<<"Compression level: "<<file->GetCompressionLevel()<<endl;
    //cout<<"Compression algroithm: "<<file->GetCompressionAlgorithm()<<endl;
    // Work with the TFile pointer
    file->cd();
    // Perform any operations on the file
    hCount->Write();
    file->Write();
    file->Close();
  }

  if(_data){
    //Save Passing and failing Json
    string skimfname = string(_SkimFileName);
    size_t pos = skimfname.find(".root");
    if (pos != std::string::npos) {
      skimfname = skimfname.substr(0, pos); // Trim the ".root"
    }
    std::ofstream json_file1(skimfname+"_passing_json.json");
    //cout<<"Passing JSON: "<<passing_json.dump(4)<<endl;
    json_file1 << passing_json.dump(4);
    json_file1.close();
    std::ofstream json_file2(skimfname+"_failing_json.json");
    //cout<<"Failing JSON: "<<failing_json.dump(4)<<endl;
    json_file2 << failing_json.dump(4);
    json_file2.close();
  }
  //Output to screen
  cout<<" "<<endl;
  cout<<"##### SUCCESSFUL JOBS | JOBS REPORT #####"<<endl;
  cout<<"-----------------------------------------"<<endl;
  cout<<"          Total events = "<<nEvtTotal<<endl;
  cout<<"      Total events ran = "<<nEvtRan<<endl;
  cout<<"Total Triggered events = "<< nEvtTrigger <<endl;
  cout<<"          genEventsumw = "<<genEventsumw<<endl;
  cout<<"              sampleID = "<<_sample<<endl;
  cout<<"                  Year = "<<_year<<endl;
  cout<<"                   Era = "<<_era<<endl;
  cout<<"           Analysis WP = "<<AnalysisWP<<endl;
  cout<<"                         "<<endl;
  cout<<"-----------------------------------------------"<<endl;
  cout<<"               Event Selection                 "<<endl;
  cout<<"-----------------------------------------------"<<endl;
  cout<<"No of 4L   Events     = "<<n_4L<<endl;
  cout<<"No of 3L1T Events     = "<<n_3L1T<<endl;
  cout<<"No of 3L   Events     = "<<n_3L<<endl;
  cout<<"No of 2L2T Events     = "<<n_2L2T<<endl;
  cout<<"No of 2L1T Events     = "<<n_2L1T<<endl;
  cout<<"No of 1L3T Events     = "<<n_1L3T<<endl;
  cout<<"No of 1L2T Events     = "<<n_1L2T<<endl;
  cout<<"No of 2L   Events     = "<<n_2L<<endl;  
  cout<<"No of 1L2J Events     = "<<n_1L2J<<endl;
  cout<<"__________________________________________________"<<endl;
  
  cout<<"No of muJJ Events pass global selection   = "<<n_muJJ_passglobalsel<<endl;
  cout<<"No of muJJ Events pass skimtree selection = "<<n_muJJ_passtreesel<<endl;
  cout<<"___________________________________________________"<<endl;
  
  //Open the text output file
  ofstream fout(_SumFileName);
  //Put text output in the summary file.
  fout<<"          Total events = "<<nEvtTotal<<endl;
  fout<<"      Total events ran = "<<nEvtRan<<endl;
  fout<<"Total Triggered events = "<< nEvtTrigger <<endl;
  fout<<"          genEventsumw = "<<genEventsumw<<endl;
  fout<<"              sampleID = "<<_sample<<endl;
  fout<<"                  Year = "<<_year<<endl;
  fout<<"                   Era = "<<_era<<endl;
  fout<<"           Analysis WP = "<<AnalysisWP<<endl;
  fout<<"                         "<<endl;
  fout<<"-----------------------------------------------"<<endl;
  fout<<"               Event Selection                 "<<endl;
  fout<<"-----------------------------------------------"<<endl;
  fout<<"No of 4L   Events     = "<<n_4L<<endl;
  fout<<"No of 3L1T Events     = "<<n_3L1T<<endl;
  fout<<"No of 3L   Events     = "<<n_3L<<endl;
  fout<<"No of 2L2T Events     = "<<n_2L2T<<endl;
  fout<<"No of 2L1T Events     = "<<n_2L1T<<endl;
  fout<<"No of 1L3T Events     = "<<n_1L3T<<endl;
  fout<<"No of 1L2T Events     = "<<n_1L2T<<endl;
  fout<<"No of 2L   Events     = "<<n_2L<<endl;  
  fout<<"No of 1L2J Events     = "<<n_1L2J<<endl;
  fout<<"__________________________________________________"<<endl;
  
  fout<<"No of muJJ Events pass global selection   = "<<n_muJJ_passglobalsel<<endl;
  fout<<"No of muJJ Events pass skimtree selection = "<<n_muJJ_passtreesel<<endl;
  fout<<"___________________________________________________"<<endl;
  
  time(&end);
  
  
  double time_taken = double(end-start);
  cout<<"Time taken by the programe is= "<<fixed<<time_taken<<setprecision(5);
  cout<<"sec"<<endl;
  
  
}


void VLLAna::Terminate()
{
  //   cout<<"Inside Terminate()"<<endl;
}

Bool_t VLLAna::Process(Long64_t entry)
{
  fReader_common.SetLocalEntry(entry);
  if(_year<2020)      fReader_Run2.SetLocalEntry(entry);
  else if(_year>2020) fReader_Run3.SetLocalEntry(entry);

  if(_data == 0){
    fReader_commonMC.SetLocalEntry(entry);
    if(_year<2020)      fReader_Run2MC.SetLocalEntry(entry);
    else if(_year>2020) fReader_Run3MC.SetLocalEntry(entry);
  }
  if(_year==2016)fReader_2016.SetLocalEntry(entry);
  if(_year==2017)fReader_2017.SetLocalEntry(entry);
  if(_year==2018)fReader_2018.SetLocalEntry(entry);
  if(_year==2017 || _year==2018)fReader_2017_2018.SetLocalEntry(entry);
  if(_year==2018 || _year==2022)fReader_2018_2022.SetLocalEntry(entry);  
  if(((_sample < 9) || (_sample>19)) && !_data)fReader_commonMCSpecial.SetLocalEntry(entry);  
  
  //Verbosity determines the number of processed events after which the root prompt is supposed to display a status update.
  if(_verbosity==0 && nEvtTotal%10000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;
  else if(_verbosity>0 && nEvtTotal%10000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;
  
  nEvtTotal++;
  h.nevt->Fill(0);
  nEvtTree=1;
  nEvt = nEvtTree;

  //Debugging for one event
  bool PRINT_DEBUG=(nEvtTotal==1);
  
  if(PRINT_DEBUG)
    cout<<"##### CHECKPOST: fReader successfully set..processing start...####"<<endl;
  
  //ERA settings
  string ERA = "";
  if(_year==2018)ERA="2018";
  else if(_year==2017)ERA="2017";
  else if(_year==2016 && _era=="preVFP")ERA="2016preVFP";
  else if(_year==2016 && _era=="postVFP")ERA="2016postVFP";
  else if(_year==2022)ERA="2022";
  else if(_year==2022 && _era=="EE")ERA="2022EE";
  

  //GenEventSumW
  if(_data==0){
    genEventsumw+= *Generator_weight;
    genEventSumw = *Generator_weight;
    //cout<<"genEventsumw"<<genEventsumw<<"/"<<*Generator_weight<<endl;
  }
  else genEventSumw=1.0;
  
  //check valid LS/RUNNo (important for data)
  bool validRunAndLS=true; //always true by default
  if(_data){
    int runno = *run;
    int lsno  = *luminosityBlock;  
    validRunAndLS=checkJson(_data,runno,lsno);
    CreatePassFailOutputJsonFiles(validRunAndLS,runno,lsno);
  }
  if(validRunAndLS)nEvtValidRunLS++;

  //GoodEvent Flags(MetFilters etc)
  GoodEvt = PassEventFlags(_data,_year);
  GoodEvt = GoodEvt && validRunAndLS;
  
  if(PRINT_DEBUG)cout<<"##### CHECKPOST: Flag applied for selecting good events #####"<<endl;   
  
  if(GoodEvt){
    
    nEvtRan++; //only good events
    h.nevt->Fill(1);

    //Trigger-blocks
    passmuTrigger  = true;
    passeleTrigger = true;
    passmuTrigger  = PassLeptonTrigger(_year,1);
    passeleTrigger = PassLeptonTrigger(_year,0);
    if(_data) triggerRes =(((_lep==1) & passmuTrigger)||((_lep==0) & passeleTrigger));
    else      triggerRes = passmuTrigger || passeleTrigger;
    
    if(PRINT_DEBUG)cout<<"##### CHECKPOST: Trigger condition applied...#####"<<endl;
    
    if(triggerRes){
      nEvtTrigger++; //only triggered events
      h.nevt->Fill(2);
      
      
      //-----------------------------------------------------------------------------------------------------------------------
      //                                              GEN OBJECTS BLOCK BEGINS                                                |
      //-----------------------------------------------------------------------------------------------------------------------
      
      // Filled up with code when Gen studies are needed
      //if(_data==0){
      //float Wmass=get_WMass();
      //}
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
      unsigned int iterator_muon = (unsigned int)*nMuon;
      for(unsigned int i=0; i< (iterator_muon); i++){
	Lepton temp; temp.v.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],0.105); 
	temp.id = -13*Muon_charge[i]; temp.ind = i;  temp.charge = Muon_charge[i];
	temp.iso=Muon_pfRelIso04_all[i];
	temp.sip3d=Muon_sip3d[i];
	temp.deepjet=Jet_btagDeepFlavB[Muon_jetIdx[i]];
	h.ptlep[0]->Fill(temp.v.Pt());
	
	//Selections
	bool cut_pt=(temp.v.Pt()>10);
	bool cut_eta=(fabs(temp.v.Eta())<2.4);
	bool cut_dxy=(fabs(Muon_dxy[i])<0.05);
	bool cut_dz =(fabs(Muon_dz[i])<0.1);
	bool cut_ID =Muon_mediumId[i];
	bool cut_LooseIso=(Muon_pfRelIso04_all[i]<1.0);
	bool cut_TightIso=(Muon_pfRelIso04_all[i]<0.15);
	bool cut_customID=LeptonCustomIDCut(_year,temp.sip3d,temp.deepjet);
	
	bool passcut= temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4 && Muon_mediumId[i]&& Muon_pfRelIso04_all[i]<1.0;
	bool is_promptmuon= fabs(Muon_dxy[i])<0.05 && fabs(Muon_dz[i])<0.1;	
	
	bool analysis_loosemuon= cut_pt & cut_eta & cut_dxy & cut_dz & cut_ID & cut_customID & cut_LooseIso;
	bool analysis_tightmuon= analysis_loosemuon & cut_TightIso;
	
	bool analysis_muon_selection=analysis_tightmuon;//default
	if     (AnalysisWP=="TIGHT") analysis_muon_selection = analysis_tightmuon;
	else if(AnalysisWP=="LOOSE") analysis_muon_selection = analysis_loosemuon;
	
	//Define Loose Muon
	if(analysis_loosemuon){LooseLep.push_back(temp),loosemuon.push_back(temp);}	  
	
	//Define Analysis Level Muon
	if(analysis_muon_selection){
	  nmu++;
	  h.ptlep[1]->Fill(temp.v.Pt());		
	  Muon.push_back(temp);
	  llep.push_back(temp);
	}
      }
      
      //RecoElectron
      h.nlep[2]->Fill(*nElectron);
      unsigned int iterator_ele = (unsigned int)*nElectron;
      for(unsigned int i=0; i< (iterator_ele); i++){
	Lepton temp; temp.v.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],0.000511); 
	temp.id = -11*Electron_charge[i]; temp.ind = i; temp.charge = Electron_charge[i];
	temp.iso=Electron_pfRelIso03_chg[i];
	temp.sip3d=Electron_sip3d[i];
	temp.deepjet=Jet_btagDeepFlavB[Electron_jetIdx[i]];
	temp.muoncleaning=MuonCleaning(temp.v,0);//Cleaning against selected muon: 0=loosemuon,1=analysis muon
	h.ptlep[2]->Fill(temp.v.Pt());
	
	//Object quality criteria
	bool cut_pt = temp.v.Pt()>10;
	bool cut_eta= fabs(temp.v.Eta())<2.4;
	bool is_barrel=fabs(temp.v.Eta())<=1.479;
	bool is_endcap=fabs(temp.v.Eta())>1.479;
	bool cut_dxy= ((is_barrel &(fabs(Electron_dxy[i])<0.05)) ||(is_endcap &(fabs(Electron_dxy[i])<0.1)));  
	bool cut_dz = ((is_barrel &(fabs(Electron_dz[i])<0.1)) ||(is_endcap &(fabs(Electron_dxy[i])<0.2)));
	bool cut_eleMediumID = electronCustomID(Electron_vidNestedWPBitmap[i],3,-1);
	bool cut_eleMediumIDNoIso  = electronCustomID(Electron_vidNestedWPBitmap[i],3,7);
	bool cut_LooseIso = Electron_pfRelIso03_chg[i]<1.0;
	bool cut_customID = LeptonCustomIDCut(_year,temp.sip3d,temp.deepjet);

	//Selection
	bool analysis_looseele=cut_pt & cut_eta & cut_dxy & cut_dz & cut_customID & cut_eleMediumIDNoIso & cut_LooseIso & temp.muoncleaning;
	bool analysis_tightele = analysis_looseele && cut_eleMediumID; //inbuilt iso as cutbased mediumID
	
	//Loose Electron array
	if(analysis_looseele)LooseLep.push_back(temp);
	
	bool analysis_ele_selection=analysis_tightele; //default
	if     (AnalysisWP=="TIGHT") analysis_ele_selection = analysis_tightele;
	else if(AnalysisWP=="LOOSE") analysis_ele_selection = analysis_looseele;
	
	if(analysis_ele_selection){ //We may fill for loose ele condition since we want to predict QCD bkg from Isolation sidebands
	  nel++;
	  h.ptlep[3]->Fill(temp.v.Pt());
	  Electron.push_back(temp);
	  llep.push_back(temp);
	}
      }
      
      //RecoTau
      h.nlep[4]->Fill(*nTau);
      unsigned int iterator_Tau = (unsigned int)*nTau;
      for(unsigned int i=0; i< (iterator_Tau); i++){
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
      //DeepJetWPThreshold (https://btv-wiki.docs.cern.ch/ScaleFactors/#useful-links)
      //  Year          L           M          T
      //  2018         0.0490     0.2783     0.7100
      //  2017         0.0532     0.3040     0.7476
      //  2016preVFP   0.0508     0.2598     0.6502
      //  2016postVFP  0.0480     0.2489     0.6377
      //  2022         0.0583     0.3086     0.7183
      //  2022EE       0.0614     0.3196     0.73
      //  2023         0.0479     0.2431     0.6553
      //  2023BPix     0.048      0.2435     0.6563
      //-----------------------------------------------//
      unsigned int iterator_Jet = (unsigned int)*nJet;
      for(unsigned int i=0; i< (iterator_Jet); i++){
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
	    //nominal JER correction to jets (Only on MC)
	    if(_data==0){
	      //int genind = Jet_genJetIdx[temp.ind];
	      //float jer_nom = jetRF(temp.v.Pt(),temp.v.Eta(),temp.v.Phi(),genind,GenJet_pt[genind],GenJet_eta[genind],GenJet_phi[genind],*fixedGridRhoFastjetAll,ERA,"nom");
	      float jer_nom=1.0;
	      temp.v *=jer_nom;
	    }
	    
	    jets.push_back(temp);
	    njet++;
	    h.ptlep[7]->Fill(temp.v.Pt());
	    
	    //b-tagging	    
	    if(_year ==2022 ? Jet_btagDeepFlavB[i]>=0.3086: (_year ==2018) ? Jet_btagDeepFlavB[i]>=0.2783: (_year == 2017) ? Jet_btagDeepFlavB[i]>=0.3040 :(_year==2016 && _era=="preVFP") ? Jet_btagDeepFlavB[i]>=0.2598:(_year==2016 && _era=="postVFP") ? Jet_btagDeepFlavB[i]>=0.2489:1.0)
	      //B tag discriminator medium working point(float variable)//2016preVFP: 0.2598 and 2016postVFP:2489
	      {
		BTaggedJet.push_back(temp);
		nbjet++;
		h.ptlep[8]->Fill(temp.v.Pt());
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
      Sort(LooseLep);
      Sort(taus);
      Sort(jets);
      Sort(BTaggedJet);
      
      //MET
      metpt = *MET_pt;
      metphi = *MET_phi;
      
      h.etmiss[0]->Fill(metpt);//fill a histogram with the missing Et of the event.
      h.etmiss[1]->Fill(metphi);
      if(_year<2020){
	//XY-correctedMET
	pair<float,float> met_xycorr = METXYCorr_Met_MetPhi(metpt,metphi,*run,ERA,_data,*PV_npvs);
	metpt  = met_xycorr.first;
	metphi = met_xycorr.second;
	h.etmiss[2]->Fill(metpt);
	h.etmiss[3]->Fill(metphi);
      }
      //2018HEMVeto
      //https://indico.cern.ch/event/920157/contributions/3866209/attachments/2055030/3445708/HEM_issue_forMETXMeeting.pdf
      if(_year==2018){
	bool HEMregion=false;
	for(unsigned int i=0;i<jets.size();i++){
	  bool etaregion = jets[i].v.Eta()>-3.0 && jets[i].v.Eta()<-1.3;
	  bool phiregion = jets[i].v.Phi()>-1.57 && jets[i].v.Phi()<-0.87;
	  if(etaregion && phiregion)HEMregion=true;
	  if(HEMregion){
	    //cout<<"Njet="<<jets.size()<<"/jet "<<i<<" /eta/phi/HEMRegion ="<<jets[i].v.Eta()<<"/"<<jets[i].v.Phi()<<"/"<<HEMregion<<endl;
	    HEMVeto=false;break;
	  }
	  else HEMVeto=1.0;
	}

      }
      else HEMVeto=1.0;
      
      
      //  cout<<"jet & met stored"<<endl;

      // Object Selection is completed with different set of criterias on object
      //-----------------------------------------------------------------------------------------------------------------------
      //                                              RECO OBJECTS BLOCK ENDS                                                 |
      //-----------------------------------------------------------------------------------------------------------------------
      
      
      
      
      //-----------------------------------------------------------------------------------------------------------------------
      //                                        EVENT SELECTION                                                               |
      //-----------------------------------------------------------------------------------------------------------------------
      
      //Lepton pt per year
      float trigmupt =0;float trigelept =0;
      if(_year==2016)trigmupt=26,trigelept=30;
      if(_year==2017)trigmupt=29,trigelept=35;
      if(_year==2018)trigmupt=26,trigelept=35;
      if(_year==2022)trigmupt=26,trigelept=35;
      
      bool keepThisEvent=false;
      
      bool is_4L_event   = false;
      bool is_3L1T_event = false;
      bool is_3L_event   = false;
      bool is_2L2T_event = false;
      bool is_2L1T_event = false;
      bool is_1L3T_event = false;
      bool is_1L2T_event = false;
      bool is_2L_event   = false;  //Current VLL analysis
      bool is_1L2J_event = false;  //Current VLL analysis
      
      if((int)llep.size()>3){ //>=4L
	bool passTrigger = false;
	if(abs(llep.at(0).id)==13)
	  passTrigger = llep.at(0).v.Pt()>trigmupt;
	else if(abs(llep.at(0).id)==11)
	  passTrigger = llep.at(0).v.Pt()>trigelept;
	
	if(passTrigger)
	  is_4L_event=true,n_4L++;
      }
      else if((int)llep.size()==3){ //3L
	bool passTrigger = false;
	if(abs(llep.at(0).id)==13)
	  passTrigger = llep.at(0).v.Pt()>trigmupt;
	else if(abs(llep.at(0).id)==11)
	  passTrigger = llep.at(0).v.Pt()>trigelept;
	
	if(passTrigger){
	  if((int)taus.size()>0)is_3L1T_event=true,n_3L1T++; //3L1T
	  else is_3L_event=true,n_3L++; //3L	  
	}
      }
      else if((int)llep.size()==2){ //2L
	bool passTrigger = false;
	if(abs(llep.at(0).id)==13)
	  passTrigger = llep.at(0).v.Pt()>trigmupt;
	else if(abs(llep.at(0).id)==11)
	  passTrigger = llep.at(0).v.Pt()>trigelept;
	
	if(passTrigger){
	  if((int)taus.size()>1)is_2L2T_event=true,n_2L2T++;      //2L2T
	  else if((int)taus.size()>0)is_2L1T_event=true,n_2L1T++; //2L1T
	  
	  //Count2L(any no of taus is allowed etc)
	  is_2L_event=true,n_2L++; //2L(Exactly 2L)
	}
      }
      else if((int)llep.size()==1){ //1L
	bool passTrigger = false;
	if(abs(llep.at(0).id)==13)
	  passTrigger = llep.at(0).v.Pt()>trigmupt;
	else if(abs(llep.at(0).id)==11)
	  passTrigger = llep.at(0).v.Pt()>trigelept;
	
	if(passTrigger){
	  if((int)taus.size()>2)is_1L3T_event=true,n_1L3T++;      //1L3T
	  else if((int)taus.size()>1)is_1L2T_event=true,n_1L2T++; //1L2T
	  
	  //count 1L2J
	  if((int)jets.size()>1)is_1L2J_event=true,n_1L2J++; //1L2J(any no of taus)
	}
      }
      
      
      //-----------------------------------------------------------------------------------------------------------------------
      //                                          Test Scale Factors                                                          |
      //-----------------------------------------------------------------------------------------------------------------------
      if((int)llep.size()>0 && (int)jets.size()>0){
	int lepflav = abs(llep[0].id);
	Lepton lep0 =llep[0];
	string era  = ERA;
	float lepIDSF_nom  = 1.0;float lepIDSF_up  = 1.0;float lepIDSF_down  = 1.0;
	float lepISOSF_nom = 1.0;float lepISOSF_up = 1.0;float lepISOSF_down = 1.0;
	
	lepIDSF_nom    = LeptonIDSF(lepflav,lep0.v.Pt(),lep0.v.Eta(),lep0.v.Phi(),era,"nom");
	lepIDSF_up     = LeptonIDSF(lepflav,lep0.v.Pt(),lep0.v.Eta(),lep0.v.Phi(),era,"up");
	lepIDSF_down   = LeptonIDSF(lepflav,lep0.v.Pt(),lep0.v.Eta(),lep0.v.Phi(),era,"down");
	
	lepISOSF_nom   = LeptonISOSF(lepflav,lep0.v.Pt(),lep0.v.Eta(),era,"nom");
	lepISOSF_up    = LeptonISOSF(lepflav,lep0.v.Pt(),lep0.v.Eta(),era,"up");
	lepISOSF_down  = LeptonISOSF(lepflav,lep0.v.Pt(),lep0.v.Eta(),era,"down");
	
	//print	
	//cout<<"lepIDSF: flavor/nom/up/down: "<<lepflav<<"/"<<lepIDSF_nom<<"/"<<lepIDSF_up<<"/"<<lepIDSF_down<<endl;
	//cout<<"lepISOSF: flavor/nom/up/down: "<<lepflav<<"/"<<lepISOSF_nom<<"/"<<lepISOSF_up<<"/"<<lepISOSF_down<<endl;
	
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
    unsigned int iterator_GenPart = (unsigned int)*nGenPart;
    for(unsigned int i=0; i< (iterator_GenPart); i++){
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
  if      (year==2016 && sip3d<10 && deepjet<0.6)passCut=true;
  else if (year==2017 && sip3d<12 && deepjet<0.4)passCut=true;
  else if (year==2018 && sip3d<9 && deepjet<0.3)passCut=true;
  else if (year==2022 && sip3d<9 && deepjet<0.3)passCut=true;  
  return passCut;
}

bool VLLAna:: PassEventFlags(bool isData, int year){
  bool passFlags=true;  
  bool passFlags_MC= *Flag_goodVertices & *Flag_globalSuperTightHalo2016Filter & *Flag_HBHENoiseFilter & *Flag_HBHENoiseIsoFilter & *Flag_EcalDeadCellTriggerPrimitiveFilter & *Flag_BadPFMuonFilter & *Flag_BadPFMuonDzFilter & *Flag_ecalBadCalibFilter;
  bool passFlags_Data= *Flag_goodVertices & *Flag_globalSuperTightHalo2016Filter & *Flag_HBHENoiseFilter & *Flag_HBHENoiseIsoFilter & *Flag_EcalDeadCellTriggerPrimitiveFilter & *Flag_BadPFMuonFilter & *Flag_BadPFMuonDzFilter & *Flag_ecalBadCalibFilter & *Flag_eeBadScFilter;
  
  if(isData)passFlags=passFlags_Data;
  else passFlags=passFlags_MC;

  return passFlags;
}
bool VLLAna::PassLeptonTrigger(int year,int LeptonType){
  bool passTrigger=false;
  
  //SingleMuonTrigger  
  if(LeptonType==1){
    if     (year>=2018)passTrigger= *HLT_IsoMu24;
    else if(year==2017)passTrigger= *HLT_IsoMu27;
    else if(year==2016)passTrigger= (*HLT_IsoMu24==1 || *HLT_IsoTkMu24);
  }
  if(LeptonType==0){
    if     (year>=2018)passTrigger = *HLT_Ele32_WPTight_Gsf;
    else if(year==2017)passTrigger = *HLT_Ele32_WPTight_Gsf_L1DoubleEG;
    else if(year==2016)passTrigger = *HLT_Ele27_WPTight_Gsf;
  }

  return passTrigger;
}

void VLLAna::CreatePassFailOutputJsonFiles(bool validRunLSFlag, int runno, int lsno){
  string strrunno = std::to_string(runno);
  if (validRunLSFlag) {
    // If the run number already exists, append the LS to the existing array
    if (passing_json.contains(strrunno)) {
      if (std::find(passing_json[strrunno].begin(), passing_json[strrunno].end(), lsno) == passing_json[strrunno].end()) {
	passing_json[strrunno].push_back(lsno);  // Append if LS is not already present
      }
    } else {
      // If run number doesn't exist, create a new entry with the LS
      passing_json[strrunno] = json::array();
      passing_json[strrunno].push_back(lsno);
    }
  }
  else {
    // Same logic for failing_json
    if (failing_json.contains(strrunno)) {
      if (std::find(failing_json[strrunno].begin(), failing_json[strrunno].end(), lsno) == failing_json[strrunno].end()) {
	failing_json[strrunno].push_back(lsno);  // Append if LS is not already present
      }
    } else {
      // If run number doesn't exist, create a new entry with the LS
      failing_json[strrunno] = json::array();
      failing_json[strrunno].push_back(lsno);
    }
  }
}
std::pair<float,float> VLLAna::METXYCorr_Met_MetPhi(float uncormet, float uncormet_phi,int runnb, string year, bool isData,int npv){
  
  bool isUL=true;
  bool isMC=!isData;
  
  //Initialize
  std::pair<float,float>  TheXYCorr_Met_MetPhi(uncormet,uncormet_phi);
  if(npv>100) npv=100;
  
  string runera = "";
  
  if     (isMC && year == "2016preVFP" && isUL) runera = "yUL2016MCAPV";
  else if(isMC && year == "2016postVFP" && isUL) runera = "yUL2016MCnonAPV";
  else if(isMC && year == "2017" && isUL) runera = "yUL2017MC";
  else if(isMC && year == "2018" && isUL) runera = "yUL2018MC";
  
  else if(!isMC && runnb >=315252 && runnb <=316995 && isUL) runera = "yUL2018A";
  else if(!isMC && runnb >=316998 && runnb <=319312 && isUL) runera = "yUL2018B";
  else if(!isMC && runnb >=319313 && runnb <=320393 && isUL) runera = "yUL2018C";
  else if(!isMC && runnb >=320394 && runnb <=325273 && isUL) runera = "yUL2018D";

  else if(!isMC && runnb >=297020 && runnb <=299329 && isUL){ runera = "yUL2017B";}
  else if(!isMC && runnb >=299337 && runnb <=302029 && isUL){ runera = "yUL2017C";}
  else if(!isMC && runnb >=302030 && runnb <=303434 && isUL){ runera = "yUL2017D";}
  else if(!isMC && runnb >=303435 && runnb <=304826 && isUL){ runera = "yUL2017E";}
  else if(!isMC && runnb >=304911 && runnb <=306462 && isUL){ runera = "yUL2017F";}

  else if(!isMC && runnb >=272007 && runnb <=275376 && isUL) runera = "yUL2016B";
  else if(!isMC && runnb >=275657 && runnb <=276283 && isUL) runera = "yUL2016C";
  else if(!isMC && runnb >=276315 && runnb <=276811 && isUL) runera = "yUL2016D";
  else if(!isMC && runnb >=276831 && runnb <=277420 && isUL) runera = "yUL2016E";
  else if(!isMC && ((runnb >=277772 && runnb <=278768) || runnb==278770) && isUL) runera = "yUL2016F";
  else if(!isMC && ((runnb >=278801 && runnb <=278808) || runnb==278769) && isUL) runera = "yUL2016Flate";
  else if(!isMC && runnb >=278820 && runnb <=280385 && isUL) runera = "yUL2016G";
  else if(!isMC && runnb >=280919 && runnb <=284044 && isUL) runera = "yUL2016H";


  double METxcorr(0.),METycorr(0.);

  //UL2016
  if(runera=="yUL2016B") METxcorr = -(-0.0214894*npv +-0.188255);
  if(runera=="yUL2016B") METycorr = -(0.0876624*npv +0.812885);
  if(runera=="yUL2016C") METxcorr = -(-0.032209*npv +0.067288);
  if(runera=="yUL2016C") METycorr = -(0.113917*npv +0.743906);
  if(runera=="yUL2016D") METxcorr = -(-0.0293663*npv +0.21106);
  if(runera=="yUL2016D") METycorr = -(0.11331*npv +0.815787);
  if(runera=="yUL2016E") METxcorr = -(-0.0132046*npv +0.20073);
  if(runera=="yUL2016E") METycorr = -(0.134809*npv +0.679068);
  if(runera=="yUL2016F") METxcorr = -(-0.0543566*npv +0.816597);
  if(runera=="yUL2016F") METycorr = -(0.114225*npv +1.17266);
  if(runera=="yUL2016Flate") METxcorr = -(0.134616*npv +-0.89965);
  if(runera=="yUL2016Flate") METycorr = -(0.0397736*npv +1.0385);
  if(runera=="yUL2016G") METxcorr = -(0.121809*npv +-0.584893);
  if(runera=="yUL2016G") METycorr = -(0.0558974*npv +0.891234);
  if(runera=="yUL2016H") METxcorr = -(0.0868828*npv +-0.703489);
  if(runera=="yUL2016H") METycorr = -(0.0888774*npv +0.902632);
  if(runera=="yUL2016MCnonAPV") METxcorr = -(-0.153497*npv +-0.231751);
  if(runera=="yUL2016MCnonAPV") METycorr = -(0.00731978*npv +0.243323);
  if(runera=="yUL2016MCAPV") METxcorr = -(-0.188743*npv +0.136539);
  if(runera=="yUL2016MCAPV") METycorr = -(0.0127927*npv +0.117747);
  
  //UL2017
  if(runera=="yUL2017B") METxcorr = -(-0.211161*npv +0.419333);
  if(runera=="yUL2017B") METycorr = -(0.251789*npv +-1.28089);
  if(runera=="yUL2017C") METxcorr = -(-0.185184*npv +-0.164009);
  if(runera=="yUL2017C") METycorr = -(0.200941*npv +-0.56853);
  if(runera=="yUL2017D") METxcorr = -(-0.201606*npv +0.426502);
  if(runera=="yUL2017D") METycorr = -(0.188208*npv +-0.58313);
  if(runera=="yUL2017E") METxcorr = -(-0.162472*npv +0.176329);
  if(runera=="yUL2017E") METycorr = -(0.138076*npv +-0.250239);
  if(runera=="yUL2017F") METxcorr = -(-0.210639*npv +0.72934);
  if(runera=="yUL2017F") METycorr = -(0.198626*npv +1.028);
  if(runera=="yUL2017MC") METxcorr = -(-0.300155*npv +1.90608);
  if(runera=="yUL2017MC") METycorr = -(0.300213*npv +-2.02232);
  
  //UL2018
  if(runera=="yUL2018A") METxcorr = -(0.263733*npv +-1.91115);
  if(runera=="yUL2018A") METycorr = -(0.0431304*npv +-0.112043);
  if(runera=="yUL2018B") METxcorr = -(0.400466*npv +-3.05914);
  if(runera=="yUL2018B") METycorr = -(0.146125*npv +-0.533233);
  if(runera=="yUL2018C") METxcorr = -(0.430911*npv +-1.42865);
  if(runera=="yUL2018C") METycorr = -(0.0620083*npv +-1.46021);
  if(runera=="yUL2018D") METxcorr = -(0.457327*npv +-1.56856);
  if(runera=="yUL2018D") METycorr = -(0.0684071*npv +-0.928372);
  if(runera=="yUL2018MC") METxcorr = -(0.183518*npv +0.546754);
  if(runera=="yUL2018MC") METycorr = -(0.192263*npv +-0.42121);


  float CorrectedMET_x = uncormet *cos( uncormet_phi)+METxcorr;
  float CorrectedMET_y = uncormet *sin( uncormet_phi)+METycorr;  

  float CorrectedMET = sqrt(CorrectedMET_x*CorrectedMET_x+CorrectedMET_y*CorrectedMET_y);
  float CorrectedMETPhi;
  if(CorrectedMET_x==0 && CorrectedMET_y>0) CorrectedMETPhi = TMath::Pi();
  else if(CorrectedMET_x==0 && CorrectedMET_y<0 )CorrectedMETPhi = -TMath::Pi();
  else if(CorrectedMET_x >0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x);
  else if(CorrectedMET_x <0&& CorrectedMET_y>0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) + TMath::Pi();
  else if(CorrectedMET_x <0&& CorrectedMET_y<0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) - TMath::Pi();
  else CorrectedMETPhi =0;
  
  TheXYCorr_Met_MetPhi.first= CorrectedMET;
  TheXYCorr_Met_MetPhi.second= CorrectedMETPhi;

  return TheXYCorr_Met_MetPhi;  
  
}

bool VLLAna::TriggerObjectMatching(int id,TLorentzVector t){
  bool result=false;
  float dr=999.0,drmin=999.0;
  float trigobjpt=0;
  int n_HLTid13=0;
  int n_matchedHLT=0;
  int n_HLTid13_drmatched=0;
  int trigobjindex=-1;
  unsigned int iterator_TrigObj = (unsigned int)*nTrigObj;
  for(unsigned int i=0; i< (iterator_TrigObj); i++){
    int triggerObj_id =1;
    if(_year<2020)triggerObj_id=abs(TrigObj_id[i]);
    else triggerObj_id=TrigObj_id[i];
    if(abs(id)==triggerObj_id){
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
  h.etmiss[0]= new TH1F("MET","Missing E_{T}",1000,0,1000); h.etmiss[0]->Sumw2();
  h.etmiss[1]= new TH1F("MET_Phi","Missing E_{T} phi",100,-5,5); h.etmiss[1]->Sumw2();
  h.etmiss[2]= new TH1F("corrected_MET","Missing E_{T}",1000,0,1000); h.etmiss[2]->Sumw2();
  h.etmiss[3]= new TH1F("corrected MET Phi","Missing E_{T} phi",100,-5,5); h.etmiss[3]->Sumw2();
  //Channel 1L2J
  h.channel1L2J[0] = new TH1F("channel1L2J_Count","channel 1L2J count: 1=incl,2=mu,3=ele",10,0,10);
  for(int i=0; i<1; i++) h.channel1L2J[i]->Sumw2();
  
  //GenStudy
  h.genstudy[0] = new TH1F("W_mass_gen","W mass at gen level",1000,0,1000);
  h.genstudy[1] = new TH1F("W_mass_LHE","W mass at LHE level",1000,0,1000);
  for(int i=0; i<2; i++)h.genstudy[i]->Sumw2();
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

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
//Analysis Functions//
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
void VLLAna::Analysis_LJJ(TTree* outputTree, string era, string jetsystvar){
  
  //Verify robustness
  //Applying Jer/Jec successively can apply these correction on nominal
  //pt multiple times which is incorrect. This is what we do if we call
  //this function multiple times in the same event loop. To ease the
  //workflow we call this function multiple times for 4 different variation
  //of jerup, jerdown, jecup and jecdown. To ensure that correction only applied
  //for once, we have to create a copy of struct, Lepton<vector> jets. To verify the nominal jet pt
  //let's fill histogram of nominal jet pt before the function starts doing what it is
  //supposed to do. For all iteration, the nominal pt must not change.
  
  //h.anaverify[0][1]->Fill(jets[0].v.Pt()); //leading pt is enough
  
  //create copy of Lepton Objects
  vector<Lepton> jetscopy = jets;
  
  //Best Jet Pairs coming from W
  int jet0index=999,jet1index=999;
  vector<int>jetindex=BestJetPairFinder(jetscopy); //This function returns a vector array of indices
  jet0index=jetindex[0];
  jet1index=jetindex[1];
  
  
  //create instance of the object
  Lepton lep0 = llep[0];
  Lepton jet0 = jetscopy[jet0index];
  Lepton jet1 = jetscopy[jet1index];
  int lepflav = abs(lep0.id);
  
  
  //EventWeight Calculation
  float evtwt = 1.0;
  //LeptonID_ISO
  float lepIDSF_nom  = 1.0;float lepIDSF_up  = 1.0;float lepIDSF_down  = 1.0;
  float lepISOSF_nom = 1.0;float lepISOSF_up = 1.0;float lepISOSF_down = 1.0;
  float btagSF_nom = 1.0;
  float btagSFbc_upCorrelated=1.0;
  float btagSFbc_downCorrelated=1.0;
  float btagSFbc_upUncorrelated=1.0;
  float btagSFbc_downUncorrelated=1.0;
  float btagSFlight_upCorrelated=1.0;
  float btagSFlight_downCorrelated=1.0;
  float btagSFlight_upUncorrelated=1.0;
  float btagSFlight_downUncorrelated=1.0;
  float pileupwt_nom =1.0;
  float pileupwt_up  =1.0;
  float pileupwt_down=1.0;
  float pdf_up =1.0;
  float pdf_down=1.0;
  float qcdscale_up=1.0;
  float qcdscale_down=1.0;

  bool apply_SF=false;
  
  if((_data==0) & (apply_SF)){
    //PILEUP
    //MUON ID/ISO
    //EGAMMA ID
    //BTV
    //JEC
    //JER

    //PDF and QCD scale Uncertainties(LHEPdfWeight & LHEScaleWeight)
    //https://cms-talk.web.cern.ch/t/pileup-reweighting-questions/30405
    //https://cms-talk.web.cern.ch/t/lhescaleweight-and-lhepdfweight/44715
    if((_sample<9) || (_sample>19)){ //_sample:9-19>>QCD Multijet sample. These branches are not available.
      vector<float>pdfweights(LHEPdfWeight.begin() + 1, LHEPdfWeight.end()); 
      vector<float>qcdscaleweights={LHEScaleWeight[1], LHEScaleWeight[3], LHEScaleWeight[5],LHEScaleWeight[7]};  
      
      pdf_up         = *std::max_element(pdfweights.begin(),pdfweights.end());
      pdf_down       = *std::min_element(pdfweights.begin(),pdfweights.end());
      qcdscale_up    = *std::max_element(qcdscaleweights.begin(),qcdscaleweights.end());
      qcdscale_down  = *std::min_element(qcdscaleweights.begin(),qcdscaleweights.end());
    }
    
    //PILE-UP
    pileupwt_nom   = pileupWeight(*Pileup_nTrueInt,era,"nom");
    pileupwt_up    = pileupWeight(*Pileup_nTrueInt,era,"up");
    pileupwt_down  = pileupWeight(*Pileup_nTrueInt,era,"down");
    
    //LEPTON SF
    if(lepflav==13){
      lepIDSF_nom  = muonIDSF(lep0.v.Pt(),lep0.v.Eta(),era,"nom");
      lepIDSF_up   = muonIDSF(lep0.v.Pt(),lep0.v.Eta(),era,"up");
      lepIDSF_down = muonIDSF(lep0.v.Pt(),lep0.v.Eta(),era,"down");
      
      lepISOSF_nom   = muonIsoSF(lep0.v.Pt(),lep0.v.Eta(),era,"nom");
      lepISOSF_up    = muonIsoSF(lep0.v.Pt(),lep0.v.Eta(),era,"up");
      lepISOSF_down  = muonIsoSF(lep0.v.Pt(),lep0.v.Eta(),era,"down");
    }
    else if (lepflav==11){
      lepIDSF_nom = egammaIDSF(lep0.v.Pt(),lep0.v.Eta(),era,"nom");
      lepIDSF_up = egammaIDSF(lep0.v.Pt(),lep0.v.Eta(),era,"up");
      lepIDSF_down = egammaIDSF(lep0.v.Pt(),lep0.v.Eta(),era,"down");
    }
    
    //BTV	    
    //btagSF_nom               = btagIDSF(_sample,jetscopy,era,"nom");
    //btagSF_upUncorrelated    = btagIDSF(_sample,jetscopy,era,"upUncorrelated");
    //btagSF_downUncorrelated  = btagIDSF(_sample,jetscopy,era,"downUncorrelated");
    //btagSF_upCorrelated      = btagIDSF(_sample,jetscopy,era,"upCorrelated");
    //btagSF_downCorrelated    = btagIDSF(_sample,jetscopy,era,"downCorrelated");
    
    vector<float>vec_btagSF_nom=btagIDSFv2(_sample,jetscopy,era,"bc","nom");
    vector<float>vec_btagSF_upUncorr=btagIDSFv2(_sample,jetscopy,era,"bc","upUncorrelated");
    vector<float>vec_btagSF_downUncorr=btagIDSFv2(_sample,jetscopy,era,"bc","downUncorrelated");
    vector<float>vec_btagSF_upCorr=btagIDSFv2(_sample,jetscopy,era,"bc","upCorrelated");
    vector<float>vec_btagSF_downCorr=btagIDSFv2(_sample,jetscopy,era,"bc","downCorrelated");
    
    btagSF_nom                    = vec_btagSF_nom[0]*vec_btagSF_nom[1];
    btagSFbc_upUncorrelated       = vec_btagSF_upUncorr[0];
    btagSFbc_downUncorrelated     = vec_btagSF_downUncorr[0];
    btagSFbc_upCorrelated         = vec_btagSF_upCorr[0];
    btagSFbc_downCorrelated       = vec_btagSF_downCorr[0];
    btagSFlight_upUncorrelated    = vec_btagSF_upUncorr[1];
    btagSFlight_downUncorrelated  = vec_btagSF_downUncorr[1];
    btagSFlight_upCorrelated      = vec_btagSF_upCorr[1];
    btagSFlight_downCorrelated    = vec_btagSF_downCorr[1];
    
    //cout<<"btag sf/ btagsf_bc/btagsf_light"<<btagSF_nom<<"/"<<btagSFbc_nom*btagSFlight_nom<<"/"<<btagSFbc_nom<<"/"<<btagSFlight_nom<<"/"<<endl;    
    //cout<<"btagsf(nom) ="<<btagSF_nom<<"/SAMPLEID="<<_sample<<endl;
    
    //Trigger
    float triggerSF_nom = 1.0;
    
    //JME
    for(int i=0;i<(int)jetscopy.size();i++){
      
      //JER      
      int genind = Jet_genJetIdx[jetscopy[i].ind];
      float jer_up = jetRF(jetscopy[i].v.Pt(),jetscopy[i].v.Eta(),jetscopy[i].v.Phi(),genind,GenJet_pt[genind],GenJet_eta[genind],GenJet_phi[genind],*fixedGridRhoFastjetAll,era,"up");
      float jer_down = jetRF(jetscopy[i].v.Pt(),jetscopy[i].v.Eta(),jetscopy[i].v.Phi(),genind,GenJet_pt[genind],GenJet_eta[genind],GenJet_phi[genind],*fixedGridRhoFastjetAll,era,"down");

      //JEC      
      float jec_up = jetSF(jetscopy[i].v.Pt(),jetscopy[i].v.Eta(),era,"up");
      float jec_down = jetSF(jetscopy[i].v.Pt(),jetscopy[i].v.Eta(),era,"down");
      
      if(jetsystvar=="jer_up")   jetscopy[i].v *= jer_up;
      if(jetsystvar=="jer_down") jetscopy[i].v *= jer_down;
      if(jetsystvar=="jec_up")	 jetscopy[i].v *= jec_up;
      if(jetsystvar=="jec_down") jetscopy[i].v *= jec_down;
      
      //if(i<1){
      //check one element
      //float j0pt=jetscopy[i].v.Pt();
      //cout<<"Event/JetSyst: "<<nEvtTotal<<"/"<<jetsystvar<<endl;
      //cout<<"jerup/jerdown/jecup/jecdown: "<<jer_up<<"/"<<jer_down<<"/"<<jec_up<<"/"<<jec_down<<endl;
      //cout<<"jets vector(before/after): "<<jets[i].v.Pt()<<"/"<<j0pt<<"/"<<j0pt*jer_up<<j0pt*jer_down<<"/"<<j0pt*jec_up<<"/"<<j0pt*jec_down<<endl; 
      //}
      
    }
    
  }
  
  
  
  //Muon Rochester Correction
  if(apply_SF){
    float muonRoccRcorr=1.0;
    bool isGenMatched=false;
    bool isMC = !(_data);
    float genMuPt = 0;
    if(isMC && lepflav==13){
      int matchedgenMuonIndex=Muon_genPartIdx[lep0.ind];
      if(matchedgenMuonIndex<0)isGenMatched=false;
      else isGenMatched=true,genMuPt = GenPart_pt[matchedgenMuonIndex];
      //cout<<"genmu index/pt/genpt="<<matchedgenMuonIndex<<"/"<<lep0.v.Pt()<<"/"<<genMuPt<<endl;
    }
    //muonRoccRcorr= muon_RochesterSF(lep0.charge,lep0.v.Pt(),lep0.v.Phi(),lep0.v.Eta(),isMC,isGenMatched,genMuPt,Muon_nTrackerLayers[lep0.ind],era,"nom");
    //float corrected_muon_pt=lep0.v.Pt()*muonRoccRcorr;
    
    //shift the muon energy scale
    //lep0.v.SetPtEtaPhiM(corrected_muon_pt,lep0.v.Eta(),lep0.v.Phi(),lep0.v.M());
  }
  bool DEBUG=false;
  if(DEBUG){
    cout<<"Year / JetSystVar="<<era<<"/"<<jetsystvar<<endl;
    cout<<"lepIDSF/lepIDSF_up/lepIDSF_down"<<lepIDSF_nom<<"/"<<lepIDSF_up<<"/"<<lepIDSF_down<<endl;
    cout<<"lepISOSF/lepISOSF_up/lepISOSF_down"<<lepISOSF_nom<<"/"<<lepISOSF_up<<"/"<<lepISOSF_down<<endl;
    cout<<"btagSFbc(nom/upUncorrelated/upCorrelated/downUncorrelated/downCorrelated) = "<< btagSF_nom <<"/"<<btagSFbc_upUncorrelated<<"/"<<btagSFbc_upCorrelated<<"/"<<btagSFbc_downUncorrelated<<"/"<<btagSFbc_downCorrelated<<endl;
    cout<<"btagSFlight(nom/upUncorrelated/upCorrelated/downUncorrelated/downCorrelated) = "<< btagSF_nom <<"/"<<btagSFlight_upUncorrelated<<"/"<<btagSFlight_upCorrelated<<"/"<<btagSFlight_downUncorrelated<<"/"<<btagSFlight_downCorrelated<<endl;
    cout<<"\n-------------------------"<<endl;
    cout<<"PILE-UP"<<endl;
    cout<<"pileup weight (nom/up/down) = "<<pileupwt_nom <<"/"<<pileupwt_up <<"/"<<pileupwt_down<<endl;
    cout<<"PDF uncertainties (up/down) = "<<pdf_up<<"/"<<pdf_down<<endl;
    cout<<"QCD scale uncertainties (up/down) = "<<qcdscale_up<<"/"<<qcdscale_down<<endl;
    //cout<<"muon rochester corr="<<muonRoccRcorr<<"/"<<"Before/After Pt="<<lep0.v.Pt()<<"/"<<corrected_muon_pt<<endl;        
  }
  
  //Redefine
  jet0 = jetscopy[jet0index];
  jet1 = jetscopy[jet1index];
  //check jetpt(they should be the same)
  //cout<<"Event/Jet0 instance pt/jetscopy vec pt"<<nEvtTotal<<"/"<<jet0.v.Pt()<<"/"<<jetscopy[jet0index].v.Pt()<<endl;
  
  ///analysis variables  
  //Properties 
  float btagscore=0;float avgbtagscore=0;
  float avgCvsBscore=0;float avgCvsLscore=0;
  float avgQGscore=0;int nparticle=0;int avgnparticle=0;
  float HT=0.0;
  
  for(int i=0;i<(int)jetscopy.size();i++){
    HT=HT+jetscopy.at(i).v.Pt();
    btagscore += Jet_btagDeepFlavB[jetscopy[i].ind];
    nparticle += Jet_nConstituents[jetscopy[i].ind];
    avgCvsBscore+= Jet_btagDeepFlavCvB[jetscopy[i].ind];
    avgCvsLscore+= Jet_btagDeepFlavCvL[jetscopy[i].ind];
    avgQGscore+= Jet_btagDeepFlavQG[jetscopy[i].ind];
  }
  avgbtagscore = btagscore/jetscopy.size();
  avgnparticle = nparticle/jetscopy.size();
  avgCvsBscore = avgCvsBscore/jetscopy.size();
  avgCvsLscore = avgCvsLscore/jetscopy.size();
  avgQGscore   = avgQGscore/jetscopy.size();
  
  //physics variables  
  float mtlep0       = sqrt(2.0*lep0.v.Pt()*metpt*(1-cos(delta_phi(lep0.v.Phi(),metphi))));
  float mtjet0       = sqrt(2.0*jet0.v.Pt()*metpt*(1-cos(delta_phi(jet0.v.Phi(),metphi))));
  float mtjet1       = sqrt(2.0*jet1.v.Pt()*metpt*(1-cos(delta_phi(jet1.v.Phi(),metphi))));
  float MET          = metpt;
  float ST           = lep0.v.Pt()+MET+HT;
  float dijetMass    = getInvMass(jet0.v,jet1.v);
  float drjet01      = jet0.v.DeltaR(jet1.v);
  float dphijet01    = delta_phi(jet0.v.Phi(),jet1.v.Phi());
  float dijetMT      = sqrt(2.0*(jet0.v+jet1.v).Pt()*metpt*(1-cos(delta_phi((jet0.v+jet1.v).Phi(),metphi))));
  float dijetPT      = (jet0.v+jet1.v).Pt();
  int njet           = jetscopy.size();
  int nbjet          = BTaggedJet.size();
  float dphimetjet0  = delta_phi(jet0.v.Phi(),metphi);
  float dphimetjet1  = delta_phi(jet1.v.Phi(),metphi);
  float dphimetlep0  = delta_phi(lep0.v.Phi(),metphi);
  float dphijet0lep0 = delta_phi(lep0.v.Phi(),jet0.v.Phi());
  float dphijet1lep0 = delta_phi(lep0.v.Phi(),jet1.v.Phi());
  float dphidijetlep0= delta_phi(lep0.v.Phi(),(jet0.v+jet1.v).Phi());
  float dphimetdijet = delta_phi(metphi,(jet0.v+jet1.v).Phi());
  float LeadingJetPt = jet0.v.Pt();
  float drjet0lep0   = lep0.v.DeltaR(jet0.v);
  float drjet1lep0   = lep0.v.DeltaR(jet1.v);
  float ljjsysPT       = (lep0.v +jet0.v + jet1.v).Pt();
  float dphimetljjsys  =  delta_phi((lep0.v +jet0.v + jet1.v).Phi(),metphi);
  float dphiljjsyslep0 = delta_phi((lep0.v +jet0.v + jet1.v).Phi(),lep0.v.Phi());
  float dphiljjsysjet0 = delta_phi((lep0.v +jet0.v + jet1.v).Phi(),jet0.v.Phi());
  float dphiljjsysjet1 = delta_phi((lep0.v +jet0.v + jet1.v).Phi(),jet1.v.Phi());
  float ljjsysPTMETratio = (ljjsysPT/MET);
  float ljjsysMass       = getInvMass(jet0.v+jet1.v,lep0.v);
  float ljjsysMassMETratio = (ljjsysMass/MET);

  //New Variables
  float dphidijetjet0=delta_phi(jet0.v.Phi(),(jet0.v+jet1.v).Phi());
  float dphidijetjet1=delta_phi(jet1.v.Phi(),(jet0.v+jet1.v).Phi());  
  float detadijetjet0=jet0.v.Eta()-(jet0.v+jet1.v).Eta();
  float detadijetjet1=jet1.v.Eta()-(jet0.v+jet1.v).Eta();
  // Event PT balance and Zeppenfield variable
  float ptW= sqrt(pow(lep0.v.Pt(),2)+pow(metpt,2)+2*lep0.v.Pt()*cos(delta_phi(lep0.v.Phi(),metphi)));
  float R = sqrt(pow((lep0.v+jet0.v+jet1.v).Pt(),2)+pow(metpt,2)+2*(lep0.v+jet0.v+jet1.v).Pt()*metpt*cos(delta_phi((lep0.v+jet0.v+jet1.v).Phi(),metphi)))/(jet0.v.Pt()+jet1.v.Pt()+ptW);
  
  
  //Fill the Analysis Variables
  lep0_flavor=lep0.id;
  lep0_pt=lep0.v.Pt();
  jet0_pt=jet0.v.Pt();
  jet1_pt=jet1.v.Pt();
  lep0_mt = mtlep0;
  jet0_mt =mtjet0;
  jet1_mt =mtjet1;
  dijet_mass = dijetMass;
  deltaR_jet01 = drjet01;
  deltaPhi_metjet0 = dphimetjet0;
  deltaPhi_metjet1 = dphimetjet1;
  deltaPhi_metlep0 = dphimetlep0;
  deltaPhi_jet0lep0 = dphijet0lep0;
  deltaPhi_jet1lep0 = dphijet1lep0;
  deltaPhi_dijetlep0 = dphidijetlep0;
  deltaPhi_metdijet = dphimetdijet;
  dijet_pt = dijetPT;
  dijet_mt = dijetMT;
  event_MET = MET;
  event_METPhi = metphi;
  event_HT = HT;
  event_ST = ST;
  n_Jet = njet;
  n_bJet = nbjet;
  
  //Extra for Analysis Tree
  lep0_phi=lep0.v.Phi();
  lep0_eta=lep0.v.Eta();
  jet0_phi=jet0.v.Phi();
  jet0_eta=jet0.v.Eta();
  jet1_phi=jet1.v.Phi();
  jet1_eta=jet1.v.Eta();
  
  //Flavor dependent Lepton variables

  lep0_sip3d=lep0.sip3d;
  lep0_deepjet=lep0.deepjet;
  lep0_iso=lep0.iso;
  if(lepflav==13)
    lep0_tight=lep0.iso<0.15;
  else if(lepflav==11)
    lep0_tight=Electron_cutBased[lep0.ind]>2;//MediumID Isolation cuts will be applied
  
  //Extra variables for wjets separation
  ljjsys_PT = ljjsysPT;
  ljjsys_mass= ljjsysMass;
  deltaPhi_ljjsysmet  = dphimetljjsys;
  deltaPhi_ljjsyslep0 = dphiljjsyslep0;
  deltaPhi_ljjsysjet0 = dphiljjsysjet0;
  deltaPhi_ljjsysjet1 = dphiljjsysjet1;
  deepjetQG_jet0= Jet_btagDeepFlavQG[jet0.ind];
  deepjetQG_jet1= Jet_btagDeepFlavQG[jet1.ind];
  event_avgCvsBscore=avgCvsBscore;
  event_avgCvsLscore=avgCvsLscore;
  event_avgQGscore=avgQGscore;
  event_Rpt       = R;
  
  
  //Add SF
  event_lepIDSF       = lepIDSF_nom;
  event_lepIDSF_up    = lepIDSF_up;
  event_lepIDSF_down  = lepIDSF_down;
  event_lepISOSF      = lepISOSF_nom;
  event_lepISOSF_up   = lepISOSF_up;
  event_lepISOSF_down = lepISOSF_down;
  
  event_btagsf        = btagSF_nom;
  event_btagsfCorrelated_up = btagSFbc_upCorrelated*btagSFlight_upCorrelated;
  event_btagsfCorrelated_down = btagSFbc_downCorrelated*btagSFlight_downCorrelated;
  event_btagsfUncorrelated_up = btagSFbc_upUncorrelated*btagSFlight_upUncorrelated;
  event_btagsfUncorrelated_down = btagSFbc_downUncorrelated*btagSFlight_downUncorrelated;  
  event_btagsfbcCorrelated_up = btagSFbc_upCorrelated;
  event_btagsfbcCorrelated_down = btagSFbc_downCorrelated;
  event_btagsfbcUncorrelated_up = btagSFbc_upUncorrelated;
  event_btagsfbcUncorrelated_down = btagSFbc_downUncorrelated;
  
  //TriggerInfo
  event_passmuTrigger = passmuTrigger;
  event_passeleTrigger = passeleTrigger;

  //HEM(2018)
  event_HEMVeto=HEMVeto;
  
  //OnlyMC
  if(_data==0){
    GenWeight = *Generator_weight;
    PU_N      = *Pileup_nPU;
    PU_TrueN  = *Pileup_nTrueInt;
    PileUpWt_nom = pileupwt_nom;
    PileUpWt_up  = pileupwt_up;
    PileUpWt_down  = pileupwt_down;
    //L1PreFireWt_nom  = *L1PreFiringWeight_Nom;
    //L1PreFireWt_up   = *L1PreFiringWeight_Up;
    //L1PreFireWt_down = *L1PreFiringWeight_Dn;
    
    //PDF and QCD scale uncertainties
    PDF_up=pdf_up;
    PDF_down=pdf_down;
    qcd_scale_up=qcdscale_up;
    qcd_scale_down=qcdscale_down;
  }
  
  /*-------------------------------------------------------------------------------------------
    Sep11,2022[decided to not cut on njets]
    Event can have
    - minimum 2 jets, maximum 3 jets
    - dijet Invariant mass must be >40 GeV: ignore the region where soft jets can arise
    - decided on Nov9,2022
    --------------------------------------------------------------------------------------------*/
  //L/J0/J1 Pt : >26/30/30 GeV
  //deltaR between selected object >0.4
  //Lepton Isolation<1.0
  
  bool LJJselection_global= lep0.v.Pt()>26 && jet0.v.Pt()>30 && jet1.v.Pt()>30 && drjet01>0.4 && drjet0lep0>0.4 && drjet1lep0>0.4;
  bool LJJselection_produceTree = LJJselection_global && (nbjet==0);
  
  
  //Fill Analysis Tree
  if(LJJselection_global)n_muJJ_passglobalsel++;
  if(LJJselection_produceTree)outputTree->Fill(),n_muJJ_passtreesel++;
  
}

//Functions
vector<int>VLLAna::BestJetPairFinder(vector<Lepton>jetarray){
  float invmass=-1.0;float dMmin=999999.0;
  float Wmass=80.3;float bestMass=-1.0;
  float dRjets=999.0;
  
  vector<int>jetindex;
  //jetindex.clear();
  if(jetarray.size()>2){
    for(int i=0;i<(int)jetarray.size();i++){
      for(int j=0;j<(int)jetarray.size();j++){
	if(i!=j && j>i){
	  jetindex.clear();
	  invmass=getInvMass(jetarray[i].v,jetarray[j].v);
	  //cout<<"Pair("<<i<<","<<j<<"), mass is = "<<invmass<<endl;	
	  //Find Out the closest one to W Mass
	  //bool onWpair = invmass>60 && invmass<100;
	  bool onWpair=true;
	  if(onWpair){
	    if(abs(Wmass-invmass)<=dMmin){
	      bestMass=invmass;
	      dMmin=abs(Wmass-invmass);
	      jetindex.push_back(i);
	      jetindex.push_back(j);
	    }
	  }
	}
      }
    }
  }
  else jetindex.push_back(0),jetindex.push_back(1); // If only two jets found return 0 & 1;
  
  return jetindex;
}

//void VLLAna::InitializeAnalysisTreeBranch(const std::string &outputFileName, TFile*& outputFile, TTree*& outputTree){
void VLLAna::InitializeAnalysisTreeBranch(const std::string &outputFileName,TTree*& outputTree){
  
  TFile* newFile = TFile::Open(outputFileName.c_str(), "RECREATE");
  //newFile->SetCompressionAlgorithm(1); //1:ZLIB,2:LZMA,4:LZ4,5:ZSTD
  outputTree = new TTree("Events","AnalysisTree");
  
  //Lepton
  outputTree->Branch("lep0_flavor",&lep0_flavor);                //2.Leading Lepton Flavor
  outputTree->Branch("lep0_pt",&lep0_pt);                        //3.Leading Lepton Pt
  outputTree->Branch("lep0_phi",&lep0_phi);                      //4.Leading Lepton Phi
  outputTree->Branch("lep0_eta",&lep0_eta);                      //5.Leading Lepton Eta
  outputTree->Branch("lep0_iso",&lep0_iso);                      //6.Leading Lepton Isolation
  outputTree->Branch("lep0_tight",&lep0_tight);                  //6.Leading Lepton Tight or Loose
  outputTree->Branch("lep0_mt",&lep0_mt);                        //7.Leading Lepton MT
  outputTree->Branch("lep0_sip3d",&lep0_sip3d);                  //8.Leading Lepton SIP3D
  outputTree->Branch("lep0_deepjet",&lep0_deepjet);              //9.Leading Lepton Closer Jet DeepJet Score  
  
  //FirstJet
  outputTree->Branch("jet0_pt",&jet0_pt);                        //9.Leading Jet Pt
  outputTree->Branch("jet0_phi",&jet0_phi);                      //10.Leading Jet Phi
  outputTree->Branch("jet0_eta",&jet0_eta);                      //11.Leading Jet Eta
  outputTree->Branch("jet0_mt",&jet0_mt);                        //12.Leading Jet MT
  
  //SecondJet
  outputTree->Branch("jet1_pt",&jet1_pt);                        //13.Subleading Jet Pt
  outputTree->Branch("jet1_phi",&jet1_phi);                      //14.Subleading Jet Phi
  outputTree->Branch("jet1_eta",&jet1_eta);                      //15.Subleading Jet Eta
  outputTree->Branch("jet1_mt",&jet1_mt);                        //16.Subleading Jet MT
  
  //dijetSystem
  outputTree->Branch("dijet_mass",&dijet_mass);                  //17.Dijet Invariant Mass
  outputTree->Branch("dijet_pt",&dijet_pt);                      //18.Dijet System PT
  outputTree->Branch("dijet_mt",&dijet_mt);                      //19.Dijet System MT
  
  //AngularVariable
  outputTree->Branch("deltaR_jet01",&deltaR_jet01);              //20.deltaR(Leading Jet,Subleading Jet)
  outputTree->Branch("deltaPhi_metjet0",&deltaPhi_metjet0);      //21.deltaPhi(MET,Leading Jet)
  outputTree->Branch("deltaPhi_metjet1",&deltaPhi_metjet1);      //22.deltaPhi(MET,Subleading Jet)
  outputTree->Branch("deltaPhi_metlep0",&deltaPhi_metlep0);      //23.deltaPhi(MET,Leading Lepton)
  outputTree->Branch("deltaPhi_jet0lep0",&deltaPhi_jet0lep0);    //24.deltaPhi(Leading Jet,Leading Lepton)
  outputTree->Branch("deltaPhi_jet1lep0",&deltaPhi_jet1lep0);    //25.deltaPhi(Subleading Jet,Leading Lepton)
  outputTree->Branch("deltaPhi_dijetlep0",&deltaPhi_dijetlep0);  //26.deltaPhi(Dijet System,Leading Lepton)
  outputTree->Branch("deltaPhi_metdijet",&deltaPhi_metdijet);    //27.deltaPhi(Dijet System,MET)
  
  //EventLevel
  outputTree->Branch("event_MET",&event_MET);                    //28.MET
  outputTree->Branch("event_METPhi",&event_METPhi);              //28.METPhi
  outputTree->Branch("event_HT",&event_HT);                      //29.HT
  outputTree->Branch("event_ST",&event_ST);                      //30.ST
  outputTree->Branch("n_Jet",&n_Jet);                            //31.No of Jets
  outputTree->Branch("n_bJet",&n_bJet);                          //32.No of bjets    
  
  //Added extra variables
  outputTree->Branch("ljjsys_PT",&ljjsys_PT);                    //25.ljjsys_PT
  outputTree->Branch("deltaPhi_ljjsysmet",&deltaPhi_ljjsysmet);  //25.deltaPhi(ljjsys,MET)
  outputTree->Branch("deltaPhi_ljjsyslep0",&deltaPhi_ljjsyslep0);//26.deltaPhi(ljjsys,Lepton)
  outputTree->Branch("deltaPhi_ljjsysjet0",&deltaPhi_ljjsysjet0);//27.deltaPhi(ljjsys,Ledaing Jet)
  outputTree->Branch("deltaPhi_ljjsysjet1",&deltaPhi_ljjsysjet1);//28.deltaPhi(ljjsys,Subleading Jet)
  outputTree->Branch("ljjsys_mass",&ljjsys_mass);                //29.ljj system mass
  //outputTree->Branch("dijet_deltaeta",&dijet_deltaeta);          //30.deltaEta(leading jet,subleading jet)
  outputTree->Branch("deepjetQG_jet0",&deepjetQG_jet0);          //31.deepjetQG_leading jet
  outputTree->Branch("deepjetQG_jet1",&deepjetQG_jet1);          //32.deepjetQG_subleading jet
  //outputTree->Branch("npart_jet0",&npart_jet0);                  //33.nparticle_leading jet
  //outputTree->Branch("npart_jet1",&npart_jet1);                  //34.nparticle_subleading jet
  //outputTree->Branch("avgnpart",&avgnpart);                      //35.average no of particle in jets
  //outputTree->Branch("event_avgCvsBscore",&event_avgCvsBscore);
  //outputTree->Branch("event_avgCvsLscore",&event_avgCvsLscore);
  outputTree->Branch("event_avgQGscore",&event_avgQGscore);
  outputTree->Branch("event_Rpt",&event_Rpt);
  //outputTree->Branch("event_zstar",&event_zstar);

  //SF
  outputTree->Branch("event_lepIDSF",                 &event_lepIDSF);
  outputTree->Branch("event_lepIDSF_up",              &event_lepIDSF_up);
  outputTree->Branch("event_lepIDSF_down",            &event_lepIDSF_down);
  outputTree->Branch("event_lepISOSF",                &event_lepISOSF);
  outputTree->Branch("event_lepISOSF_up",             &event_lepISOSF_up);
  outputTree->Branch("event_lepISOSF_down",           &event_lepISOSF_down);
  outputTree->Branch("event_btagsf",                  &event_btagsf);              //btag SF for the event(or weight)
  outputTree->Branch("event_btagsfCorrelated_up",     &event_btagsfCorrelated_up);
  outputTree->Branch("event_btagsfCorrelated_down",   &event_btagsfCorrelated_down);
  outputTree->Branch("event_btagsfUncorrelated_up",   &event_btagsfUncorrelated_up);
  outputTree->Branch("event_btagsfUncorrelated_down", &event_btagsfUncorrelated_down);
  outputTree->Branch("event_btagsfbcCorrelated_up",   &event_btagsfbcCorrelated_up);
  outputTree->Branch("event_btagsfbcCorrelated_down", &event_btagsfbcCorrelated_down);
  outputTree->Branch("event_btagsfbcUncorrelated_up", &event_btagsfbcUncorrelated_up);
  outputTree->Branch("event_btagsfbcUncorrelated_down",&event_btagsfbcUncorrelated_down);

  //trigger info
  outputTree->Branch("event_passmuTrigger",                 &event_passmuTrigger);
  outputTree->Branch("event_passeleTrigger",                 &event_passeleTrigger);

  //HEM(2018)
  outputTree->Branch("event_HEMVeto",&event_HEMVeto);
  
  //onlyMC
  outputTree->Branch("GenWeight",&GenWeight);
  outputTree->Branch("PU_N",&PU_N);
  outputTree->Branch("PU_TrueN",&PU_TrueN);
  outputTree->Branch("PileUpWt_nom", &PileUpWt_nom);
  outputTree->Branch("PileUpWt_up", &PileUpWt_up);
  outputTree->Branch("PileUpWt_down", &PileUpWt_down);  
  outputTree->Branch("L1PreFireWt_nom",&L1PreFireWt_nom);
  outputTree->Branch("L1PreFireWt_up",&L1PreFireWt_up);
  outputTree->Branch("L1PreFireWt_down",&L1PreFireWt_down);
  outputTree->Branch("PDF_up", &PDF_up);
  outputTree->Branch("PDF_down", &PDF_down);
  outputTree->Branch("qcd_scale_up", &qcd_scale_up);
  outputTree->Branch("qcd_scale_down", &qcd_scale_down);

  
  //push in fileVector
  outputFileVector.push_back(newFile);
  
  //print
  cout<<"Created file: "<<outputFileName<<endl;
}
		     
