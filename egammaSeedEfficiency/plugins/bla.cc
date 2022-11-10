#ifndef NTUPLIZER_H
#define NTUPLIZER_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include <TNtuple.h>
#include <TString.h>
#include <bitset>
#include <typeinfo>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Utilities/interface/InputTag.h>



#include <DataFormats/PatCandidates/interface/Electron.h>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/TriggerResults.h"



#include "tParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"


//Set this variable to decide the number of triggers that you want to check simultaneously
#define NUMBER_OF_MAXIMUM_TRIGGERS 64


class Ntuplizer : public edm::EDAnalyzer {

	public:
		explicit Ntuplizer(const edm::ParameterSet&);
		virtual ~Ntuplizer();

	private:
		//----edm control---
		virtual void beginJob() ;
		virtual void beginRun(edm::Run const&, edm::EventSetup const&);
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob();
		virtual void endRun(edm::Run const&, edm::EventSetup const&);
		void Initialize();

		TTree *_L1EGTree;			
		std::string _treeName;

		// -------------------------------------
		// variables to be filled in output tree
		ULong64_t       _indexevents;
		Int_t           _runNumber;
		Int_t           _lumi;
		unsigned long _eleProbeTriggerBits;
		unsigned long _eleTagTriggerBits;


		//L1EG
        float _l1objPt;
        float _l1objEta; 
        float _l1objPhi;
        float _matchedrechitE; 
        float _matchedrechitEta;
        float _matchedrechitPhi;
        float _matchedrechitTime;
        float _DR_rechit_L1; 
        int   _severity;
        int   _isSpike;

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          
		edm::EDGetTokenT<edm::View<reco::GsfElectron> >  _electronsTag;
		edm::EDGetTokenT<edm::ValueMap<bool> > _eleMediumIdMapTag;
		edm::EDGetTokenT<edm::ValueMap<bool> > _eleLooseIdMapTag;

		edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> _triggerObjects;
		edm::EDGetTokenT< edm::TriggerResults> _triggerBits;

        edm::EDGetTokenT< l1t::EGammaBxCollection > _L1EGTag;
        edm::EDGetTokenT< l1t::EGammaBxCollection > _L1EmuEGTag;
        edm::EDGetTokenT<std::vector<reco::Vertex>> _VtxTag;
        edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> EBRecHitsToken_; 
        edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> EERecHitsToken_; 

		bool _useHLTMatch;

		//!Contains the parameters
		tVParameterSet _parametersTag;
		tVParameterSet _parametersProbe;

		edm::InputTag _processName;
		//! Maximum
		std::bitset<NUMBER_OF_MAXIMUM_TRIGGERS> _eleProbeTriggerBitSet;
		std::bitset<NUMBER_OF_MAXIMUM_TRIGGERS> _eleTagTriggerBitSet;

		HLTConfigProvider _hltConfig;


};

// ----Constructor and Destructor -----
Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig) :
	_electronsTag       (consumes<edm::View<reco::GsfElectron> >       (iConfig.getParameter<edm::InputTag>("electrons"))),
	_eleMediumIdMapTag  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
	_eleLooseIdMapTag  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
	_triggerObjects (consumes<pat::TriggerObjectStandAloneCollection> (iConfig.getParameter<edm::InputTag>("triggerSet"))),
	_triggerBits    (consumes<edm::TriggerResults>                    (iConfig.getParameter<edm::InputTag>("triggerResultsLabel"))),
	_L1EGTag       (consumes<l1t::EGammaBxCollection>                   (iConfig.getParameter<edm::InputTag>("L1EG"))),
	_L1EmuEGTag    (consumes<l1t::EGammaBxCollection>                   (iConfig.getParameter<edm::InputTag>("L1EmuEG"))),
	_VtxTag         (consumes<std::vector<reco::Vertex>>              (iConfig.getParameter<edm::InputTag>("Vertices"))),
	EBRecHitsToken_ (consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(iConfig.getParameter<edm::InputTag>("recHitsEB"))),	 
        EERecHitsToken_ (consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(iConfig.getParameter<edm::InputTag>("recHitsEE"))),      
        CaloTowerToken_ (consumes<l1t::CaloTowerBxCollection>	(iConfig.getParameter<edm::InputTag>("calotower")))
{
	this -> _treeName = iConfig.getParameter<std::string>("treeName");
	this -> _processName = iConfig.getParameter<edm::InputTag>("triggerResultsLabel");
	this -> _useHLTMatch = iConfig.getParameter<bool>("useHLTMatch");

	this -> Initialize();
	return;
}

Ntuplizer::~Ntuplizer()
{}

void Ntuplizer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
	Bool_t changedConfig = false;

	if(!this -> _hltConfig.init(iRun, iSetup, this -> _processName.process(), changedConfig)){
		edm::LogError("HLTMatchingFilter") << "Initialization of HLTConfigProvider failed!!";
		return;
	}

	const edm::TriggerNames::Strings& triggerNames = this -> _hltConfig.triggerNames();

}

void Ntuplizer::Initialize() {

    this -> _indexevents = 0;
    this -> _runNumber = 0;
    this -> _lumi = 0;
    this -> _DR_rechit_L1      = 999; 		 
}


void Ntuplizer::beginJob()
{
    edm::Service<TFileService> fs;
    TString Spikes;
    this -> _L1EGTree = fs -> make<TTree>("L1EG", "L1EG");
    
    TString SpikesEmu;
    this -> _L1EGTreeEmu = fs -> make<TTree>("EmulatedL1EG", "EmulatedL1EG");  

    this -> _L1EGTree -> Branch("EventNumber",&_indexevents,"EventNumber/l");
    this -> _L1EGTree -> Branch("RunNumber",&_runNumber,"RunNumber/I");
    this -> _L1EGTree -> Branch("lumi",&_lumi,"lumi/I");

	return;
}


void Ntuplizer::endJob()
{
	return;
}


void Ntuplizer::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
	return;
}


void Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	
        ////////	 Define the tag electrons of the event	//////////////
	
	edm::Handle<edm::View<reco::GsfElectron> > electrons;
	edm::Handle<edm::View<reco::GenParticle> > genParticles;
	edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
	edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
	edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
	edm::Handle<edm::TriggerResults> triggerBits;
	edm::Handle<std::vector<reco::Vertex> >  vertices;

	//////////////
	iEvent.getByToken(this -> _electronsTag, electrons);
	iEvent.getByToken(this -> _eleMediumIdMapTag, medium_id_decisions);
	iEvent.getByToken(this -> _eleLooseIdMapTag, loose_id_decisions);
	
	if(this->_useHLTMatch)
		iEvent.getByToken(this -> _triggerObjects, triggerObjects);

	iEvent.getByToken(this -> _triggerBits, triggerBits);
	iEvent.getByToken(this -> _VtxTag,vertices);

	const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
	

	//ECAL rechits Handles
 	edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> EBRecHits;	//Barrel
	iEvent.getByToken(EBRecHitsToken_, EBRecHits);
        edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> EERecHits;	//Endcap
        iEvent.getByToken(EERecHitsToken_, EERecHits);

	//L1 Handles
        edm::Handle< BXVector<l1t::EGamma> >  L1EGHandle;
        iEvent.getByToken(_L1EGTag, L1EGHandle);

        edm::Handle< BXVector<l1t::EGamma> >  L1EmuEGHandle;
        try {iEvent.getByToken(_L1EmuEGTag, L1EmuEGHandle);} catch (...) {;}
	//iEvent.getByToken(_L1EmuEGTag, L1EmuEGHandle);

	//get Calorimeter geometry
	edm::ESHandle<CaloGeometry> geometry;
	iSetup.get<CaloGeometryRecord>().get(geometry);
	const CaloGeometry* geo = geometry.product();

	//get severity level
	edm::ESHandle<EcalSeverityLevelAlgo> sevlv;                     
        iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevlv);              
        const EcalSeverityLevelAlgo* sevLevel = sevlv.product(); 


	//	Initialize arrays used to identify matched unpacked/Emulated L1EGs with robe electrons	///////////

	for (int i =0; i<100; ++i){
		isL1EGMatchedWithProbe[i]    = 0;
		isL1EGEmuMatchedWithProbe[i] = 0;
	}

	////////////////	Select Unpacked and Emulated  L1EG candidates via Tag and Probe 	////////////////////////////

	for (unsigned int i = 0; i< electrons->size(); ++i){

		const auto eleTag = electrons->ptrAt(i);
		int isTagIDMedium = (*medium_id_decisions)[eleTag];
		//cout<<"Is tag medium = "<< isTagIDMedium <<endl;
		if(!isTagIDMedium || eleTag->p4().Pt()<30) continue;


		for (unsigned int j = 0; j< electrons->size(); ++j){

			if(i==j) continue; //Tag and probe not the same electron


                        //this -> Initialize();
                        _indexevents = iEvent.id().event();
                        _runNumber = iEvent.id().run();
                        _lumi = iEvent.luminosityBlock();
                                               
			const auto eleProbe = electrons->ptrAt(j);
			int isProbeIDLoose = (*loose_id_decisions)[eleProbe];
			if(!isProbeIDLoose) continue;   

			// In barrel or endcap volume
			float eleProbeEta = eleProbe->p4().Eta();
	                if((abs(eleProbeEta)>1.4442 && abs(eleProbeEta)<1.566)) continue;

			// Tag and probe electron opposite sign
			bool isOS = (eleTag->charge() / eleProbe->charge() < 0) ? true : false;
			if(!isOS) continue;

			// Tag and probe invariant mass consisted with Z mass
			float Mee = (eleTag->p4() + eleProbe->p4()).M();
			if(!(Mee>60 && Mee<120)) continue;
			//cout<<" Mee = "<<Mee<<" is OS = "<<isOS<<endl;

                        ////////        Unpacked L1EG matching with probe electron      ////////////////////

			float minDR = 0.3; 
			int k = 0;
			int ID = -1;
			for (l1t::EGammaBxCollection::const_iterator bx0EGIt = L1EGHandle->begin(0); bx0EGIt != L1EGHandle->end(0) ; bx0EGIt++)
			{
				const float dR = deltaR(*eleProbe, *bx0EGIt);
				const l1t::EGamma& l1tEG = *bx0EGIt;
				
				if (dR < minDR) 
				{
					minDR = dR; 
					ID = k;
				}
			++k;	//Value of constant iterator
			}
			
			if(minDR != 0.3 && ID != -1) isL1EGMatchedWithProbe[ID] = 1;
			
			//cout<<" probe pT = "<<eleProbe->pt()<<" eta = "<<eleProbe->eta()<<" phi = "<< eleProbe->phi()<<endl;
			//cout<<" It got matched with L1EG with ID "<<ID<<endl;
			//if(ID!=-1){
			//	cout<<" Has it got the same pT = "<<(L1EGHandle->at(0,ID)).pt()<<endl;
                	//        cout<<" Has it got the same eta = "<<(L1EGHandle->at(0,ID)).eta()<<endl;
                        //        cout<<" Has it got the same phi = "<<(L1EGHandle->at(0,ID)).phi()<<endl;
			//}

			////////	Emulated L1EG matching with probe electron	////////////////////

			if (L1EmuEGHandle.isValid())
			{
				minDR = 0.3;
	                        k = 0;
        	                ID = -1;

				for (l1t::EGammaBxCollection::const_iterator bx0EmuEGIt = L1EmuEGHandle->begin(0); bx0EmuEGIt != L1EmuEGHandle->end(0) ; bx0EmuEGIt++)
				{
					const float dR = deltaR(*eleProbe, *bx0EmuEGIt);
					const l1t::EGamma& l1tEmuEG = *bx0EmuEGIt;
                                          
					if (dR < minDR) 
					{
						minDR = dR; 
	                                        ID = k;
					}
       	
 		                ++k;    //Value of constant iterator
	
				}
			if(minDR != 0.3 && ID != -1) isL1EGEmuMatchedWithProbe[ID] = 1;
			}

                        //////////////////////////////////////////////////////////////////////////////////

		}
	}
	
	////	Cross check		/////
	int k = 0;
	for (l1t::EGammaBxCollection::const_iterator bx0EGIt = L1EGHandle->begin(0); bx0EGIt != L1EGHandle->end(0); bx0EGIt++)
        {
		++count;
		if(isL1EGMatchedWithProbe[k] == 1) ++countM;
		++k;
	}
	k = 0;
        for (l1t::EGammaBxCollection::const_iterator bx0EmuEGIt = L1EmuEGHandle->begin(0); bx0EmuEGIt != L1EmuEGHandle->end(0); bx0EmuEGIt++ )
        {
        	++countEmu;
		if(isL1EGEmuMatchedWithProbe[k] == 1) ++countEmuM;
		++k;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////		L1EG and RecHit matching Algo 		//////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//	First Step :  Match the Rechits identified as Spikes with a L1EG candidate

	//Here the matching is performed by looping over rechits that are ECAL spikes and finding the closest unpacked/emulated L1EG candidate

        float min_DR;
	int ID = -1;
	int j;

	for(int j= 0;j<100;++j){
		l1[j]    = 0;
		l1emu[j] = 0;
	}
	 

	EcalRecHitCollection::const_iterator start;
        EcalRecHitCollection::const_iterator end;
                                        
        for (int i=0; i<2; ++i){
                 if(i==0) {
			// Barrel recHits
                 	start = EBRecHits->begin();
                        end = EBRecHits->end();
                 }
                                        
                 if(i==1) {
			// EndCap recHits
                        start = EERecHits->begin();
                        end = EERecHits->end();
                 }
                         
                                
                 for (EcalRecHitCollection::const_iterator hit = start; hit != end; ++hit) {

			_severity = sevlv->severityLevel(*hit);		//Severity level 3 OR 4 corresponds to Spikes
			if((_severity == 3)||(_severity == 4)){

				if(hit->energy()<4.) continue;		// Cut off RecHits with energy smaller than 4 GeV

        		/////////////////       L1EG unpacked Data              ////////////////////////////////////////

 				min_DR = 999;
				j = -1;
                                for (l1t::EGammaBxCollection::const_iterator bx0EGIt = L1EGHandle->begin(0); bx0EGIt != L1EGHandle->end(0); bx0EGIt++)
                                {
					++j;
					if(isL1EGMatchedWithProbe[j]!=1)	continue;
										
			
                                        const l1t::EGamma& l1tEG = *bx0EGIt;
        				//cout<<"Size = "<<L1EGHandle->size(0)<<endl;                                
                                        float dphi = (l1tEG.phi())-(geo->getPosition(hit->detid()).phi()); //Takes into account the 2PI periodicity
                                        float DR = sqrt(pow(l1tEG.eta()-geo->getPosition(hit->detid()).eta(),2)+pow(dphi,2));

                                        if ((DR < min_DR) && (DR < 0.1)){
                                                min_DR = DR;
                                                this -> _matchedrechitE = hit->energy();
                                                this -> _DR_rechit_L1 = min_DR;
                                                this -> _matchedrechitEta = geo->getPosition(hit->detid()).eta();
                                                this -> _matchedrechitPhi = geo->getPosition(hit->detid()).phi();
                                                this -> _severity   = sevlv->severityLevel(*hit);
                                                this -> _matchedrechitTime = hit->time();
                                                this -> _l1objPt    = l1tEG.pt();
                                                this -> _l1objEta   = l1tEG.eta();
                                                this -> _l1objPhi   = l1tEG.phi();
						this -> _Nvtx 	    = vertices->size();
						this -> _isSpike = 1;
						ID = j;
						//cout<<" I am matched ! "<<endl;
						//cout<<"my pT is = "<<l1tEG.pt()<<" and my eta is= "<<l1tEG.eta()<< " and ID  "<<ID<<endl;
                        		}
                                }        
				
				if(min_DR != 999){
					//cout<<"This L1EG got mathced with an ECAL Spike."<<endl;
					//cout<<"Its pT is = "<<(L1EGHandle->at(0,ID)).pt()<<" and eta = "<<(L1EGHandle->at(0,ID)).eta()<<endl;
        				this -> _L1EGTree -> Fill();
					l1[ID] = 1;	// In this array we keep the ID of the L1EG candidate that got mathched with a Spike
				}

	
			

        /////////////////       L1EG Emulated Data              ////////////////////////////////////////

		        	min_DR = 999;
				j = -1;
                		for (l1t::EGammaBxCollection::const_iterator bx0EmuEGIt = L1EmuEGHandle->begin(0); bx0EmuEGIt != L1EmuEGHandle->end(0); bx0EmuEGIt++ )
                		{
					++j;

					if(isL1EGEmuMatchedWithProbe[j]!=1)	continue;
					
                			const l1t::EGamma& l1tEmuEG = *bx0EmuEGIt;

	 				float dphi = (l1tEmuEG.phi())-(geo->getPosition(hit->detid()).phi());
                                	float DR = sqrt(pow(l1tEmuEG.eta()-geo->getPosition(hit->detid()).eta(),2)+pow(dphi,2));

                                	if ((DR < min_DR) && (DR < 0.1)){
                                        	min_DR = DR;
                                        	this -> _DR_Emu_rechit_L1       = min_DR;
                                        	this -> _matchedrechitEmuE      = hit->energy();
                                        	this -> _severityEmu            = sevlv->severityLevel(*hit);
                                        	this -> _matchedrechitEmuTime   = hit->time();
                                        	this -> _l1EmuobjPt             = l1tEmuEG.pt();
                                        	this -> _l1EmuobjEta            = l1tEmuEG.eta();
                                        	this -> _l1EmuobjPhi            = l1tEmuEG.phi();
                                        	this -> _l1EmuobjnTT            = l1tEmuEG.nTT();
                                        	this -> _l1EmuobjShape          = l1tEmuEG.shape();
                                        	this -> _l1EmuobjQual           = l1tEmuEG.hwQual();
                                        	this -> _l1EmuobjIso            = l1tEmuEG.hwIso();
                                        	this -> _l1EmuobjHoE            = l1tEmuEG.towerHoE();
						this -> _NvtxEmu 		= vertices->size();
						this -> _isSpikeEmu = 1;
						ID = j;
                                 	}
				
                           	}

				if(min_DR != 999){
					this -> _L1EGTreeEmu -> Fill();
					l1emu[ID] = 1;
				}
	////////////////////////////////////////////////////////////////////////////////////////////////////

			}

		}

	}


//	Second Step :  Match the remaining L1EG candidates with their closest RecHit (that is not a spike)

//	Here the matching is performed by looping over the remaining unpacked/emulated L1EGs and finding the closest RecHit that is not an ECAL spike ---> We do this in order to compare the rechit distributions between L1EG candidates that are or are not spikes. 

	
        /////////////////       L1EG unpacked Data              ////////////////////////////////////////

	j = -1;
        for (l1t::EGammaBxCollection::const_iterator bx0EGIt = L1EGHandle->begin(0); bx0EGIt != L1EGHandle->end(0); bx0EGIt++)
        {
        	const l1t::EGamma& l1tEG = *bx0EGIt;
		++j;
		//if(l1[j]==1){
		//	cout<<" I found a matched L1EG with an ID = "<< j <<endl;
		//	cout<<" Has it got the same pT = "<<(L1EGHandle->at(0,j)).pt()<<endl;
		//	cout<<" Has it got the same eta = "<<(L1EGHandle->at(0,j)).eta()<<endl;
		//}

		min_DR = 999;

		if( (l1[j]==1) || (isL1EGMatchedWithProbe[j]!=1) ) continue;	// If the L1EG candidate has already been matched with an ECAL spike OR has not been matched with a probe electron, do not consider it


                for (int i=0; i<2; ++i){
                	if(i==0) {
				// Barrel recHits
                 		start = EBRecHits->begin();
                        	end = EBRecHits->end();
                 	}
                                        
                 	if(i==1) {
				// EndCap recHits
                        	start = EERecHits->begin();
                        	end = EERecHits->end();
                 	}
                          
                 	for (EcalRecHitCollection::const_iterator hit = start; hit != end; ++hit) {

				_severity = sevlv->severityLevel(*hit);		
				if((_severity != 3)&&(_severity != 4)){		//Consider Rechits that are not ECAL spikes
			
	                                if(hit->energy()<4.) continue;          // Cut off RecHits with energy smaller than 4 GeV

		        		float dphi = (l1tEG.phi())-(geo->getPosition(hit->detid()).phi());
		                        float DR = sqrt(pow(l1tEG.eta()-geo->getPosition(hit->detid()).eta(),2)+pow(dphi,2));
			
		                        if ((DR < min_DR)&&(DR < 0.1)){
		                        	min_DR = DR;
		                                this -> _matchedrechitE = hit->energy();
		                                this -> _DR_rechit_L1 = min_DR;
		                                this -> _matchedrechitEta = geo->getPosition(hit->detid()).eta();
		                                this -> _matchedrechitPhi = geo->getPosition(hit->detid()).phi();
		                                this -> _severity   = sevlv->severityLevel(*hit);
		                                this -> _matchedrechitTime = hit->time();
		                                this -> _l1objPt    	= l1tEG.pt();
		                                this -> _l1objEta   	= l1tEG.eta();
		                                this -> _l1objPhi   	= l1tEG.phi();
						this -> _Nvtx 		= vertices->size();
						this -> _isSpike 	= 0;
		                        }
				}
			}
		}
	if ( min_DR != 999 ) this -> _L1EGTree -> Fill();
	}

        /////////////////       L1EG Emulated Data              ////////////////////////////////////////

	j = -1;
        for (l1t::EGammaBxCollection::const_iterator bx0EmuEGIt = L1EmuEGHandle->begin(0); bx0EmuEGIt != L1EmuEGHandle->end(0); bx0EmuEGIt++ )
        {

		min_DR = 999;

        	const l1t::EGamma& l1tEmuEG = *bx0EmuEGIt;
		++j;
		if( l1emu[j]==1 || isL1EGEmuMatchedWithProbe[j]!=1 ) continue;	// If the Emulated L1EG candidate has already been matched with an ECAL spike OR has not been matched with a probe electron, do not consider it
		

                for (int i=0; i<2; ++i){
                	if(i==0) {
				// Barrel recHits
                 		start = EBRecHits->begin();
                        	end = EBRecHits->end();
                 	}
                                        
                 	if(i==1) {
				// EndCap recHits
                        	start = EERecHits->begin();
                        	end = EERecHits->end();
                 	}
                          
                 	for (EcalRecHitCollection::const_iterator hit = start; hit != end; ++hit) {

				_severity = sevlv->severityLevel(*hit);		
				if((_severity != 3)&&(_severity != 4)){		//Consider Rechits that are not ECAL spikes

	                                if(hit->energy()<4.) continue;          // Cut off RecHits with energy smaller than 4 GeV
		
		 			float dphi = (l1tEmuEG.phi())-(geo->getPosition(hit->detid()).phi());
		                        float DR = sqrt(pow(l1tEmuEG.eta()-geo->getPosition(hit->detid()).eta(),2)+pow(dphi,2));                		
					if ((DR < min_DR) && (DR<0.1)){
		                                min_DR = DR;
		                                this -> _DR_Emu_rechit_L1       = min_DR;
		                                this -> _matchedrechitEmuE      = hit->energy();
		                                this -> _severityEmu            = sevlv->severityLevel(*hit);
		                                this -> _matchedrechitEmuTime   = hit->time();
		                                this -> _l1EmuobjPt             = l1tEmuEG.pt();
		                                this -> _l1EmuobjEta            = l1tEmuEG.eta();
		                                this -> _l1EmuobjPhi            = l1tEmuEG.phi();
		                                this -> _l1EmuobjnTT            = l1tEmuEG.nTT();
		                                this -> _l1EmuobjShape          = l1tEmuEG.shape();
		                                this -> _l1EmuobjQual           = l1tEmuEG.hwQual();
		                                this -> _l1EmuobjIso            = l1tEmuEG.hwIso();
		                                this -> _l1EmuobjHoE            = l1tEmuEG.towerHoE();
						this -> _NvtxEmu		= vertices->size();
						this -> _isSpikeEmu 		= 0;
		                       }
				}
			}
		}
	if ( min_DR != 999 ) this -> _L1EGTreeEmu -> Fill();
	}



  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

cout<< " Total L1EGs are "<<count<< endl;  
cout<< " Total Emulated L1EGs are "<<countEmu<<endl; 
cout<< " Total L1EGs that got matched with a probe electron are "<<countM<<endl;
cout<< " Total Emulated L1EGs that got matched with a probe electron are "<<countEmuM<<endl;


}




#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(Ntuplizer);

#endif //NTUPLIZER_H
