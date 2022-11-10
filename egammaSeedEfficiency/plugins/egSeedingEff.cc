#ifndef EGSEEDINGEFF_H
#define EGSEEDINGEFF_H

// Standard C++ includes
#include <memory>
#include <vector>
#include <iostream>

// ROOT includes
#include "TTree.h"
#include "TLorentzVector.h"

// CMSSW framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 

// CMSSW data formats
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
//#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/TrackReco/interface/TrackBase.h" 
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


class egSeedingEff : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
	public:
		explicit egSeedingEff(const edm::ParameterSet&);
		~egSeedingEff();
		
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
	private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
        virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
        
        const edm::EDGetTokenT<std::vector<reco::Electron>>  electronToken;
		const edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken;

        TTree* tree;

		std::vector<float>  electron_pt;
		std::vector<float>  electron_eta;
		std::vector<float>  electron_phi;

		int run_, lumi_, event_;

};

//Constructor
egSeedingEff::egSeedingEff(const edm::ParameterSet& iConfig): 
				electronToken  (consumes<std::vector<reco::Electron> >(iConfig.getParameter<edm::InputTag>("electron"))),
				genParticlesToken  (consumes<reco::GenParticleCollection> (iConfig.getParameter<edm::InputTag>("genParticles")))

{
	usesResource("TFileService");	
}

//Destructor
egSeedingEff::~egSeedingEff() {}

void egSeedingEff::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}

void egSeedingEff::beginJob() 
{

	// Access the TFileService
	edm::Service<TFileService> fs; 
  
    // Create the TTree
    tree = fs->make<TTree>("tree"  , "tree");

    tree->Branch("event"  ,&event_ , "event/I");   
    tree->Branch("lumi"   ,&lumi_  , "lumi/I");  
    tree->Branch("run"    ,&run_   , "run/I");    

    tree->Branch("electron_pt" , "std::vector<float>", &electron_pt  , 32000, 0);
    tree->Branch("electron_eta", "std::vector<float>", &electron_eta , 32000, 0);
    tree->Branch("electron_phi", "std::vector<float>", &electron_phi , 32000, 0);
  
    return ;
}

void egSeedingEff::endJob() {}
void egSeedingEff::endRun(edm::Run const&, edm::EventSetup const&) {}

void egSeedingEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    using namespace edm;
    using namespace std;
    
	run_ = 0;
	lumi_= 0;
	event_ = 0;

	electron_pt.clear();
	electron_eta.clear();
	electron_phi.clear();


    edm::Handle<vector<reco::Electron>> electronH;
    iEvent.getByToken(electronToken, electronH);

    edm::Handle<reco::GenParticleCollection> genParticlesH;
    iEvent.getByToken(genParticlesToken, genParticlesH);
 

	//-------------- Get Event Info -----------------------------------

	run_    = iEvent.id().run();
	event_  = iEvent.id().event();
	lumi_   = iEvent.id().luminosityBlock();


	for (auto genItr = genParticlesH->begin(); genItr != genParticlesH->end(); ++genItr) 
	{
		std::cout<<" Testing gen variables "<< genItr->pdgId() <<std::endl;
	}

	for (auto eleItr = electronH->begin(); eleItr != electronH->end(); ++eleItr) 
	{	
		electron_pt.push_back( eleItr->pt() );
		electron_eta.push_back( eleItr->eta() );
		electron_phi.push_back( eleItr->phi() );
	} 

    tree->Fill();	
}


void egSeedingEff::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {}
void egSeedingEff::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}
void egSeedingEff::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(egSeedingEff);

#endif 
