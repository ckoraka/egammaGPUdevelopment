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

// CMSSW DataFormats
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/TrackBase.h" 
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

// CMSSW sim DataFormats
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"

// Other includes 
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetType.h"

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
		virtual void initialize();

		const edm::EDGetTokenT<reco::ElectronCollection>  electronToken;
		const edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken;
		const edm::EDGetTokenT<std::vector<SimTrack>> simtracksToken;
		const edm::EDGetTokenT<edm::View<TrackingParticle>> trackingParticlesToken;		
		const edm::EDGetTokenT<edm::PSimHitContainer> tokPixelBarrelHitsLowTof_;
		const edm::EDGetTokenT<edm::PSimHitContainer> tokPixelBarrelHitsHighTof_;
		const edm::EDGetTokenT<edm::PSimHitContainer> tokPixelEndcapHitsLowTof_;
		const edm::EDGetTokenT<edm::PSimHitContainer> tokPixelEndcapHitsHighTof_;
		const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
		const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;

		TTree* tree;

		std::vector<float>  electron_pt;
		std::vector<float>  electron_eta;
		std::vector<float>  electron_phi;

		int run_, lumi_, event_;
		bool verbose_;

};

//Constructor
egSeedingEff::egSeedingEff(const edm::ParameterSet& iConfig): 
							electronToken  (consumes<reco::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electron"))),
							genParticlesToken  (consumes<reco::GenParticleCollection> (iConfig.getParameter<edm::InputTag>("genParticles"))),
							trackingParticlesToken (consumes<edm::View<TrackingParticle>>(iConfig.getParameter<edm::InputTag>("trackingParticles"))),
							tokPixelBarrelHitsLowTof_(consumes<edm::PSimHitContainer>(edm::InputTag("g4SimHits", "TrackerHitsPixelBarrelLowTof"))),
							tokPixelBarrelHitsHighTof_(consumes<edm::PSimHitContainer>(edm::InputTag("g4SimHits", "TrackerHitsPixelBarrelHighTof"))),
							tokPixelEndcapHitsLowTof_(consumes<edm::PSimHitContainer>(edm::InputTag("g4SimHits", "TrackerHitsPixelEndcapLowTof"))),
							tokPixelEndcapHitsHighTof_(consumes<edm::PSimHitContainer>(edm::InputTag("g4SimHits", "TrackerHitsPixelEndcapHighTof"))),
							topoToken_(esConsumes()),
						    geomToken_(esConsumes()),
							verbose_(iConfig.getParameter<bool>("verbose"))
{
	initialize();
	usesResource("TFileService");	
}

//Destructor
egSeedingEff::~egSeedingEff() {}

void egSeedingEff::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}

void egSeedingEff::initialize() {
	this -> run_ = 0;
	this -> lumi_= 0;
	this -> event_ = 0;

	this -> electron_pt.clear();
	this -> electron_eta.clear();
	this -> electron_phi.clear();
}

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

	const TrackerTopology* tTopo = &iSetup.getData(topoToken_);
	const TrackerGeometry* tGeom = &iSetup.getData(geomToken_);

	edm::Handle<reco::ElectronCollection> electronH;
	iEvent.getByToken(electronToken, electronH);
	edm::Handle<reco::GenParticleCollection> genParticlesH;
	iEvent.getByToken(genParticlesToken, genParticlesH);
	edm::Handle<edm::View<TrackingParticle>> TrackingParticleH;
	iEvent.getByToken(trackingParticlesToken, TrackingParticleH);
	const edm::View<TrackingParticle>& trackingParticles = *TrackingParticleH;

	//-------------- Event Info -----------------------------------
	run_    = iEvent.id().run();
	event_  = iEvent.id().event();
	lumi_   = iEvent.id().luminosityBlock();


	//-------------- Gen particle info -----------------------------------
	for (auto genItr = genParticlesH->begin(); genItr != genParticlesH->end(); ++genItr) 
	{
		// Needs refinement - idea is to find the electrons originating fromt the Z boson
		if(abs(genItr->pdgId())==11 && genItr->status()==1){
			auto daughters = genItr->daughterRefVector();
			auto mothers = genItr->motherRefVector();
			std::cout << " Particle mother pdgid & status : "<< mothers.at(0)->pdgId() <<" , "<<mothers.at(0)->status() <<std::endl;
			//if (daughters.size()==1 && abs(daughters.at(0)->pdgId())==11)	
		}
	}

	//-------------- Sim track info -----------------------------------
	// Unfortunately no simHit attached to sim track - cannot extract info on how may pixel layers have hits :/
	// At least this is my current understanding 
	// https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h
	for (unsigned long ntrackingparticle = 0; ntrackingparticle < trackingParticles.size(); ntrackingparticle++) 
	{
		const auto& tp = trackingParticles.at(ntrackingparticle);
		bool isElectron = (std::abs(tp.pdgId()) == 11);

      	if(!(tp.eventId().bunchCrossing() == 0 && tp.eventId().event() == 0)){continue;}
		if(tp.charge() == 0 || tp.status() != 1) {continue;}

		if(isElectron)
			std::cout << " Number of layers  "<< tp.numberOfTrackerLayers() << std::endl;

    	//for (std::vector<SimTrack>::const_iterator itrk = tp.g4Track_begin();itrk != tp.g4Track_end();++itrk) 
		//{
		//	std::cout << "  id " << itrk->trackId() <<std::endl;
		//	std::cout << " Track Momentum  " << itrk->momentum() << std::endl;
		//	std::cout << " Track ID & type " << itrk->trackId() << " " << itrk->type() << std::endl;
		//}
		
		if(verbose_)
			std::cout<<" TrackingParticle with PdgId: "<<tp.pdgId()<<" 4-momentum :("<<tp.p4().pt() <<","<<tp.p4().eta()<<","<<tp.p4().phi() <<","<< tp.p4().e()<<")"<<std::endl;
	}

	//-------------- simhit info ---------------------------------------

	edm::Handle<edm::PSimHitContainer> PixelBarrelHitsLowTof;
	iEvent.getByToken(tokPixelBarrelHitsLowTof_, PixelBarrelHitsLowTof);
	edm::Handle<edm::PSimHitContainer> PixelBarrelHitsHighTof;
	iEvent.getByToken(tokPixelBarrelHitsHighTof_, PixelBarrelHitsHighTof);
	edm::Handle<edm::PSimHitContainer> PixelEndcapHitsLowTof;
	iEvent.getByToken(tokPixelEndcapHitsLowTof_, PixelEndcapHitsLowTof);
	edm::Handle<edm::PSimHitContainer> PixelEndcapHitsHighTof;
	iEvent.getByToken(tokPixelEndcapHitsHighTof_, PixelEndcapHitsHighTof);


	for (unsigned long simHitCounter = 0; simHitCounter < PixelBarrelHitsLowTof->size(); ++simHitCounter) 
	{
    	const PSimHit& simHit = (*PixelBarrelHitsLowTof)[simHitCounter];
		std::cout<< " simHit.trackId() "<<simHit.trackId()<<std::endl;
		std::cout<< " simHit.detUnitId() "<<simHit.detUnitId()<<std::endl;
		std::cout<< " simHit.processType() "<<simHit.processType()<<std::endl;
		std::cout<< " simHit.particleType() "<<simHit.particleType()<<std::endl;
	}


	//-------------- hltGsfElectrons -----------------------------------
	if(electronH.isValid()) 
	{
		for (auto eleItr = electronH->begin(); eleItr != electronH->end(); ++eleItr) 
		{	
			auto gsfTrack = eleItr->gsfTrack();
			auto seed = gsfTrack->seedRef();
			auto rhits = seed->recHits();

			// Check number of hits 
            std::cout<<" Number of hits per seed : "<< seed->nHits() <<std::endl;	

			// Check hits per Pixel layer
			for(auto const& rhit: rhits){
				if(rhit.isValid() && rhit.det() != nullptr)
				{
					if(verbose_)
						std::cout<<" rechits "<< rhit.isValid()<<std::endl;

                    DetId det = rhit.geographicalId();					
					auto subdet = rhit.det()->geographicalId().subdetId();
                    auto trkSubDet = tGeom->geomDetSubDetector(subdet);

					if(GeomDetEnumerators::isTrackerPixel(trkSubDet))
					{
                        if(subdet == PixelSubdetector::PixelBarrel)
						{
                            std::cout << "Barrel Layer: " << tTopo->pxbLayer(det) << std::endl;
                            //bpixlayer->push_back(ttopo.pxbLayer(det));
                        }
                        if(subdet == PixelSubdetector::PixelEndcap)
						{
                            std::cout << "Endcap disk: " << tTopo->pxfDisk(det) << std::endl;
                            //fpixlayer->push_back(tTopo.pxfDisk(det));
                        }
                    }
				}
			}
			electron_pt.push_back( eleItr->pt() );
			electron_eta.push_back( eleItr->eta() );
			electron_phi.push_back( eleItr->phi() );
		} 
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
