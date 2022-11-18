#ifndef EGSEEDINGEFF_H
#define EGSEEDINGEFF_H

/**\class egSeedingEff egSeedingEff.cc egammaGPUdevelopment/egammaSeedEfficiency/plugins/egSeedingEff.cc
 Description: Look into the number of electron SimHits in the Pixel tracker 
 Implementation:
     [Notes on implementation]
*/

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

		std::vector<float>  genEle_pt;
		std::vector<float>  genEle_eta;
		std::vector<float>  genEle_phi;
		std::vector<float>  genEle_E;

		std::vector<float>  recoEle_pt;
		std::vector<float>  recoEle_eta;
		std::vector<float>  recoEle_phi;
		std::vector<float>  recoEle_E;

		std::vector<float>  simEle_pt;
		std::vector<float>  simEle_eta;
		std::vector<float>  simEle_phi;
		std::vector<float>  simEle_E;

		std::vector<int>   nSimHitLayers;
		std::vector<int>  nSimHitsLayer1;
		std::vector<int>  nSimHitsLayer2;
		std::vector<int>  nSimHitsLayer3;
		std::vector<int>  nSimHitsLayer4;

		std::vector<int>  nRecoHitLayers;

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

	this -> genEle_pt.clear();
	this -> genEle_eta.clear();
	this -> genEle_phi.clear();
	this -> genEle_E.clear();

	this -> recoEle_pt.clear();
	this -> recoEle_eta.clear();
	this -> recoEle_phi.clear();
	this -> recoEle_E.clear();

	this -> simEle_pt.clear();
	this -> simEle_eta.clear();
	this -> simEle_phi.clear();
	this -> simEle_E.clear();

	this -> nSimHitLayers.clear();
	this -> nSimHitsLayer1.clear();
	this -> nSimHitsLayer2.clear();
	this -> nSimHitsLayer3.clear();
	this -> nSimHitsLayer4.clear();

	this -> nRecoHitLayers.clear();

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

	tree->Branch("simEle_pt" , "std::vector<float>", &simEle_pt );
	tree->Branch("simEle_eta", "std::vector<float>", &simEle_eta );
	tree->Branch("simEle_phi", "std::vector<float>", &simEle_phi );
	tree->Branch("simEle_E", "std::vector<float>", &simEle_E );
	tree->Branch("nSimHitLayers", "std::vector<int>", &nSimHitLayers  );
	tree->Branch("nSimHitsLayer1", "std::vector<int>", &nSimHitsLayer1 );
	tree->Branch("nSimHitsLayer2", "std::vector<int>", &nSimHitsLayer2 );
	tree->Branch("nSimHitsLayer3", "std::vector<int>", &nSimHitsLayer3 );
	tree->Branch("nSimHitsLayer4", "std::vector<int>", &nSimHitsLayer4 );

	tree->Branch("recoEle_pt" , "std::vector<float>", &recoEle_pt );
	tree->Branch("recoEle_eta", "std::vector<float>", &recoEle_eta );
	tree->Branch("recoEle_phi", "std::vector<float>", &recoEle_phi );
	tree->Branch("recoEle_E", "std::vector<float>", &recoEle_E );
	tree->Branch("nRecoHitLayers", "std::vector<int>", &nRecoHitLayers );

	return ;
}

void egSeedingEff::endJob() {}
void egSeedingEff::endRun(edm::Run const&, edm::EventSetup const&) {}

reco::GenParticle get_lastcopy_prefsrSN(reco::GenParticle part);
reco::GenParticle get_lastcopySN(reco::GenParticle part);

void egSeedingEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

	using namespace edm;
	using namespace std;

	const TrackerTopology* tTopo = &iSetup.getData(topoToken_);
	const TrackerGeometry* tGeom = &iSetup.getData(geomToken_);

	edm::Handle<reco::ElectronCollection> electronH;
	iEvent.getByToken(electronToken, electronH);
	edm::Handle<reco::GenParticleCollection> genParticlesH;
	iEvent.getByToken(genParticlesToken, genParticlesH);

	//-------------- Event Info -----------------------------------
	run_    = iEvent.id().run();
	event_  = iEvent.id().event();
	lumi_   = iEvent.id().luminosityBlock();

	//-------------- Gen particle info -----------------------------------
	// https://cmssdt.cern.ch/dxr/CMSSW/source/DataFormats/HepMCCandidate/interface/GenParticle.h

	std::vector<reco::GenParticle> genElectrons;

	for (auto genItr = genParticlesH->begin(); genItr != genParticlesH->end(); ++genItr) 
	{
		const auto& genPart = *genItr;
		if(abs(genItr->pdgId())==11)
		{
			if(genPart.mother()->pdgId() == 23 && genPart.isPromptFinalState())
			{
				if(genPart.status()==1)
					genElectrons.push_back(genPart);
				else if(genPart.status()==2)	
					genElectrons.push_back(get_lastcopy_prefsrSN(genPart));
				else if(genPart.status()==3)	
					genElectrons.push_back(get_lastcopySN(genPart));
			}	
			else
			{
				if(genPart.numberOfMothers() == 0)
					genElectrons.push_back(genPart);
			}
		}
	}

	if(verbose_){
		std::cout<<" --- Collecion of gen electrons ---- "<<std::endl;
		for(size_t i = 0; i != genElectrons.size(); ++i) {
			std::cout<<" Gen electron "<< i <<" 4-momentum :("<< genElectrons.at(i).pt() <<","<<genElectrons.at(i).eta()<<","<<genElectrons.at(i).phi() <<","<< genElectrons.at(i).energy()<<")"<<std::endl;	
		}
		std::cout<<" ------------------------------------"<<std::endl;
	}

	//-------------- Sim track info -----------------------------------
	// Unfortunately no simHit attached to sim track - cannot extract info on how may pixel layers have hits :/
	// At least this is my understanding 
	// https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h
	// Using the collection of genElectrons identified previously & a DR matching the collection of electron track IDs is stored 

	edm::Handle<edm::View<TrackingParticle>> TrackingParticleH;
	iEvent.getByToken(trackingParticlesToken, TrackingParticleH);
	const edm::View<TrackingParticle>& trackingParticles = *TrackingParticleH;

	std::vector<unsigned int> eleTrackIds;
	auto sq = [](float x) { return x * x; };
	auto dr = [sq](float x1,float x2,float y1,float y2) {return std::sqrt(sq(x1 - x2) + sq(y1 - y2));};

	for (unsigned long ntrackingparticle = 0; ntrackingparticle < trackingParticles.size(); ntrackingparticle++) 
	{
		const auto& tp = trackingParticles.at(ntrackingparticle);
		bool isElectron = (std::abs(tp.pdgId()) == 11);

		if(!(tp.eventId().bunchCrossing() == 0 && tp.eventId().event() == 0)){continue;}
		if(tp.charge() == 0 || tp.status() != 1) {continue;}
		if(!isElectron || tp.g4Tracks().empty()){continue;}
		for(size_t genPartCounter = 0; genPartCounter < genElectrons.size(); ++genPartCounter)
		{
			auto genElec = genElectrons.at(genPartCounter);
			if(dr(genElec.eta(),tp.p4().eta(),genElec.phi(),tp.p4().phi()) < 0.1)
			{
				eleTrackIds.push_back(tp.g4Tracks().at(0).trackId());
				if(verbose_)
				{
					std::cout<< " Gen electron kinematics from TrackingParticle vs GenPartoicle collections "<<std::endl;
					std::cout<<"Pt : "<< genElec.pt() << "," << tp.p4().pt()<<", eta : "<< genElec.eta() << "," << tp.p4().eta() <<", phi : "<< genElec.phi() << "," << tp.p4().phi()<<std::endl;
					std::cout<<" From TrackingParticle "<<std::endl;
					std::cout << " Number of layers  "<< tp.numberOfTrackerLayers() << std::endl;
					std::cout << " Track ID & type " << tp.g4Tracks().at(0).trackId() << " " << tp.g4Tracks().at(0).type() << std::endl;
					std::cout<< " Event ID : " << tp.g4Tracks().at(0).eventId().event() << " and size " << tp.g4Tracks().size() << std::endl;
					std::cout<<" TrackingParticle with PdgId: "<<tp.pdgId()<<" and 4-momentum :("<<tp.p4().pt() <<","<<tp.p4().eta()<<","<<tp.p4().phi() <<","<< tp.p4().e()<<")"<<std::endl;
				}
			}
		}
	}	

	if(verbose_){
		std::cout<<" --- Collecion of gen electron track Ids ---- "<<std::endl;
		for (size_t i = 0; i != eleTrackIds.size(); ++i) {
			std::cout<<" Gen electron "<< i <<" track ID : "<< eleTrackIds.at(i) <<std::endl;	
		}
		std::cout<<" ---------------------------------------------"<<std::endl;
	}

	//-------------- Sim Hit info ---------------------------------------
	// https://cmssdt.cern.ch/dxr/CMSSW/source/SimDataFormats/TrackingHit/interface/PSimHit.h
	// Here we use the simHit parent sim trackID & the track Id determined previously from the TrackingParticle collection
	// to associate sim hits with the gen electrons

	edm::Handle<edm::PSimHitContainer> PixelBarrelHitsLowTof;
	iEvent.getByToken(tokPixelBarrelHitsLowTof_, PixelBarrelHitsLowTof);
	edm::Handle<edm::PSimHitContainer> PixelBarrelHitsHighTof;
	iEvent.getByToken(tokPixelBarrelHitsHighTof_, PixelBarrelHitsHighTof);
	edm::Handle<edm::PSimHitContainer> PixelEndcapHitsLowTof;
	iEvent.getByToken(tokPixelEndcapHitsLowTof_, PixelEndcapHitsLowTof);
	edm::Handle<edm::PSimHitContainer> PixelEndcapHitsHighTof;
	iEvent.getByToken(tokPixelEndcapHitsHighTof_, PixelEndcapHitsHighTof);

	for(size_t eleTrackIdsCounter = 0; eleTrackIdsCounter < eleTrackIds.size(); ++eleTrackIdsCounter)
	{

		for (unsigned long simHitCounter = 0; simHitCounter < PixelBarrelHitsLowTof->size(); ++simHitCounter) 
		{
			const PSimHit& simHit = (*PixelBarrelHitsLowTof)[simHitCounter];

			if(!(simHit.eventId().bunchCrossing() == 0 && simHit.eventId().event() == 0)){
				continue;
			}
			if(abs(simHit.particleType())!=11)
				continue;

			if(eleTrackIds.at(eleTrackIdsCounter)!=simHit.trackId())
				continue;

			if(verbose_)
			{
				std::cout<< " simHit.trackId() "<<simHit.trackId()<<std::endl;
				std::cout<< " simHit.detUnitId() "<<simHit.detUnitId()<<std::endl;
				std::cout<< " simHit.processType() "<<simHit.processType()<<std::endl;
				std::cout<< " simHit.particleType() "<<simHit.particleType()<<std::endl;
				std::cout<< " simHit.eventId() "<<simHit.eventId().event()<<std::endl;
			}

			DetId theDetUnitId = simHit.detUnitId();					
			auto theDet = tGeom->idToDet(theDetUnitId);
			if(theDetUnitId.subdetId() == PixelSubdetector::PixelBarrel){
				if(verbose_)
					std::cout << "Barrel Layer: " << tTopo->pxbLayer(theDetUnitId) << std::endl;
				//if(tTopo->pxbLayer(theDetUnitId)==1)
			}
			if(theDetUnitId.subdetId() == PixelSubdetector::PixelEndcap){
				if(verbose_)
					std::cout << "Endcap Layer: " << tTopo->pxbLayer(theDetUnitId) << std::endl;
			}
		}	
	}
	//-------------- hltGsfElectrons -----------------------------------
	// https://cmssdt.cern.ch/dxr/CMSSW/source/DataFormats/EgammaCandidates/interface/GsfElectron.h

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
			for(auto const& rhit: rhits)
			{
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
						}
						if(subdet == PixelSubdetector::PixelEndcap)
						{
							std::cout << "Endcap disk: " << tTopo->pxfDisk(det) << std::endl;
						}
					}
				}
			}
			recoEle_pt.push_back( eleItr->pt() );
			recoEle_eta.push_back( eleItr->eta() );
			recoEle_phi.push_back( eleItr->phi() );
		} 
	}
	tree->Fill();	
}

// Taken from here : https://github.com/waredjeb/ElectronPixelMatching/blob/add_gen_particles/SeedFromTrack/SeedFromTrackAnalyzer/plugins/ElectronMatchSeed.cc#L206

reco::GenParticle get_lastcopy_prefsrSN(reco::GenParticle part) {
	auto daughters = part.daughterRefVector();
	if (daughters.size() == 1 && daughters.at(0)->pdgId() == part.pdgId()) {
		return get_lastcopy_prefsrSN(*(daughters.at(0)));
	} else {
		return part;
	}
}

reco::GenParticle get_lastcopySN(reco::GenParticle part) {
	auto daughters = part.daughterRefVector();
	for (size_t p = 0; p != daughters.size(); ++p) {
		if (daughters.at(p)->pdgId() == part.pdgId()) {
			return get_lastcopySN(*(daughters.at(p)));
		}
	}
	return part;
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
