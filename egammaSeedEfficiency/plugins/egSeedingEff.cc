#ifndef EGSEEDINGEFF_H
#define EGSEEDINGEFF_H

/**\class egSeedingEff egSeedingEff.cc egammaGPUdevelopment/egammaSeedEfficiency/plugins/egSeedingEff.cc
 Description: Macro that looks into the pixel tracker hit reconstruction efficiency of electrons. Goal is to
	  		  understand why we loose efficiency when doubles are removed from the electron seeding.
 Implementation:
    - Find the TrackingParticles that are electrons
	- Find the simHits & RecHits associated with these tracking particles
	- Look into the number & pattern of reconstructed pixel hits for each TP
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
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

// CMSSW sim DataFormats
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/UniqueSimTrackId.h"
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

#include "SimTracker/TrackerHitAssociation/interface/ClusterTPAssociation.h"

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
		const edm::EDGetTokenT<TrackingParticleCollection> trackingParticlesToken;		
		const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
		const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
		std::vector<edm::EDGetTokenT<edm::PSimHitContainer>> simHit_;
		const edm::EDGetTokenT<ClusterTPAssociation>  clusterTPAssocToken_;
		const edm::EDGetTokenT<SiPixelRecHitCollection>  pixelRecHitToken_;

		TTree* tree;

		std::vector<float>  recoEle_pt;
		std::vector<float>  recoEle_eta;
		std::vector<float>  recoEle_phi;
		std::vector<float>  recoEle_E;

		std::vector<float>  simEle_pt;
		std::vector<float>  simEle_eta;
		std::vector<float>  simEle_phi;
		std::vector<float>  simEle_E;

		std::vector<float>  dr_RecoSim;
		std::vector<int>    isMatched;

		std::vector<int>   nSimHitLayersBPIX;
		std::vector<int>   nSimHitLayersFPIX;
		std::vector<int>  nSimHitsLayer1_BPIX;
		std::vector<int>  nSimHitsLayer2_BPIX;
		std::vector<int>  nSimHitsLayer3_BPIX;
		std::vector<int>  nSimHitsLayer4_BPIX;
		std::vector<int>  nSimHitsLayer1_FPIX_pos;
		std::vector<int>  nSimHitsLayer2_FPIX_pos;
		std::vector<int>  nSimHitsLayer3_FPIX_pos;
		std::vector<int>  nSimHitsLayer1_FPIX_neg;
		std::vector<int>  nSimHitsLayer2_FPIX_neg;
		std::vector<int>  nSimHitsLayer3_FPIX_neg;

		std::vector<int>  nRecoHitLayersBPIX;
		std::vector<int>  nRecoHitLayersFPIX;
		std::vector<int>  nRecoHitsLayer1_BPIX;
		std::vector<int>  nRecoHitsLayer2_BPIX;
		std::vector<int>  nRecoHitsLayer3_BPIX;
		std::vector<int>  nRecoHitsLayer4_BPIX;
		std::vector<int>  nRecoHitsLayer1_FPIX_pos;
		std::vector<int>  nRecoHitsLayer2_FPIX_pos;
		std::vector<int>  nRecoHitsLayer3_FPIX_pos;
		std::vector<int>  nRecoHitsLayer1_FPIX_neg;
		std::vector<int>  nRecoHitsLayer2_FPIX_neg;
		std::vector<int>  nRecoHitsLayer3_FPIX_neg;

		// Hits per module split into BPIX layers and FPIX layers (in -z & +z disks)
		// Vectors (differebt from all other variables)
		std::vector<std::vector<int>>  nSimHitsPerModuleLayer1_BPIX;
		std::vector<std::vector<int>>  nSimHitsPerModuleLayer2_BPIX;
		std::vector<std::vector<int>>  nSimHitsPerModuleLayer3_BPIX;
		std::vector<std::vector<int>>  nSimHitsPerModuleLayer4_BPIX;
		std::vector<std::vector<int>>  nSimHitsPerModuleLayer1_FPIX_pos;
		std::vector<std::vector<int>>  nSimHitsPerModuleLayer2_FPIX_pos;
		std::vector<std::vector<int>>  nSimHitsPerModuleLayer3_FPIX_pos;
		std::vector<std::vector<int>>  nSimHitsPerModuleLayer1_FPIX_neg;
		std::vector<std::vector<int>>  nSimHitsPerModuleLayer2_FPIX_neg;
		std::vector<std::vector<int>>  nSimHitsPerModuleLayer3_FPIX_neg;

		std::vector<std::vector<int>>  nRecoHitsPerModuleLayer1_BPIX;
		std::vector<std::vector<int>>  nRecoHitsPerModuleLayer2_BPIX;
		std::vector<std::vector<int>>  nRecoHitsPerModuleLayer3_BPIX;
		std::vector<std::vector<int>>  nRecoHitsPerModuleLayer4_BPIX;
		std::vector<std::vector<int>>  nRecoHitsPerModuleLayer1_FPIX_pos;
		std::vector<std::vector<int>>  nRecoHitsPerModuleLayer2_FPIX_pos;
		std::vector<std::vector<int>>  nRecoHitsPerModuleLayer3_FPIX_pos;
		std::vector<std::vector<int>>  nRecoHitsPerModuleLayer1_FPIX_neg;
		std::vector<std::vector<int>>  nRecoHitsPerModuleLayer2_FPIX_neg;
		std::vector<std::vector<int>>  nRecoHitsPerModuleLayer3_FPIX_neg;

		std::vector<std::vector<float>>  recHit_x;
		std::vector<std::vector<float>>  recHit_y;
		std::vector<std::vector<float>>  recHit_z;
		std::vector<std::vector<int>>    recHit_l;

		std::vector<std::vector<float>>  simHit_x;
		std::vector<std::vector<float>>  simHit_y;
		std::vector<std::vector<float>>  simHit_z;
		std::vector<std::vector<int>>    simHit_l;

		int run_, lumi_, event_;
		bool verbose_;
		double DeltaR_;
		const bool useGsfElectrons = false;
		bool isBarrel_;
};

//Constructor
egSeedingEff::egSeedingEff(const edm::ParameterSet& iConfig): 
							electronToken  (consumes<reco::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electron"))),
							genParticlesToken  (consumes<reco::GenParticleCollection> (iConfig.getParameter<edm::InputTag>("genParticles"))),
							trackingParticlesToken (consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("trackingParticles"))),
							topoToken_(esConsumes()),
							geomToken_(esConsumes()),
							simHit_(),
							clusterTPAssocToken_(consumes<ClusterTPAssociation>(iConfig.getParameter<edm::InputTag>("cluster2TPSrc"))),
							pixelRecHitToken_(consumes<SiPixelRecHitCollection>(iConfig.getParameter<edm::InputTag>("pixelRecHits"))),							
							verbose_(iConfig.getParameter<bool>("verbose")),
							DeltaR_(iConfig.getParameter<double>("deltaR")),
							isBarrel_(iConfig.getParameter<bool>("isBarrel"))
{
	initialize();
	usesResource("TFileService");	

	std::vector<edm::InputTag> tags = iConfig.getParameter<std::vector<edm::InputTag>>("simHitSrc");
	simHit_.reserve(tags.size());
	for (auto const &tag : tags) 
	{
		simHit_.emplace_back(consumes<edm::PSimHitContainer>(tag));
	}	
}

//Destructor
egSeedingEff::~egSeedingEff() {}

void egSeedingEff::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}

void egSeedingEff::initialize() {

	this -> run_ = 0;
	this -> lumi_= 0;
	this -> event_ = 0;

	this -> recoEle_pt.clear();
	this -> recoEle_eta.clear();
	this -> recoEle_phi.clear();
	this -> recoEle_E.clear();

	this -> simEle_pt.clear();
	this -> simEle_eta.clear();
	this -> simEle_phi.clear();
	this -> simEle_E.clear();

	this-> dr_RecoSim.clear();
	this-> isMatched.clear();

	this -> nSimHitLayersBPIX.clear();
	this -> nSimHitLayersFPIX.clear();
	this -> nSimHitsLayer1_BPIX.clear();
	this -> nSimHitsLayer2_BPIX.clear();
	this -> nSimHitsLayer3_BPIX.clear();
	this -> nSimHitsLayer4_BPIX.clear();
	this -> nSimHitsLayer1_FPIX_pos.clear();
	this -> nSimHitsLayer2_FPIX_pos.clear();
	this -> nSimHitsLayer3_FPIX_pos.clear();
	this -> nSimHitsLayer1_FPIX_neg.clear();
	this -> nSimHitsLayer2_FPIX_neg.clear();
	this -> nSimHitsLayer3_FPIX_neg.clear();

	this -> nSimHitsPerModuleLayer1_BPIX.clear();
	this -> nSimHitsPerModuleLayer2_BPIX.clear();
	this -> nSimHitsPerModuleLayer3_BPIX.clear();
	this -> nSimHitsPerModuleLayer4_BPIX.clear();
	this -> nSimHitsPerModuleLayer1_FPIX_pos.clear();
	this -> nSimHitsPerModuleLayer2_FPIX_pos.clear();
	this -> nSimHitsPerModuleLayer3_FPIX_pos.clear();
	this -> nSimHitsPerModuleLayer1_FPIX_neg.clear();
	this -> nSimHitsPerModuleLayer2_FPIX_neg.clear();
	this -> nSimHitsPerModuleLayer3_FPIX_neg.clear();

	this -> nRecoHitLayersBPIX.clear();
	this -> nRecoHitLayersFPIX.clear();
	this -> nRecoHitsLayer1_BPIX.clear();
	this -> nRecoHitsLayer2_BPIX.clear();
	this -> nRecoHitsLayer3_BPIX.clear();
	this -> nRecoHitsLayer4_BPIX.clear();
	this -> nRecoHitsLayer1_FPIX_pos.clear();
	this -> nRecoHitsLayer2_FPIX_pos.clear();
	this -> nRecoHitsLayer3_FPIX_pos.clear();
	this -> nRecoHitsLayer1_FPIX_neg.clear();
	this -> nRecoHitsLayer2_FPIX_neg.clear();
	this -> nRecoHitsLayer3_FPIX_neg.clear();

	this -> nRecoHitsPerModuleLayer1_BPIX.clear();
	this -> nRecoHitsPerModuleLayer2_BPIX.clear();
	this -> nRecoHitsPerModuleLayer3_BPIX.clear();
	this -> nRecoHitsPerModuleLayer4_BPIX.clear();
	this -> nRecoHitsPerModuleLayer1_FPIX_pos.clear();
	this -> nRecoHitsPerModuleLayer2_FPIX_pos.clear();
	this -> nRecoHitsPerModuleLayer3_FPIX_pos.clear();
	this -> nRecoHitsPerModuleLayer1_FPIX_neg.clear();
	this -> nRecoHitsPerModuleLayer2_FPIX_neg.clear();
	this -> nRecoHitsPerModuleLayer3_FPIX_neg.clear();

	this -> recHit_x.clear();
	this -> recHit_y.clear();
	this -> recHit_z.clear();
	this -> recHit_l.clear();

	this -> simHit_x.clear();
	this -> simHit_y.clear();
	this -> simHit_z.clear();
	this -> simHit_l.clear();
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

	tree->Branch("dr_RecoSim","std::vector<float>",&dr_RecoSim);
	tree->Branch("isMatched","std::vector<int>",&isMatched);

	tree->Branch("nSimHitLayersBPIX", "std::vector<int>", &nSimHitLayersBPIX  );
	tree->Branch("nSimHitLayersFPIX", "std::vector<int>", &nSimHitLayersFPIX  );
	tree->Branch("nSimHitsLayer1_BPIX", "std::vector<int>", &nSimHitsLayer1_BPIX );
	tree->Branch("nSimHitsLayer2_BPIX", "std::vector<int>", &nSimHitsLayer2_BPIX );
	tree->Branch("nSimHitsLayer3_BPIX", "std::vector<int>", &nSimHitsLayer3_BPIX );
	tree->Branch("nSimHitsLayer4_BPIX", "std::vector<int>", &nSimHitsLayer4_BPIX );
	tree->Branch("nSimHitsLayer1_FPIX_pos", "std::vector<int>", &nSimHitsLayer1_FPIX_pos );
	tree->Branch("nSimHitsLayer2_FPIX_pos", "std::vector<int>", &nSimHitsLayer2_FPIX_pos );
	tree->Branch("nSimHitsLayer3_FPIX_pos", "std::vector<int>", &nSimHitsLayer3_FPIX_pos );
	tree->Branch("nSimHitsLayer1_FPIX_neg", "std::vector<int>", &nSimHitsLayer1_FPIX_neg );
	tree->Branch("nSimHitsLayer2_FPIX_neg", "std::vector<int>", &nSimHitsLayer2_FPIX_neg );
	tree->Branch("nSimHitsLayer3_FPIX_neg", "std::vector<int>", &nSimHitsLayer3_FPIX_neg );

	tree->Branch("nSimHitsPerModuleLayer1_BPIX", "std::vector<std::vector<int>>", &nSimHitsPerModuleLayer1_BPIX );
	tree->Branch("nSimHitsPerModuleLayer2_BPIX", "std::vector<std::vector<int>>", &nSimHitsPerModuleLayer2_BPIX );
	tree->Branch("nSimHitsPerModuleLayer3_BPIX", "std::vector<std::vector<int>>", &nSimHitsPerModuleLayer3_BPIX );
	tree->Branch("nSimHitsPerModuleLayer4_BPIX", "std::vector<std::vector<int>>", &nSimHitsPerModuleLayer4_BPIX );
	tree->Branch("nSimHitsPerModuleLayer1_FPIX_pos", "std::vector<std::vector<int>>", &nSimHitsPerModuleLayer1_FPIX_pos );
	tree->Branch("nSimHitsPerModuleLayer2_FPIX_pos", "std::vector<std::vector<int>>", &nSimHitsPerModuleLayer2_FPIX_pos );
	tree->Branch("nSimHitsPerModuleLayer3_FPIX_pos", "std::vector<std::vector<int>>", &nSimHitsPerModuleLayer3_FPIX_pos );
	tree->Branch("nSimHitsPerModuleLayer1_FPIX_neg", "std::vector<std::vector<int>>", &nSimHitsPerModuleLayer1_FPIX_neg );
	tree->Branch("nSimHitsPerModuleLayer2_FPIX_neg", "std::vector<std::vector<int>>", &nSimHitsPerModuleLayer2_FPIX_neg );
	tree->Branch("nSimHitsPerModuleLayer3_FPIX_neg", "std::vector<std::vector<int>>", &nSimHitsPerModuleLayer3_FPIX_neg );

	tree->Branch("nRecoHitLayersBPIX", "std::vector<int>", &nRecoHitLayersBPIX );
	tree->Branch("nRecoHitLayersFPIX", "std::vector<int>", &nRecoHitLayersFPIX );
	tree->Branch("nRecoHitsLayer1_BPIX", "std::vector<int>", &nRecoHitsLayer1_BPIX );
	tree->Branch("nRecoHitsLayer2_BPIX", "std::vector<int>", &nRecoHitsLayer2_BPIX );
	tree->Branch("nRecoHitsLayer3_BPIX", "std::vector<int>", &nRecoHitsLayer3_BPIX );
	tree->Branch("nRecoHitsLayer4_BPIX", "std::vector<int>", &nRecoHitsLayer4_BPIX );
	tree->Branch("nRecoHitsLayer1_FPIX_pos", "std::vector<int>", &nRecoHitsLayer1_FPIX_pos );
	tree->Branch("nRecoHitsLayer2_FPIX_pos", "std::vector<int>", &nRecoHitsLayer2_FPIX_pos );
	tree->Branch("nRecoHitsLayer3_FPIX_pos", "std::vector<int>", &nRecoHitsLayer3_FPIX_pos );
	tree->Branch("nRecoHitsLayer1_FPIX_neg", "std::vector<int>", &nRecoHitsLayer1_FPIX_neg );
	tree->Branch("nRecoHitsLayer2_FPIX_neg", "std::vector<int>", &nRecoHitsLayer2_FPIX_neg );
	tree->Branch("nRecoHitsLayer3_FPIX_neg", "std::vector<int>", &nRecoHitsLayer3_FPIX_neg );

	tree->Branch("nRecoHitsPerModuleLayer1_BPIX", "std::vector<std::vector<int>>", &nRecoHitsPerModuleLayer1_BPIX );
	tree->Branch("nRecoHitsPerModuleLayer2_BPIX", "std::vector<std::vector<int>>", &nRecoHitsPerModuleLayer2_BPIX );
	tree->Branch("nRecoHitsPerModuleLayer3_BPIX", "std::vector<std::vector<int>>", &nRecoHitsPerModuleLayer3_BPIX );
	tree->Branch("nRecoHitsPerModuleLayer4_BPIX", "std::vector<std::vector<int>>", &nRecoHitsPerModuleLayer4_BPIX );
	tree->Branch("nRecoHitsPerModuleLayer1_FPIX_pos", "std::vector<std::vector<int>>", &nRecoHitsPerModuleLayer1_FPIX_pos );
	tree->Branch("nRecoHitsPerModuleLayer2_FPIX_pos", "std::vector<std::vector<int>>", &nRecoHitsPerModuleLayer2_FPIX_pos );
	tree->Branch("nRecoHitsPerModuleLayer3_FPIX_pos", "std::vector<std::vector<int>>", &nRecoHitsPerModuleLayer3_FPIX_pos );
	tree->Branch("nRecoHitsPerModuleLayer1_FPIX_neg", "std::vector<std::vector<int>>", &nRecoHitsPerModuleLayer1_FPIX_neg );
	tree->Branch("nRecoHitsPerModuleLayer2_FPIX_neg", "std::vector<std::vector<int>>", &nRecoHitsPerModuleLayer2_FPIX_neg );
	tree->Branch("nRecoHitsPerModuleLayer3_FPIX_neg", "std::vector<std::vector<int>>", &nRecoHitsPerModuleLayer3_FPIX_neg );

	tree->Branch("recHit_x" , "std::vector<std::vector<float>>", &recHit_x );
	tree->Branch("recHit_y", "std::vector<std::vector<float>>", &recHit_y );
	tree->Branch("recHit_z", "std::vector<std::vector<float>>", &recHit_z );
	tree->Branch("recHit_l", "std::vector<std::vector<int>>", &recHit_l );

	tree->Branch("simHit_x" , "std::vector<std::vector<float>>", &simHit_x );
	tree->Branch("simHit_y", "std::vector<std::vector<float>>", &simHit_y );
	tree->Branch("simHit_z", "std::vector<std::vector<float>>", &simHit_z );
	tree->Branch("simHit_l", "std::vector<std::vector<int>>", &simHit_l );

	tree->Branch("recoEle_pt" , "std::vector<float>", &recoEle_pt );
	tree->Branch("recoEle_eta", "std::vector<float>", &recoEle_eta );
	tree->Branch("recoEle_phi", "std::vector<float>", &recoEle_phi );
	tree->Branch("recoEle_E", "std::vector<float>", &recoEle_E );

	return ;
}

void egSeedingEff::endJob() {}
void egSeedingEff::endRun(edm::Run const&, edm::EventSetup const&) {}

reco::GenParticle get_lastcopy_prefsrSN(reco::GenParticle part);
reco::GenParticle get_lastcopySN(reco::GenParticle part);
using P = std::pair<OmniClusterRef, TrackingParticleRef>;
bool compare(const P &i, const P &j) { return i.second.index() > j.second.index(); }

void egSeedingEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

	using namespace edm;
	using namespace std;

	initialize();

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

	edm::Handle<TrackingParticleCollection> TrackingParticleH;
	iEvent.getByToken(trackingParticlesToken, TrackingParticleH);
	const TrackingParticleCollection& trackingParticles = *TrackingParticleH;
	TrackingParticleRefVector trackingPartRefs;
	std::vector<unsigned int> eleTrackIds;
	std::vector<std::pair<uint32_t, EncodedEventId>> TrackIds;
	
	auto sq = [](float x) { return x * x; };
	auto dr = [sq](float x1,float x2,float y1,float y2) {return std::sqrt(sq(x1 - x2) + sq(y1 - y2));};

	for (unsigned long ntrackingparticle = 0; ntrackingparticle < trackingParticles.size(); ntrackingparticle++) 
	{
		const auto& tp = trackingParticles.at(ntrackingparticle);
		bool isElectron = (std::abs(tp.pdgId()) == 11);

		TrackingParticleRef trackingParticleRef(TrackingParticleH, ntrackingparticle);

		if(!(tp.eventId().bunchCrossing() == 0 && tp.eventId().event() == 0)){continue;}
		if(tp.charge() == 0 || tp.status() != 1) {continue;}
		if(!isElectron || tp.g4Tracks().empty()){continue;}

		for(size_t genPartCounter = 0; genPartCounter < genElectrons.size(); ++genPartCounter)
		{
			auto genElec = genElectrons.at(genPartCounter);

			bool pass = false;
			if(isBarrel_ && (abs(tp.p4().eta()) < 1.45))
				pass = true;
			if(!isBarrel_ && (abs(tp.p4().eta()) > 1.55))
				pass = true;
			if(!pass)
				continue;
				
			if(dr(genElec.eta(),tp.p4().eta(),genElec.phi(),tp.p4().phi()) < 0.01)
			{
				eleTrackIds.push_back(tp.g4Tracks().at(0).trackId());
				simEle_pt.push_back(tp.p4().pt());
				simEle_eta.push_back(tp.p4().eta());
				simEle_phi.push_back(tp.p4().phi());
				simEle_E.push_back(tp.p4().energy());
				isMatched.push_back(0);
				trackingPartRefs.push_back(trackingParticleRef);

				for (auto const &trk : tp.g4Tracks()) 
				{
					EncodedEventId eid(tp.eventId());
					UniqueSimTrackId trkid(trk.trackId(), eid);
					TrackIds.push_back(trkid);
				}

				if(verbose_)
				{
					std::cout<< " Gen electron kinematics from TrackingParticle vs GenParticle collections : "<<std::endl;
					std::cout<<"Pt : "<< genElec.pt() << "," << tp.p4().pt()<<", eta : "<< genElec.eta() << "," << tp.p4().eta() <<", phi : "<< genElec.phi() << "," << tp.p4().phi()<<std::endl;
					std::cout<<" From TrackingParticle : "<<std::endl;
					std::cout << " Number of layers  "<< tp.numberOfTrackerLayers() << std::endl;
					std::cout << " Track ID & type " << tp.g4Tracks().at(0).trackId() << " " << tp.g4Tracks().at(0).type() << std::endl;
					std::cout<< " Event ID : " << tp.g4Tracks().at(0).eventId().event() << " and size : " << tp.g4Tracks().size() << "and raw ID " << tp.g4Tracks().at(0).eventId().rawId() << std::endl;
					std::cout<<" TrackingParticle with PdgId: "<<tp.pdgId()<<" and 4-momentum :("<<tp.p4().pt() <<","<<tp.p4().eta()<<","<<tp.p4().phi() <<","<< tp.p4().e()<<")"<<std::endl;
				}
			}
		}
	}	

	if(verbose_){
		std::cout<<" --- Collecion of track IDs from eleTrackIds ---- "<<std::endl;
		for (size_t i = 0; i != eleTrackIds.size(); ++i) {
			std::cout<<" Gen electron "<< i <<" track ID : "<< eleTrackIds.at(i) <<std::endl;	
		}
		std::cout<<" --- Collecion of all track ids from TrackIds  ---- "<<std::endl;
		for (size_t i = 0; i != TrackIds.size(); ++i) {
			std::cout<<" Gen electron "<< i <<" track ID : "<< TrackIds.at(i).first << std::endl;	
		}
		std::cout<<" ---------------------------------------------"<<std::endl;
	}

	//-------------- Sim Hit info ---------------------------------------
	// https://cmssdt.cern.ch/dxr/CMSSW/source/SimDataFormats/TrackingHit/interface/PSimHit.h
	// Here we use the simHit parent sim trackID & the track Id determined previously from the TrackingParticle collection
	// to associate sim hits with the gen electrons
	// Inspired from here : https://cmssdt.cern.ch/dxr/CMSSW/source/SimGeneral/TrackingAnalysis/plugins/SimHitTPAssociationProducer.cc

	for(size_t eleTrackIdsCounter = 0; eleTrackIdsCounter < eleTrackIds.size(); ++eleTrackIdsCounter)
	{
		unsigned int nLayersFPIX = 0;
		unsigned int nLayersBPIX = 0;	
		unsigned int nHitsPerLayer_BPIX[4] = {0};				
		unsigned int nHitsPerLayer_FPIX[6] = {0};
		std::vector<int> pxbL1(0),pxbL2(0),pxbL3(0),pxbL4(0),pxfL1p(0),pxfL2p(0),pxfL3p(0),pxfL1n(0),pxfL2n(0),pxfL3n(0);
		std::vector<float> rx(0),ry(0),rz(0);
		std::vector<int> rl(0);

		for(size_t TrackIdsCounter = 0; TrackIdsCounter < TrackIds.size(); ++TrackIdsCounter)
		{
			for (auto const &psit : simHit_) 
			{
				edm::Handle<edm::PSimHitContainer> PSimHitCollectionH;
				iEvent.getByToken(psit, PSimHitCollectionH);
				auto const &pSimHitCollection = *PSimHitCollectionH;			

				for (unsigned int simHitCounter = 0, size = pSimHitCollection.size(); simHitCounter < size; ++simHitCounter) 
				{
					auto const &simHit = pSimHitCollection[simHitCounter];				

					// For associating with TrackingParticle collection
					UniqueSimTrackId simTkIds(simHit.trackId(), simHit.eventId());

					if(!(TrackIds.at(TrackIdsCounter)==simTkIds))
						continue;

					if(!(simHit.eventId().bunchCrossing() == 0 && simHit.eventId().event() == 0))
						continue;

					if(abs(simHit.particleType())!=11)
						continue;

					if(verbose_)
					{
						std::cout<< " simHit.trackId() "<<simHit.trackId()<<std::endl;
						std::cout<< " simHit.detUnitId() "<<simHit.detUnitId()<<std::endl;
						std::cout<< " simHit.processType() "<<simHit.processType()<<std::endl;
						std::cout<< " simHit.particleType() "<<simHit.particleType()<<std::endl;
						std::cout<< " simHit.eventId() "<<simHit.eventId().event()<<std::endl;
						std::cout<< " simHit.eventId().rawId() "<< simHit.eventId().rawId()<<std::endl;
					}

					DetId theDetUnitId = simHit.detUnitId();	
					GlobalPoint gp = tGeom->idToDet(theDetUnitId)->surface().toGlobal(simHit.localPosition());
					rx.push_back(gp.x());
					ry.push_back(gp.y());
					rz.push_back(gp.z());
					
					if(theDetUnitId.subdetId() == PixelSubdetector::PixelBarrel)
					{						
						if(verbose_)
							std::cout<<" Barrel Layer: "<< tTopo->pxbLayer(theDetUnitId) << "  and module " << tTopo->pxbModule(theDetUnitId) <<std::endl;

						if(tTopo->pxbLayer(theDetUnitId)==1){
							++nHitsPerLayer_BPIX[0];
							pxbL1.push_back(tTopo->pxbModule(theDetUnitId));
							rl.push_back(1);
						}							
						if(tTopo->pxbLayer(theDetUnitId)==2){
							++nHitsPerLayer_BPIX[1];
							pxbL2.push_back(tTopo->pxbModule(theDetUnitId));
							rl.push_back(2);
						}	
						if(tTopo->pxbLayer(theDetUnitId)==3){
							++nHitsPerLayer_BPIX[2];
							pxbL3.push_back(tTopo->pxbModule(theDetUnitId));
							rl.push_back(3);
						}	
						if(tTopo->pxbLayer(theDetUnitId)==4){
							++nHitsPerLayer_BPIX[3];
							pxbL4.push_back(tTopo->pxbModule(theDetUnitId));
							rl.push_back(4);
						}
					}

					if(theDetUnitId.subdetId() == PixelSubdetector::PixelEndcap)
					{
						if(verbose_)
							std::cout<< "Endcap Layer : "<< tTopo->pxfDisk(theDetUnitId) << " module : " << tTopo->pxfModule(theDetUnitId) << " and side : " << tTopo->pxfSide(theDetUnitId) << std::endl; //side=1 for -z, 2 for +z

						if(tTopo->pxfDisk(theDetUnitId)==1 && tTopo->pxfSide(theDetUnitId)==2){
							++nHitsPerLayer_FPIX[0];
							pxfL1p.push_back(tTopo->pxfModule(theDetUnitId));
							rl.push_back(5);
						}
						if(tTopo->pxfDisk(theDetUnitId)==1 && tTopo->pxfSide(theDetUnitId)==1){
							++nHitsPerLayer_FPIX[3];							
							pxfL1n.push_back(tTopo->pxfModule(theDetUnitId));
							rl.push_back(6);
						}
						if(tTopo->pxfDisk(theDetUnitId)==2 && tTopo->pxfSide(theDetUnitId)==2){
							++nHitsPerLayer_FPIX[1];
							pxfL2p.push_back(tTopo->pxfModule(theDetUnitId));
							rl.push_back(7);
						}
						if(tTopo->pxfDisk(theDetUnitId)==2 && tTopo->pxfSide(theDetUnitId)==1){
							++nHitsPerLayer_FPIX[4];	
							pxfL2n.push_back(tTopo->pxfModule(theDetUnitId));
							rl.push_back(8);
						}
						if(tTopo->pxfDisk(theDetUnitId)==3 && tTopo->pxfSide(theDetUnitId)==2){
							++nHitsPerLayer_FPIX[2];
							pxfL3p.push_back(tTopo->pxfModule(theDetUnitId));
							rl.push_back(9);
						}
						if(tTopo->pxfDisk(theDetUnitId)==3 && tTopo->pxfSide(theDetUnitId)==1){
							++nHitsPerLayer_FPIX[5];
							pxfL3n.push_back(tTopo->pxfModule(theDetUnitId));
							rl.push_back(10);
						}
					}
				}
			}
		}

		for(unsigned int i=0;i<4;++i)
		{
			if(nHitsPerLayer_BPIX[i]>0)
				++nLayersBPIX;
		}		
		for(unsigned int i=0;i<3;++i)
		{
			if(nHitsPerLayer_FPIX[i]>0 || nHitsPerLayer_FPIX[i+3]>0)
				++nLayersFPIX;
		}

		if(verbose_)
		{
			std::cout<< " Sim electron "<<std::endl;
			std::cout<< " nLayers with hits BPIX : "<<nLayersBPIX << " and FPIX " << nLayersFPIX << std::endl;
			std::cout<< " Layer 1 hits BPIX : "<<nHitsPerLayer_BPIX[0] << " and FPIX pos " << nHitsPerLayer_FPIX[0] << " and FPIX neg "<< nHitsPerLayer_FPIX[3]<< std::endl;
			std::cout<< " Layer 2 hits BPIX : "<<nHitsPerLayer_BPIX[1] << " and FPIX pos " << nHitsPerLayer_FPIX[1] << " and FPIX neg "<< nHitsPerLayer_FPIX[4]<< std::endl;
			std::cout<< " Layer 3 hits BPIX : "<<nHitsPerLayer_BPIX[2] << " and FPIX pos " << nHitsPerLayer_FPIX[2] << " and FPIX neg "<< nHitsPerLayer_FPIX[5]<< std::endl;
			std::cout<< " Layer 4 hits BPIX : "<<nHitsPerLayer_BPIX[3] << std::endl;
		}

		nSimHitLayersBPIX.push_back(nLayersBPIX);
		nSimHitLayersFPIX.push_back(nLayersFPIX);
		nSimHitsLayer1_BPIX.push_back(nHitsPerLayer_BPIX[0]);
		nSimHitsLayer2_BPIX.push_back(nHitsPerLayer_BPIX[1]);
		nSimHitsLayer3_BPIX.push_back(nHitsPerLayer_BPIX[2]);
		nSimHitsLayer4_BPIX.push_back(nHitsPerLayer_BPIX[3]);
		nSimHitsLayer1_FPIX_pos.push_back(nHitsPerLayer_FPIX[0]);
		nSimHitsLayer2_FPIX_pos.push_back(nHitsPerLayer_FPIX[1]);
		nSimHitsLayer3_FPIX_pos.push_back(nHitsPerLayer_FPIX[2]);
		nSimHitsLayer1_FPIX_neg.push_back(nHitsPerLayer_FPIX[3]);
		nSimHitsLayer2_FPIX_neg.push_back(nHitsPerLayer_FPIX[4]);
		nSimHitsLayer3_FPIX_neg.push_back(nHitsPerLayer_FPIX[5]);
		nSimHitsPerModuleLayer1_BPIX.push_back(pxbL1);
		nSimHitsPerModuleLayer2_BPIX.push_back(pxbL2);
		nSimHitsPerModuleLayer3_BPIX.push_back(pxbL3);
		nSimHitsPerModuleLayer4_BPIX.push_back(pxbL4);
		nSimHitsPerModuleLayer1_FPIX_pos.push_back(pxfL1p);
		nSimHitsPerModuleLayer2_FPIX_pos.push_back(pxfL2p);
		nSimHitsPerModuleLayer3_FPIX_pos.push_back(pxfL3p);
		nSimHitsPerModuleLayer1_FPIX_neg.push_back(pxfL1n);
		nSimHitsPerModuleLayer2_FPIX_neg.push_back(pxfL2n);
		nSimHitsPerModuleLayer3_FPIX_neg.push_back(pxfL3n);
		simHit_x.push_back(rx);
		simHit_y.push_back(ry);
		simHit_z.push_back(rz);
		simHit_l.push_back(rl);
	}

	// ---------------- TP to Cluster Association ------------------------------------------------------
	// Access the Tracker RecHits & associate them with the corresponding TrackingParticle by using the ClusterTPAssociationProducer
	// https://cmssdt.cern.ch/dxr/CMSSW/source/SimTracker/TrackerHitAssociation/plugins/ClusterTPAssociationProducer.cc

	edm::Handle<ClusterTPAssociation> clusterTPAssocH;
	iEvent.getByToken(clusterTPAssocToken_,clusterTPAssocH);
	const ClusterTPAssociation& clusterToTPMap = *clusterTPAssocH;
	edm::Handle<SiPixelRecHitCollection> recHitCollH;
	iEvent.getByToken(pixelRecHitToken_, recHitCollH);
	// ---------------- RecHits -----------------------------------------------------------------------
	const SiPixelRecHitCollection* rechits;
	rechits = recHitCollH.product();
	// --------- Look for equal_range of a given TrackingParticle --------------------------------------
	// Use the collection of TrackingParticleRefs stored in the previous step after being matched with a 
	// gen electron coming from a Z boson decay

	auto clusterTPmap = clusterToTPMap.map();
	std::sort(clusterTPmap.begin(), clusterTPmap.end(), compare);

	if(!useGsfElectrons)
	{
		for (RefVector<TrackingParticleCollection>::const_iterator it = trackingPartRefs.begin(); it != trackingPartRefs.end();it++) 
		{
			if(verbose_)
				std::cout<<"Tracking Particle pT : "<< (*it)->pt() << std::endl;

			unsigned int nLayersFPIX = 0;
			unsigned int nLayersBPIX = 0;
			unsigned int nHitsPerLayer_BPIX[4] = {0};				
			unsigned int nHitsPerLayer_FPIX[6] = {0};	
			std::vector<int> pxbL1(0),pxbL2(0),pxbL3(0),pxbL4(0),pxfL1p(0),pxfL2p(0),pxfL3p(0),pxfL1n(0),pxfL2n(0),pxfL3n(0);
			std::vector<float> rx(0),ry(0),rz(0);
			std::vector<int> rl(0);

			auto clusterRange = std::equal_range(clusterTPmap.begin(), clusterTPmap.end(), std::make_pair(OmniClusterRef(), (*it)), compare);
			if (clusterRange.first != clusterRange.second) 
			{
				for (auto ip = clusterRange.first; ip != clusterRange.second; ++ip) 
				{
					const OmniClusterRef &cluster = ip->first;
					const TrackingParticleRef  &trackingPart = ip->second;

					if (!(cluster.isPixel() && cluster.isValid())) 
						continue;

					if(verbose_)	
						std::cout<< " TrackingPart pT: "<< trackingPart->p4() <<std::endl;

					for (SiPixelRecHitCollection::const_iterator itr = rechits->begin(); itr != rechits->end(); itr++) 
					{
						SiPixelRecHitCollection::DetSet hits = *itr;
						DetId detId = DetId(hits.detId());
						SiPixelRecHitCollection::const_iterator recHitMatch = rechits->find(detId);
						const SiPixelRecHitCollection::DetSet recHitRange = *recHitMatch;
						for (SiPixelRecHitCollection::DetSet::const_iterator recHitIterator = recHitRange.begin();recHitIterator != recHitRange.end();++recHitIterator) 
						{
							const SiPixelRecHit* recHit = &(*recHitIterator);
							int clusterSize = recHit->cluster()->size();

							// rechit position in global coordinates
							GlobalPoint gp = tGeom->idToDet(detId)->surface().toGlobal(recHit->localPosition());
					
							if(OmniClusterRef(recHit->cluster()) == cluster )
							{
								rx.push_back(gp.x());
								ry.push_back(gp.y());
								rz.push_back(gp.z());

								if (detId.det() == DetId::Tracker && detId.subdetId() == PixelSubdetector::PixelBarrel)
								{
									if(verbose_)
										std::cout<< "ClusterSize : "<< clusterSize << " & detector layer : "<< tTopo->layer(detId) << std::endl;
									if(tTopo->pxbLayer(detId)==1){
										++nHitsPerLayer_BPIX[0];
										pxbL1.push_back(tTopo->pxbModule(detId));
										rl.push_back(1);
									}							
									if(tTopo->pxbLayer(detId)==2){
										++nHitsPerLayer_BPIX[1];
										pxbL2.push_back(tTopo->pxbModule(detId));
										rl.push_back(2);
									}	
									if(tTopo->pxbLayer(detId)==3){
										++nHitsPerLayer_BPIX[2];
										pxbL3.push_back(tTopo->pxbModule(detId));
										rl.push_back(3);
									}	
									if(tTopo->pxbLayer(detId)==4){
										++nHitsPerLayer_BPIX[3];
										pxbL4.push_back(tTopo->pxbModule(detId));
										rl.push_back(4);
									}
								}
								if (detId.det() == DetId::Tracker && detId.subdetId() == PixelSubdetector::PixelEndcap)
								{
									if(verbose_)
										std::cout<< "ClusterSize : "<< clusterSize << " and Endcap Layer : "<< tTopo->pxfDisk(detId) << " module : " << tTopo->pxfModule(detId) << " and side : " << tTopo->pxfSide(detId) << std::endl; //side=1 for -z, 2 for +z

									if(tTopo->pxfDisk(detId)==1 && tTopo->pxfSide(detId)==2){
										++nHitsPerLayer_FPIX[0];
										pxfL1p.push_back(tTopo->pxfModule(detId));
										rl.push_back(5);
									}
									if(tTopo->pxfDisk(detId)==1 && tTopo->pxfSide(detId)==1){
										++nHitsPerLayer_FPIX[3];							
										pxfL1n.push_back(tTopo->pxfModule(detId));
										rl.push_back(6);
									}
									if(tTopo->pxfDisk(detId)==2 && tTopo->pxfSide(detId)==2){
										++nHitsPerLayer_FPIX[1];
										pxfL2p.push_back(tTopo->pxfModule(detId));
										rl.push_back(7);
									}
									if(tTopo->pxfDisk(detId)==2 && tTopo->pxfSide(detId)==1){
										++nHitsPerLayer_FPIX[4];	
										pxfL2n.push_back(tTopo->pxfModule(detId));
										rl.push_back(8);
									}
									if(tTopo->pxfDisk(detId)==3 && tTopo->pxfSide(detId)==2){
										++nHitsPerLayer_FPIX[2];
										pxfL3p.push_back(tTopo->pxfModule(detId));
										rl.push_back(9);
									}
									if(tTopo->pxfDisk(detId)==3 && tTopo->pxfSide(detId)==1){
										++nHitsPerLayer_FPIX[5];
										pxfL3n.push_back(tTopo->pxfModule(detId));
										rl.push_back(10);
									}
								}
							}
						}
					}
				}
			}	
	
			for(unsigned int i=0;i<4;++i)
			{
				if(nHitsPerLayer_BPIX[i]>0)
					++nLayersBPIX;
			}		
			for(unsigned int i=0;i<3;++i)
			{
				if(nHitsPerLayer_FPIX[i]>0 || nHitsPerLayer_FPIX[i+3]>0)
					++nLayersFPIX;
			}

			if(verbose_)
			{
				std::cout<< " Matched RecHits Per TrackingParticle "<<std::endl;
				std::cout<< " nLayers with hits BPIX : "<<nLayersBPIX << " and FPIX " << nLayersFPIX << std::endl;
				std::cout<< " Layer 1 hits BPIX : "<<nHitsPerLayer_BPIX[0] << " and FPIX pos " << nHitsPerLayer_FPIX[0] << " and FPIX neg "<< nHitsPerLayer_FPIX[3]<< std::endl;
				std::cout<< " Layer 2 hits BPIX : "<<nHitsPerLayer_BPIX[1] << " and FPIX pos " << nHitsPerLayer_FPIX[1] << " and FPIX neg "<< nHitsPerLayer_FPIX[4]<< std::endl;
				std::cout<< " Layer 3 hits BPIX : "<<nHitsPerLayer_BPIX[2] << " and FPIX pos " << nHitsPerLayer_FPIX[2] << " and FPIX neg "<< nHitsPerLayer_FPIX[5]<< std::endl;
				std::cout<< " Layer 4 hits BPIX : "<<nHitsPerLayer_BPIX[3] << std::endl;
			}

			nRecoHitsLayer1_BPIX.push_back(nHitsPerLayer_BPIX[0]);
			nRecoHitsLayer2_BPIX.push_back(nHitsPerLayer_BPIX[1]);
			nRecoHitsLayer3_BPIX.push_back(nHitsPerLayer_BPIX[2]);
			nRecoHitsLayer4_BPIX.push_back(nHitsPerLayer_BPIX[3]);
			nRecoHitsLayer1_FPIX_pos.push_back(nHitsPerLayer_FPIX[0]);
			nRecoHitsLayer2_FPIX_pos.push_back(nHitsPerLayer_FPIX[1]);
			nRecoHitsLayer3_FPIX_pos.push_back(nHitsPerLayer_FPIX[2]);
			nRecoHitsLayer1_FPIX_neg.push_back(nHitsPerLayer_FPIX[3]);
			nRecoHitsLayer2_FPIX_neg.push_back(nHitsPerLayer_FPIX[4]);
			nRecoHitsLayer3_FPIX_neg.push_back(nHitsPerLayer_FPIX[5]);

			nRecoHitLayersBPIX.push_back(nLayersBPIX);
			nRecoHitLayersFPIX.push_back(nLayersFPIX);

			nRecoHitsPerModuleLayer1_BPIX.push_back(pxbL1);
			nRecoHitsPerModuleLayer2_BPIX.push_back(pxbL2);
			nRecoHitsPerModuleLayer3_BPIX.push_back(pxbL3);
			nRecoHitsPerModuleLayer4_BPIX.push_back(pxbL4);
			nRecoHitsPerModuleLayer1_FPIX_pos.push_back(pxfL1p);
			nRecoHitsPerModuleLayer2_FPIX_pos.push_back(pxfL2p);
			nRecoHitsPerModuleLayer3_FPIX_pos.push_back(pxfL3p);
			nRecoHitsPerModuleLayer1_FPIX_neg.push_back(pxfL1n);
			nRecoHitsPerModuleLayer2_FPIX_neg.push_back(pxfL2n);
			nRecoHitsPerModuleLayer3_FPIX_neg.push_back(pxfL3n);
			recHit_x.push_back(rx);
			recHit_y.push_back(ry);
			recHit_z.push_back(rz);
			recHit_l.push_back(rl);
		}
	}
	
	//-------------- hltGsfElectrons -----------------------------------~
	// DR matching with gen electron identified previously
	// https://cmssdt.cern.ch/dxr/CMSSW/source/DataFormats/EgammaCandidates/interface/GsfElectron.h

	if(electronH.isValid()) 
	{
		for(size_t genElectronCounter = 0; genElectronCounter < simEle_pt.size(); ++genElectronCounter)
		{
			for (auto eleItr = electronH->begin(); eleItr != electronH->end(); ++eleItr) 
			{	
				bool matched = false;

				float deltaR = dr(simEle_eta.at(genElectronCounter),eleItr->eta(),simEle_phi.at(genElectronCounter),eleItr->phi());
			
				if( (deltaR < DeltaR_) && !matched)
				{
					matched = true;
					isMatched.at(genElectronCounter) = 1;
					dr_RecoSim.push_back(deltaR);

					unsigned int nLayersFPIX = 0;
					unsigned int nLayersBPIX = 0;
					unsigned int nHitsPerLayer_BPIX[4] = {0};				
					unsigned int nHitsPerLayer_FPIX[6] = {0};		

					auto gsfTrack = eleItr->gsfTrack();
					auto seed = gsfTrack->seedRef();
					auto rhits = seed->recHits();

					// Check number of hits 
					if(verbose_)
						std::cout<<" Number of hits per seed : "<< seed->nHits() <<std::endl;	

					// Check hits per Pixel layer
					for(auto const& rhit: rhits)
					{
						if(!useGsfElectrons)
							continue;

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
									if(verbose_)	
										std::cout << "Barrel Layer: " << tTopo->pxbLayer(det) << std::endl;
									if(tTopo->pxbLayer(det)==1)
										++nHitsPerLayer_BPIX[0];
									if(tTopo->pxbLayer(det)==2)
										++nHitsPerLayer_BPIX[1];
									if(tTopo->pxbLayer(det)==3)
										++nHitsPerLayer_BPIX[2];
									if(tTopo->pxbLayer(det)==4)
										++nHitsPerLayer_BPIX[3];
								}
								if(subdet == PixelSubdetector::PixelEndcap)
								{
									if(verbose_)	
										std::cout << "Endcap disk: " << tTopo->pxfDisk(det) << std::endl;
									if(tTopo->pxfDisk(det)==1 && tTopo->pxfSide(det)==2)
										++nHitsPerLayer_FPIX[0];
									if(tTopo->pxfDisk(det)==1 && tTopo->pxfSide(det)==1)
										++nHitsPerLayer_FPIX[3];							
									if(tTopo->pxfDisk(det)==2 && tTopo->pxfSide(det)==2)
										++nHitsPerLayer_FPIX[1];
									if(tTopo->pxfDisk(det)==2 && tTopo->pxfSide(det)==1)
										++nHitsPerLayer_FPIX[4];							
									if(tTopo->pxfDisk(det)==3 && tTopo->pxfSide(det)==2)
										++nHitsPerLayer_FPIX[2];
									if(tTopo->pxfDisk(det)==3 && tTopo->pxfSide(det)==1)
										++nHitsPerLayer_FPIX[5];										
								}
							}
						}
					}

					for(unsigned int i=0;i<4;++i)
					{
						if(nHitsPerLayer_BPIX[i]>0)
							++nLayersBPIX;
					}		
					for(unsigned int i=0;i<3;++i)
					{
						if(nHitsPerLayer_FPIX[i]>0 || nHitsPerLayer_FPIX[i+3]>0)
							++nLayersFPIX;
					}	
					
					if(verbose_)
					{
						std::cout<<" Matched reco electron "<<std::endl;
						std::cout<< " nLayers with hits BPIX : "<<nLayersBPIX << " and FPIX " << nLayersFPIX << std::endl;
						std::cout<< " Layer 1 hits BPIX : "<<nHitsPerLayer_BPIX[0] << " and FPIX pos " << nHitsPerLayer_FPIX[0] << " and FPIX neg "<< nHitsPerLayer_FPIX[3]<< std::endl;
						std::cout<< " Layer 2 hits BPIX : "<<nHitsPerLayer_BPIX[1] << " and FPIX pos " << nHitsPerLayer_FPIX[1] << " and FPIX neg "<< nHitsPerLayer_FPIX[4]<< std::endl;
						std::cout<< " Layer 3 hits BPIX : "<<nHitsPerLayer_BPIX[2] << " and FPIX pos " << nHitsPerLayer_FPIX[2] << " and FPIX neg "<< nHitsPerLayer_FPIX[5]<< std::endl;
						std::cout<< " Layer 4 hits BPIX : "<<nHitsPerLayer_BPIX[3] << std::endl;
					}

					if(useGsfElectrons)
					{
						nRecoHitsLayer1_BPIX.push_back(nHitsPerLayer_BPIX[0]);
						nRecoHitsLayer2_BPIX.push_back(nHitsPerLayer_BPIX[1]);
						nRecoHitsLayer3_BPIX.push_back(nHitsPerLayer_BPIX[2]);
						nRecoHitsLayer4_BPIX.push_back(nHitsPerLayer_BPIX[3]);
						nRecoHitsLayer1_FPIX_pos.push_back(nHitsPerLayer_FPIX[0]);
						nRecoHitsLayer2_FPIX_pos.push_back(nHitsPerLayer_FPIX[1]);
						nRecoHitsLayer3_FPIX_pos.push_back(nHitsPerLayer_FPIX[2]);
						nRecoHitsLayer1_FPIX_neg.push_back(nHitsPerLayer_FPIX[3]);
						nRecoHitsLayer2_FPIX_neg.push_back(nHitsPerLayer_FPIX[4]);
						nRecoHitsLayer3_FPIX_neg.push_back(nHitsPerLayer_FPIX[5]);

						nRecoHitLayersBPIX.push_back(nLayersBPIX);
						nRecoHitLayersFPIX.push_back(nLayersFPIX);
					}

					recoEle_pt.push_back( eleItr->pt() );
					recoEle_eta.push_back( eleItr->eta() );
					recoEle_phi.push_back( eleItr->phi() );
					recoEle_E.push_back( eleItr->energy() );

				}
			} 
		}
	}

	// Fill if at least one sim electron is found
	if(simEle_pt.size()>0)
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
