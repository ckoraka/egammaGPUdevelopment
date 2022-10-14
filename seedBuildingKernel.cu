// Kernel that does the seed building / 
// Everything should be done in here since kernel launches are costly
// Maybe use device only kernels for functions?

/* Info to be stored in SoA seed struct */
//nhits
// Collection of rechits per seed
//rechits.isValid
//rechits.globalPosition
//rechits.geographicalID

__global__  void processSeed(seedInfo, ecalSCinfo,trajState,matchedSeeds*){


// Easy to pass  - in the Ecal SC struct
eta =  ecalSCinfo.eta
et = ecalSCinfo.energy * sin(ecalSCinfo.theta)

// this is 1 or -1
trajState.charge

// Parameter from config - should be passed as kernel argument
const auto nCuts = cfg_.matchingCuts.size(); 
const auto enableHitSkipping  = cfg_.enableHitSkipping
//These I think are the same for all seeds 
//////


//=== Loop over seeds 

// Create a vector of stucts or cuda friendly object to store matches
matches [1]
//Free trajStare object
firstMatchFreeTraj
//And two Global point structures
prevHitPosition
vertex 


//*** Loop over seed hits
for (size_t iHit = 0;matches.size() < nCuts && iHit < seed.nHits() && (cfg_.enableHitSkipping || iHit == matches.size());iHit++) {

    //access rechit --> Should be part of seed struct
    auto const& recHit = *(seed.recHits().begin() + iHit);

    //check if rechit is valid
    if (!recHit.isValid())
        continue;
    // Add var that check if this is the first attempt to perfrom a matching
    const bool doFirstMatch = matches.empty();

    // Check if the first match should be made 
    auto const& trajState = doFirstMatch? getTrajStateFromVtx(recHit, initialTrajState, backwardPropagator_)  : getTrajStateFromPoint(recHit, firstMatchFreeTraj, prevHitPos, forwardPropagator_);

    // If true -> getTrajStateFromVtx
    // If false -> getTrajStateFromPoint
    // 
    
    *Propagator* https://cmssdt.cern.ch/dxr/CMSSW/source/TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h#25
 https://cmssdt.cern.ch/dxr/CMSSW/source/RecoEgamma/EgammaElectronAlgos/src/TrajSeedMatcher.cc#115
    enum PropagationDirection { oppositeToMomentum, alongMomentum, anyDirection, invalidDirection };
    https://cmssdt.cern.ch/dxr/CMSSW/source/MagneticField/Engine/interface/MagneticField.h#19
    https://cmssdt.cern.ch/dxr/CMSSW/source/TrackingTools/MaterialEffects/src/PropagatorWithMaterial.cc





//***


//===


}


From initial trajectory :
trajState(ecalSCinfo,charge)

charge 



[1] 
  struct SCHitMatch {

    const DetId detId = 0;

    const GlobalPoint hitPos; <-- this we can make a struct 

    const float dRZ = std::numeric_limits<float>::max();
    const float dPhi = std::numeric_limits<float>::max();

    const TrackingRecHit& hit; <--- Check this

    const float et = 0.f;
    const float eta = 0.f;
    const float phi = 0.f;
    const int charge = 0;
    const int nrClus = 0;
  };

