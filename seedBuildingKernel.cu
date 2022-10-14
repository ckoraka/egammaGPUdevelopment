

// Kernel that does the seed building / 
// Everything should be done in here since kernel launches are costly
// Mayve use device only kernels for functions?


seedInfo 
nhits
rechits.isValid
rechits.globalPosition
rechits.geographicalID


__global__  void processSeed(seedInfo, ecalSCinfo,trajState,matchedSeeds*){


// Easy to pass  - in the Ecal SC struct
eta =  ecalSCinfo.eta
et = ecalSCinfo.energy * sin(ecalSCinfo.theta)
trajState.charge
// Parameter from config - should be passed as kernel argument
const auto nCuts = cfg_.matchingCuts.size(); 
const auto enableHitSkipping  = cfg_.enableHitSkipping
//These I think are the same for all seeds 
//////


//=== Loop over seeds 

// Create a ?vector? or cuda friendly object to store matches
matches
//Free trajStare object
firstMatchFreeTraj
//And two Global point structures
prevHitPosition
vertex 


//*** Loop over seed hits



//***


//===


}


From initial trajectory :
trajState(ecalSCinfo,charge)

charge 
