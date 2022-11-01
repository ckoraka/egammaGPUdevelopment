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


//Need to replace
/*
TsosWP with cuda compatable data structure I think tsos will be replaced by a data structure  
*/

//Information needed for PropogateWithPath
/*
Class Inputs: propagation direction, mass and magnetic field are sepcified
Defualt parameters: maxDphi = 1.6, useRungeKutta = false, ptMin = -1, useOldGeoPropLogic = true)
Objects created by the class: AnalyticalPropagator()
Arguements for propagate with path: fts and surface (either plane or cylinder) these come from the initial state and the surface type of hit.det()


*/

//AnalyticalPropagator
/*
Class inputs: Magnetic field (passed in from MaterialPropagator), PropagationDir dir (Also passed in from MaterailPropagator), maxDPhi = 1.6, isOld = True.
PropogateWithPath Method inputs: fts, and surface (either cylinder or plane)


*/

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
    
    //vtx vs point check (difference between bacwardPropagator_ and forwardPropagator_)
    if (doFirstMatch){
      auto propDir = oppositeToMomentum;
    }else{
      auto propDir = alongMomentum;
    }
    
    //Attempt to reproduce PropagatorWithMaterial::propagateWithPath methods
    //
    // add material at starting surface, if requested
    //
    TsosWP newTsosWP(tsos, 0.);
    if (materialAtSource()) {
      bool updateOk = theMEUpdator->updateStateInPlace(newTsosWP.first, propDir); //MEupdator.updateStateInPlace needs to be implemented
      // MEUpdator functionality reproduced
       
      if UNLIKELY (!updateOk)
        newTsosWP.first = TrajectoryStateOnSurface();
    }
    if UNLIKELY (!newTsosWP.first.isValid())
      return newTsosWP; //Break?
    //
    // geometrical propagation (Analytic propagators)
    //
    if (dynamic_cast<Cylinder*>(recHit.det()->surface()) == nullptr){ //Planar Geometry
      // check curvature
      float rho = fts.transverseCurvature();
      // propagate parameters
      GlobalPoint x;
      GlobalVector p;
      double s;
      // check if already on plane
      if LIKELY (plane.localZclamped(fts.position()) != 0) {
        // propagate
        bool parametersOK = this->propagateParametersOnPlane(fts, plane, x, p, s); //This method needs to be implemented
        // check status and deltaPhi limit
        float dphi2 = float(s) * rho;
        dphi2 = dphi2 * dphi2 * fts.momentum().perp2();
        if UNLIKELY (!parametersOK || dphi2 > theMaxDPhi2 * fts.momentum().mag2())
          newTsosWP = TsosWP(TrajectoryStateOnSurface(), 0.);
        } else {
          LogDebug("AnalyticalPropagator") << "not going anywhere. Already on surface.\n"
                                           << "plane.localZ(fts.position()): " << plane.localZ(fts.position()) << "\n"
                                           << "plane.position().mag(): " << plane.position().mag() << "\n"
                                           << "plane.posPrec: " << plane.posPrec(); //Is the log statement realy needed?

          x = fts.position();
          p = fts.momentum();
          s = 0;
        }
        //
        // Compute propagated state and check change in curvature
        //
        GlobalTrajectoryParameters gtp(x, p, fts.charge(), theField);
        if UNLIKELY (std::abs(gtp.transverseCurvature() - rho) > theMaxDBzRatio * std::abs(rho))
          newTsosWP = TsosWP(TrajectoryStateOnSurface(), 0.);
        newTsosWP = propagatedStateWithPath(fts, plane, gtp, s);
    }
    if (dynamic_cast<Plane*>(recHit.det()->surface()) == nullptr){ //Cylindrical Geometry
      // check curvature
      auto rho = fts.transverseCurvature();
      // propagate parameters
      GlobalPoint x;
      GlobalVector p;
      double s = 0;
      bool parametersOK = this->propagateParametersOnCylinder(fts, cylinder, x, p, s); //This method needs to be implemented
      // check status and deltaPhi limit
      float dphi2 = s * rho;
      dphi2 = dphi2 * dphi2 * fts.momentum().perp2();

      if UNLIKELY (!parametersOK || dphi2 > theMaxDPhi2 * fts.momentum().mag2())
        newTsosWP = TsosWP(TrajectoryStateOnSurface(), 0.);

      GlobalTrajectoryParameters gtp(x, p, fts.charge(), theField);

      if UNLIKELY (std::abs(gtp.transverseCurvature() - rho) > theMaxDBzRatio * std::abs(rho))
        newTsosWP = TsosWP(TrajectoryStateOnSurface(), 0.);

      ConstReferenceCountingPointer<TangentPlane> plane(
          cylinder.tangentPlane(x));  // need to be here until tsos is created!
      newTsosWP = propagatedStateWithPath(fts, *plane, gtp, s); //This needs to be implemented
      }

      if UNLIKELY (!(newTsosWP.first).isValid() || materialAtSource())
        return newTsosWP;
      //
      // add material at destination surface, if requested
      //
      bool updateOk = theMEUpdator->updateStateInPlace(
          newTsosWP.first, PropagationDirectionFromPath()(newTsosWP.second, propDir));
      if UNLIKELY (!updateOk)
        newTsosWP.first = TrajectoryStateOnSurface(); 
      return newTsosWP;

    


 
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




  std::pair<TrajectoryStateOnSurface, double> AnalyticalPropagator::propagatedStateWithPath(
        const FreeTrajectoryState& fts,
        const Surface& surface,
        const GlobalTrajectoryParameters& gtp,
        const double& s) const {
      //
      // for forward propagation: state is before surface,
      // for backward propagation: state is after surface
      //
      SurfaceSide side =
          PropagationDirectionFromPath()(s, propagationDirection()) == alongMomentum ? beforeSurface : afterSurface;
      //
      //
      // error propagation (if needed) and conversion to a TrajectoryStateOnSurface
      //
      if (fts.hasError()) {
        //
        // compute jacobian
        //
        AnalyticalCurvilinearJacobian analyticalJacobian(fts.parameters(), gtp.position(), gtp.momentum(), s);
        const AlgebraicMatrix55& jacobian = analyticalJacobian.jacobian();
        // CurvilinearTrajectoryError cte(ROOT::Math::Similarity(jacobian, fts.curvilinearError().matrix()));
        return TsosWP(
            TrajectoryStateOnSurface(gtp, ROOT::Math::Similarity(jacobian, fts.curvilinearError().matrix()), surface, side),
            s);
      } else {
        //
        // return state without errors
        //
        return TsosWP(TrajectoryStateOnSurface(gtp, surface, side), s);
      }
    }

