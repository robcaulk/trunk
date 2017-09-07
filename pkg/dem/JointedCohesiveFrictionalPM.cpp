/* LucScholtes2010  */

#include"JointedCohesiveFrictionalPM.hpp"
#include<core/Scene.hpp>
#include<pkg/dem/ScGeom.hpp>
#include<core/Omega.hpp>
#include<pkg/common/Sphere.hpp>

YADE_PLUGIN((JCFpmMat)(JCFpmState)(JCFpmPhys)(Ip2_JCFpmMat_JCFpmMat_JCFpmPhys)(Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM));


/********************** Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM ****************************/
CREATE_LOGGER(Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM);

static boost::mutex nearbyInts_mutex;
static boost::mutex clusterInts_mutex;

bool Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact){

	const int &id1 = contact->getId1();
	const int &id2 = contact->getId2();
	ScGeom* geom = static_cast<ScGeom*>(ig.get()); 
	JCFpmPhys* phys = static_cast<JCFpmPhys*>(ip.get());
	
	Body* b1 = Body::byId(id1,scene).get();
	Body* b2 = Body::byId(id2,scene).get();

	Real Dtensile=phys->FnMax/phys->kn;
	
	string fileCracks = "cracks_"+Key+".txt";
	string fileMoments = "moments_"+Key+".txt";
	/// Defines the interparticular distance used for computation
	Real D = 0;

	/*this is for setting the equilibrium distance between all cohesive elements in the first contact detection*/
	if ( contact->isFresh(scene) ) { 
	  phys->normalForce = Vector3r::Zero(); 
	  phys->shearForce = Vector3r::Zero();
	  if ((smoothJoint) && (phys->isOnJoint)) {
	    phys->jointNormal = geom->normal.dot(phys->jointNormal)*phys->jointNormal; //to set the joint normal colinear with the interaction normal
	    phys->jointNormal.normalize();
	    phys->initD = std::abs((b1->state->pos - b2->state->pos).dot(phys->jointNormal)); // to set the initial gap as the equilibrium gap
	  } else { 
	    phys->initD = geom->penetrationDepth; 
	  }
	}
	
	if ( smoothJoint && phys->isOnJoint ) {
	  if ( phys->more || ( phys-> jointCumulativeSliding > (2*min(geom->radius1,geom->radius2)) ) ) { 
	    if (!neverErase) return false; 
	    else {
	      phys->shearForce = Vector3r::Zero();
	      phys->normalForce = Vector3r::Zero();
	      phys->isCohesive =0;
	      phys->FnMax = 0;
	      phys->FsMax = 0;
	      return true; // do we need this? -> yes if it ends the loop (avoid the following calculations)
	      }
	  } else { 
	    D = phys->initD - std::abs((b1->state->pos - b2->state->pos).dot(phys->jointNormal)); 
	  }
	} else { 
	  D = geom->penetrationDepth - phys->initD; 
	}
	
	phys->crackJointAperture = D<0? -D : 0.; // for DFNFlow
	phys->separation = D; // track the separation between particles
	if (!phys->momentBroken) phys->strainEnergy = 0.5*((pow(phys->normalForce.norm(),2)/phys->kn) + (pow(phys->shearForce.norm(),2)/phys->ks));

//for clustered events:
	if (phys->momentBroken && recordMoments && !phys->momentCalculated){
		//if (!phys->checkedForCluster) checkForCluster(phys, geom, b1, b2, contact);
		//if (!phys->interactionsAdded) clusterInteractions(phys, contact);
		if (phys->originalClusterEvent && !phys->computedCentroid) computeCentroid(phys);
		if (phys->originalClusterEvent) computeClusteredMoment(phys);
		
		if (phys->momentCalculated){
			std::ofstream file (fileMoments.c_str(), !momentsFileExist ? std::ios::trunc : std::ios::app);
			if(file.tellp()==0){ file <<"i p0 p1 p2 moment numInts eventNum"<<endl; }
			file << boost::lexical_cast<string> ( scene->iter )<<" "<< boost::lexical_cast<string> ( phys->momentCentroid[0] ) <<" "<< boost::lexical_cast<string> ( phys->momentCentroid[1] ) <<" "<< boost::lexical_cast<string> ( phys->momentCentroid[2] ) <<" "<< boost::lexical_cast<string> ( phys->momentMagnitude ) << " " << boost::lexical_cast<string> ( phys->clusterInts.size() ) << " " << boost::lexical_cast<string> ( phys->eventNumber ) << endl;
			momentsFileExist=true;
		}
		
	}

	/* Determination of interaction */
	if (D < 0) { //spheres do not touch 
	  if ( !phys->isCohesive) {
	    if (!neverErase) return false;
	    else {
	      phys->shearForce = Vector3r::Zero();
	      phys->normalForce = Vector3r::Zero();
	      phys->isCohesive =0;
	      phys->FnMax = 0;
	      phys->FsMax = 0;
	      return true; // do we need this? not sure -> yes, it ends the loop (avoid the following calculations)
	    }
	  }
	  
	  if ( phys->isCohesive && (phys->FnMax>0) && (std::abs(D)>Dtensile) ) {
	    
            nbTensCracks++;

	    // update body state with the number of broken bonds
	    JCFpmState* st1=dynamic_cast<JCFpmState*>(b1->state.get());
	    JCFpmState* st2=dynamic_cast<JCFpmState*>(b2->state.get());
            phys->breakType = 0;
	    phys->breakOccurred = true;
	    st1->tensBreak+=1;
	    st2->tensBreak+=1;
	    st1->tensBreakRel+=1.0/st1->noIniLinks;
	    st2->tensBreakRel+=1.0/st2->noIniLinks;
            phys->momentBroken = true;

            Real scalarNF=phys->normalForce.norm();
	    Real scalarSF=phys->shearForce.norm();
            totalTensCracksE+=0.5*( ((scalarNF*scalarNF)/phys->kn) + ((scalarSF*scalarSF)/phys->ks) );
	    
    	    // create a text file to record properties of the broken bond (iteration, position, type (tensile), cross section and contact normal orientation)
	    if (recordCracks){
	      std::ofstream file (fileCracks.c_str(), !cracksFileExist ? std::ios::trunc : std::ios::app);
	      if(file.tellp()==0){ file <<"i p0 p1 p2 t s norm0 norm1 norm2 onFrac nrg"<<endl; }
	      Vector3r crackNormal=Vector3r::Zero();
	      if ((smoothJoint) && (phys->isOnJoint)) { crackNormal=phys->jointNormal; } else {crackNormal=geom->normal;}
	      file << boost::lexical_cast<string> ( scene->iter )<<" "<< boost::lexical_cast<string> ( geom->contactPoint[0] ) <<" "<< boost::lexical_cast<string> ( geom->contactPoint[1] ) <<" "<< boost::lexical_cast<string> ( geom->contactPoint[2] ) <<" "<< 0 <<" "<< boost::lexical_cast<string> ( 0.5*(geom->radius1+geom->radius2) ) <<" "<< boost::lexical_cast<string> ( crackNormal[0] ) <<" "<< boost::lexical_cast<string> ( crackNormal[1] ) <<" "<< boost::lexical_cast<string> ( crackNormal[2] ) << " " << boost::lexical_cast<string> ( phys->onFracture) << " " <<boost::lexical_cast<string> ( 0.5*( ((scalarNF*scalarNF)/phys->kn) + ((scalarSF*scalarSF)/phys->ks) ) ) << endl;
	    
            }
            cracksFileExist=true;
	    /// Timos
	    if (!neverErase) return false; 
	    else {
	      phys->shearForce = Vector3r::Zero();
	      phys->normalForce = Vector3r::Zero();
	      phys->isCohesive =0;
	      phys->FnMax = 0;
	      phys->FsMax = 0;
	      phys->isBroken = true;
	      return true; // do we need this? not sure -> yes, it ends the loop (avoid the following calculations)
	    }

            if (recordMoments && !phys->momentCalculated){
                checkForCluster(phys, geom, b1, b2, contact);
		clusterInteractions(phys, contact);
            }

	  }
	}
	
	/* NormalForce */
	Real Fn = 0;
	Fn = phys->kn*D; 

	/* ShearForce */
	Vector3r& shearForce = phys->shearForce; 
	Real jointSliding=0;

	if ((smoothJoint) && (phys->isOnJoint)) {
	  
	  /// incremental formulation (OK?)
	  Vector3r relativeVelocity = (b2->state->vel - b1->state->vel); // angVel are not taken into account as particles on joint don't rotate ????
	  Vector3r slidingVelocity = relativeVelocity - phys->jointNormal.dot(relativeVelocity)*phys->jointNormal; 
	  Vector3r incrementalSliding = slidingVelocity*scene->dt;
	  shearForce -= phys->ks*incrementalSliding;
	  
	  jointSliding = incrementalSliding.norm();
	  phys->jointCumulativeSliding += jointSliding;
  
	} else {

	  shearForce = geom->rotate(phys->shearForce);
	  const Vector3r& incrementalShear = geom->shearIncrement();
	  shearForce -= phys->ks*incrementalShear;
	  
	}
	
	/* Mohr-Coulomb criterion */
	Real maxFs = phys->FsMax + Fn*phys->tanFrictionAngle;
	Real scalarShearForce = shearForce.norm();	  
	if (scalarShearForce > maxFs) {	
	  if (scalarShearForce != 0)
	    shearForce*=maxFs/scalarShearForce;
	  else
	    shearForce=Vector3r::Zero();
	  if ((smoothJoint) && (phys->isOnJoint)) {phys->dilation=phys->jointCumulativeSliding*phys->tanDilationAngle-D; phys->initD+=(jointSliding*phys->tanDilationAngle);}
	  // take into account shear cracking -> are those lines critical? -> TODO testing with and without
	  if ( phys->isCohesive ) { 
            nbShearCracks++; 

	    // update body state with the number of broken bonds
	    JCFpmState* st1=dynamic_cast<JCFpmState*>(b1->state.get());
	    JCFpmState* st2=dynamic_cast<JCFpmState*>(b2->state.get());
            phys->breakOccurred = true;
            phys->breakType=1; // for DFNFlow
	    st1->shearBreak+=1;
	    st2->shearBreak+=1;
	    st1->shearBreakRel+=1.0/st1->noIniLinks;
	    st2->shearBreakRel+=1.0/st2->noIniLinks;
            phys->momentBroken = true;
	    Real scalarNF=phys->normalForce.norm();
	    Real scalarSF=phys->shearForce.norm();
            totalShearCracksE+=0.5*( ((scalarNF*scalarNF)/phys->kn) + ((scalarSF*scalarSF)/phys->ks) );

	    // extend smooth joint model
	    if (phys->breakType==1 && extendSmoothJoint && !phys->isOnJoint) orientJointNormal(phys, geom, b1, b2, contact);
	    // create a text file to record properties of the broken bond (iteration, position, type (shear), cross section and contact normal orientation)
	    if (recordCracks){
	      std::ofstream file (fileCracks.c_str(), !cracksFileExist ? std::ios::trunc : std::ios::app);
	      if(file.tellp()==0){ file <<"i p0 p1 p2 t s norm0 norm1 norm2 onFrac nrg"<<endl; }
	      Vector3r crackNormal=Vector3r::Zero();
	      if ((smoothJoint) && (phys->isOnJoint)) { crackNormal=phys->jointNormal; } else {crackNormal=geom->normal;}
	      file << boost::lexical_cast<string> ( scene->iter )<<" "<< boost::lexical_cast<string> ( geom->contactPoint[0] ) <<" "<< boost::lexical_cast<string> ( geom->contactPoint[1] ) <<" "<< boost::lexical_cast<string> ( geom->contactPoint[2] ) <<" "<< 1 <<" "<< boost::lexical_cast<string> ( 0.5*(geom->radius1+geom->radius2) ) <<" "<< boost::lexical_cast<string> ( crackNormal[0] ) <<" "<< boost::lexical_cast<string> ( crackNormal[1] ) <<" "<< boost::lexical_cast<string> ( crackNormal[2] ) << " " << boost::lexical_cast<string> ( phys->onFracture) << " " <<boost::lexical_cast<string> ( 0.5*( ((scalarNF*scalarNF)/phys->kn) + ((scalarSF*scalarSF)/phys->ks) ) ) << endl;
	    	    
            }
	    cracksFileExist=true;
	    // set the contact properties to friction if in compression, delete contact if in tension
	    phys->isBroken = true;
	    phys->isCohesive = 0;
	    phys->FnMax = 0;
	    phys->FsMax = 0;
           
// 	    shearForce *= Fn*phys->tanFrictionAngle/scalarShearForce; // now or at the next timestep?
	    if ( D < 0 ) { // spheres do not touch
	      if (!neverErase) return false;
	      else {
		phys->shearForce = Vector3r::Zero();
		phys->normalForce = Vector3r::Zero();
		return true; // do we need this? not sure -> yes, it ends the loop (avoid the following calculations)
	      }

	    }
            if (recordMoments && !phys->momentCalculated){
                checkForCluster(phys, geom, b1, b2, contact);
		clusterInteractions(phys, contact);
            }
	  }
	}

	/* Apply forces */
	if ((smoothJoint) && (phys->isOnJoint)) { phys->normalForce = Fn*phys->jointNormal; } else { phys->normalForce = Fn*geom->normal; }
	
	Vector3r f = phys->normalForce + shearForce;
	
	/// applyForceAtContactPoint computes torque also and, for now, we don't want rotation for particles on joint (some errors in calculation due to specific geometry) 
 	//applyForceAtContactPoint(f, geom->contactPoint, I->getId2(), b2->state->pos, I->getId1(), b1->state->pos, scene);
	scene->forces.addForce (id1,-f);
	scene->forces.addForce (id2, f);
	
	// simple solution to avoid torque computation for particles interacting on a smooth joint 
	if ( (phys->isOnJoint)&&(smoothJoint) ) return true;
	
	/// those lines are needed if rootBody->forces.addForce and rootBody->forces.addMoment are used instead of applyForceAtContactPoint -> NOTE need to check for accuracy!!!
	scene->forces.addTorque(id1,(geom->radius1-0.5*geom->penetrationDepth)* geom->normal.cross(-f));
	scene->forces.addTorque(id2,(geom->radius2-0.5*geom->penetrationDepth)* geom->normal.cross(-f));
	return true;
	
}


// function used to parse through new interactions and only add unique ints to cluster list	
void Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM::addUniqueIntsToList(JCFpmPhys* phys, JCFpmPhys* nearbyPhys){
	unsigned int size = phys->nearbyInts.size();
	for (unsigned int i=0; i<nearbyPhys->nearbyInts.size(); i++){
		if (!nearbyPhys->nearbyInts[i]) continue;
		bool pushBack = true;
		for (unsigned int j=0; j<size; j++){
			if (!phys->nearbyInts[j]) continue;
			if (phys->nearbyInts[j] == nearbyPhys->nearbyInts[i]) {
				pushBack = false;
				break;
			}
		}
		boost::mutex::scoped_lock lock(nearbyInts_mutex);
		if (pushBack && nearbyPhys->nearbyInts[i]){
			phys->nearbyInts.push_back(nearbyPhys->nearbyInts[i]);
			JCFpmPhys* strainPhys = YADE_CAST<JCFpmPhys*> (nearbyPhys->nearbyInts[i]->phys.get()); 
			phys->momentStrainEnergy += strainPhys->strainEnergy;  // while we are here we update the reference strain by adding the strain of the newly added ints
		}
			
	}
}


// function used to check if the newly broken bond belongs in a cluster or not, if so attach to proper cluster and set appropriate flags
void Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM::checkForCluster(JCFpmPhys* phys, ScGeom* geom, Body* b1, Body* b2, Interaction* contact){

	const shared_ptr<Shape>& sphere1 = b1->shape;
	const shared_ptr<Shape>& sphere2 = b2->shape;
	const Real sphereRadius1 = static_cast<Sphere*>(sphere1.get())->radius;
	const Real sphereRadius2 = static_cast<Sphere*>(sphere2.get())->radius;
	const Real avgDiameter = (sphereRadius1+sphereRadius2);
	Vector3r& brokenInteractionLocation = geom->contactPoint;
	phys->nearbyFound=0;

	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
	//#endif
		JCFpmPhys* nearbyPhys;
		const ScGeom* nearbyGeom;
		if (!I || !I->geom.get() || !I->phys.get()) continue;
		if (I && I->isReal() && JCFpmPhys::getClassIndexStatic()==I->phys->getClassIndex()){
			nearbyPhys = YADE_CAST<JCFpmPhys*>(I->phys.get());
			if (!nearbyPhys) continue;
			if (I->geom.get() /*&& !nearbyPhys->momentBroken*/){
				nearbyGeom = YADE_CAST<ScGeom*> (I->geom.get());
                    		if (!nearbyGeom) continue;
				Vector3r nearbyInteractionLocation = nearbyGeom->contactPoint;
				Vector3r proximityVector = nearbyInteractionLocation-brokenInteractionLocation;
				Real proximity = proximityVector.norm();
				
				// this logic is finding interactions within a radius of the broken one, and identifiying if it is an original event of if it belongs in a cluster
				if (proximity < avgDiameter*momentRadiusFactor && proximity != 0){
					if (nearbyPhys->originalClusterEvent && !nearbyPhys->momentCalculated && !phys->clusteredEvent && clusterMoments) {
						phys->eventNumber = nearbyPhys->eventNumber; 
						phys->clusteredEvent = true;
						phys->originalEvent = I.get();
					} else if (nearbyPhys->clusteredEvent && !phys->clusteredEvent && clusterMoments){
						JCFpmPhys* originalPhys = YADE_CAST<JCFpmPhys*>(nearbyPhys->originalEvent->phys.get());
						if (!originalPhys->momentCalculated){
							phys->eventNumber = nearbyPhys->eventNumber;
							phys->clusteredEvent = true;
							phys->originalEvent = nearbyPhys->originalEvent;
						}
					} 

					if (nearbyPhys->momentBroken) continue;
					phys->nearbyInts.push_back(I.get());
				}
			}
		}
	}
	if (!phys->clusteredEvent) {
		phys->originalClusterEvent = true; // if not clustered, its an original event. We use this interaction as the master for the cluster. Its list of nearbyInts will expand if other ints break nearby. 
		phys->originalEvent = contact;
		eventNumber += 1;
		phys->eventNumber = eventNumber;
	}
	phys->checkedForCluster = true;
}

// function used for clustering broken bonds and nearby interactions
void Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM::clusterInteractions(JCFpmPhys* phys, Interaction* contact){
	JCFpmPhys* originalPhys = YADE_CAST<JCFpmPhys*>(phys->originalEvent->phys.get());
	addUniqueIntsToList(originalPhys, phys);  //NEED TO PUSHBACK ONLY UNIQUE INTS. we don't want a list with duplicate events (also updates reference strain)
	phys->interactionsAdded = true;
	originalPhys->elapsedIter = 1;  // reset the temporal window? do we want this?
	//originalPhys->firstMomentCalc=true; // do we need a new reference strain energy for the calculation of strain energy change?
	phys->momentMagnitude = 0; // dirty way to avoid recording these clustered events twice? maybe dont need this if proper recording is applied
	originalPhys->computedCentroid=false;  // set flag to compute a new centroid since we added more ints
	boost::mutex::scoped_lock lock(clusterInts_mutex);
	originalPhys->clusterInts.push_back(contact);  // add this broken interaction to list of broken bonds in this event
} 


// function cycles through the list of interactions associated with a cluster and computes moment
void Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM::computeClusteredMoment(JCFpmPhys* phys){
	Real totalMomentStrainEnergy = 0;
	Real momentStrainEnergyChange = 0;
	for (unsigned int i=0; i<phys->nearbyInts.size(); i++){
		const JCFpmPhys* nearbyPhys;
		if (!phys->nearbyInts[i] || !phys->nearbyInts[i]->geom.get() || !phys->nearbyInts[i]->phys.get()) continue;
		nearbyPhys = YADE_CAST<JCFpmPhys*>(phys->nearbyInts[i]->phys.get());
		if (!nearbyPhys) continue;
		totalMomentStrainEnergy += nearbyPhys->strainEnergy;
	}
	if(phys->firstMomentCalc){
		phys->momentStrainEnergy = totalMomentStrainEnergy;
		phys->firstMomentCalc = false;
	}
	momentStrainEnergyChange = totalMomentStrainEnergy - phys->momentStrainEnergy;
	phys->elapsedIter += 1;
	if (momentStrainEnergyChange > phys->momentStrainEnergyChange) phys->momentStrainEnergyChange = momentStrainEnergyChange;
	if (phys->elapsedIter >= elapsedIterFactor*momentRadiusFactor){ // the elapsed time should reflect 20*particlediameters radius Hazzard and Damjanac 2013
		phys->momentCalculated=true;
		phys->originalClusterEvent=false; // this event no longer exists, so we need to allow other new events to occur nearby.    
		//cout << "MomentEnergyChange" << momentStrainEnergyChance << endl;
		if(phys->momentStrainEnergyChange!=0) phys->momentMagnitude = (2./3.)*log(phys->momentStrainEnergyChange*momentStrainFudgeFactor)-3.2;  //empirical equation for energy magnitude (Hazzard and Damjanac 2013) 
		//if(phys->momentStrainEnergyChange==0) cout<<"avgDiameter " << avgDiameter << " found nearby interaciton? " << phys->nearbyFound << "over " << phys->elapsedIter << " iterations" <<endl;  // debugging. It appears some strain energy searches yield decreases of strain energy in neighbor hood. We are handling these by assining a 0 magnitude...but that seems wrong. Turns out these are just very very small events. There is actually no change of strain around them. Due to weibull dist int areas?
	}
					
}


// function computes the centroid of a cluster
void Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM::computeCentroid(JCFpmPhys* phys){
	JCFpmPhys* originalPhys = YADE_CAST<JCFpmPhys*>(phys->originalEvent->phys.get());
	Vector3r summedLocations = Vector3r::Zero();
	for (unsigned int i=0; i<originalPhys->clusterInts.size(); i++){
		ScGeom* nearbyGeom;
		if (!originalPhys->clusterInts[i]) continue;
		if (originalPhys->clusterInts[i]->geom.get()){
			nearbyGeom = YADE_CAST<ScGeom*> (originalPhys->clusterInts[i]->geom.get());
			Vector3r nearbyInteractionLocation = nearbyGeom->contactPoint;
			summedLocations += nearbyInteractionLocation;
		}
	}
	originalPhys->momentCentroid = summedLocations/originalPhys->clusterInts.size(); // new location of event is average of all clustered events
	originalPhys->computedCentroid = true;
	
}

// function for computing moments based on single breaks
//void Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM::computeMoment(JCFpmPhys* phys, ScGeom* geom, Body* b1, Body* b2, Interaction* contact){

//	const shared_ptr<Shape>& sphere1 = b1->shape;
//	const shared_ptr<Shape>& sphere2 = b2->shape;
//	const Real sphereRadius1 = static_cast<Sphere*>(sphere1.get())->radius;
//	const Real sphereRadius2 = static_cast<Sphere*>(sphere2.get())->radius;
//	const Real avgDiameter = (sphereRadius1+sphereRadius2);
//	Vector3r& brokenInteractionLocation = geom->contactPoint;
//	Real totalMomentStrainEnergy = 0;
//	Real momentStrainEnergyChange = 0;
//	phys->nearbyFound=0;
//	

//	// search through all interactions to find nearby interactions for energy calc
//	//#ifdef YADE_OPENMP
//	//int nIntr=(int)scene->interactions->size();
//	//#pragma omp parallel for
//	//for (int j=0; j<nIntr; j++){
//	//	const shared_ptr<Interaction>& I = (*scene->interactions)[j];
//	//#else 
//	if(phys->firstMomentCalc){
//	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
//	//#endif
//		const JCFpmPhys* nearbyPhys;
//		const ScGeom* nearbyGeom;
//		if (!I || !I->geom.get() || !I->phys.get()) continue;
//		if (I && I->isReal() && JCFpmPhys::getClassIndexStatic()==I->phys->getClassIndex()){
//			nearbyPhys = YADE_CAST<JCFpmPhys*>(I->phys.get());
//			if (!nearbyPhys) continue;
//			if (I->geom.get() && !nearbyPhys->momentBroken){		//nearbyPhys->isOnJoint && && !nearbyPhys->onFracture	
//				nearbyGeom = YADE_CAST<ScGeom*> (I->geom.get());
//				//cout<<"nearbyPGeom declared"<<endl;
//                    		if (!nearbyGeom) continue;
//				Vector3r nearbyInteractionLocation = nearbyGeom->contactPoint;
//				//cout<<"nearbyInteraction location"<< nearbyInteractionLocation << endl;
//				Vector3r proximityVector = nearbyInteractionLocation-brokenInteractionLocation;
//				Real proximity = proximityVector.norm();
//				//cout<<"proximity computed"<< proximity << endl;
//				//cout << "allowable radius" << avgDiameter*momentRadiusFactor << endl;
//				
//				if (proximity < avgDiameter*momentRadiusFactor && proximity != 0){
//					totalMomentStrainEnergy += nearbyPhys->strainEnergy;
//					phys->nearbyInts.push_back(I);
//					//cout << "nearbyInts set in vector" << endl;					
//				}
//			}
//		}
//	}
//	} else {
//		for (unsigned int i=0; i<phys->nearbyInts.size(); i++){
//			const JCFpmPhys* nearbyPhys;
//			if (!phys->nearbyInts[i] || !phys->nearbyInts[i]->geom.get() || !phys->nearbyInts[i]->phys.get()) continue;
//			nearbyPhys = YADE_CAST<JCFpmPhys*>(phys->nearbyInts[i]->phys.get());
//			if (!nearbyPhys) continue; 
//			totalMomentStrainEnergy += nearbyPhys->strainEnergy;
//			//cout<< "nearbyInts indexed" << endl;
//		}
//	}

//	if(phys->firstMomentCalc){
//		phys->momentStrainEnergy = totalMomentStrainEnergy;
//		phys->firstMomentCalc = false;
//	}
//	momentStrainEnergyChange = totalMomentStrainEnergy - phys->momentStrainEnergy;
//	//cout << "Moment Strain Energy Change" << momentStrainEnergyChange << endl; 
//	//cout << "Moment strain energy" << phys->momentStrainEnergy << endl; 
//	//cout << "Total moment strain energy" << totalMomentStrainEnergy << endl;
//	//phys->momentStrainEnergy = totalMomentStrainEnergy;
//	phys->elapsedIter += 1;
//	//cout << "Elapsed iter"<< phys->elapsedIter << endl;
//	if (momentStrainEnergyChange > phys->momentStrainEnergyChange) phys->momentStrainEnergyChange = momentStrainEnergyChange;
//	if (phys->elapsedIter >= elapsedIterFactor*momentRadiusFactor){ // the elapsed time should reflect 20*particlediameters radius Hazzard and Damjanac 2013
//		phys->momentCalculated=true;   
//		//cout << "MomentEnergyChange" << momentStrainEnergyChance << endl;
//		if(phys->momentStrainEnergyChange!=0) phys->momentMagnitude = (2./3.)*log(phys->momentStrainEnergyChange)-3.2;  //empirical equation for energy magnitude (Hazzard and Damjanac 2013) 
//		//if(phys->momentStrainEnergyChange==0) cout<<"avgDiameter " << avgDiameter << " found nearby interaciton? " << phys->nearbyFound << "over " << phys->elapsedIter << " iterations" <<endl;  // debugging. It appears some strain energy searches yield decreases of strain energy in neighbor hood. We are handling these by assining a 0 magnitude...but that seems wrong. 
//	}
//					
//}

// function used for extending smooth joint to newly broken interactions (function works, but results are mixed)
void Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM::orientJointNormal(JCFpmPhys* phys, ScGeom* geom, Body* b1, Body* b2, Interaction* contact){
		
		const shared_ptr<Shape>& sphere1 = b1->shape;
		const shared_ptr<Shape>& sphere2 = b2->shape;
		const Real sphereRadius1 = static_cast<Sphere*>(sphere1.get())->radius;
		const Real sphereRadius2 = static_cast<Sphere*>(sphere2.get())->radius;
		Real maxAllowableDistance = ((sphereRadius1 + sphereRadius2)/2) * fracProximityFactor;
		bool foundClosest = false;
		//bool foundNextClosest = false;
		
		
//			cout<<"Extending smooth joint" <<endl;			
			Vector3r& brokenInteractionLocation = geom->contactPoint;
//			cout<<"contactPoint found" << brokenInteractionLocation << endl;
			Real closestProximity = 10000;  //big number
			//Real nextClosestProximity = 10000;  //big number
			Real planeAngle=0;
//			shared_ptr<Interaction> closestInteraction;
//			shared_ptr<Interaction> nextClosestInteraction;
			Vector3r interactionPlane1 = Vector3r::Zero();
			Vector3r interactionPlane2 = Vector3r::Zero();
			Vector3r closestProximityVector = Vector3r::Zero();

			FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
                    const JCFpmPhys* nearbyPhys;
                    const ScGeom* nearbyGeom;
			if (!I || !I->geom.get() || !I->phys.get()) continue;
			if (I && I->isReal() && JCFpmPhys::getClassIndexStatic()==I->phys->getClassIndex()) {
//				const int &id1 = I->getId1();
//				const int &id2 = I->getId2();
//				Body* b1 = Body::byId(id1,scene).get();
//				Body* b2 = Body::byId(id2,scene).get();
//				JCFpmState* st1=dynamic_cast<JCFpmState*>(Body::byId(I->getId1(),scene)->state.get());
//				JCFpmState* st2=dynamic_cast<JCFpmState*>(Body::byId(I->getId2(),scene)->state.get());

//				Real particleProximity = (st1->pos_get()-st2->pos_get()).norm();
//				if (particleProximity> b1->bound->sweepLength) continue;
//				cout << "particleProximity " << particleProximity << " Sweep length "<< b1->bound->sweepLength << endl;
                nearbyPhys = YADE_CAST<JCFpmPhys*>(I->phys.get());
                if (!nearbyPhys) continue;
//				cout << "nearbyPhys assigned" << endl;				
				if (nearbyPhys->isOnJoint && I->geom.get() && !nearbyPhys->onFracture){
					//cout<<"nearbyPhys is broken and onFracture"<<endl;				

					nearbyGeom = YADE_CAST<ScGeom*> (I->geom.get());
					//cout<<"nearbyPGeom declared"<<endl;
                    			if (!nearbyGeom) continue;
					Vector3r nearbyInteractionLocation = nearbyGeom->contactPoint;
					//cout<<"nearbyInteraction location"<< nearbyInteractionLocation << endl;
					Vector3r proximityVector = nearbyInteractionLocation-brokenInteractionLocation;
					//cout<<"proximityVector computed"<<endl;
					Real proximity = proximityVector.norm();

					//cout<<"proximity computed "<< proximity << "proximityVector " << proximityVector << endl;
					if (proximity < closestProximity && proximity < maxAllowableDistance && proximity!=0) {
						closestProximity = proximity;
						closestProximityVector = proximityVector;
					//	cout<<"closestProximity set"<< closestProximity<<endl;
						//closestInteraction = I;
						
						interactionPlane1 = proximityVector/proximity; // the normalized interaction plane
					//	cout<<"interactionPlane1 computed "<<interactionPlane1<< endl;
						foundClosest = true;
					}

					if (foundClosest && proximity < maxAllowableDistance){ 
						interactionPlane2 = proximityVector/proximity;
						Real angle = acos((interactionPlane1.dot(interactionPlane2))/(closestProximity*proximity));	
						if (angle>planeAngle) planeAngle = angle;
					}
							
			
					
					//if (proximity!=0 && angle>planeAngle) { 
					//	nextClosestProximity = proximity;
					//	cout<<"nextclosestProximity set "<< nextClosestProximity<< endl;
						//nextClosestInteraction = I;
					//	interactionPlane2 = proximityVector/proximity; // normalized interaction plane 
					//	cout<<"interactionPlane2 computed "<<interactionPlane2<<endl;
					//	Vector3r forceFracturePlane = proximityVector - closestProximityVector;
					//	planeAngle = acos((interactionPlane1.dot(interactionPlane2))/(closestProximity*nextClosestProximity))
					//	if (forceFracturePlane.norm() > 0.5*maxAllowableDistance) foundNextClosest = true;
					//}
				}
				}
			}
		if (foundClosest) {
			const shared_ptr<JCFpmMat>& yade1 = YADE_PTR_CAST<JCFpmMat>(b1->material);
			const shared_ptr<JCFpmMat>& yade2 = YADE_PTR_CAST<JCFpmMat>(b2->material);
			JCFpmState* st1=dynamic_cast<JCFpmState*>(Body::byId(contact->getId1(),scene)->state.get());
			JCFpmState* st2=dynamic_cast<JCFpmState*>(Body::byId(contact->getId2(),scene)->state.get());
			st1->onJoint=st2->onJoint=true;
		

			Real jf1 	= yade1->jointFrictionAngle;
			Real jf2 	= yade2->jointFrictionAngle;
			Real jkn1 	= yade1->jointNormalStiffness;
			Real jkn2 	= yade2->jointNormalStiffness;
			Real jks1 	= yade1->jointShearStiffness;
			Real jks2 	= yade2->jointShearStiffness;
			Real jdil1 	= yade1->jointDilationAngle;
			Real jdil2 	= yade2->jointDilationAngle;
//			Real jcoh1 	= yade1->jointCohesion;
//			Real jcoh2 	= yade2->jointCohesion;
//			Real jSigT1	= yade1->jointTensileStrength;
//			Real jSigT2	= yade2->jointTensileStrength;
		
			phys->tanFrictionAngle = std::tan(std::min(jf1,jf2));
			
			//contactPhysics->kn = jointNormalStiffness*2.*R1*R2/(R1+R2); // very first expression from Luc
			//contactPhysics->kn = (jkn1+jkn2)/2.0*2.*R1*R2/(R1+R2); // after putting jointNormalStiffness in material
			phys->kn = ( jkn1 + jkn2 ) /2.0 * phys->crossSection; // for a size independant expression
			phys->ks = ( jks1 + jks2 ) /2.0 * phys->crossSection; // for a size independant expression
			
			phys->tanDilationAngle = std::tan(std::min(jdil1,jdil2));
			phys->isOnJoint=phys->onFracture=1;
			cout<<"set the current interaction to onJoint"<<endl;
			phys->jointNormal=interactionPlane1.cross(interactionPlane2);	
			cout<<"set the interaction jointNormal "<<phys->jointNormal<<endl;
		}
		
}

CREATE_LOGGER(Ip2_JCFpmMat_JCFpmMat_JCFpmPhys);

void Ip2_JCFpmMat_JCFpmMat_JCFpmPhys::go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction){

	/* avoid updates of interaction if it already exists */
	if( interaction->phys ) return; 

	ScGeom* geom=dynamic_cast<ScGeom*>(interaction->geom.get());
	assert(geom);

	const shared_ptr<JCFpmMat>& yade1 = YADE_PTR_CAST<JCFpmMat>(b1);
	const shared_ptr<JCFpmMat>& yade2 = YADE_PTR_CAST<JCFpmMat>(b2);
	JCFpmState* st1=dynamic_cast<JCFpmState*>(Body::byId(interaction->getId1(),scene)->state.get());
	JCFpmState* st2=dynamic_cast<JCFpmState*>(Body::byId(interaction->getId2(),scene)->state.get());
	
	shared_ptr<JCFpmPhys> contactPhysics(new JCFpmPhys()); 
	
	/* From material properties */
	Real E1 	= yade1->young;
	Real E2 	= yade2->young;
	Real v1 	= yade1->poisson;
	Real v2 	= yade2->poisson;
	Real f1 	= yade1->frictionAngle;
	Real f2 	= yade2->frictionAngle;
	Real SigT1	= yade1->tensileStrength;
	Real SigT2	= yade2->tensileStrength;
	Real Coh1	= yade1->cohesion;
	Real Coh2	= yade2->cohesion;

	/* From interaction geometry */
	Real R1= geom->radius1;
	Real R2= geom->radius2;
	
	// control the radius used for cross-sectional area computation
	if (useAvgRadius){	
	contactPhysics->crossSection = Mathr::PI*pow((R1+R2)/2,2);  // use the average radius instead of the minimum
	} else if (totalAvgRadius > 0) {
	contactPhysics->crossSection = Mathr::PI*pow(totalAvgRadius, 2);  // use total average radius
	} else if(xSectionWeibullShapeParameter>0 && xSectionWeibullScaleParameter>0) distributeCrossSectionsWeibull(contactPhysics, R1, R2);	
	else {
	contactPhysics->crossSection = Mathr::PI*pow(min(R1,R2),2); 
	}
	/* Pass values to JCFpmPhys. In case of a "jointed" interaction, the following values will be replaced by other ones later (in few if(){} blocks)*/
	
	// frictional properties
	contactPhysics->kn = 2.*E1*R1*E2*R2/(E1*R1+E2*R2);
	if ( (v1==0)&&(v2==0) )
	  contactPhysics->ks = 0;
	else
	  contactPhysics->ks = 2.*E1*R1*v1*E2*R2*v2/(E1*R1*v1+E2*R2*v2);
	contactPhysics->tanFrictionAngle = std::tan(std::min(f1,f2));

	if (stiffnessWeibullShapeParameter!=0) distributeStiffnesses(contactPhysics);
	// cohesive properties
	///to set if the contact is cohesive or not
	if ( ((cohesiveTresholdIteration < 0) || (scene->iter < cohesiveTresholdIteration)) && (std::min(SigT1,SigT2)>0 || std::min(Coh1,Coh2)>0) && (yade1->type == yade2->type)){ 
	  contactPhysics->isCohesive=true;
	  st1->noIniLinks++;
	  st2->noIniLinks++;
	}
	
	
	if ( contactPhysics->isCohesive ) {
		contactPhysics->FnMax = std::min(SigT1,SigT2)*contactPhysics->crossSection;
	  	contactPhysics->FsMax = std::min(Coh1,Coh2)*contactPhysics->crossSection;
		// stochastic bond strength assignment
		if (yade1->tensileStrengthDeviation>0 && yade1->cohStrengthDeviation>0) distributeStrengthsNormal(contactPhysics, yade1, yade2);
		if (strengthWeibullShapeParameter!=0) distributeStrengthsWeibull(contactPhysics, yade1, yade2);
		if (yade1->tensileStrengthDeviation>0 && strengthWeibullShapeParameter!=0) cerr << "You tried to distribute the strengths according to gaussian and weibull. Using weibull." << endl;
	}

	/// +++ Jointed interactions ->NOTE: geom->normal is oriented from 1 to 2 / jointNormal from plane to sphere 
	if ( st1->onJoint && st2->onJoint )
	{
		if ( (((st1->jointNormal1.cross(st2->jointNormal1)).norm()<0.1) && (st1->jointNormal1.dot(st2->jointNormal1)<0)) || (((st1->jointNormal1.cross(st2->jointNormal2)).norm()<0.1) && (st1->jointNormal1.dot(st2->jointNormal2)<0)) || (((st1->jointNormal1.cross(st2->jointNormal3)).norm()<0.1) && (st1->jointNormal1.dot(st2->jointNormal3)<0)) )
		{
		  contactPhysics->isOnJoint = true;
		  contactPhysics->jointNormal = st1->jointNormal1;
		}
		else if ( (((st1->jointNormal2.cross(st2->jointNormal1)).norm()<0.1) && (st1->jointNormal2.dot(st2->jointNormal1)<0)) || (((st1->jointNormal2.cross(st2->jointNormal2)).norm()<0.1) && (st1->jointNormal2.dot(st2->jointNormal2)<0)) || (((st1->jointNormal2.cross(st2->jointNormal3)).norm()<0.1) && (st1->jointNormal2.dot(st2->jointNormal3)<0)) )
		{
		  contactPhysics->isOnJoint = true;
		  contactPhysics->jointNormal = st1->jointNormal2;
		}
		else if ( (((st1->jointNormal3.cross(st2->jointNormal1)).norm()<0.1) && (st1->jointNormal3.dot(st2->jointNormal1)<0)) || (((st1->jointNormal3.cross(st2->jointNormal2)).norm()<0.1) && (st1->jointNormal3.dot(st2->jointNormal2)<0)) || (((st1->jointNormal3.cross(st2->jointNormal3)).norm()<0.1) && (st1->jointNormal3.dot(st2->jointNormal3)<0)) )
		{
		  contactPhysics->isOnJoint = true;
		  contactPhysics->jointNormal = st1->jointNormal3;
		}
		else if ( (st1->joint>3 || st2->joint>3) && ( ( ((st1->jointNormal1.cross(st2->jointNormal1)).norm()>0.1) && ((st1->jointNormal1.cross(st2->jointNormal2)).norm()>0.1) && ((st1->jointNormal1.cross(st2->jointNormal3)).norm()>0.1) ) || ( ((st1->jointNormal2.cross(st2->jointNormal1)).norm()>0.1) && ((st1->jointNormal2.cross(st2->jointNormal2)).norm()>0.1) && ((st1->jointNormal2.cross(st2->jointNormal3)).norm()>0.1) ) || ( ((st1->jointNormal3.cross(st2->jointNormal1)).norm()>0.1) && ((st1->jointNormal3.cross(st2->jointNormal2)).norm()>0.1) && ((st1->jointNormal3.cross(st2->jointNormal3)).norm()>0.1) ) ) )  {  contactPhysics->isOnJoint = true; contactPhysics->more = true; contactPhysics->jointNormal = geom->normal; }
	}
	
	///to specify joint properties 
	if ( contactPhysics->isOnJoint ) {
			Real jf1 	= yade1->jointFrictionAngle;
			Real jf2 	= yade2->jointFrictionAngle;
			Real jkn1 	= yade1->jointNormalStiffness;
			Real jkn2 	= yade2->jointNormalStiffness;
			Real jks1 	= yade1->jointShearStiffness;
			Real jks2 	= yade2->jointShearStiffness;
			Real jdil1 	= yade1->jointDilationAngle;
			Real jdil2 	= yade2->jointDilationAngle;
			Real jcoh1 	= yade1->jointCohesion;
			Real jcoh2 	= yade2->jointCohesion;
			Real jSigT1	= yade1->jointTensileStrength;
			Real jSigT2	= yade2->jointTensileStrength;
			
			contactPhysics->tanFrictionAngle = std::tan(std::min(jf1,jf2));
			
			//contactPhysics->kn = jointNormalStiffness*2.*R1*R2/(R1+R2); // very first expression from Luc
			//contactPhysics->kn = (jkn1+jkn2)/2.0*2.*R1*R2/(R1+R2); // after putting jointNormalStiffness in material
			contactPhysics->kn = ( jkn1 + jkn2 ) /2.0 * contactPhysics->crossSection; // for a size independant expression
			contactPhysics->ks = ( jks1 + jks2 ) /2.0 * contactPhysics->crossSection; // for a size independant expression
			
			contactPhysics->tanDilationAngle = std::tan(std::min(jdil1,jdil2));
		  
			///to set if the contact is cohesive or not
			if ( ((cohesiveTresholdIteration < 0) || (scene->iter < cohesiveTresholdIteration)) && (std::min(jcoh1,jcoh2)>0 || std::min(jSigT1,jSigT2)>0) ) {
			  contactPhysics->isCohesive=true;
			  st1->noIniLinks++;
			  st2->noIniLinks++;
			} 
			else { contactPhysics->isCohesive=false; contactPhysics->FnMax=0; contactPhysics->FsMax=0; }
		  
			if ( contactPhysics->isCohesive ) {
				contactPhysics->FnMax = std::min(jSigT1,jSigT2)*contactPhysics->crossSection;
				contactPhysics->FsMax = std::min(jcoh1,jcoh2)*contactPhysics->crossSection;
			}
	}
	interaction->phys = contactPhysics;
}


void Ip2_JCFpmMat_JCFpmMat_JCFpmPhys::distributeStiffnesses(shared_ptr<JCFpmPhys> contactPhysics){
	std::random_device rd;
	std::mt19937 e2(rd());
	std::weibull_distribution<Real> weibullDistribution(stiffnessWeibullShapeParameter, 1.);
	Real cStiff = weibullDistribution(e2); 
	contactPhysics->kn = cStiff*contactPhysics->kn;
	contactPhysics->ks = cStiff*contactPhysics->ks;	
}

void Ip2_JCFpmMat_JCFpmMat_JCFpmPhys::distributeStrengthsNormal(shared_ptr<JCFpmPhys> contactPhysics,const shared_ptr<JCFpmMat>& yade1,const shared_ptr<JCFpmMat>& yade2 ){
	std::random_device rd;
	std::mt19937 e2(rd());
	std::normal_distribution<Real> tensileDistribution(min(yade1->tensileStrength,yade2->tensileStrength), yade1->tensileStrengthDeviation);
	std::normal_distribution<Real> cohDistribution(min(yade1->cohesion,yade2->cohesion), yade1->cohStrengthDeviation);
	Real SigT = tensileDistribution(e2); 
	Real Coh = cohDistribution(e2);
	contactPhysics->FnMax = SigT*contactPhysics->crossSection;
	contactPhysics->FsMax = Coh*contactPhysics->crossSection;
}

void Ip2_JCFpmMat_JCFpmMat_JCFpmPhys::distributeStrengthsWeibull(shared_ptr<JCFpmPhys> contactPhysics,const shared_ptr<JCFpmMat>& yade1,const shared_ptr<JCFpmMat>& yade2 ){
	std::random_device rd;
	std::mt19937 e2(rd());
	std::weibull_distribution<Real> weibullDistribution(strengthWeibullShapeParameter, 1.);
	Real cStrength = weibullDistribution(e2);
	contactPhysics->FnMax = cStrength*min(yade1->tensileStrength,yade2->tensileStrength)*contactPhysics->crossSection;
	contactPhysics->FsMax = cStrength*min(yade1->cohesion,yade2->cohesion)*contactPhysics->crossSection;
}


void Ip2_JCFpmMat_JCFpmMat_JCFpmPhys::distributeCrossSectionsWeibull(shared_ptr<JCFpmPhys> contactPhysics, Real R1, Real R2){
	std::random_device rd;
	std::mt19937 e2(rd());
	std::weibull_distribution<Real> weibullDistribution(xSectionWeibullShapeParameter, /*xSectionWeibullScaleParameter*/ 1.);
	//Real interactingRadius = scaleFactor*weibullDistribution(e2)/2.;
	Real correction = weibullDistribution(e2);
	if (correction < weibullCutOff) correction = weibullCutOff;
	//cout << "correction " << correction << " min R1 R2 "<< min(R1, R2) << endl;
	Real interactingRadius = correction*min(R1, R2);  // correcting radius to account for grain interactions
	contactPhysics->crossSection = Mathr::PI*pow(interactingRadius,2);
}
	


JCFpmPhys::~JCFpmPhys(){}
