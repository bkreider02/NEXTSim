#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4OpticalPhoton.hh"
#include "G4UnitsTable.hh"

#include "nDetSteppingAction.hh"
#include "nDetConstruction.hh"
#include "nDetRunAction.hh"


nDetSteppingAction::nDetSteppingAction(nDetRunAction* runAct) : runAction(runAct) {
	neutronTrack = false;
}

nDetSteppingAction::~nDetSteppingAction(){ }

void nDetSteppingAction::UserSteppingAction(const G4Step* aStep){
	G4Track *track = aStep->GetTrack();
	if(track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){ // Check for detected optical photons.
		if(aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary && aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName().find("psSiPM") != std::string::npos) {
			runAction->AddDetectedPhoton(aStep);
		}
	}
	else if(track->GetTrackStatus() != fAlive) return;
	else if(neutronTrack){ // Normal scattering event.
		if(track->GetTrackID() != 1)
			neutronTrack = false;
		else if(aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() == "expHall_physV"){ // Exiting the detector.
			runAction->finalizeNeutron(aStep);
			neutronTrack = false;
			return;
		}
		else
			runAction->scatterNeutron(aStep);
	}
	else if(track->GetTrackID() == 1){ // This is a primary particle scatter event.
		runAction->initializeNeutron(aStep);
		neutronTrack = true;
	}


	/*
	// check for secondary electrons and add their energy to the running total
	auto secondary = aStep->GetSecondaryInCurrentStep();
	size_t size_secondary = (*secondary).size();
	if (size_secondary){
		double energy;
		for (size_t i=0; i<(size_secondary);i++){
			auto secstep = (*secondary)[i];

			G4String secondaryName = secstep->GetDefinition()->GetParticleName();

			if (secondaryName == "e-") {
				energy = secstep->GetKineticEnergy();

				runAction->AddDetectedSecondary(aStep,energy);
			}
		}
	}
	*/


}

void nDetSteppingAction::Reset(){
	neutronTrack = false;
}
