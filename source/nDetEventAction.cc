#include "nDetEventAction.hh"
#include "nDetRunAction.hh"

#include "G4Event.hh"
 
nDetEventAction::nDetEventAction(nDetRunAction* run) : runAct(run) { }

nDetEventAction::~nDetEventAction(){ 
}

void nDetEventAction::BeginOfEventAction(const G4Event* evt){
	runAct->setEventNumber(evt->GetEventID());

	// Get initial starting positions and initial particle kinetic energy
	runAct->getDebugData()->nStartPosX = -(evt->GetPrimaryVertex(0)->GetZ0());
	runAct->getDebugData()->nStartPosY = (evt->GetPrimaryVertex(0)->GetY0());
	runAct->getDebugData()->nStartPosZ = evt->GetPrimaryVertex(0)->GetX0();
	runAct->getDebugData()->initialE = evt->GetPrimaryVertex(0)->GetPrimary()->GetKineticEnergy();
}

void nDetEventAction::EndOfEventAction(const G4Event*){
	// Process primary scatters. 
	runAct->process();
}
