#include "time.h"

#include "Randomize.hh"
#include "G4Timer.hh"
#include "G4Run.hh"

#include "nDetDataPack.hh"
#include "nDetThreadContainer.hh"
#include "nDetRunAction.hh"
#include "nDetEventAction.hh"
#include "nDetStackingAction.hh"
#include "nDetTrackingAction.hh"
#include "nDetSteppingAction.hh"
#include "nDetParticleSource.hh"
#include "nDetMasterOutputFile.hh"
#include "nDetDetector.hh"
#include "termColors.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include <random>

const double KINETIC_ENERGY_THRESHOLD = 0.03; // MeV

std::default_random_engine generator;
std::normal_distribution<double> distribution(6.551,0.0);

/// Returns the dot product of two vectors i.e. v1 * v2
double dotProduct(const G4ThreeVector &v1, const G4ThreeVector &v2){
	return (v1.getX()*v2.getX() + v1.getY()*v2.getY() + v1.getZ()*v2.getZ());
}

double calcRecMass(G4double E1, double p1, G4double E2, double p2, G4double Edel, double pdel){
	double v1 = pow(2*E1/p1,2);
	double v2 = pow(2*E2/p2,2);
	double vdel = pow(2*Edel/pdel,2);
	
	G4double mx = (v1-v2)/vdel;
	return mx;
}

primaryTrackInfo::primaryTrackInfo(const G4Step *step){
	this->setValues(step->GetTrack());
	if(step->GetPreStepPoint()->GetPhysicalVolume()->GetName().find("Scint") != std::string::npos) // Scatter event occured inside a scintillator.
		inScint = true;
	//dkE = -1*(step->GetPostStepPoint()->GetKineticEnergy()-step->GetPreStepPoint()->GetKineticEnergy());
	dkE = -1*step->GetDeltaEnergy();
	atomicMass = calcRecMass(step->GetPreStepPoint()->GetKineticEnergy(),step->GetPreStepPoint()->GetMomentum().mag(),
							 step->GetPostStepPoint()->GetKineticEnergy(),step->GetPostStepPoint()->GetMomentum().mag(),
							 dkE, step->GetDeltaMomentum().mag());

}

primaryTrackInfo::primaryTrackInfo(const G4Track *track){
	this->setValues(track);
}

bool primaryTrackInfo::compare(const G4Track *rhs){
	return (this->pos == rhs->GetPosition());
}

double primaryTrackInfo::getAngle(const G4Track *rhs){
	if(compare(rhs))
		angle = getAngle(rhs->GetMomentumDirection());
	else
		angle = -1;
	return angle;
}

double primaryTrackInfo::getAngle(const G4ThreeVector &dir2){
	angle = std::acos(dotProduct(dir, dir2)/(dir.mag()*dir2.mag())) * 180/3.14159;
	return angle;
}

double primaryTrackInfo::getPathLength(const G4ThreeVector &rhs){
	return (plength = (pos-rhs).mag());
}

void primaryTrackInfo::print(){
	std::cout << "name=" << part->GetParticleName() << ", kE=" << kE << ", dkE=" << dkE << ", pos=(" << pos.getX() << ", " << pos.getY() << ", " << pos.getZ() << "), dir=(" << dir.getX() << ", " << dir.getY() << ", " << dir.getZ() << ")\n";
}

void primaryTrackInfo::setValues(const G4Track *track){
	pos = track->GetPosition();
	dir = track->GetMomentumDirection();
	kE = track->GetKineticEnergy();
	dkE = -1;
	gtime = track->GetGlobalTime();
	plength = track->GetTrackLength();
	part = track->GetParticleDefinition();
	
	copyNum = track->GetTouchable()->GetCopyNumber();
	trackID = track->GetTrackID();
	//atomicMass = part->GetBaryonNumber();
	inScint = false;
}

nDetRunAction::nDetRunAction(){
	timer = new G4Timer;
	
	outputTraces = false;
	outputDebug = false;
	outputMultiDebug = false;
	verbose = false;
	
	eventAction = NULL;
	stacking = NULL;
	tracking = NULL;
	stepping = NULL;
	
	// Set the particle source.
	source = &nDetParticleSource::getInstance();

	numPhotonsTotal = 0;
	numPhotonsDetTotal = 0;
	
	// Pointer to the start detector (if available)
	startDetector = NULL;

	// Get a pointer to the detector.
	detector = &nDetConstruction::getInstance(); // The detector builder is a singleton class.

	// Set data structure addresses.
	data.setDataAddresses(&evtData, &outData, &outImplantData, &multData, &debugData, &traceData);
}

nDetRunAction::~nDetRunAction(){
	delete timer;
}

void nDetRunAction::BeginOfRunAction(const G4Run* aRun)
{
	numPhotonsTotal = 0;
	numPhotonsDetTotal = 0;

	// Get RunId and threadID
	evtData.runNb = aRun->GetRunID();
	evtData.threadID = G4Threading::G4GetThreadId();

	if(G4Threading::G4GetThreadId() >= 0) return; // Master thread only.

	// Update the source. Only need to do this once since it's a singleton
	source->UpdateAll();
	
	// initialize all pmt responses with correct output information
	initializePmtResponses();

	G4cout << "nDetRunAction::BeginOfRunAction()->"<< G4endl;
	G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl; 
	
	timer->Start();

	// Get a pointer to the output file
	nDetMasterOutputFile *outputFile = &nDetMasterOutputFile::getInstance();

	nDetThreadContainer *container = &nDetThreadContainer::getInstance();
	if(userDetectors.size() > 1){
		if(outputFile->getOutputDebug()){
			Display::WarningPrint("Main debug output is not supported for more than one detector! Using reduced debug readout.", "nDetRunAction");
			outputFile->setOutputDebug((outputDebug = false));
			outputMultiDebug = true;
		}
		if(outputFile->getOutputTraces()){
			Display::WarningPrint("Trace output is not supported for more than one detector!", "nDetRunAction");
			outputFile->setOutputTraces((outputTraces = false));	
		}
		for(size_t index = 0; index < container->size(); index++){ // Set options per thread
			container->getActionManager(index)->getRunAction()->setOutputDebug(false);
			container->getActionManager(index)->getRunAction()->setOutputTraces(false);
			container->getActionManager(index)->getRunAction()->setOutputMultiDebug(outputMultiDebug);
		}
		outputFile->setMultiDetectorMode(true);
	}
	else if(userImplants.size() > 1){
		if(outputFile->getOutputDebug()){
			Display::WarningPrint("Main debug output is not supported for more than one implant! Using reduced debug readout.", "nDetRunAction");
			outputFile->setOutputDebug((outputDebug = false));
			outputMultiDebug = true;
		}
		if(outputFile->getOutputTraces()){
			Display::WarningPrint("Trace output is not supported for more than one implant!", "nDetRunAction");
			outputFile->setOutputTraces((outputTraces = false));	
		}
		for(size_t index = 0; index < container->size(); index++){ // Set options per thread
			container->getActionManager(index)->getRunAction()->setOutputDebug(false);
			container->getActionManager(index)->getRunAction()->setOutputTraces(false);
			container->getActionManager(index)->getRunAction()->setOutputMultiDebug(outputMultiDebug);
		}
		outputFile->setMultiDetectorMode(true);
	}
	else{ // Single detector mode
		for(size_t index = 0; index < container->size(); index++){ // Set options per thread
			container->getActionManager(index)->getRunAction()->setOutputDebug(outputFile->getOutputDebug());
			container->getActionManager(index)->getRunAction()->setOutputTraces(outputFile->getOutputTraces());
			container->getActionManager(index)->getRunAction()->initializePmtResponses();
		}
		outputFile->setMultiDetectorMode(false);
	}

	// Open a root file.
	outputFile->openRootFile(aRun);

	// Set the total number of events
	outputFile->setTotalEvents(aRun->GetNumberOfEventToBeProcessed());
	numberOfEvents = aRun->GetNumberOfEventToBeProcessed();
}

void nDetRunAction::EndOfRunAction(const G4Run* aRun)
{
	if(!IsMaster()) return; // Master thread only.
	
	timer->Stop();
	G4cout << "number of event = " << aRun->GetNumberOfEvent() << " " << *timer << G4endl;
}

void nDetRunAction::initializePmtResponses() {

	// FOR DETECTORS!!!!! VERIFY THIS WORKS!!!
	// (SOMETHING ABOUT TWO ARRAYS IN THE ROOT TREE (?))
	for (auto& det : userDetectors) {
		centerOfMass *cmL = det.getCenterOfMassL();
		centerOfMass *cmR = det.getCenterOfMassR();
		pmtResponse *pmtL = cmL->getPmtResponse();
		pmtResponse *pmtR = cmR->getPmtResponse();
		copyResponseParameters(pmtL,detector->GetCenterOfMassL()->getPmtResponse());
		copyResponseParameters(pmtR,detector->GetCenterOfMassR()->getPmtResponse());

		std::vector<pmtResponse> *anodeResponseL = cmL->getAnodeResponse();
		std::vector<pmtResponse> *anodeResponseR = cmR->getAnodeResponse();

		std::vector<pmtResponse> *tempAL = detector->GetCenterOfMassL()->getAnodeResponse();
		std::vector<pmtResponse> *tempAR = detector->GetCenterOfMassR()->getAnodeResponse();

		for(size_t i = 0; i < 4; i++){
			copyResponseParameters(&anodeResponseL->at(i),&tempAL->at(i));
			copyResponseParameters(&anodeResponseR->at(i),&tempAR->at(i));
		}

		std::vector<pmtResponse> *pixelResponseL = cmL->getPixelResponse();
		std::vector<pmtResponse> *pixelResponseR = cmR->getPixelResponse();

		std::vector<pmtResponse> *tempPL = detector->GetCenterOfMassL()->getPixelResponse();
		std::vector<pmtResponse> *tempPR = detector->GetCenterOfMassR()->getPixelResponse();

		if (detector->GetDetectorParameters().GetNumPmtRows() == 8 && detector->GetDetectorParameters().GetNumPmtColumns() == 8) { // only record pixel response for 8x8 pmts
			for (size_t i = 0; i < 8; i++) {
				for (size_t j = 0; j < 8; j++) {
					copyResponseParameters(&pixelResponseL->at(8*i+j),&tempPL->at(8*i+j));
					copyResponseParameters(&pixelResponseR->at(8*i+j),&tempPR->at(8*i+j));
				}
			}
		}
	}


	for (auto& imp : userImplants) {
		centerOfMass *cmI = imp.getCenterOfMass();
		pmtResponse *pmtI = cmI->getPmtResponse();
		copyResponseParameters(pmtI,detector->GetCenterOfMass()->getPmtResponse());

		std::vector<pmtResponse> *anodeResponseI = cmI->getAnodeResponse();
		std::vector<pmtResponse> *tempA = detector->GetCenterOfMass()->getAnodeResponse();
		for(size_t i = 0; i < 4; i++){
			copyResponseParameters(&anodeResponseI->at(i),&tempA->at(i));
		}

		std::vector<pmtResponse> *pixelResponseI = cmI->getPixelResponse();
		std::vector<pmtResponse> *tempP = detector->GetCenterOfMass()->getPixelResponse();
		if (detector->GetDetectorParameters().GetNumPmtRows() == 8 && detector->GetDetectorParameters().GetNumPmtColumns() == 8) { // only record pixel response for 8x8 pmts
			for (size_t i = 0; i < 8; i++) {
				for (size_t j = 0; j < 8; j++) {
					copyResponseParameters(&pixelResponseI->at(8*i+j),&tempP->at(8*i+j));
				}
			}
		}
	}
}

void nDetRunAction::updateDetector(nDetConstruction *construction){
	// Clear all currently defined detectors
	userDetectors.clear();
	userImplants.clear();

	// Copy the list of detectors
	construction->GetCopiesOfDetectors(userDetectors);
	construction->GetCopiesOfImplants(userImplants);	
	
	// Search for a start detector. Currently only one start is supported, break after finding the first one
	startDetector = NULL;
	for(std::vector<nDetDetector>::iterator iter = userDetectors.begin(); iter != userDetectors.end(); iter++){
		if(iter->GetIsStart()){
			startDetector = &(*iter);
			 break;
		}
		else
			endDetector=true;
	}

	startImplant = NULL;
	for(std::vector<nDetImplant>::iterator iter = userImplants.begin(); iter != userImplants.end(); iter++){
		if(iter->GetIsStart()){
			startImplant = &(*iter);
			 break;
		}
		else
			endImplant=true;
	}
}

G4int nDetRunAction::checkCopyNumber(const G4int &num) const {
	for(std::vector<nDetDetector>::const_iterator iter = userDetectors.begin(); iter != userDetectors.end(); iter++){
		if(iter->checkCopyNumber(num))
			return iter->getParentCopyNumber();
	}

	for(std::vector<nDetImplant>::const_iterator iter = userImplants.begin(); iter != userImplants.end(); iter++){
		if(iter->checkCopyNumber(num))
			return iter->getParentCopyNumber();
	}
	return -1;
}

bool nDetRunAction::getSegmentFromCopyNum(const G4int &copyNum, G4int &col, G4int &row) const {
	for(std::vector<nDetDetector>::const_iterator iter = userDetectors.begin(); iter != userDetectors.end(); iter++){
		if(iter->getSegmentFromCopyNum(copyNum, col, row))
			return true;
	}

	for(std::vector<nDetImplant>::const_iterator iter = userImplants.begin(); iter != userImplants.end(); iter++){
		if(iter->getSegmentFromCopyNum(copyNum, col, row))
			return true;
	}
	return false;
}

bool nDetRunAction::processDetector(nDetDetector* det){
	if(!det || det->empty()) // No detected photons, do not process
		return false;

	// Get the time offset due to straggling in the target (if applicable)
	double targetTimeOffset = source->GetTargetTimeOffset();
	bool Ftrigger = false;
	// Get pointers to the CoM calculators
	centerOfMass *cmL = det->getCenterOfMassL();
	centerOfMass *cmR = det->getCenterOfMassR();

	debugData.nPhotons[0] += cmL->getNumDetected();
	debugData.nPhotons[1] += cmR->getNumDetected();
	
	// Compute the total number of detected photons
	outData.nPhotonsDet += cmL->getNumDetected() + cmR->getNumDetected();
		
	// Check for valid bar detection
	if(cmL->getNumDetected() > 0 && cmR->getNumDetected() > 0)
		evtData.goodEvent = true;
		
	// Compute the photon detection efficiency
	outData.nPhotonsTot = stacking->GetNumPhotonsProduced();
	if(outData.nPhotonsTot > 0)
		outData.photonDetEff = outData.nPhotonsDet/(double)outData.nPhotonsTot;
	else
		outData.photonDetEff = -1;		

	// Get the photon center-of-mass positions
	G4ThreeVector centerL = cmL->getCenter();
	G4ThreeVector centerR = cmR->getCenter();
	debugData.photonDetComX[0] = centerL.getX(); debugData.photonDetComX[1] = centerR.getX(); 
	debugData.photonDetComY[0] = centerL.getY(); debugData.photonDetComY[1] = centerR.getY(); 
	//debugData.photonDetComZ[0] = centerL.getZ(); debugData.photonDetComZ[1] = centerR.getZ(); 

	// Get photon arrival times at the PMTs
	debugData.photonMinTime[0] = cmL->getMinArrivalTime();
	debugData.photonAvgTime[0] = cmL->getAvgArrivalTime();

	debugData.photonMinTime[1] = cmR->getMinArrivalTime();
	debugData.photonAvgTime[1] = cmR->getAvgArrivalTime();
	
	pmtResponse *pmtL = cmL->getPmtResponse();
	pmtResponse *pmtR = cmR->getPmtResponse();
	//std::cout<<"!!! Process !!! "<<detector->GetCenterOfMass()->getNumColumns()<<" "<<cmL->getNumColumns()<<std::endl;

	// The ADC clock is running continuously so the PMT signals may arrive at any time relative
	// to the clock. To remove any inherent bias from this, we latch the ADCs at a random time 
	// in the range [-adcTick/2, +adcTick/2]. The pmtResponse class will automatically remove 
	// the latching time after running the CFD.
	double latch = G4UniformRand()-0.5; // [-0.5, 0.5]
	pmtL->setAdcLatchTicks(latch);
	pmtR->setAdcLatchTicks(latch);
	
	// "Digitize" the light pulses.
	pmtL->digitize();
	pmtR->digitize();
	
	// Check for saturated pulse.
	if(pmtL->getPulseIsSaturated() || pmtR->getPulseIsSaturated()){
		Display::WarningPrint("One or more PMT's traces have saturated! Recommend lowering the gain.", "nDetRunAction");
	}
	
	// Copy the trace into the trace vector.
	if(outputTraces){
		pmtL->copyTrace(traceData.left);
		pmtR->copyTrace(traceData.right);
		traceData.mult++;
	}		
	
	// Do some light pulse analysis
	debugData.pulsePhase[0] = pmtL->analyzePolyCFD() + targetTimeOffset;
	debugData.pulseQDC[0] = pmtL->integratePulseFromMaximum();
	debugData.pulseMax[0] = pmtL->getMaximum();
	debugData.pulseMaxTime[0] = pmtL->getMaximumTime();
	debugData.pulseArrival[0] = pmtL->getWeightedPhotonArrivalTime();

	debugData.pulsePhase[1] = pmtR->analyzePolyCFD() + targetTimeOffset;
	debugData.pulseQDC[1] = pmtR->integratePulseFromMaximum();
	debugData.pulseMax[1] = pmtR->getMaximum();
	debugData.pulseMaxTime[1] = pmtR->getMaximumTime();
	debugData.pulseArrival[1] = pmtR->getWeightedPhotonArrivalTime();		

	// Set Trigger boolean for events that would register in DAQ
	if(pmtR->getTrigger() && pmtL->getTrigger() && abs(pmtL->getMaximumTime()-pmtR->getMaximumTime()) < 5) Ftrigger = true;

	// Print the digitized traces.
	if(pmtL->getPrintTrace() || pmtR->getPrintTrace()){
		size_t traceLength = pmtL->getPulseLength();
		unsigned short *traceL = pmtL->getDigitizedPulse();
		unsigned short *traceR = pmtR->getDigitizedPulse();
		std::cout << "***********************************************************\n";
		std::cout << "* PhotonsTot     : " << outData.nPhotonsTot << std::endl;
		std::cout << "* PhotonsDet     : " << outData.nPhotonsDet << std::endl;
		std::cout << "* MaxIndex       : " << pmtL->getMaximumIndex() << "\t" << pmtR->getMaximumIndex() << std::endl;
		std::cout << "* Baseline       : " << pmtL->getBaseline() << "\t" << pmtR->getBaseline() << std::endl;	
		std::cout << "* Maximum        : " << pmtL->getMaximum() << "\t" << pmtR->getMaximum() << std::endl;
		std::cout << "* MaxTime        : " << pmtL->getMaximumTime() << "\t" << pmtR->getMaximumTime() << std::endl;
		std::cout << "* WeightedArrival: " << pmtL->getWeightedPhotonArrivalTime() << "\t" << pmtR->getWeightedPhotonArrivalTime() << std::endl;
		std::cout << "* MinimumArrival : " << pmtL->getMinimumPhotonArrivalTime() << "\t" << pmtR->getMinimumPhotonArrivalTime() << std::endl;
		std::cout << "***********************************************************\n";

		int adcClockTick = pmtL->getAdcClockInNanoseconds();
		for(size_t i = 0; i < traceLength; i++){
			std::cout << i*adcClockTick << "\t" << traceL[i] << "\t" << traceR[i] << std::endl;
		}
	}
	
	// Get the digitizer response of the anodes.
	std::vector<pmtResponse> *anodeResponseL = cmL->getAnodeResponse();
	std::vector<pmtResponse> *anodeResponseR = cmR->getAnodeResponse();

	// Digitize anode waveforms and integrate.
	double anodeQDC[2][4];
	for(size_t i = 0; i < 4; i++){
		anodeResponseL->at(i).digitize();
		anodeResponseR->at(i).digitize();
		debugData.anodeQDC[0][i] = anodeResponseL->at(i).integratePulseFromMaximum();
		debugData.anodeQDC[1][i] = anodeResponseR->at(i).integratePulseFromMaximum();
		anodeQDC[0][i] = debugData.anodeQDC[0][i];
		anodeQDC[1][i] = debugData.anodeQDC[1][i];
	}	
	
	// Compute the anode positions.
	for(size_t i = 0; i < 2; i++){
		debugData.reconDetComX[i] = -((anodeQDC[i][0]+anodeQDC[i][1])-(anodeQDC[i][2]+anodeQDC[i][3]))/(anodeQDC[i][0]+anodeQDC[i][1]+anodeQDC[i][2]+anodeQDC[i][3]);
		debugData.reconDetComY[i] = ((anodeQDC[i][1]+anodeQDC[i][2])-(anodeQDC[i][3]+anodeQDC[i][0]))/(anodeQDC[i][0]+anodeQDC[i][1]+anodeQDC[i][2]+anodeQDC[i][3]);
	}
	outData.reconComX = (debugData.reconDetComX[0] + debugData.reconDetComX[1]) / 2;
	outData.reconComY = (debugData.reconDetComY[0] + debugData.reconDetComY[1]) / 2;
	
	if(outputDebug || outputMultiDebug){
		// Perform CFD on digitized anode waveforms.
		/*for(size_t i = 0; i < 4; i++){
			debugData.anodePhase[0][i] = anodeResponseL[i].analyzePolyCFD() + targetTimeOffset; // left
			debugData.anodePhase[1][i] = anodeResponseR[i].analyzePolyCFD() + targetTimeOffset; // right
		}*/

		G4ThreeVector nCenterMass(debugData.nComX, debugData.nComY, debugData.nComZ);
		G4ThreeVector nIncidentPos(debugData.nEnterPosX, debugData.nEnterPosY, debugData.nEnterPosZ);
		G4ThreeVector nExitPos(debugData.nExitPosX, debugData.nExitPosY, debugData.nExitPosZ);

		// Compute the neutron scatter center-of-mass.
		nCenterMass = (1/debugData.neutronWeight)*nCenterMass;
	
		// Convert the neutron incident/exit positions to the frame of the detector.
		nIncidentPos = nIncidentPos;
		nExitPos = nExitPos;

		// Now in the rotated frame of the detector.
		debugData.nComX = nCenterMass.getX();
		debugData.nComY = nCenterMass.getY();
		debugData.nComZ = nCenterMass.getZ();
		debugData.nEnterPosX = nIncidentPos.getX();
		debugData.nEnterPosY = nIncidentPos.getY();
		debugData.nEnterPosZ = nIncidentPos.getZ();
		debugData.nExitPosX = nExitPos.getX();
		debugData.nExitPosY = nExitPos.getY();
		debugData.nExitPosZ = nExitPos.getZ();
	}

	// Compute the light balance (Z).
	outData.lightBalance = (debugData.pulseQDC[0]-debugData.pulseQDC[1])/(debugData.pulseQDC[0]+debugData.pulseQDC[1]);
	outData.tdiff = (debugData.pulsePhase[0]-debugData.pulsePhase[1]);
	outData.photonTdiff = (debugData.photonAvgTime[0] - debugData.photonAvgTime[1]);
	
	// Compute "bar" variables.
	double offset = distribution(generator);
	outData.barTOF = (debugData.pulsePhase[0]+debugData.pulsePhase[1])/2-offset;
	outData.barQDC = std::sqrt(debugData.pulseQDC[0]*debugData.pulseQDC[1]);
	outData.barMaxADC = std::sqrt(abs(debugData.pulseMax[0])*abs(debugData.pulseMax[1]));
	outData.barTrig = Ftrigger;
	outData.photonTOF = (debugData.photonAvgTime[0]+debugData.photonAvgTime[1])/2.0-offset;
	outData.photonComX = (debugData.photonDetComX[0] + debugData.photonDetComX[1]) / 2;
	outData.photonComY = (debugData.photonDetComY[0] + debugData.photonDetComY[1]) / 2;

	// Get the segment of the detector where the photon CoM occurs.
	cmL->getCenterSegment(debugData.centerOfMassColumn[0], debugData.centerOfMassRow[0]);
	cmR->getCenterSegment(debugData.centerOfMassColumn[1], debugData.centerOfMassRow[1]);		

	// Update photon statistics.
	numPhotonsTotal += outData.nPhotonsTot;
	numPhotonsDetTotal += outData.nPhotonsDet;

	return true;
}

bool nDetRunAction::processImplant(nDetImplant* imp){
	if(!imp || imp->empty()) // No detected photons, do not process
		return false;

	// Get the time offset due to straggling in the target (if applicable)
	double targetTimeOffset = source->GetTargetTimeOffset();
	bool Ftrigger = false;
	// Get pointers to the CoM calculators
	params = detector->GetDetectorParameters();
	centerOfMass *cmI = imp->getCenterOfMass();

	//std::cout<<"!!! Process !!! "<<detector->GetCenterOfMass()->getNumColumns()<<" "<<cmI->getNumColumns()<<std::endl;
	cmI->setSegmentedPmt(&params);
	if (params.GetNumPmtRows() == 16 && params.GetNumPmtColumns() == 16) {
		cmI->loadGainMatrix("hamamatsuH13700.dat");
	}
	else {
		cmI->loadGainMatrix("hamamatsuH12700A_LA0967.dat");
	}

	debugData.nPhotons[0] += cmI->getNumDetected();
	
	// Compute the total number of detected photons
	outImplantData.nPhotonsDet += cmI->getNumDetected();
		
	// Check for valid bar detection
	if(cmI->getNumDetected() > 0)
		evtData.goodEvent = true;
		
	// Compute the photon detection efficiency
	outImplantData.nPhotonsTot = stacking->GetNumPhotonsProduced();
	if(outImplantData.nPhotonsTot > 0)
		outImplantData.photonDetEff = outImplantData.nPhotonsDet/(double)outImplantData.nPhotonsTot;
	else
		outImplantData.photonDetEff = -1;			

	// Get the photon center-of-mass positions
	G4ThreeVector centerI = cmI->getCenter();
	debugData.photonDetComX[0] = centerI.getX(); 
	debugData.photonDetComY[0] = centerI.getY();
	//debugData.photonDetComZ[0] = centerI.getZ();

	// Get photon arrival times at the PMTs
	debugData.photonMinTime[0] = cmI->getMinArrivalTime();
	debugData.photonAvgTime[0] = cmI->getAvgArrivalTime();
	
	pmtResponse *pmtI = cmI->getPmtResponse();

	// The ADC clock is running continuously so the PMT signals may arrive at any time relative
	// to the clock. To remove any inherent bias from this, we latch the ADCs at a random time 
	// in the range [-adcTick/2, +adcTick/2]. The pmtResponse class will automatically remove 
	// the latching time after running the CFD.
	double latch = G4UniformRand()-0.5; // [-0.5, 0.5]
	pmtI->setAdcLatchTicks(latch);
	
	// "Digitize" the light pulses.
	pmtI->digitize();
	
	// Check for saturated pulse.
	if(pmtI->getPulseIsSaturated()){
		Display::WarningPrint("Implant PMT's traces have saturated! Recommend lowering the gain.", "nDetRunAction");
	}
	
	// Copy the trace into the trace vector.
	if(outputTraces){
		pmtI->copyTrace(traceData.left); //This might be wrong 
		traceData.mult++;
	}		
	
	// Do some light pulse analysis
	debugData.pulsePhase[0] = pmtI->analyzePolyCFD() + targetTimeOffset;
	debugData.pulseQDC[0] = pmtI->integratePulseFromMaximum();


	debugData.pulseMax[0] = pmtI->getMaximum();
	debugData.pulseMaxTime[0] = pmtI->getMaximumTime();
	debugData.pulseArrival[0] = pmtI->getWeightedPhotonArrivalTime();	

	// Set Trigger boolean for events that would register in DAQ
	if(pmtI->getTrigger() && abs(pmtI->getMaximumTime()) < 5) Ftrigger = true;

	// Print the digitized traces.
	if(pmtI->getPrintTrace()){
		size_t traceLength = pmtI->getPulseLength();
		unsigned short *traceI = pmtI->getDigitizedPulse();
		std::cout << "***********************************************************\n";
		std::cout << "* PhotonsTot     : " << outImplantData.nPhotonsTot << std::endl;
		std::cout << "* PhotonsDet     : " << outImplantData.nPhotonsDet << std::endl;
		std::cout << "* MaxIndex       : " << pmtI->getMaximumIndex() << std::endl;
		std::cout << "* Baseline       : " << pmtI->getBaseline() << std::endl;	
		std::cout << "* Maximum        : " << pmtI->getMaximum() <<  std::endl;
		std::cout << "* MaxTime        : " << pmtI->getMaximumTime() << std::endl;
		std::cout << "* WeightedArrival: " << pmtI->getWeightedPhotonArrivalTime() << std::endl;
		std::cout << "* MinimumArrival : " << pmtI->getMinimumPhotonArrivalTime() << std::endl;
		std::cout << "***********************************************************\n";

		int adcClockTick = pmtI->getAdcClockInNanoseconds();
		for(size_t i = 0; i < traceLength; i++){
			std::cout << i*adcClockTick << "\t" << traceI[i] << std::endl;
		}
	}
	
	// Get the digitizer response of the anodes.
	std::vector<pmtResponse> *anodeResponseI = cmI->getAnodeResponse();

	// Digitize anode waveforms and integrate.
	double anodeQDC[4]={0};
	for(size_t i = 0; i < 4; i++){
		anodeResponseI->at(i).digitize();
		debugData.anodeQDC[0][i] = anodeResponseI->at(i).integratePulseFromMaximum();
		anodeQDC[i] = debugData.anodeQDC[0][i];
	}

	std::vector<pmtResponse> *pixelResponseI = cmI->getPixelResponse();

	for (size_t i = 0; i < 8; i++)
		for (size_t j = 0; j < 8; j++) {
			debugData.photonsPerPixel[i][j] = 0;
			debugData.pixelQDC[i][j] = 0;
		}

	/*

	double pixelQDC[8][8];
	double photonsPerPixel[8][8];

	// initialize to zero
	for (size_t i = 0; i < 8; i++)
		for (size_t j = 0; j < 8; j++) {
			pixelQDC[i][j] = 0;
			photonsPerPixel[i][j] = 0;
		}

	if (detector->GetDetectorParameters().GetNumPmtRows() == 8 && detector->GetDetectorParameters().GetNumPmtColumns() == 8) { // only record pixel response for 8x8 pmts
		for (size_t i = 0; i < 8; i++) {
			for (size_t j = 0; j < 8; j++) {
				if (pixelResponseI->at(8*i+j).getPhotonCount() != 0) { // don't bother with pixels that aren't hit
					
					// get number of photons
					debugData.photonsPerPixel[i][j] = pixelResponseI->at(8*i+j).getPhotonCount();

					pixelResponseI->at(8*i+j).digitize();
					debugData.pixelQDC[i][j] = pixelResponseI->at(8*i+j).integratePulseFromMaximum();

					pixelQDC[i][j] = debugData.pixelQDC[i][j];
					photonsPerPixel[i][j] = debugData.photonsPerPixel[i][j];
					
				}
			}
		}
	}
	*/
	
	// Compute the anode positions.
	//debugData.reconDetComX[0] = -((anodeQDC[0]+anodeQDC[1])-(anodeQDC[2]+anodeQDC[3]))/(anodeQDC[0]+anodeQDC[1]+anodeQDC[2]+anodeQDC[3]);
	//debugData.reconDetComY[0] = ((anodeQDC[1]+anodeQDC[2])-(anodeQDC[3]+anodeQDC[0]))/(anodeQDC[0]+anodeQDC[1]+anodeQDC[2]+anodeQDC[3]);
	if(params.GetReconLogic()){
		debugData.reconDetComX[0] = cmI->getSipmX();
		debugData.reconDetComY[0] = cmI->getSipmY();
	}
	else{
		debugData.reconDetComX[0] = cmI->getReconstructedX();
		debugData.reconDetComY[0] = cmI->getReconstructedY();
	}
	//std::cout<<"recon "<<debugData.reconDetComX[0]<<std::endl;
	outImplantData.reconComX = debugData.reconDetComX[0];
	outImplantData.reconComY = debugData.reconDetComY[0];
	
	if(outputDebug || outputMultiDebug){
		// Perform CFD on digitized anode waveforms.
		/*for(size_t i = 0; i < 4; i++){
			debugData.anodePhase[0][i] = anodeResponseL[i].analyzePolyCFD() + targetTimeOffset; // left
			debugData.anodePhase[1][i] = anodeResponseR[i].analyzePolyCFD() + targetTimeOffset; // right
		}*/

		G4ThreeVector nCenterMass(debugData.nComX, debugData.nComY, debugData.nComZ);
		G4ThreeVector nIncidentPos(debugData.nEnterPosX, debugData.nEnterPosY, debugData.nEnterPosZ);
		G4ThreeVector nExitPos(debugData.nExitPosX, debugData.nExitPosY, debugData.nExitPosZ);

		// Compute the neutron scatter center-of-mass.
		nCenterMass = (1/debugData.neutronWeight)*nCenterMass;
	
		// Convert the neutron incident/exit positions to the frame of the detector.
		nIncidentPos = nIncidentPos;
		nExitPos = nExitPos;

		// Now in the rotated frame of the detector.
		debugData.nComX = nCenterMass.getX();
		debugData.nComY = nCenterMass.getY();
		debugData.nComZ = nCenterMass.getZ();
		debugData.nEnterPosX = nIncidentPos.getX();
		debugData.nEnterPosY = nIncidentPos.getY();
		debugData.nEnterPosZ = nIncidentPos.getZ();
		debugData.nExitPosX = nExitPos.getX();
		debugData.nExitPosY = nExitPos.getY();
		debugData.nExitPosZ = nExitPos.getZ();
	}

	// Compute the light balance (Z).
	outImplantData.lightBalance = (debugData.pulseQDC[0]-debugData.pulseQDC[1])/(debugData.pulseQDC[0]+debugData.pulseQDC[1]);
	outImplantData.tdiff = (debugData.pulsePhase[0]-debugData.pulsePhase[1]);
	outImplantData.photonTdiff = (debugData.photonAvgTime[0]);
	
	// Compute "bar" variables.
	double offset = distribution(generator);
	outImplantData.barTOF = (debugData.pulsePhase[0])-offset;
	outImplantData.barQDC = (debugData.pulseQDC[0]);
	outImplantData.barMaxADC = (abs(debugData.pulseMax[0]));
	outImplantData.barTrig = Ftrigger;
	outImplantData.photonTOF = (debugData.photonAvgTime[0])-offset;
	outImplantData.photonComX = (debugData.photonDetComX[0]);
	outImplantData.photonComY = (debugData.photonDetComY[0]);
	// Get the segment of the detector where the photon CoM occurs.
	cmI->getCenterSegment(debugData.centerOfMassColumn[0], debugData.centerOfMassRow[0]);	

	// THIS MAY NEED TO BE CHANGED FOR PMTs WITH GAP
	// calculate the location of contact within the pmt segment
	double pmtWidth = params.GetPmtWidth();
	double pmtHeight = params.GetPmtHeight();
	double pixelWidth = cmI->getPmtPixelWidth();
	double pixelHeight = cmI->getPmtPixelHeight();
	int col = debugData.centerOfMassColumn[0];
	int row = debugData.centerOfMassRow[0];
	debugData.pmtX[0] = 0.5*pmtWidth+debugData.photonDetComX[0]-(col+0.5)*pixelWidth;
	debugData.pmtY[0] = 0.5*pmtHeight+debugData.photonDetComY[0]-(row+0.5)*pixelHeight;

	// Update photon statistics.
	numPhotonsTotal += outImplantData.nPhotonsTot;
	numPhotonsDetTotal += outImplantData.nPhotonsDet;
	
	debugData.secondaryEnergy = cmI->getSecondaryEnergy();
	
	return true;
}


void nDetRunAction::copyResponseParameters(pmtResponse *dest, pmtResponse *src) {
	dest->setRisetime(src->getRisetime());
	dest->setFalltime(src->getFalltime());
	dest->setTransitTimeSpread(src->getTimeSpread());
	dest->setGain(src->getGain());
	dest->setBaselineFraction(src->getBaselineFraction());
	dest->setBaselineJitterFraction(src->getBaselineJitterFraction());
	dest->setTraceDelay(src->getTraceDelay());
	dest->setPulseLength(src->getPulseLength());
	dest->setAdcBins(src->getBitRange());
	dest->setPolyCfdFraction(src->getCfdFraction());
	dest->setPulseIntegralLow(src->getPulseIntegralLow());
	dest->setPulseIntegralHigh(src->getPulseIntegralHigh());
	dest->setFunctionType(src->getFunctionType());
}


bool nDetRunAction::processStartDetector(nDetDetector* det, double &startTime){
	if(!processDetector(det)) // No detected photons, do not process
		return false;
		
	// Return the time-of-flight from the start detector
	startTime = outData.barTOF;
		
	return true;
}

bool nDetRunAction::processStartImplant(nDetImplant* imp, double &startTime){
	//std::cout<<"!!! Start !!! "<<detector->GetCenterOfMass()->getNumColumns()<<" "<<imp->getCenterOfMass()->getNumColumns()<<std::endl;
	if(!processImplant(imp)){ // No detected photons, do not process
		//std::cout<<"No photons"<<std::endl;
		return false;
	}
	//std::cout<<"Photons"<<std::endl;
	// Return the time-of-flight from the start detector
	startTime = outImplantData.barTOF;

	return true;
}

void nDetRunAction::process(){
	while(this->scatterEvent()){
	}
	
	if(counter){ 
		debugData.photonsProd.clear(); // A little messy :(
		for(short i = evtData.nScatters+1; i >= 2; i--) // The first recoil particle ID is equal to 2
			debugData.photonsProd.push_back(counter->getPhotonCount(i));
	}

	short detID = 0;
	double startTime;
	if(processStartDetector(startDetector, startTime)){ // Check for valid start signal
			for(std::vector<nDetDetector>::iterator iter = userDetectors.begin(); iter != userDetectors.end(); iter++){
				// Skip the start detector because we already processed it
				if(&(*iter) != startDetector && !processDetector(&(*iter))){ // Skip events with no detected photons
					detID++;
					continue;
				}

				// Update the time-of-flight of the event
				outData.barTOF = outData.barTOF - startTime;
				
				// Push data onto the output branch for multiple detectors
				if(userDetectors.size() > 1)
					multData.Append(outData, detID++);
			}
		}
	else if(endImplant){ // Un-triggered mode (default)
		for(std::vector<nDetImplant>::iterator iter = userImplants.begin(); iter != userImplants.end(); iter++){
	//std::cout<<"!!! End !!! "<<detector->GetCenterOfMass()->getNumColumns()<<" "<<iter->getCenterOfMass()->getNumColumns()<<std::endl;
			if(!processImplant(&(*iter))){ // Skip the start detector because we already processed it
				detID++;
				continue;
			}
			
			// Push data onto the output branch for multiple detectors
			if(userImplants.size() > 1)
				multData.Append(outImplantData, detID++);
		}
	}
	else if(endDetector){ // Un-triggered mode (default)
		for(std::vector<nDetDetector>::iterator iter = userDetectors.begin(); iter != userDetectors.end(); iter++){
			if(!processDetector(&(*iter))){ // Skip the start detector because we already processed it
				detID++;
				continue;
			}
			
			// Push data onto the output branch for multiple detectors
			if(userDetectors.size() > 1 || userImplants.size()>1)
				multData.Append(outData, detID++);
		}
	}
	
	else{ // Start triggered mode
		if(processStartImplant(startImplant, startTime)){ // Check for valid start signal
			for(std::vector<nDetImplant>::iterator iter = userImplants.begin(); iter != userImplants.end(); iter++){
	//std::cout<<"!!! Detector !!! "<<detector->GetCenterOfMass()->getNumColumns()<<" "<<iter->getCenterOfMass()->getNumColumns()<<std::endl;
				// Skip the start detector because we already processed it
				if(&(*iter) != startImplant && !processImplant(&(*iter))){ // Skip events with no detected photons
					detID++;
					continue;
				}

				// Update the time-of-flight of the event
				outImplantData.barTOF = outImplantData.barTOF - startTime;
				
				// Push data onto the output branch for multiple detectors
				if(userImplants.size() > 1)
					multData.Append(outImplantData, detID++);
			}
		}
	}
	
	
	if(outputMultiDebug) 
		multData.Append(debugData, evtData.nScatters);
		
	
	// Write the data (mutex protected, thread safe).
	nDetMasterOutputFile::getInstance().fillBranch(data); // The master output file is a singleton class.

	// Clear all data structures.
	data.clear();

	// Clear all statistics.
	for(std::vector<nDetDetector>::iterator iter = userDetectors.begin(); iter != userDetectors.end(); iter++)
		iter->clear();

	for(std::vector<nDetImplant>::iterator iter = userImplants.begin(); iter != userImplants.end(); iter++)
		iter->clear();
	
	if(stacking) stacking->Reset();
	if(tracking) tracking->Reset();
	if(stepping) stepping->Reset();
}

void nDetRunAction::setActions(userActionManager *action){
	eventAction = action->getEventAction();
	stacking = action->getStackingAction();
	tracking = action->getTrackingAction();
	stepping = action->getSteppingAction();

	// Get the photon counter
	counter = stacking->GetCounter();
}

bool nDetRunAction::AddDetectedPhoton(const G4Step *step, const double &mass/*=1*/){
	if(!step->GetPostStepPoint()){
		Display::WarningPrint("INVALID POST POINT!", "nDetRunAction");
		return false;
	}

	// Find which detector this optical photon is inside.
	G4int copyNum = step->GetPostStepPoint()->GetTouchable()->GetCopyNumber();
	bool foundMatch=false, isLeft=0, isImplant=0;
	G4ThreeVector *detPos;
	G4RotationMatrix *detRot;
	centerOfMass *hitDetPmtL;
	centerOfMass *hitDetPmtR;
	centerOfMass *hitDetPmtI;

	for(std::vector<nDetDetector>::iterator iter = userDetectors.begin(); iter != userDetectors.end(); iter++){
		if(iter->checkPmtCopyNumber(copyNum, isLeft)){
			foundMatch = true;
			detPos = iter->getPosition();
			detRot = iter->getRotation();
			hitDetPmtL = iter->getCenterOfMassL();
			hitDetPmtR = iter->getCenterOfMassR();
			hitDetPmtL->setNumColumns(iter->GetNumPmtColumns());
			hitDetPmtL->setNumRows(iter->GetNumPmtRows());
			hitDetPmtR->setNumColumns(iter->GetNumPmtColumns());
			hitDetPmtR->setNumRows(iter->GetNumPmtRows());
			break;
		}
	}
	for(std::vector<nDetImplant>::iterator iter = userImplants.begin(); iter != userImplants.end(); iter++){
		if(iter->checkPmtCopyNumber(copyNum, isImplant)){
			foundMatch = true;
			detPos = iter->getPosition();
			detRot = iter->getRotation();
			hitDetPmtI = iter->getCenterOfMass();
			hitDetPmtI->setNumColumns(iter->GetNumPmtColumns());
			hitDetPmtI->setNumRows(iter->GetNumPmtRows());
			break;
		}
	}
	
	if(!foundMatch){
		Display::WarningPrint("Failed to find matching detector for detected photon?", "nDetRunAction");
		return false;
	}
	
	double energy = step->GetTrack()->GetTotalEnergy();
	double time = step->GetPostStepPoint()->GetGlobalTime();
	G4ThreeVector position = step->GetPostStepPoint()->GetPosition() - (*detPos);
	position = (*detRot)*position;
	
	if(isLeft){
		if(hitDetPmtL->addPoint(energy, time, position, mass))
			return true;
	}
	else if(isImplant){
		if(hitDetPmtI->addPoint(energy, time, position, mass)){
			return true;
		}
	}
	else{
		if(hitDetPmtR->addPoint(energy, time, position, mass))
			return true;
	}
	return false;
}


bool nDetRunAction::AddDetectedSecondary(const G4Step *step, const double &energy){
	if(!step->GetPostStepPoint()){
		Display::WarningPrint("INVALID POST POINT!", "nDetRunAction");
		return false;
	}

	/*
	// Find which detector this secondary is inside.
	G4int copyNum = step->GetPostStepPoint()->GetTouchable()->GetCopyNumber();
	bool foundMatch=false, isLeft=0, isImplant=0;
	G4ThreeVector *detPos;
	G4RotationMatrix *detRot;
	centerOfMass *hitDetPmtL;
	centerOfMass *hitDetPmtR;
	centerOfMass *hitDetPmtI;

	for(std::vector<nDetDetector>::iterator iter = userDetectors.begin(); iter != userDetectors.end(); iter++){
		if(iter->checkPmtCopyNumber(copyNum, isLeft)){
			foundMatch = true;
			//detPos = iter->getPosition();
			//detRot = iter->getRotation();
			hitDetPmtL = iter->getCenterOfMassL();
			hitDetPmtR = iter->getCenterOfMassR();
			//hitDetPmtL->setNumColumns(iter->GetNumPmtColumns());
			//hitDetPmtL->setNumRows(iter->GetNumPmtRows());
			//hitDetPmtR->setNumColumns(iter->GetNumPmtColumns());
			//hitDetPmtR->setNumRows(iter->GetNumPmtRows());
			break;
		}
	}
	for(std::vector<nDetImplant>::iterator iter = userImplants.begin(); iter != userImplants.end(); iter++){
		if(iter->checkPmtCopyNumber(copyNum, isImplant)){
			foundMatch = true;
			//detPos = iter->getPosition();
			//detRot = iter->getRotation();
			hitDetPmtI = iter->getCenterOfMass();
			//hitDetPmtI->setNumColumns(iter->GetNumPmtColumns());
			//hitDetPmtI->setNumRows(iter->GetNumPmtRows());
			break;
		}
	}
	
	if(!foundMatch){
		Display::WarningPrint("Failed to find matching detector for detected secondary?", "nDetRunAction");
		std::cout << "Energy is " << energy << std::endl;
		// Not sure what it means when this fails...?
		return false;
	}
	*/
	centerOfMass *pmtI = startImplant->getCenterOfMass();

	pmtI->addSecondaryElectron(energy);
	return true;
}


bool nDetRunAction::scatterEvent(){
	if(primaryTracks.size() <= 1)
		return false;

	// Get a primary particle step off the stack.
	primaryTrackInfo *priTrack = &primaryTracks.back();

	if(evtData.nScatters++ == 0){ // First scatter event (CAREFUL! These are in reverse order).
		if(primaryTracks.size() < 2){
			if(verbose){
				Display::WarningPrint("Expected at least two events in primaryTracks.", "nDetRunAction");
				priTrack->print();
			}
			return false;
		}
		if(outputDebug || outputMultiDebug){
			G4ThreeVector firstScatter = primaryTracks.at(0).pos + primaryTracks.at(1).pos;
			debugData.nFirstScatterTime = primaryTracks.at(1).gtime;
			debugData.nFirstScatterLen = firstScatter.mag();
		}
	}
	
	// Add the deposited energy to the total.
	evtData.nDepEnergy += priTrack->dkE;

	if(verbose){
		G4ThreeVector vertex = priTrack->pos;
		std::cout << " id=" << priTrack->trackID << " "; priTrack->print();	
		std::cout << "	pos=(" << vertex.getX() << ", " << vertex.getY() << ", " << vertex.getZ() << ")\n";
		std::cout << "	dE=" << priTrack->dkE << ", dE2=" << priTrack->kE << ", size=" << primaryTracks.size() << std::endl;
	}

	if(outputDebug || outputMultiDebug){
		G4int segCol, segRow;
		G4ThreeVector vertex = priTrack->pos;
		getSegmentFromCopyNum(priTrack->copyNum, segCol, segRow);
		debugData.Append(vertex.getX(), vertex.getY(), vertex.getZ(), priTrack->angle, priTrack->plength, (priTrack->gtime - debugData.nEnterTime), 
		                 priTrack->dkE, segCol, segRow, 0, priTrack->atomicMass, priTrack->inScint);
		
		// Update the neutron scatter center-of-mass
		debugData.nComX += priTrack->dkE*vertex.getX();
		debugData.nComY += priTrack->dkE*vertex.getY();
		debugData.nComZ += priTrack->dkE*vertex.getZ();
		debugData.neutronWeight += priTrack->dkE;
	}
	
	primaryTracks.pop_back();
	
	return true;
}

void nDetRunAction::initializeNeutron(const G4Step *step){
	evtData.nInitEnergy = step->GetPreStepPoint()->GetKineticEnergy();
	if(outputDebug || outputMultiDebug){
		debugData.nEnterTime = step->GetTrack()->GetGlobalTime();
		debugData.nEnterPosX = step->GetPostStepPoint()->GetPosition().getX();
		debugData.nEnterPosY = step->GetPostStepPoint()->GetPosition().getY();
		debugData.nEnterPosZ = step->GetPostStepPoint()->GetPosition().getZ();
		//debugData.nTimeInMat = 0;
		//debugData.nTimeInMat = track->GetGlobalTime() - debugData.nEnterTime;
		//debugData.nTimeInMat = track->GetGlobalTime(); //assuming particle enters at t=0
		debugData.nExitPosX = 0;
		debugData.nExitPosY = 0;
		debugData.nExitPosZ = 0;
		primaryTracks.clear();
		primaryTracks.push_back(step);
		prevDirection = primaryTracks.back().dir;
		prevPosition = primaryTracks.back().pos;
	}
	if(verbose){ 
		std::cout << "IN: "; primaryTracks.back().print();
	}
}

bool nDetRunAction::scatterNeutron(const G4Step *step){
	double dE = -1*(step->GetPostStepPoint()->GetKineticEnergy()-step->GetPreStepPoint()->GetKineticEnergy());
	if(dE > 0){
		G4Track *track = step->GetTrack();
		primaryTracks.push_back(step);
		if(outputDebug || outputMultiDebug){
			primaryTracks.back().getAngle(prevDirection);
			primaryTracks.back().getPathLength(prevPosition);
			prevDirection = primaryTracks.back().dir;
			prevPosition = primaryTracks.back().pos;
		}
		if(verbose){
			std::cout << "SC: "; primaryTracks.back().print();
		}
		if(track->GetKineticEnergy() < KINETIC_ENERGY_THRESHOLD){ // Check particle kinematic energy.
			evtData.nAbsorbed = true;
			track->SetTrackStatus(fStopAndKill); // Kill the neutron.
			if(verbose) 
				std::cout << "OT: final kE=0 (absorbed)" << std::endl;
		}

		return true;
	}
	return false;
}

void nDetRunAction::finalizeNeutron(const G4Step *step){
	evtData.nAbsorbed = false;
	G4Track *track = step->GetTrack();
	if(outputDebug || outputMultiDebug){
		debugData.nTimeInMat = track->GetGlobalTime() - debugData.nEnterTime;
		debugData.nExitPosX = step->GetPreStepPoint()->GetPosition().getX();
		debugData.nExitPosY = step->GetPreStepPoint()->GetPosition().getY();
		debugData.nExitPosZ = step->GetPreStepPoint()->GetPosition().getZ();
	}
	if(verbose) 
		std::cout << "OT: final kE=" << track->GetKineticEnergy() << std::endl;
}
