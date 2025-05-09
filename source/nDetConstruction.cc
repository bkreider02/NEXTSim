#include "G4LogicalVolume.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4VPhysicalVolume.hh"

#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SubtractionSolid.hh"
#include "G4Box.hh"

#include "G4GeometryManager.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"

#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"

#include "nDetConstruction.hh"
#include "nDetConstructionMessenger.hh"
#include "nDetThreadContainer.hh"
#include "nDetParticleSource.hh"
#include "nDetWorld.hh"
#include "termColors.hh"
#include "optionHandler.hh" // split_str

///////////////////////////////////////////////////////////////////////////////
// class nDetConstruction
///////////////////////////////////////////////////////////////////////////////

nDetConstruction &nDetConstruction::getInstance(){
	// The only instance
	// Guaranteed to be lazy initialized
	// Guaranteed that it will be destroyed correctly
	static nDetConstruction instance;
	return instance;
}

nDetConstruction::nDetConstruction(){
	currentDetector = NULL;
	currentImplant = NULL;

	fDetectorMessenger = new nDetConstructionMessenger(this);

	fCheckOverlaps = false;
	
	// Initialize the detector parameter messenger
	params.InitializeMessenger();

	// Link the materials handler to the detector parameter object
	params.SetMaterials(&materials);
	
	// Setup the experimental setup area
	expHall = new nDetWorld();
}

nDetConstruction::~nDetConstruction(){
}

G4VPhysicalVolume* nDetConstruction::Construct(){
	if(!expHall->getPhysicalVolume()){
		this->ConstructDetector();
		this->ConstructImplant();
	}
	
	return expHall->getPhysicalVolume();
}
	
G4VPhysicalVolume* nDetConstruction::ConstructDetector(){
	if(!materials.materialsAreDefined())
		materials.initialize();

	// Build experiment hall.
	expHall->buildExpHall(&materials);

	for(auto det : userDetectors){
		currentDetector = det;
		// Build the detector

		det->construct();
		// Place all detectors.
		det->placeDetector(expHall->getLogicalVolume());
	}

	//for(auto imp : userImplants){
	//	currentImplant = imp;
	//
	//	imp->construct();
	//	imp->placeImplant(expHall->getLogicalVolume());
	//}
	
	return expHall->getPhysicalVolume();
}

G4VPhysicalVolume* nDetConstruction::ConstructImplant(){
	//if(!materials.materialsAreDefined())
	//	materials.initialize();

	// Build experiment hall.
	//expHall->buildExpHall(&materials);

	// Place all detectors.
	for(auto imp : userImplants){
		currentImplant = imp;

		// Build the detector
		imp->construct();
		// Place the detector into the world.
		imp->placeImplant(expHall->getLogicalVolume());
	}
	//for(auto det : userDetectors){
	//	currentDetector = det;
	//	// Build the detector
	//	det->construct();
	//	// Place all detectors.
	//	det->placeDetector(expHall->getLogicalVolume());
	//}

	return expHall->getPhysicalVolume();
}

void nDetConstruction::ClearGeometry(){
	// Clean-up previous geometry
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();
    G4LogicalSkinSurface::CleanSurfaceTable();
    G4LogicalBorderSurface::CleanSurfaceTable();
	G4SolidStore::GetInstance()->Clean();
	G4LogicalVolumeStore::GetInstance()->Clean();
	G4PhysicalVolumeStore::GetInstance()->Clean();
	
	// Reset the world volume. Why is this needed? CRT
	expHall->reset();
	
	// Clear previous construction.
	for(auto det : userDetectors)
		delete det;
	userDetectors.clear();
	for(auto imp : userImplants)
		delete imp;
	userImplants.clear();
	
	// Reset the scintillator copy number.
	params.SetScintillatorCopyNumber(1);
}

void nDetConstruction::ClearImplantGeometry(){
	// Clean-up previous geometry
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();
    G4LogicalSkinSurface::CleanSurfaceTable();
    G4LogicalBorderSurface::CleanSurfaceTable();
	G4SolidStore::GetInstance()->Clean();
	G4LogicalVolumeStore::GetInstance()->Clean();
	G4PhysicalVolumeStore::GetInstance()->Clean();
	
	// Reset the world volume. Why is this needed? CRT
	expHall->reset();
	
	// Clear previous construction.
	for(auto det : userDetectors)
		delete det;
	userDetectors.clear();
	for(auto imp : userImplants)
		delete imp;
	userImplants.clear();
	
	// Reset the scintillator copy number.
	params.SetScintillatorCopyNumber(1);
}

void nDetConstruction::UpdateGeometry(){

	// Define new one
	G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());
	G4RunManager::GetRunManager()->DefineWorldVolume(ConstructImplant());
	G4RunManager::GetRunManager()->GeometryHasBeenModified();
	G4RunManager::GetRunManager()->ReinitializeGeometry();

	// Update the particle source
	if(currentDetector)
		nDetParticleSource::getInstance().SetDetector(currentDetector);
	else if(currentImplant)
		nDetParticleSource::getInstance().SetImplant(currentImplant);

	// Update the detector lists of all user run actions
	nDetThreadContainer *container = &nDetThreadContainer::getInstance();
	if(container->getMultithreadingMode())
		container->getMasterRunAction()->updateDetector(this);
	for(size_t index = 0; index < container->size(); index++){
		container->getActionManager(index)->getRunAction()->updateDetector(this);
	}
}

bool nDetConstruction::AddGeometry(const G4String &geom){
	// Define a new detector of the specified type
	nDetDetector *newDetector = nDetDetectorTypes::getDetectorType(geom, this, &materials);
	
	if(!newDetector){ // Invalid detector type
		std::cout << Display::ErrorStr("nDetConstruction") << "User specified un-recognized detector type (" << geom << ")!" << Display::ResetStr() << std::endl;
		return false;
	}
	
	// Set the current detector to the new one
	currentDetector = newDetector;

	// Segment the PMT photo-sensitive surface
	if(params.PmtIsSegmented()){
		setSegmentedPmt();
		if(!gainMatrixFilename.empty())
			loadPmtGainMatrix();
	}
	
	// Load the anode quantum efficiency
	if(!spectralResponseFilename.empty())
		loadPmtSpectralResponse();

	// Copy the center-of-mass calculators to the new detector
	centerOfMass *cmL = currentDetector->getCenterOfMassL();
	centerOfMass *cmR = currentDetector->getCenterOfMassR();
	currentDetector->copyCenterOfMass(center[0], center[1]);

	// Segment the PMT anodes
	if(params.PmtIsSegmented()){
		cmL->setSegmentedPmt(&params);
		cmR->setSegmentedPmt(&params);	

		// Copy the PMT anode gain matrix
		cmL->copyGainMatrix(&center[0]);
		cmR->copyGainMatrix(&center[1]);
	}

	// Copy the PMT anode quantum efficiency curve
	if(center[0].getPmtResponse()->getSpectralResponseEnabled())
		cmL->copySpectralResponse(&center[0]);
	if(center[1].getPmtResponse()->getSpectralResponseEnabled())
		cmR->copySpectralResponse(&center[1]);

	// Add the new detector assembly to the vector of detectors
	userDetectors.push_back(currentDetector);
	
	// Enable/disable overlap checking
	currentDetector->setCheckOverlaps(fCheckOverlaps);

	// Update the detector's copy numbers.
	currentDetector->setParentCopyNumber(userDetectors.size()-1);

	// Set true isotropic source mode for multiple detectors
	if(userDetectors.size() > 1)
		nDetParticleSource::getInstance().SetRealIsotropicMode(true);

	// Set true isotropic source mode for multiple detectors
	if(userImplants.size() > 1)
		nDetParticleSource::getInstance().SetRealIsotropicMode(true);
	
	return true;
}

bool nDetConstruction::AddImplantGeometry(const G4String &geom){
	// Define a new detector of the specified type
	nDetImplant *newImplant = nDetDetectorTypes::getImplantType(geom, this, &materials);
	
	if(!newImplant){ // Invalid detector type
		std::cout << Display::ErrorStr("nDetConstruction") << "User specified un-recognized detector type (" << geom << ")!" << Display::ResetStr() << std::endl;
		return false;
	}
	
	// Set the current detector to the new one
	currentImplant = newImplant;

	// Segment the PMT photo-sensitive surface
	if(params.PmtIsSegmented()){
		setSegmentedPmt();
		if(!gainMatrixFilename.empty()){
			loadPmtGainMatrix();
		}
	}
	
	// Load the anode quantum efficiency
	if(!spectralResponseFilename.empty()){
		loadPmtSpectralResponse();
	}

	// Copy the center-of-mass calculators to the new detector
	centerOfMass *cmI = currentImplant->getCenterOfMass();
	currentImplant->copyCenterOfMass(center[0]);

	// Segment the PMT anodes
	if(params.PmtIsSegmented()){
		cmI->setSegmentedPmt(&params);

		// Copy the PMT anode gain matrix
		cmI->copyGainMatrix(&center[0]);
	}

	// Copy the PMT anode quantum efficiency curve
	if(center[0].getPmtResponse()->getSpectralResponseEnabled())
		cmI->copySpectralResponse(&center[0]);

	// Add the new detector assembly to the vector of detectors
	userImplants.push_back(currentImplant);
	
	// Enable/disable overlap checking
	currentImplant->setCheckOverlaps(fCheckOverlaps);

	// Update the detector's copy numbers.
	currentImplant->setParentCopyNumber(userImplants.size()-1);

	// Set true isotropic source mode for multiple detectors
	if(userImplants.size() > 1)
		nDetParticleSource::getInstance().SetRealIsotropicMode(true);
	
	return true;
}

void nDetConstruction::setSegmentedPmt(){
	center[0].setSegmentedPmt(&params);
	center[1].setSegmentedPmt(&params);
	std::cout << " nDetConstruction: Set segmented PMTs with WxH=(" << params.GetPmtWidth() << " x " << params.GetPmtHeight() << ") and " << params.GetNumPmtColumns() << " columns and " << params.GetNumPmtRows() << " rows.\n";
}

bool nDetConstruction::loadPmtSpectralResponse(){
	if(!(center[0].loadSpectralResponse(spectralResponseFilename.c_str()) && center[1].loadSpectralResponse(spectralResponseFilename.c_str()))){
		Display::ErrorPrint("Failed to load PMT spectral response from file!", "nDetConstruction");
		return false;
	}
	std::cout << " nDetConstruction: Successfully loaded PMT spectral response function\n";
	spectralResponseFilename = ""; // Unset the filename so that it will not be loaded again later
	return true;
}

bool nDetConstruction::loadPmtGainMatrix(){
	if(!(center[1].loadGainMatrix(gainMatrixFilename.c_str()) && center[0].loadGainMatrix(gainMatrixFilename.c_str()))){
		Display::ErrorPrint("Failed to load PMT anode gain matrix from file!", "nDetConstruction");
		return false;
	}
	std::cout << " nDetConstruction: Successfully loaded PMT anode gain matrix\n";
	gainMatrixFilename = ""; // Unset the filename so that it will not be loaded again later
	return true;
}

void nDetConstruction::AddLightGuideGDML(const G4String &input){
	if(currentDetector)
		currentDetector->addLightGuideGDML(input);
	else if(currentImplant)
		currentImplant->addLightGuideGDML(input);
	else
		Display::ErrorPrint("Cannot add GDML light-guide before a detector is defined!", "nDetConstruction");
}

void nDetConstruction::AddGrease(const G4String &input){
	if(currentDetector)
		currentDetector->addGreaseLayer(input);
	else if(currentImplant)
		currentImplant->addGreaseLayer(input);
	else
		Display::ErrorPrint("Cannot add grease layer before a detector is defined!", "nDetConstruction");
}

void nDetConstruction::AddDiffuser(const G4String &input){
	if(currentDetector)
		currentDetector->addDiffuserLayer(input);
	else if(currentImplant)
		currentImplant->addDiffuserLayer(input);
	else
		Display::ErrorPrint("Cannot add diffuser layer before a detector is defined!", "nDetConstruction");
}

void nDetConstruction::AddLightGuide(const G4String &input){
	if(currentDetector)
		currentDetector->addLightGuideLayer(input);
	else if(currentImplant)
		currentImplant->addLightGuideLayer(input);
	else
		Display::ErrorPrint("Cannot add light-guide before a detector is defined!", "nDetConstruction");
}

void nDetConstruction::AddSegmentedLightGuide(const G4String &input){
	if(currentImplant)
		currentImplant->addSegmentedLightGuideLayer(input);
	else
		Display::ErrorPrint("Cannot add segmented light-guide before a detector is defined!", "nDetConstruction");
}

void nDetConstruction::AddPhoswich(const G4String & input) {
	if (currentImplant)
		currentImplant->addPhoswichLayer(input);
	else
		Display::ErrorPrint("Cannot add phoswich scintillator before a detector is defined!", "nDetConstruction");
}

void nDetConstruction::AddBox(const G4String &input) {
	if (currentImplant)
		currentImplant->addBox(input);
	else
		Display::ErrorPrint("Cannot add a box before a detector is defined!", "nDetConstruction");
}


void nDetConstruction::SetDomeParameters(const G4String &input) {
	if (currentImplant)
		currentImplant->setDomeParameters(input);
	else
		Display::ErrorPrint("Cannot set dome parameters before a detector is defined!", "nDetConstruction");
}


void nDetConstruction::AddDetectorArray(const G4String &input){
	// Expects a space-delimited string of the form:
	//  "addDiffuserLayer width(mm) height(mm) thickness(mm) material"
	std::vector<std::string> args;
	unsigned int Nargs = split_str(input, args);
	if(Nargs < 3){
		std::cout << " nDetConstruction: Invalid number of arguments given to ::AddDetectorArray(). Expected 5, received " << Nargs << ".\n";
		std::cout << " nDetConstruction:  SYNTAX: addArray <geom> <r0> <startTheta> <stopTheta> <Ndet>\n";
		return;
	}
	double r0 = strtod(args.at(1).c_str(), NULL)*cm;
	double startTheta = strtod(args.at(2).c_str(), NULL);
	double stopTheta = strtod(args.at(3).c_str(), NULL);
	int Ndet = strtol(args.at(4).c_str(), NULL, 10);
	
	double dTheta = 0;
	if(Ndet > 1)
		dTheta = (stopTheta-startTheta)/(Ndet-1);
	for(int i = 0; i < Ndet; i++){
		//std::cout << " nDetConstruction: Adding detector (type=" << args.at(0) << ", r=" << r0 << ", theta=" << startTheta+dTheta*i << ")\n";
		params.SetPositionCylindrical(G4ThreeVector(r0, startTheta+dTheta*i, 0));
		params.SetRotation(G4ThreeVector(90, 0, startTheta+dTheta*i));
		AddGeometry(args.at(0));
	}
}

void nDetConstruction::AddImplantArray(const G4String &input){
	// Expects a space-delimited string of the form:
	//  "addDiffuserLayer width(mm) height(mm) thickness(mm) material"
	std::vector<std::string> args;
	unsigned int Nargs = split_str(input, args);
	if(Nargs < 3){
		std::cout << " nDetConstruction: Invalid number of arguments given to ::AddDetectorArray(). Expected 5, received " << Nargs << ".\n";
		std::cout << " nDetConstruction:  SYNTAX: addArray <geom> <r0> <startTheta> <stopTheta> <Ndet>\n";
		return;
	}
	double r0 = strtod(args.at(1).c_str(), NULL);
	double startTheta = strtod(args.at(2).c_str(), NULL);
	double stopTheta = strtod(args.at(3).c_str(), NULL);
	int Ndet = strtol(args.at(4).c_str(), NULL, 10);
	
	if(args.at(0)=="hagrid"){
		// in this case r0 will be used as the half-width of a square around the point
		if(Ndet!=8 && Ndet!=16){
			std::cout<<"!!!!!! ----- !!!!!! \n !!!!!! For preliminary testing please use Ndet=8 or 16 !!!!! \n !!!!!! ----- !!!!!!"<<std::endl;
		}
		// above or below implant
		G4double yp[8] = {0,1,1,1,0,-1,-1,-1};
		// left or right of implant
		G4double zp[8] = {1,1,0,-1,-1,-1,0,1};
		for(int i=0;i<8;i++){
			// shift along beam axis
			G4double xp = 2.88*cm;
			params.SetPosition(G4ThreeVector(xp,yp[i]*r0*cm,zp[i]*r0*cm));
			params.SetRotation(G4ThreeVector(0,270,0));
			AddImplantGeometry(args.at(0));
			if(Ndet==16){
				params.SetPosition(G4ThreeVector(-1*xp,yp[i]*r0*cm,zp[i]*r0*cm));
				params.SetRotation(G4ThreeVector(0,90,0));
				AddImplantGeometry(args.at(0));
			}
		}
	}
	else{
		double dTheta = 0;
		if(Ndet > 1)
			dTheta = (stopTheta-startTheta)/(Ndet-1);
		for(int i = 0; i < Ndet; i++){
			//std::cout << " nDetConstruction: Adding detector (type=" << args.at(0) << ", r=" << r0 << ", theta=" << startTheta+dTheta*i << ")\n";
			params.SetPositionCylindrical(G4ThreeVector(r0, startTheta+dTheta*i, 0));
			params.SetRotation(G4ThreeVector(90, 0, startTheta+dTheta*i));
			AddImplantGeometry(args.at(0));
		}
	}
}

void nDetConstruction::BuildExp(std::string expName_){
	expHall->SetExp(expName_);
}

void nDetConstruction::SetLightYieldMultiplier(const G4double &yield){ 
	materials.setLightYield(yield);
}

void nDetConstruction::PrintAllDetectors() const {
	int detCount = 0;
	for(auto det : userDetectors){
		std::cout << "***********************************************************\n";
		std::cout << " Detector ID          = " << detCount++ << std::endl;
		det->Print();
	}
	std::cout << "***********************************************************\n";
}

void nDetConstruction::GetCopiesOfDetectors(std::vector<nDetDetector> &detectors) const {
	for(auto det : userDetectors)
		detectors.push_back(det->clone());
}

void nDetConstruction::GetCopiesOfImplants(std::vector<nDetImplant> &implants) const {
	for(auto imp : userImplants)
		implants.push_back(imp->clone());
}
