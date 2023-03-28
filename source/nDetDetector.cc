
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4OpticalSurface.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4Trap.hh"

#include "nDetDetector.hh"
#include "nDetDetectorLayer.hh"
#include "nDetConstruction.hh"
#include "nDetDetectorMessenger.hh"
#include "optionHandler.hh" // split_str

///////////////////////////////////////////////////////////////////////////////
// class nDetDetectorParams
///////////////////////////////////////////////////////////////////////////////

void nDetDetectorParams::InitializeMessenger(){
	if(fMessenger) 
		return;
	fMessenger = new nDetDetectorMessenger(this);
}

void nDetDetectorParams::SetDetectorWidth(const G4double &val){ 
	fDetectorWidth = val; 
	constantWidth = true;
}

void nDetDetectorParams::SetDetectorHeight(const G4double &val){ 
	fDetectorHeight = val; 
	constantHeight = true;
}

void nDetDetectorParams::SetSegmentWidth(const G4double &val){
	fSegmentWidth = val;
	constantWidth = false;
}

void nDetDetectorParams::SetSegmentHeight(const G4double &val){
	fSegmentHeight = val;
	constantHeight = false;
}

void nDetDetectorParams::SetPmtDimension(const G4String &input){
	// Expects a space-delimited string of the form:
	//  "setPmtDimensions <sizeX> [sizeY]"
	std::vector<std::string> args;
	unsigned int Nargs = split_str(input, args);
	if(Nargs < 1){
		std::cout << " nDetDetectorParams: Invalid number of arguments given to ::SetPmtDimension(). Expected 1, received " << Nargs << ".\n";
		std::cout << " nDetDetectorParams:  SYNTAX: setPmtDimensions <sizeX> [sizeY]\n";
		return;
	}
	G4double sizeX = strtod(args.at(0).c_str(), NULL);
	G4double sizeY = (Nargs >= 2 ? strtod(args.at(1).c_str(), NULL) : -1);
	SetPmtDimension(sizeX, sizeY);
}

void nDetDetectorParams::SetPmtDimension(const G4double &width, const G4double &height/*=-1*/){
	pmtWidth = width;
	pmtHeight = (height > 0 ? height : width);
}

void nDetDetectorParams::SetPmtGapThickness(const G4double &gapSize)
{
	pmtGapThickness = gapSize;
	pmtHasGaps = true;
	pmtPixelHeight = (pmtHeight - (fNumRowsPmt - 1) * pmtGapThickness)/fNumRowsPmt;
	pmtPixelWidth = (pmtWidth - (fNumColumnsPmt - 1) * pmtGapThickness)/fNumColumnsPmt;
}

void nDetDetectorParams::SetPositionCylindrical(const G4ThreeVector &position){ 
	double x = position.getX()*std::cos(position.getY()*deg);
	double z = position.getX()*std::sin(position.getY()*deg);
	detectorPosition = G4ThreeVector(x*cm, position.getZ()*cm, z*cm);
}

void nDetDetectorParams::SetPositionSpherical(const G4ThreeVector &position){ 
	double x = position.getX()*std::sin(position.getY()*deg)*std::cos(position.getZ()*deg); 
	double y = position.getX()*std::sin(position.getY()*deg)*std::sin(position.getZ()*deg); 
	double z = position.getX()*std::cos(position.getY()*deg);
	detectorPosition = G4ThreeVector(x*cm, y*cm, z*cm);
}

void nDetDetectorParams::SetRotation(const G4ThreeVector &rotation){
	detectorRotation = G4RotationMatrix();
	detectorRotation.rotateX(rotation.getX()*deg);
	detectorRotation.rotateY(rotation.getY()*deg); 
	detectorRotation.rotateZ(rotation.getZ()*deg);  
}

void nDetDetectorParams::SetDetectorMaterial(const G4String &material){ 
	detectorMaterialName = material;
}

void nDetDetectorParams::SetWrappingMaterial(const G4String &material){ 
	wrappingMaterialName = material; 
}

void nDetDetectorParams::SetUnsegmented(){
	fSegmentWidth = 0;
	fSegmentHeight = 0;
	fNumColumns = 1;
	fNumRows = 1;
	constantWidth = true;
	constantHeight = true;
}

void nDetDetectorParams::Print() const {
	G4ThreeVector rowX = detectorRotation.rowX();
	G4ThreeVector rowY = detectorRotation.rowY();
	G4ThreeVector rowZ = detectorRotation.rowZ();
	std::cout << " Address              = " << this << std::endl;
	std::cout << " Detector Length (Z)  = " << fDetectorLength << " mm\n";
	std::cout << " Detector Size (X, Y) = " << fDetectorWidth << " x " << fDetectorHeight << " mm^2\n";
	std::cout << " Pmt Size      (X, Y) = " << pmtWidth << " x " << pmtHeight << " mm^2\n";
	std::cout << " Grease Thickness     = " << fGreaseThickness << " mm\n";
	std::cout << " Window Thickness     = " << fWindowThickness << " mm\n";
	std::cout << " Sensitive Thickness  = " << fSensitiveThickness << " mm\n";
	std::cout << " Wrapping Thickness   = " << fWrappingThickness << " mm\n";
	std::cout << " Trapezoid Length     = " << fTrapezoidLength << " mm\n";
	std::cout << " Diffuser Length      = " << fDiffuserLength << " mm\n";
	std::cout << " Wrapping Material    = \"" << wrappingMaterialName << "\"\n";
	std::cout << " Detector Material    = \"" << detectorMaterialName << "\"\n";
	std::cout << " Geometry Type        = \"" << geomType << "\"\n";
	std::cout << " CopyNum  = " << scintCopyNum << std::endl;
	std::cout << " Position = (x=" << detectorPosition.getX() << ", y=" << detectorPosition.getY() << ", z=" << detectorPosition.getZ() << ")\n";
	std::cout << " Unit X   = (" << rowX.getX() << ", " << rowX.getY() << ", " << rowX.getZ() << ")\n";
	std::cout << " Unit Y   = (" << rowY.getX() << ", " << rowY.getY() << ", " << rowY.getZ() << ")\n";
	std::cout << " Unit Z   = (" << rowZ.getX() << ", " << rowZ.getY() << ", " << rowZ.getZ() << ")\n";
	std::cout << " Polished    : " << (fPolishedInterface ? "YES" : "NO") << std::endl;
	std::cout << " Square      : " << (fSquarePMTs ? "YES" : "NO") << std::endl;
	std::cout << " Start       : " << (isStart ? "YES" : "NO") << std::endl;
	if(IsSegmented() || PmtIsSegmented())
		std::cout << " Segmentation:\n";
	if(IsSegmented()){
		std::cout << "  Detector Segments (X, Y) = " << fNumColumns << " x " << fNumRows << std::endl;
		std::cout << "  Cell Size         (X, Y) = " << fSegmentWidth << " x " << fSegmentHeight << " mm^2\n";
	}
	if(PmtIsSegmented()){
		std::cout << "  Pmt Segments      (X, Y) = " << fNumColumnsPmt << " x " << fNumRowsPmt;
		if (pmtHasGaps)
			std::cout << " (with " << pmtGapThickness << " mm gap between anodes)";
		std::cout << std::endl;
	}
}

void nDetDetectorParams::UpdateSize(){
	if(constantWidth)
		fSegmentWidth = (fDetectorWidth-(fNumColumns-1)*fWrappingThickness)/fNumColumns;
	else
		fDetectorWidth = fSegmentWidth*fNumColumns + (fNumColumns-1)*fWrappingThickness;

	if(constantHeight)
		fSegmentHeight = (fDetectorHeight-(fNumRows-1)*fWrappingThickness)/fNumRows;
	else
		fDetectorHeight = fSegmentHeight*fNumRows + (fNumRows-1)*fWrappingThickness;
}

///////////////////////////////////////////////////////////////////////////////
// class nDetDetector
///////////////////////////////////////////////////////////////////////////////

nDetDetector::nDetDetector(nDetConstruction *detector, nDetMaterials *matptr) : nDetDetectorParams(detector->GetDetectorParameters()),
					assembly_logV(NULL), assembly_physV(NULL),
					layerSizeX(0), layerSizeY(0), offsetZ(0),
					parentCopyNum(0), firstSegmentCopyNum(0), lastSegmentCopyNum(0), 
					checkOverlaps(false)
{
	copyCenterOfMass(*detector->GetCenterOfMassL(), *detector->GetCenterOfMassR());
	materials = matptr;
}

nDetDetector::~nDetDetector(){
	/*for(std::vector<nDetWorldObject*>::iterator iter = userLayers.begin(); iter != userLayers.end(); iter++){
		delete (*iter); // This causes a seg-fault
	}
	userLayers.clear();*/
}

nDetDetector nDetDetector::clone() const {
	// Clone this detector's parameters
	nDetDetector retval(*this);
	
	// Copy the center of mass calculators explicitly
	retval.cmL = cmL.clone();
	retval.cmR = cmR.clone();

	return retval;
}

void nDetDetector::getCurrentOffset(G4double &x_, G4double &y_, G4double &z_){
	x_ = layerSizeX;
	y_ = layerSizeY;
	z_ = offsetZ;
}

void nDetDetector::setPositionAndRotation(const G4ThreeVector &pos, const G4RotationMatrix &rot){
	detectorPosition = pos;
	detectorRotation = rot;
}

void nDetDetector::setCurrentOffset(const G4double &x_, const G4double &y_, const G4double &z_){
	layerSizeX = x_;
	layerSizeY = y_;
	offsetZ = z_;
}

void nDetDetector::buildAllLayers(){
	for(std::vector<nDetWorldObject*>::iterator iter = userLayers.begin(); iter != userLayers.end(); iter++){
		(*iter)->construct(this);
	}	
}

void nDetDetector::placeDetector(G4LogicalVolume *parent){
	assembly_physV = new G4PVPlacement(&detectorRotation, detectorPosition, assembly_logV, "Assembly", parent, 0, 0, false);
	assembly_physV->SetCopyNo(parentCopyNum);
	
	// Do some post-placement processing
	afterPlacement();
}

void nDetDetector::clear(){
	cmL.clear();
	cmR.clear();
}

void nDetDetector::copyCenterOfMass(const centerOfMass &left, const centerOfMass &right){
	cmL = left.clone();
	cmR = right.clone();
}

bool nDetDetector::checkPmtCopyNumber(const G4int &num, bool &isLeft) const { 
	isLeft = (num % 2 == 0);
	return (num == 2*parentCopyNum || num == 2*parentCopyNum+1);
}

bool nDetDetector::getSegmentFromCopyNum(const G4int &copyNum, G4int &col, G4int &row) const {
	if(!this->checkCopyNumber(copyNum)) return false;
	col = (copyNum-firstSegmentCopyNum) / fNumRows;
	row = (copyNum-firstSegmentCopyNum) % fNumRows;
	return true;
}

void nDetDetector::addLightGuideGDML(const G4String &input){
	addLayer(new gdmlLightGuideLayer(input));
}

void nDetDetector::addGreaseLayer(const G4String &input){
	addLayer(new greaseLayer(input));
}

void nDetDetector::addDiffuserLayer(const G4String &input){
	addLayer(new diffuserLayer(input));
}

void nDetDetector::addLightGuideLayer(const G4String &input){
	addLayer(new lightGuideLayer(input));
}

void nDetDetector::construct(){
	// Update the size of the assembly in the event it has changed
	UpdateSize(); 

	// Prepare to build the detector and compute the maximum size of the detector volume
	prepareToBuild();

	// Build the assembly volume
	constructAssembly();

	// Get materials and visual attributes
	scintMaterial = materials->getUserDetectorMaterial(detectorMaterialName);
	scintVisAtt = materials->getUserVisAttributes(detectorMaterialName);
	wrappingMaterial = materials->getUserSurfaceMaterial(wrappingMaterialName);
	outerMylar = materials->fMylar;
	wrappingVisAtt = materials->getUserVisAttributes(wrappingMaterialName);
	wrappingOpSurf = materials->getUserOpticalSurface(wrappingMaterialName);

	// Build the geometry
	buildDetector();

	// Generate all user-defined layers.
	buildAllLayers();

	// Attach PMTs.
	constructPSPmts();
}

G4LogicalVolume *nDetDetector::constructAssembly(){
	// Calculate the dimensions of the detector
	assemblyWidth = maxBodySize.getX() + 2*fWrappingThickness;
	assemblyHeight = maxBodySize.getY() + 2*fWrappingThickness;
	assemblyLength = maxBodySize.getZ();

	// Account for the size of the PSPMT
	assemblyWidth = std::max(assemblyWidth, pmtWidth);
	assemblyHeight = std::max(assemblyHeight, pmtHeight);
	assemblyLength += 2*(fGreaseThickness+fWindowThickness+fSensitiveThickness);

	// Account for the additional component layers
	for(std::vector<nDetWorldObject*>::iterator iter = userLayers.begin(); iter != userLayers.end(); iter++){
		if(!(*iter)->decodeString()){
			std::cout << " nDetDetector: Invalid number of arguments given to ::decodeString(). Expected " << (*iter)->getNumRequiredArgs() << " but received " << (*iter)->getNumSuppliedArgs() << ".\n";
			std::cout << " nDetDetector:  SYNTAX: " << (*iter)->syntaxStr() << std::endl;
			continue;
		}
		assemblyWidth = std::max(assemblyWidth, (*iter)->getSizeX());
		assemblyHeight = std::max(assemblyHeight, (*iter)->getSizeY());
		assemblyLength += 2*(*iter)->getSizeZ();
	}

	// Build the assembly box
	G4Box *assembly = new G4Box("assembly", assemblyWidth/2, assemblyHeight/2, assemblyLength/2);
	assembly_logV = new G4LogicalVolume(assembly, materials->fAir, "assembly_logV");
	assembly_logV->SetVisAttributes(materials->visAssembly);

	return assembly_logV;
}

G4PVPlacement *nDetDetector::addToDetectorBody(G4LogicalVolume *volume, const G4String &name/*=""*/, const G4ThreeVector &pos/*=G4ThreeVector(0,0,0)*/, G4RotationMatrix *rot/*=NULL*/){
	return (new G4PVPlacement(rot, pos, volume, (name.empty() ? volume->GetName() : name), assembly_logV, false, 0, checkOverlaps));
}

G4PVPlacement *nDetDetector::addSegmentToBody(G4LogicalVolume *volume, const G4String &name/*=""*/, const G4ThreeVector &pos/*=G4ThreeVector(0,0,0)*/, G4RotationMatrix *rot/*=NULL*/){
	return (new G4PVPlacement(rot, pos, volume, (name.empty() ? volume->GetName() : name), assembly_logV, false, lastSegmentCopyNum++, checkOverlaps));
}

G4PVPlacement *nDetDetector::addLeftComponent(G4LogicalVolume *volume, const G4double &offset, const G4String &name/*=""*/, G4RotationMatrix *rot/*=NULL*/){
	return (new G4PVPlacement(rot, G4ThreeVector(0, 0, offset), volume, (name.empty() ? volume->GetName() : name), assembly_logV, false, 2*parentCopyNum, checkOverlaps));
}

G4PVPlacement *nDetDetector::addRightComponent(G4LogicalVolume *volume, const G4double &offset, const G4String &name/*=""*/, G4RotationMatrix *rot/*=NULL*/){
	return (new G4PVPlacement(rot, G4ThreeVector(0, 0, -offset), volume, (name.empty() ? volume->GetName() : name), assembly_logV, false, 2*parentCopyNum+1, checkOverlaps));
}

void nDetDetector::addMirroredComponents(G4LogicalVolume *volume, const G4double &offset, const G4String &name/*=""*/, G4RotationMatrix *rot/*=NULL*/){
	addLeftComponent(volume, offset, name, rot);
	addRightComponent(volume, offset, name, rot);
}
	
void nDetDetector::addMirroredComponents(G4PVPlacement* &phys1, G4PVPlacement* &phys2, G4LogicalVolume *volume, const G4double &offset, const G4String &name/*=""*/, G4RotationMatrix *rot/*=NULL*/){
	phys1 = addLeftComponent(volume, offset, name, rot);
	phys2 = addRightComponent(volume, offset, name, rot);
}

void nDetDetector::prepareToBuild(){ 
	maxBodySize = G4ThreeVector(fDetectorWidth, fDetectorHeight, fDetectorLength); 
}

void nDetDetector::constructPSPmts(){
	// Build the sensitive PMT surfaces.
	const G4String name = "psSiPM";

	G4double sensitiveZ = offsetZ + fGreaseThickness + fWindowThickness + fSensitiveThickness/2;
	G4double wrappingThickness = fGreaseThickness + fWindowThickness;

	// The optical grease layer.
	G4PVPlacement *grease_physV[2] = {NULL, NULL};
	if(fGreaseThickness > 0){
		G4double greaseZ = offsetZ + fGreaseThickness/2;
	
		G4CSGSolid *grease_solidV = getVolume("window_solidV", pmtWidth, pmtHeight, fGreaseThickness);
		G4LogicalVolume *grease_logV = new G4LogicalVolume(grease_solidV, materials->fGrease, "grease_logV");
		
		grease_logV->SetVisAttributes(materials->visGrease);

		addMirroredComponents(grease_physV[0], grease_physV[1], grease_logV, greaseZ, "Grease");
		
		if(!fPolishedInterface){
			for(std::vector<G4PVPlacement*>::iterator iter = scintBody_physV.begin(); iter != scintBody_physV.end(); iter++){
				new G4LogicalBorderSurface("GreaseInterface", (*iter), grease_physV[0], materials->fGreaseOpSurf);
				new G4LogicalBorderSurface("GreaseInterface", (*iter), grease_physV[1], materials->fGreaseOpSurf);
			}
		}
		
		// Clear all scintillator placements.
		scintBody_physV.clear();
	}

	G4PVPlacement *window_physV[2] = {NULL, NULL};
	if(fWindowThickness > 0){ // The quartz window
		G4double windowZ = offsetZ + fGreaseThickness + fWindowThickness/2;
	
		G4CSGSolid *window_solidV = getVolume("window_solidV", pmtWidth, pmtHeight, fWindowThickness);
		G4LogicalVolume *window_logV = new G4LogicalVolume(window_solidV, materials->fSiO2, "window_logV");
		
		window_logV->SetVisAttributes(materials->visWindow);

		addMirroredComponents(window_physV[0], window_physV[1], window_logV, windowZ, "Quartz");
	}

	// Build the wrapping.
	if(WrappingEnabled() && wrappingThickness > 0){
		G4CSGSolid *boundingBox = getVolume("", pmtWidth, pmtHeight, wrappingThickness);
		G4CSGSolid *wrappingBox = getVolume("", pmtWidth + 2*fWrappingThickness, pmtHeight + 2*fWrappingThickness, wrappingThickness);
		
		G4SubtractionSolid *greaseWrapping = new G4SubtractionSolid("", wrappingBox, boundingBox);
		G4LogicalVolume *greaseWrapping_logV = new G4LogicalVolume(greaseWrapping, wrappingMaterial, "greaseWrapping_logV");
		greaseWrapping_logV->SetVisAttributes(wrappingVisAtt);
		
		G4double wrappingZ = offsetZ + fGreaseThickness/2 + fWindowThickness/2;

		// Place the wrapping around the scintillator.
		G4PVPlacement *greaseWrapping_physV[2];
		
		addMirroredComponents(greaseWrapping_physV[0], greaseWrapping_physV[1], greaseWrapping_logV, wrappingZ, "Wrapping");
		
		if(grease_physV[0] && grease_physV[1]){
			new G4LogicalBorderSurface("Wrapping", grease_physV[0], greaseWrapping_physV[0], wrappingOpSurf);
			new G4LogicalBorderSurface("Wrapping", grease_physV[1], greaseWrapping_physV[1], wrappingOpSurf);
		}
		if(window_physV[0] && window_physV[1]){
			new G4LogicalBorderSurface("Wrapping", window_physV[0], greaseWrapping_physV[0], wrappingOpSurf);
			new G4LogicalBorderSurface("Wrapping", window_physV[1], greaseWrapping_physV[1], wrappingOpSurf);
		}
	}
	
    // The photon sensitive surface
    G4CSGSolid *sensitive_solidV = getVolume(name+"_solidV", pmtWidth, pmtHeight, fSensitiveThickness);
    G4LogicalVolume *sensitive_logV = new G4LogicalVolume(sensitive_solidV, materials->fSilicon, name+"_logV");
    sensitive_logV->SetVisAttributes(materials->visSensitive);
    
    // Logical skin surface.
    new G4LogicalSkinSurface(name, sensitive_logV, materials->fSiliconOpSurf);    

	addMirroredComponents(sensitive_logV, sensitiveZ, name);

    // Move the current offset past the PMT
	layerSizeX = pmtWidth;
	layerSizeY = pmtHeight;
    offsetZ += fGreaseThickness + fWindowThickness + fSensitiveThickness;
}

G4CSGSolid *nDetDetector::getVolume(const G4String &name, const G4double &width, const G4double &height, const G4double &length){
	G4CSGSolid *retval;
	if(fSquarePMTs)
		retval = new G4Box(name, width/2, height/2, length/2);
	else
		retval = new G4Tubs(name, 0, width/2, length/2, 0, 2*CLHEP::pi);
	return retval;
}

G4CSGSolid *nDetDetector::getLightGuideVolume(const G4String &name, const G4double &w1, const double &w2, const double &h1, const double &h2, const G4double &length){
	G4CSGSolid *retval;
	if(fSquarePMTs)
		retval = new G4Trd(name, w1/2, w2/2, h1/2, h2/2, length/2);
	else
		retval = new G4Cons(name, 0, w1/2, 0, w2/2, length/2, 0, 2*CLHEP::pi);
	return retval;
}

void nDetDetector::applyGreaseLayer(){
	this->applyGreaseLayer(layerSizeX, layerSizeY);
}

void nDetDetector::applyGreaseLayer(const G4double &x, const G4double &y, double thickness/*=0*/){
	if(thickness <= 0)
		thickness = fGreaseThickness;
	if(thickness > 0){
		G4CSGSolid *grease_solidV = getVolume("grease", x, y, thickness);
		G4LogicalVolume *grease_logV = new G4LogicalVolume(grease_solidV, materials->fGrease, "grease_logV");
		grease_logV->SetVisAttributes(materials->visGrease);

		// Add the optical grease to the assembly
		addMirroredComponents(grease_logV, offsetZ+thickness/2, "Grease");

		// Offset the other layers to account for the layer of optical grease
		layerSizeX = x;
		layerSizeY = y;
		offsetZ += thickness;
	}
}

void nDetDetector::applyDiffuserLayer(){
	this->applyDiffuserLayer(layerSizeX, layerSizeY, fDiffuserLength);
}

void nDetDetector::applyDiffuserLayer(const G4double &x, const G4double &y, const double &thickness){
	if(thickness > 0){ // Build the light diffusers (if needed)
		G4CSGSolid *lightDiffuser = getVolume("lightDiffuser", x, y, thickness);
		G4LogicalVolume *lightDiffuserLog = new G4LogicalVolume(lightDiffuser, materials->fSiO2, "lightDiffuser_logV");

		// Add the optical grease to the assembly
		addMirroredComponents(lightDiffuserLog, offsetZ+thickness/2, "Diffuser");

		// Offset the other layers to account for the light-diffuser
		layerSizeX = x;
		layerSizeY = y;
		offsetZ += thickness;
	}
}

void nDetDetector::applyLightGuide(){
	this->applyLightGuide(layerSizeX, pmtWidth, layerSizeY, pmtHeight, fTrapezoidLength);
}

void nDetDetector::applyLightGuide(const G4double &x2, const G4double &y2){
	this->applyLightGuide(layerSizeX, x2, layerSizeY, y2, fTrapezoidLength);
}

void nDetDetector::applyLightGuide(const G4double &x1, const G4double &x2, const G4double &y1, const G4double &y2, const double &thickness){
    if(thickness > 0){ // Build the light guides (if needed)
        const G4double trapAngleXZ = std::atan2(2*thickness, x1-x2);
        const G4double trapAngleYZ = std::atan2(2*thickness, y1-y2);
        
        const G4double deltaX = fWrappingThickness/std::sin(trapAngleXZ);
        const G4double deltaY = fWrappingThickness/std::sin(trapAngleYZ);
        
        G4double trapezoidZ = offsetZ + thickness/2;
        
        std::string trapName = "Acrylic";
        
        // Build the light-guide.
        G4CSGSolid *lightGuide = getLightGuideVolume("lightGuide", x1, x2, y1, y2, thickness);
	    G4LogicalVolume *lightGuideLog = new G4LogicalVolume(lightGuide, materials->fSiO2, "lightGuide_logV");	

		G4RotationMatrix *rightRotation = new G4RotationMatrix();
		rightRotation->rotateX(CLHEP::pi);

		// Place the light guides.
		G4PVPlacement *trapPhysicalL = addLeftComponent(lightGuideLog, trapezoidZ, trapName);
		G4PVPlacement *trapPhysicalR = addRightComponent(lightGuideLog, trapezoidZ, trapName, rightRotation);

		// Build the wrapping.
		if(WrappingEnabled()){
			G4CSGSolid *wrappingSolid = getLightGuideVolume("wrapping", x1+2*deltaX, x2+2*deltaY, y1+2*deltaX, y2+2*deltaY, thickness);
			G4SubtractionSolid *wrapping = new G4SubtractionSolid("wrapping", wrappingSolid, lightGuide);
			G4LogicalVolume *wrapping_logV = new G4LogicalVolume(wrapping, wrappingMaterial, "wrapping_logV");		
			wrapping_logV->SetVisAttributes(wrappingVisAtt);

			// Place the wrapping around the light guides.
			G4PVPlacement *trapWrappingL = addLeftComponent(wrapping_logV, trapezoidZ, trapName);
			G4PVPlacement *trapWrappingR = addRightComponent(wrapping_logV, trapezoidZ, trapName, rightRotation);
		
			// Reflective wrapping.
			new G4LogicalBorderSurface("Wrapping", trapPhysicalL, trapWrappingL, wrappingOpSurf);
			new G4LogicalBorderSurface("Wrapping", trapPhysicalR, trapWrappingR, wrappingOpSurf);
		}
		
		// Offset the other layers to account for the light-guide
		layerSizeX = x2;
		layerSizeY = y2;
    	offsetZ += thickness;
	}
}

void nDetDetector::loadGDML(gdmlSolid *solid){
	if(!solid || !solid->isLoaded()) 
		return;

	std::cout << " nDetDetector: Loaded GDML model (name=" << solid->getName() << ") with size x=" << solid->getWidth() << " mm, y=" << solid->getThickness() << " mm, z=" << solid->getLength() << " mm\n";
	
	// Place loaded model into the assembly.
	solid->placeSolid(assembly_logV, checkOverlaps);
}

void nDetDetector::loadLightGuide(gdmlSolid *solid, const G4ThreeVector &rotation){
	if(!solid || !solid->isLoaded()) 
		return;

	// Set internal reflectors
	solid->setLogicalBorders("InnerWrapping", materials->fEsrOpSurf);

	G4double trapezoidZ = offsetZ + solid->getLength()/2;
	std::cout << " nDetConstruction: Loaded GDML model (name=" << solid->getName() << ") with size x=" << solid->getWidth() << " mm, y=" << solid->getThickness() << " mm, z=" << solid->getLength() << " mm\n";

	// Place loaded model into the assembly.
	// Place the light-guide on the positive z side.
	solid->setPosition(G4ThreeVector(0, 0, trapezoidZ));
	solid->placeSolid(assembly_logV, checkOverlaps);
	
	// And on the negative z side.
	G4RotationMatrix *trapRot = new G4RotationMatrix();
	trapRot->rotateX(rotation.getX()-CLHEP::pi);
	solid->placeSolid(trapRot, G4ThreeVector(0, 0, -trapezoidZ), assembly_logV, checkOverlaps);

	layerSizeX = solid->getWidth();
	layerSizeY = solid->getThickness();
	offsetZ += solid->getLength();
	fTrapezoidLength = solid->getLength()*mm;
}

///////////////////////////////////////////////////////////////////////////////
// class nDetImplant
///////////////////////////////////////////////////////////////////////////////

nDetImplant::nDetImplant(nDetConstruction *detector, nDetMaterials *matptr) : nDetDetectorParams(detector->GetDetectorParameters()),
				assembly_logV(NULL), assembly_physV(NULL),
				layerSizeX(0), layerSizeY(0), offsetZ(0),
				parentCopyNum(0), firstSegmentCopyNum(0), lastSegmentCopyNum(0), 
				checkOverlaps(true)
{
	copyCenterOfMass(*detector->GetCenterOfMass()); //This has been adjusted
	materials = matptr;
	sipmLogic = detector->GetDetectorParameters().GetReconLogic();
}

nDetImplant::~nDetImplant(){
	/*for(std::vector<nDetWorldObject*>::iterator iter = userLayers.begin(); iter != userLayers.end(); iter++){
		delete (*iter); // This causes a seg-fault
	}
	userLayers.clear();*/
}

nDetImplant nDetImplant::clone() const {
	// Clone this detector's parameters
	nDetImplant retval(*this);
	
	// Copy the center of mass calculators explicitly
	retval.cmI = cmI.clone();

	return retval;
}

void nDetImplant::getCurrentOffset(G4double &x_, G4double &y_, G4double &z_){
	x_ = layerSizeX;
	y_ = layerSizeY;
	z_ = offsetZ;
}

void nDetImplant::setPositionAndRotation(const G4ThreeVector &pos, const G4RotationMatrix &rot){
	detectorPosition = pos;
	detectorRotation = rot;
}

void nDetImplant::setCurrentOffset(const G4double &x_, const G4double &y_, const G4double &z_){
	layerSizeX = x_;
	layerSizeY = y_;
	offsetZ = z_;
}

void nDetImplant::buildAllLayers(){
	for(std::vector<nDetWorldObject*>::iterator iter = userLayers.begin(); iter != userLayers.end(); iter++){
		(*iter)->construct(this);
	}	
}

void nDetImplant::buildBox() {
	G4Box *outerEdge = new G4Box("outerEdge", assemblyWidth/2, assemblyHeight/2, fDetectorLength/2+boxGap+boxThickness);
	G4Box *innerEdge = new G4Box("innerEdge", assemblyWidth/2-boxThickness, assemblyHeight/2-boxThickness, fDetectorLength/2+boxGap);

	G4SubtractionSolid *boxBody = new G4SubtractionSolid("box",outerEdge,innerEdge);
	G4LogicalVolume *box_logV = new G4LogicalVolume(boxBody,boxMaterial,"box_logV");

	/*
	// upper and lower edges of box
	G4Box *upperSide = new G4Box("upperSide", assemblyWidth/2, boxThickness/2, assemblyLength/2);
	G4Box *lowerSide = new G4Box("lowerSide", assemblyWidth/2, boxThickness/2, assemblyLength/2);
	G4LogicalVolume *boxUp_logV = new G4LogicalVolume(upperSide, boxMaterial, "boxUp_logV");
	G4LogicalVolume *boxLow_logV = new G4LogicalVolume(lowerSide, boxMaterial, "boxLow_logV");
	boxUp_logV->SetVisAttributes(materials->visWrapping);
	boxLow_logV->SetVisAttributes(materials->visWrapping);
	G4VPhysicalVolume *boxUp_phys = new G4PVPlacement(0,G4ThreeVector(0,(assemblyHeight-boxThickness)/2,0),boxUp_logV,"upper box edge",assembly_logV,false, 0, checkOverlaps);
	G4VPhysicalVolume *boxLow_phys = new G4PVPlacement(0,G4ThreeVector(0,-(assemblyHeight-boxThickness)/2,0),boxLow_logV,"lower box edge",assembly_logV,false, 0, checkOverlaps);

	// sides of box
	G4Box *leftSide = new G4Box("leftSide", boxThickness/2, assemblyHeight/2-boxThickness, assemblyLength/2);
	G4Box *rightSide = new G4Box("rightSide", boxThickness/2, assemblyHeight/2-boxThickness, assemblyLength/2);
	G4LogicalVolume *boxLeft_logV = new G4LogicalVolume(leftSide, boxMaterial, "boxLeft_logV");
	G4LogicalVolume *boxRight_logV = new G4LogicalVolume(rightSide, boxMaterial, "boxRight_logV");
	boxLeft_logV->SetVisAttributes(materials->visWrapping);
	boxRight_logV->SetVisAttributes(materials->visWrapping);
	G4VPhysicalVolume *boxLeft_phys = new G4PVPlacement(0,G4ThreeVector(-(assemblyWidth-boxThickness)/2,0,0),boxLeft_logV,"left box edge",assembly_logV,false, 0, checkOverlaps);
	G4VPhysicalVolume *boxRight_phys = new G4PVPlacement(0,G4ThreeVector((assemblyWidth-boxThickness)/2,0,0),boxRight_logV,"right box edge",assembly_logV,false, 0, checkOverlaps);
	*/

	addToDetectorBody(box_logV,"implant box");
}

void nDetImplant::placeImplant(G4LogicalVolume *parent){
	assembly_physV = new G4PVPlacement(&detectorRotation, detectorPosition, assembly_logV, "Assembly", parent, 0, 0, false);
	assembly_physV->SetCopyNo(parentCopyNum);
	
	// Do some post-placement processing
	afterPlacement();
}

void nDetImplant::clear(){
	cmI.clear();
}

void nDetImplant::copyCenterOfMass(const centerOfMass &implant){
	cmI = implant.clone();
}

bool nDetImplant::checkPmtCopyNumber(const G4int &num, bool &isImp) const { 
	isImp = 1;
	return (num == 2*parentCopyNum);
}

bool nDetImplant::getSegmentFromCopyNum(const G4int &copyNum, G4int &col, G4int &row) const {
	if(!this->checkCopyNumber(copyNum)) return false;
	col = (copyNum-firstSegmentCopyNum) / fNumRows;
	row = (copyNum-firstSegmentCopyNum) % fNumRows;
	return true;
}

void nDetImplant::addLightGuideGDML(const G4String &input){
	addLayer(new gdmlLightGuideLayer(input));
}

void nDetImplant::addGreaseLayer(const G4String &input){
	addLayer(new greaseLayer(input));
}

void nDetImplant::addDiffuserLayer(const G4String &input){
	addLayer(new diffuserLayer(input));
}

void nDetImplant::addLightGuideLayer(const G4String &input){
	addLayer(new lightGuideLayer(input));
}

void nDetImplant::addSegmentedLightGuideLayer(const G4String &input){
	addLayer(new segLightGuideLayer(input));
}

void nDetImplant::addPhoswichLayer(const G4String &input) {
	addLayer(new phoswichLayer(input));
}

// MAY ALSO NEED TO ADD FRONT AND BACK
void nDetImplant::addBox(const G4String &input) {
	// Expects a space-delimited string of the form:
	//  "addBox <material name> <thickness> <gap size>"
	// (dimensions are in mm)
	std::vector<std::string> args;
	unsigned int Nargs = split_str(input, args);
	boxMaterial = materials->getMaterial(args.at(0));
	boxThickness = strtod(args.at(1).c_str(),NULL)*mm;
	boxGap = strtod(args.at(2).c_str(),NULL)*mm;

	boxAdded = true;
}

void nDetImplant::construct(){
	// Update the size of the assembly in the event it has changed
	UpdateSize(); 

	// Prepare to build the detector and compute the maximum size of the detector volume
	prepareToBuild();

	// Build the assembly volume
	constructAssembly();

	// Get materials and visual attributes
	scintMaterial = materials->getUserDetectorMaterial(detectorMaterialName);
	scintVisAtt = materials->getUserVisAttributes(detectorMaterialName);
	wrappingMaterial = materials->getUserSurfaceMaterial(wrappingMaterialName);
	outerMylar = materials->fMylar;
	wrappingVisAtt = materials->getUserVisAttributes(wrappingMaterialName);
	windowVisAtt = materials->visWindow;
	wrappingOpSurf = materials->getUserOpticalSurface(wrappingMaterialName);
	wrappingOuterOpSurf = materials->fMylarOpSurf;

	// build box first if one has been added
	if (boxAdded) {
		buildBox();
	}

	// Build the geometry
	buildDetector();

	// Generate all user-defined layers.
	buildAllLayers();

	// Attach PMT 
	constructPSPmt();
}

G4LogicalVolume *nDetImplant::constructAssembly(){

	// Calculate the dimensions of the detector
	assemblyWidth = maxBodySize.getX() + 2*fWrappingThickness;
	assemblyHeight = maxBodySize.getY() + 2*fWrappingThickness;
	assemblyLength = maxBodySize.getZ();

	// Account for the size of the PSPMT
	assemblyWidth = std::max(assemblyWidth, pmtWidth);
	assemblyHeight = std::max(assemblyHeight, pmtHeight);
	assemblyLength += (2*fGreaseThickness+fWindowThickness+fSensitiveThickness+fWrappingThickness*3/2);
	std::cout<<"assemblyLength:: "<<assemblyLength<<std::endl;
	assemblyLength += (fWrappingThickness>0.1*cm) ? 0 : assemblyLength+2*dzThick;
	// assemblyLength += 2*(2*fGreaseThickness+fWindowThickness+fSensitiveThickness+dzThick);

	// Account for the additional component layers
	for(std::vector<nDetWorldObject*>::iterator iter = userLayers.begin(); iter != userLayers.end(); iter++){
		if(!(*iter)->decodeString()){
			std::cout << " nDetImplant: Invalid number of arguments given to ::decodeString(). Expected " << (*iter)->getNumRequiredArgs() << " but received " << (*iter)->getNumSuppliedArgs() << ".\n";
			std::cout << " nDetImplant:  SYNTAX: " << (*iter)->syntaxStr() << std::endl;
			continue;
		}
		assemblyWidth = std::max(assemblyWidth, (*iter)->getSizeX());
		assemblyHeight = std::max(assemblyHeight, (*iter)->getSizeY());
		assemblyLength += (*iter)->getSizeZ();
	}

	// Account for the a box around the implant if one has been added
	if (boxAdded) {
		assemblyWidth = std::max(assemblyWidth,fDetectorWidth+2*(fWrappingThickness+boxGap+boxThickness));
		assemblyHeight = std::max(assemblyHeight,fDetectorHeight+2*(fWrappingThickness+boxGap+boxThickness));
	}

	// Build the assembly box
	G4Box *assembly = new G4Box("assembly", assemblyWidth/2, assemblyHeight/2, assemblyLength/2);
	assembly_logV = new G4LogicalVolume(assembly, materials->fAir, "assembly_logV");
	assembly_logV->SetVisAttributes(materials->visAssembly);

	return assembly_logV;
	
}

G4PVPlacement *nDetImplant::addToDetectorBody(G4LogicalVolume *volume, const G4String &name/*=""*/, const G4ThreeVector &pos/*=G4ThreeVector(0,0,0)*/, G4RotationMatrix *rot/*=NULL*/){
	return (new G4PVPlacement(rot, pos, volume, (name.empty() ? volume->GetName() : name), assembly_logV, false, 0, checkOverlaps));
}

G4PVPlacement *nDetImplant::addSegmentToBody(G4LogicalVolume *volume, const G4String &name/*=""*/, const G4ThreeVector &pos/*=G4ThreeVector(0,0,0)*/, G4RotationMatrix *rot/*=NULL*/){
	return (new G4PVPlacement(rot, pos, volume, (name.empty() ? volume->GetName() : name), assembly_logV, false, lastSegmentCopyNum++, checkOverlaps));
}

G4PVPlacement *nDetImplant::addBackComponent(G4LogicalVolume *volume, const G4double &offset, const G4String &name/*=""*/, G4RotationMatrix *rot/*=NULL*/){ 
	return (new G4PVPlacement(rot, G4ThreeVector(0, 0, offset), volume, (name.empty() ? volume->GetName() : name), assembly_logV, false, 2*parentCopyNum, checkOverlaps));
}

void nDetImplant::addPmtWithGaps(const G4double &offset, const G4String &name/*=""*/, G4RotationMatrix *rot/*=NULL*/){ 

	std::vector< std::vector<G4LogicalVolume*> > volumes;
	std::vector< std::vector<G4PVPlacement*> > placements;

	std::vector<G4LogicalVolume*> initializer1(fNumColumnsPmt, NULL);
	std::vector<G4PVPlacement*> initializer2(fNumColumnsPmt, NULL);

	G4CSGSolid* sensitive_solidV = getVolume(name+"_solidV", pmtPixelWidth, pmtPixelHeight, fSensitiveThickness);

	// create logical volumes
	for (G4int i = 0; i < fNumRowsPmt; i++) {
		volumes.push_back(initializer1);
		placements.push_back(initializer2);
		for (G4int j = 0; j < fNumColumnsPmt; j++) {
		    volumes[i][j] = new G4LogicalVolume(sensitive_solidV, materials->fSilicon, name+"_logV");
			volumes[i][j]->SetVisAttributes(materials->visSensitive);

			new G4LogicalSkinSurface(name, volumes[i][j], materials->fSiliconOpSurf);
		}
	}
	
	G4double segmentX, segmentY;
	segmentY = -pmtHeight/2 + pmtPixelHeight/2;

	// create placements for each of the Pmt segments
	for (G4int i = 0; i < fNumRowsPmt; i++) {
		segmentX = -pmtWidth/2 + pmtPixelWidth/2;

		for (G4int j = 0; j < fNumColumnsPmt; j++) {
			placements[i][j] =
				new G4PVPlacement(rot, G4ThreeVector(segmentX, segmentY, offset), volumes[i][j], (name.empty() ? volumes[i][j]->GetName() : name), assembly_logV, false, 2*parentCopyNum, checkOverlaps);
			
			segmentX = segmentX + pmtPixelWidth + pmtGapThickness;
		}

		segmentY = segmentY + pmtPixelHeight + pmtGapThickness;
	}
}


void nDetImplant::addBackComponent(G4PVPlacement* &phys1, G4LogicalVolume *volume, const G4double &offset, const G4String &name/*=""*/, G4RotationMatrix *rot/*=NULL*/){
	phys1 = (new G4PVPlacement(rot, G4ThreeVector(0, 0, offset), volume, (name.empty() ? volume->GetName() : name), assembly_logV, false, 2*parentCopyNum, checkOverlaps));
} 

G4PVPlacement *nDetImplant::addFrontComponent(G4LogicalVolume *volume, const G4double &offset, const G4String &name/*=""*/, G4RotationMatrix *rot/*=NULL*/){ 
	return (new G4PVPlacement(rot, G4ThreeVector(0, 0, offset), volume, (name.empty() ? volume->GetName() : name), assembly_logV, false, 2*parentCopyNum, checkOverlaps));
}

void nDetImplant::addFrontComponent(G4PVPlacement* &phys1, G4LogicalVolume *volume, const G4double &offset, const G4String &name/*=""*/, G4RotationMatrix *rot/*=NULL*/){
	phys1 = (new G4PVPlacement(rot, G4ThreeVector(0, 0, offset), volume, (name.empty() ? volume->GetName() : name), assembly_logV, false, 2*parentCopyNum, checkOverlaps));
} 

void nDetImplant::prepareToBuild(){ 
	maxBodySize = G4ThreeVector(fDetectorWidth, fDetectorHeight, fDetectorLength); 
}

void nDetImplant::constructPSPmt(){

	// Build the sensitive PMT surface.
	const G4String name = "psSiPM";

	G4double sensitiveZ = offsetZ + fGreaseThickness + fWindowThickness + fSensitiveThickness/2;
	G4double wrappingThickness = fGreaseThickness + fWindowThickness;

	// The optical grease layer.
	G4PVPlacement *greasePhys = NULL;
	if(fGreaseThickness > 0){
		G4double greaseZ = offsetZ + fGreaseThickness/2;
	
		G4CSGSolid *grease_solidV = getVolume("window_solidV", pmtWidth, pmtHeight, fGreaseThickness);
		G4LogicalVolume *grease_logV = new G4LogicalVolume(grease_solidV, materials->fGrease, "grease_logV");
		
		grease_logV->SetVisAttributes(materials->visGrease);

		greasePhys = addBackComponent(grease_logV, greaseZ, "Grease1"); 
		
		if(!fPolishedInterface){
			for(std::vector<G4PVPlacement*>::iterator iter = scintBody_physV.begin(); iter != scintBody_physV.end(); iter++){
				new G4LogicalBorderSurface("GreaseInterface", (*iter), greasePhys, materials->fGreaseOpSurf);
			}
		}
		
		// Clear all scintillator placements.
		scintBody_physV.clear();
	}

	G4PVPlacement *windowPhys = NULL;
	if(fWindowThickness > 0){ // The quartz window
		G4double windowZ = offsetZ + fGreaseThickness + fWindowThickness/2;
	
		G4CSGSolid *window_solidV = getVolume("window_solidV", pmtWidth, pmtHeight, fWindowThickness);
		G4LogicalVolume *window_logV = new G4LogicalVolume(window_solidV, materials->fSiO2, "window_logV");
		
		window_logV->SetVisAttributes(materials->visWindow);

		windowPhys = addBackComponent(window_logV, windowZ, "Quartz");
	}

	// Build the wrapping.
	if(WrappingEnabled() && wrappingThickness > 0){
		G4CSGSolid *boundingBox = getVolume("", pmtWidth, pmtHeight, wrappingThickness);
		G4CSGSolid *wrappingBox = getVolume("", pmtWidth + 2*fWrappingThickness, pmtHeight + 2*fWrappingThickness, wrappingThickness);
		
		G4SubtractionSolid *greaseWrapping = new G4SubtractionSolid("", wrappingBox, boundingBox);
		G4LogicalVolume *greaseWrapping_logV = new G4LogicalVolume(greaseWrapping, wrappingMaterial, "greaseWrapping_logV");
		greaseWrapping_logV->SetVisAttributes(wrappingVisAtt);
		
		G4double wrappingZ = offsetZ + fGreaseThickness/2 + fWindowThickness/2;

		// Place the wrapping around the scintillator.
		G4PVPlacement *greaseWrapping_physV;
		
		if(greasePhys){
			new G4LogicalBorderSurface("Wrapping", greasePhys, greaseWrapping_physV, wrappingOpSurf);
		}
		if(windowPhys){
			new G4LogicalBorderSurface("Wrapping", windowPhys, greaseWrapping_physV, wrappingOpSurf);
		}
	}

	if (pmtHasGaps) {
		addPmtWithGaps(sensitiveZ, name, NULL);

		// Move the current offset past the PMT
		layerSizeX = pmtWidth;
		layerSizeY = pmtHeight;
  		offsetZ += fGreaseThickness + fWindowThickness + fSensitiveThickness;
	}
	else {

		// The photon sensitive surface
		G4CSGSolid *sensitive_solidV = getVolume(name+"_solidV", pmtWidth, pmtHeight, fSensitiveThickness);
		G4LogicalVolume *sensitive_logV = new G4LogicalVolume(sensitive_solidV, materials->fSilicon, name+"_logV");
		sensitive_logV->SetVisAttributes(materials->visSensitive);
	
		// Logical skin surface.
		new G4LogicalSkinSurface(name, sensitive_logV, materials->fSiliconOpSurf);    

		//addMirroredComponents(sensitive_logV, sensitiveZ, name);
		addBackComponent(sensitive_logV, sensitiveZ, name);

		// Move the current offset past the PMT
		layerSizeX = pmtWidth;
		layerSizeY = pmtHeight;
  		offsetZ += fGreaseThickness + fWindowThickness + fSensitiveThickness;
	}
}

G4CSGSolid *nDetImplant::getVolume(const G4String &name, const G4double &width, const G4double &height, const G4double &length){
	G4CSGSolid *retval;
	if(fSquarePMTs)
		retval = new G4Box(name, width/2, height/2, length/2);
	else
		retval = new G4Tubs(name, 0, width/2, length/2, 0, 2*CLHEP::pi);
	return retval;
}

G4CSGSolid *nDetImplant::getLightGuideVolume(const G4String &name, const G4double &w1, const double &w2, const double &h1, const double &h2, const G4double &length){
	G4CSGSolid *retval;
	if(fSquarePMTs)
		retval = new G4Trd(name, w1/2, w2/2, h1/2, h2/2, length/2);
	else
		retval = new G4Cons(name, 0, w1/2, 0, w2/2, length/2, 0, 2*CLHEP::pi);
	return retval;
}

void nDetImplant::applyGreaseLayer(){
	std::cout<<"Layer size x "<<layerSizeX;
	this->applyGreaseLayer(layerSizeX, layerSizeY);
}

void nDetImplant::applyGreaseLayer(const G4double &x, const G4double &y, double thickness/*=0*/, const G4String& material/*=""*/){
	if(thickness <= 0)
		thickness = fGreaseThickness;
	if(thickness > 0){
		G4CSGSolid *grease_solidV = getVolume("grease", x, y, thickness);
		if (material != "") {
			// determine which grease material the user has specified; default to regular optical grease
			G4Material* greaseMaterial;
			if (material == "noa61") {
				greaseMaterial = materials->fNOA61;
			}
			else if (material == "noa68") {
				greaseMaterial = materials->fNOA68;
			}
			else if (material == "noa136") {
				greaseMaterial = materials->fNOA136;
			}
			else if (material == "noa170") {
				greaseMaterial = materials->fNOA170;
			}
			else {
				std::cout << "\"" << material <<"\n\" is not a valid grease option; defaulting to optical grease\n";
				greaseMaterial = materials->fGrease;
			}
			G4LogicalVolume *grease_logV = new G4LogicalVolume(grease_solidV, greaseMaterial, "grease_logV");
			grease_logV->SetVisAttributes(materials->visGrease);

			// Add the optical grease to the assembly
			//addMirroredComponents(grease_logV, offsetZ+thickness/2, "Grease"); 
			G4PVPlacement *greaseLogical = addBackComponent(grease_logV, offsetZ+thickness/2, "Grease0"); 
		}
		else {
			// use optical grease if material not specified
			G4LogicalVolume *grease_logV = new G4LogicalVolume(grease_solidV, materials->fGrease, "grease_logV");
			grease_logV->SetVisAttributes(materials->visGrease);

			// Add the optical grease to the assembly
			//addMirroredComponents(grease_logV, offsetZ+thickness/2, "Grease"); 
			G4PVPlacement *greaseLogical = addBackComponent(grease_logV, offsetZ+thickness/2, "Grease0"); 
		}

		// Offset the other layers to account for the layer of optical grease
		layerSizeX = x;
		layerSizeY = y;
		offsetZ += thickness;
	}
}

void nDetImplant::applyDiffuserLayer(){
	this->applyDiffuserLayer(layerSizeX, layerSizeY, fDiffuserLength);
}

void nDetImplant::applyDiffuserLayer(const G4double &x, const G4double &y, const double &thickness){
	if(thickness > 0){ // Build the light diffusers (if needed)
		G4CSGSolid *lightDiffuser = getVolume("lightDiffuser", x, y, thickness);
		G4LogicalVolume *lightDiffuserLog = new G4LogicalVolume(lightDiffuser, materials->fSiO2, "lightDiffuser_logV");

		// Add the optical grease to the assembly
		//addMirroredComponents(lightDiffuserLog, offsetZ+thickness/2, "Diffuser");
		G4PVPlacement *lightDiffuserLogical = addBackComponent(lightDiffuserLog, offsetZ+thickness/2, "Diffuser"); 

		// Offset the other layers to account for the light-diffuser
		layerSizeX = x;
		layerSizeY = y;
		offsetZ += thickness;
	}
}

void nDetImplant::applyLightGuide(){
	this->applyLightGuide(layerSizeX, pmtWidth, layerSizeY, pmtHeight, fTrapezoidLength);
}

void nDetImplant::applyLightGuide(const G4double &x2, const G4double &y2){
	this->applyLightGuide(layerSizeX, x2, layerSizeY, y2, fTrapezoidLength);
}

void nDetImplant::applyLightGuide(const G4double &x1, const G4double &x2, const G4double &y1, const G4double &y2, const double &thickness){
    if(thickness > 0){ // Build the light guides (if needed)
        const G4double trapAngleXZ = std::atan2(2*thickness, x1-x2);
        const G4double trapAngleYZ = std::atan2(2*thickness, y1-y2);
        
        const G4double deltaX = fWrappingThickness/std::sin(trapAngleXZ);
        const G4double deltaY = fWrappingThickness/std::sin(trapAngleYZ);
        
        G4double trapezoidZ = offsetZ + thickness/2;
        
        std::string trapName = "Acrylic";
        
        // Build the light-guide.
        G4CSGSolid *lightGuide = getLightGuideVolume("lightGuide", x1, x2, y1, y2, thickness);
	    G4LogicalVolume *lightGuideLog = new G4LogicalVolume(lightGuide, materials->fSiO2, "lightGuide_logV");	

		G4RotationMatrix *rightRotation = new G4RotationMatrix();
		rightRotation->rotateX(CLHEP::pi);

		// Place the light guides.
		G4PVPlacement *trapPhysical = addBackComponent(lightGuideLog, trapezoidZ, trapName);

		// Build the wrapping.
		if(WrappingEnabled()){
			G4CSGSolid *wrappingSolid = getLightGuideVolume("wrapping", x1+2*deltaX, x2+2*deltaY, y1+2*deltaX, y2+2*deltaY, thickness);
			G4SubtractionSolid *wrapping = new G4SubtractionSolid("wrapping", wrappingSolid, lightGuide);
			G4LogicalVolume *wrapping_logV = new G4LogicalVolume(wrapping, wrappingMaterial, "wrapping_logV");		
			wrapping_logV->SetVisAttributes(wrappingVisAtt);

			// Place the wrapping around the light guides.
			G4PVPlacement *trapWrapping = addBackComponent(wrapping_logV, trapezoidZ, trapName);
		
			// Reflective wrapping.
			new G4LogicalBorderSurface("Wrapping", trapPhysical, trapWrapping, wrappingOpSurf);
		}
		
		// Offset the other layers to account for the light-guide
		layerSizeX = x2;
		layerSizeY = y2;
    offsetZ += thickness;
	}
}

void nDetImplant::applySegmentedLightGuide(G4int Xseg, G4int Yseg, G4double spacing, G4double topWidth, G4double topThick, G4double botWidth, G4double botThick, G4double zThick){
	G4int numSegs = Xseg*Yseg;
	std::string name = "segLightGuide";
	auto lg_mat = materials->fAcrylic;
	auto seg_solid = new G4Trap("trap");
	//auto seg_solid = new G4Box("seg_solid",2*mm,2*mm,10*mm);
	auto lg_box = new G4Box("lg_box",topWidth/2.,topThick/2.,zThick/2.);
	//auto lg_box = new G4Box("lg_box",fDetectorWidth/2.,fDetectorHeight/2.,zThick/2.);
	auto lg_log = new G4LogicalVolume(lg_box,materials->fAir,"lg_log");

	//Attempt to create a solid for the wrapping
	G4double cellXt = ( topWidth-  fWrappingThickness*( Xseg-1))/ Xseg;
	G4double cellXb = ( botWidth-  fWrappingThickness*( Xseg-1))/ Xseg;
	G4double cellYt = ( topThick-  fWrappingThickness*( Yseg-1))/ Yseg;
	G4double cellYb = ( botThick-  fWrappingThickness*( Yseg-1))/ Yseg;
	std::vector<G4SubtractionSolid*> wrappings_vec;
	for(int copyNo=0;copyNo<numSegs;copyNo++){
		G4double xposb = (copyNo% Xseg- Xseg/2)*(cellXb);
		G4double yposb = (copyNo/ Yseg- Yseg/2)*(cellYb);
		G4double xpost = (copyNo% Xseg- Xseg/2)*(cellXt);
		G4double ypost = (copyNo/ Yseg- Yseg/2)*(cellYt);
		G4double zpos = 0.;
		G4double pTheta;
		if(xposb<0)
			pTheta = -1*atan((sqrt(pow(xposb,2)+pow(yposb,2))-sqrt(pow(xpost,2)+pow(ypost,2)))/( zThick));
		else
			pTheta = atan((sqrt(pow(xposb,2)+pow(yposb,2))-sqrt(pow(xpost,2)+pow(ypost,2)))/( zThick));
		G4double pPhi;
		if(xposb==0)
			pPhi = -1*atan((yposb-ypost+0.0000001)/(xposb-xpost));
		else
			pPhi = atan((yposb-ypost+0.0000001)/(xposb-xpost));
		auto seg_solid = new G4Trap("trap",
					zThick/2., // pDz Half z length
					pTheta, //pTheta Polar angle of the line joining the centers of the faces at -/+ pDz
					pPhi, //pPhi Azimuthal angle of the line joining the center of the face at -pDz to the center of the face at +pDz
					cellYt/2., //pDy1 Half y length at -pDz
					cellXt/2., //pDx1 Half x length of the side at y=-pDy1 of the face at -pDz
					cellXt/2., //pDx2 Half x length of the side at y=+pDy1 of the face at -pDz
					0*degree, //pAlp1 Angle with respect to the y-axis from the center of the side (lower endcap)
					cellYb/2., //pDy2 Half y length at +pDz
					cellXb/2., //pDx3 Half x length of the side at y=-pDy2 of the face at -pDz
					cellXb/2., //pDx4 Half x length of the side at y=+pDy2 of the face at -pDz
					0*degree //pAlp2 Angle with respect to the y-axis from the center of the side (upper endcap)
					);
		G4double xpos = -0.25*(botWidth+topWidth-(cellXb+cellXt))+(copyNo%Xseg)*(0.5*(cellXb+cellXt)+ fWrappingThickness);
		G4double ypos = -0.25*(botThick+topThick-(cellYb+cellYt))+(copyNo/Yseg)*(0.5*(cellYb+cellYt)+ fWrappingThickness);
		/*auto box_solid = new G4Box("box_solid",2*mm,2*mm,10*mm);
		auto box_log = new G4LogicalVolume(box_solid,lg_mat,"box_log");*/
		auto seg_log = new G4LogicalVolume(seg_solid,lg_mat,"seg_log");
		G4PVPlacement *seg_place = new G4PVPlacement(
				nullptr,
				G4ThreeVector(xpos,ypos,0),
				seg_log,
				"seg_phys",
				lg_log,
				false,
				0,
				true
			);
		if(WrappingEnabled()){
			auto outer = new G4Trap("outer",zThick/2.,pTheta,pPhi,cellYt/2.+fWrappingThickness/2.,cellXt/2.+fWrappingThickness/2.,cellXt/2.+fWrappingThickness/2.,0*degree,cellYb/2.+fWrappingThickness/2.,cellXb/2.+fWrappingThickness/2.,cellXb/2.+fWrappingThickness/2.,0*degree);
			auto wrapping = new G4SubtractionSolid("wrapping",outer,seg_solid);
			//G4LogicalVolume *wrapping_logV = new G4LogicalVolume(wrapping, wrappingMaterial, "wrapping_logV");
			G4LogicalVolume *wrapping_logV = new G4LogicalVolume(wrapping,materials->fMylar, "wrapping_logV");		
			wrapping_logV->SetVisAttributes(wrappingVisAtt);
			G4PVPlacement *wrap_place = new G4PVPlacement(nullptr,G4ThreeVector(xpos,ypos,0),wrapping_logV,"wrap_phys",lg_log,false,0,true);
			new G4LogicalBorderSurface("Wrapping", seg_place, wrap_place, wrappingOuterOpSurf);
		}
	}
	G4PVPlacement *lightGuideBox = addBackComponent(lg_log,offsetZ+zThick/2.,name);
	offsetZ+=zThick;
	return;
}

void nDetImplant::applyPhoswich(G4int fXseg, G4int fYseg, G4double fThick, G4double fWidth, G4double fHeight, G4double fWrapping, G4String fSegMaterial, G4String fWrapMat) {

	// RECALCULATE THE SEGMENT DIMENSIONS
	const G4double cellWidth = (fWidth-(fXseg-1)*fWrapping)/fXseg;
	const G4double cellHeight = (fHeight-(fYseg-1)*fWrapping)/fYseg;
	
	int Ncol = fXseg;
	int Nrow = fYseg;

    G4Material* phoswichMaterial = materials->getUserDetectorMaterial(fSegMaterial);
	G4VisAttributes* phoswichVisAtt = materials->getUserVisAttributes(fSegMaterial);
	outerMylar = materials->fMylar;

    G4Material* wrapMat;
	G4VisAttributes* wrapVisAtt;
	G4OpticalSurface* wrapOpSurf;

	if (fWrapMat != "") {
		wrapMat = materials->getUserSurfaceMaterial(fWrapMat);
		wrapVisAtt = materials->getUserVisAttributes(fWrapMat);
		wrapOpSurf = materials->getUserOpticalSurface(fWrapMat);
	}
	else {
		wrapMat = materials->getUserSurfaceMaterial(wrappingMaterialName);
		wrappingVisAtt = materials->getUserVisAttributes(wrappingMaterialName);
		wrappingOpSurf = materials->getUserOpticalSurface(wrappingMaterialName);
	}


    // Construct the scintillator cell
    G4Box *cellScint = new G4Box("phoswich", cellWidth/2, cellHeight/2, fThick/2);
    G4LogicalVolume *cellScint_logV = new G4LogicalVolume(cellScint, phoswichMaterial, "phoswich_log");
    cellScint_logV->SetVisAttributes(phoswichVisAtt);

	G4Box *mylarVertLayer = NULL;
	G4Box *mylarHorizLayer = NULL;
	
	G4LogicalVolume *mylarVertLayer_logV = NULL;
	G4LogicalVolume *mylarHorizLayer_logV = NULL;

	// Build the wrapping.
	G4PVPlacement *wrapping_physV = NULL;
	G4PVPlacement* wrappingFacePhys = NULL;
	if(fWrapping > 0){
		// Construct the outer wrapping.
		G4Box *wrappingBox = new G4Box("wrappingBox", fWidth/2+fWrapping, fHeight/2+fWrapping, fThick/2);
		G4Box *scintBox = new G4Box("phoswichBox", fWidth/2, fHeight/2, fThick/2);
		
		G4SubtractionSolid *wrappingBody = new G4SubtractionSolid("wrapping", wrappingBox, scintBox);
		G4LogicalVolume *wrapping_logV = new G4LogicalVolume(wrappingBody, outerMylar, "wrapping_logV");
		wrapping_logV->SetVisAttributes(wrappingVisAtt);
	
		// Place the outer wrapping into the assembly.
		G4ThreeVector wrappingPos(0,0,(offsetZ + fThick/2));
		wrapping_physV = addToDetectorBody(wrapping_logV, "Wrapping", wrappingPos);

		// DON'T NEED FRONT WRAPPING FOR PHOSWICH
		//auto wrappingFace = new G4Box("wrappingFace",fDetectorWidth/2+fWrappingThickness, fDetectorHeight/2+fWrappingThickness,fWrappingThickness/4.); 
		//auto wrappingFace_log = new G4LogicalVolume(wrappingFace,outerMylar,"wrappingFace_log");
		//wrappingFace_log->SetVisAttributes(windowVisAtt);	
		//wrappingFacePhys = addFrontComponent(wrappingFace_log,-1*fDetectorLength/2-fWrappingThickness/4.,"wrappingFacePhys");
		
		// Construct vertical and horizontal reflector layers for later use.
		mylarVertLayer = new G4Box("mylarVertLayer", fWrapping/2, fHeight/2, fThick/2);
		mylarHorizLayer = new G4Box("mylarHorizLayer", cellWidth/2, fWrapping/2, fThick/2);

		mylarVertLayer_logV = new G4LogicalVolume(mylarVertLayer, wrapMat, "mylarVertLayer_logV");
		mylarHorizLayer_logV = new G4LogicalVolume(mylarHorizLayer, wrapMat, "mylarHorizLayer_logV");
		
		mylarVertLayer_logV->SetVisAttributes(wrappingVisAtt);
		mylarHorizLayer_logV->SetVisAttributes(wrappingVisAtt);
	}

	// Place the scintillator segments into the assembly.
	std::vector<G4PVPlacement*> mylarVertLayer_physV(Ncol, NULL);
	std::vector<std::vector<G4PVPlacement*> > mylarHorizLayer_physV(Ncol, std::vector<G4PVPlacement*>(Nrow, NULL));
	std::vector<std::vector<G4PVPlacement*> > cellScint_physV(Ncol, std::vector<G4PVPlacement*>(Nrow, NULL));	
	for(int col = 0; col < Ncol; col++){
		for(int row = 0; row < Nrow; row++){
			G4ThreeVector cellCenter(-fWidth/2 + col*fWrapping + (col+0.5)*cellWidth, -fHeight/2 + row*fWrapping + (row+0.5)*cellHeight, offsetZ + fThick/2);

			// Copy numbers (segment IDs), indexed from 1
			std::stringstream stream; stream << "Phoswich-" << col << "," << row;
			cellScint_physV[col][row] = addSegmentToBody(cellScint_logV, stream.str(), cellCenter);
			scintBody_physV.push_back(cellScint_physV[col][row]);
		
			// Place vertical and horizontal reflectors.
			if(fWrapping > 0){ 
				if(row == 0 && col != Ncol-1){ // New vertical reflector layer.
					std::stringstream stream2; stream2 << "Wrapping-" << col;
					mylarVertLayer_physV[col] = addToDetectorBody(mylarVertLayer_logV, stream2.str().c_str(), G4ThreeVector(cellCenter.getX()+cellWidth/2+fWrapping/2, 0, offsetZ + fThick/2));
				}
				if(row != Nrow-1){ // New horizontal reflector layer.
					std::stringstream stream2; stream2 << "Wrapping-" << col << "," << row;
					mylarHorizLayer_physV[col][row] = addToDetectorBody(mylarHorizLayer_logV, stream2.str().c_str(), G4ThreeVector(cellCenter.getX(), cellCenter.getY()+cellHeight/2+fWrapping/2, offsetZ + fThick/2));
				}
			}
		}
	}
	offsetZ+=fThick;
	return;
}

void nDetImplant::loadGDML(gdmlSolid *solid){
	if(!solid || !solid->isLoaded()) 
		return;

	std::cout << " nDetImplant: Loaded GDML model (name=" << solid->getName() << ") with size x=" << solid->getWidth() << " mm, y=" << solid->getThickness() << " mm, z=" << solid->getLength() << " mm\n";
	
	// Place loaded model into the assembly.
	solid->placeSolid(assembly_logV, checkOverlaps);
}

void nDetImplant::loadLightGuide(gdmlSolid *solid, const G4ThreeVector &rotation){
	if(!solid || !solid->isLoaded()) 
		return;

	// Set internal reflectors
	solid->setLogicalBorders("InnerWrapping", materials->fEsrOpSurf);

	G4double trapezoidZ = offsetZ + solid->getLength()/2.0;
	std::cout << " nDetConstruction: Loaded GDML model (name=" << solid->getName() << ") with size x=" << solid->getWidth() << " mm, y=" << solid->getThickness() << " mm, z=" << solid->getLength() << " mm\n";

	// Place loaded model into the assembly.
	// Place the light-guide on the positive z side.
	solid->setPosition(G4ThreeVector(0, 0, trapezoidZ));
	solid->placeSolid(assembly_logV, checkOverlaps);
	
	// And on the negative z side.
	G4RotationMatrix *trapRot = new G4RotationMatrix();
	trapRot->rotateX(rotation.getX()-CLHEP::pi);
	solid->placeSolid(trapRot, G4ThreeVector(0, 0, -trapezoidZ), assembly_logV, checkOverlaps);

	layerSizeX = solid->getWidth();
	layerSizeY = solid->getThickness();
	offsetZ += solid->getLength();
	fTrapezoidLength = solid->getLength()*mm;
}

