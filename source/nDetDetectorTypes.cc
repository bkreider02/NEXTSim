/** @file nDetDetectorTypes.cc
  * @author Cory R. Thornsberry (cthornsb@vols.utk.edu
  * @date September 25, 2019
  * @brief Contains all pre-defined detector construction types
  */

#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"

#include "nDetDetectorTypes.hh"

using namespace nDetDetectorTypes;

nDetDetector* nDetDetectorTypes::getDetectorType(const G4String &geom){
	nDetDetector *retval = NULL;
	if(geom == "next" || geom == "module") // Next module (segmented)
		retval = new nextModuleType();
	else if(geom == "ellipse") // Elliptical detector (football)
		retval = new ellipticalType();
	else if(geom == "rectangle") // Rectangular detector
		retval = new rectangularType();
	else if(geom == "cylinder") // Cylindrical detector
		retval = new cylindricalType();
	return retval;
}

nDetImplant* nDetDetectorTypes::getImplantType(const G4String &geom){
	nDetImplant *retval = NULL;
	if(geom == "hagrid")
		retval = new hagridType();
	else if (geom == "dome") // "Dome" implant
		retval = new domeType();
	else
		retval = new implantType();
	return retval;
}

nDetDetector* nDetDetectorTypes::getDetectorType(const G4String &geom, nDetConstruction *construction, nDetMaterials *matptr){
	nDetDetector *retval = NULL;
	if(geom == "next" || geom == "module") // Next module (segmented)
		retval = new nextModuleType(construction, matptr);
	else if(geom == "ellipse") // Elliptical detector (football)
		retval = new ellipticalType(construction, matptr);
	else if(geom == "rectangle") // Rectangular detector
		retval = new rectangularType(construction, matptr);
	else if(geom == "cylinder") // Cylindrical detector
		retval = new cylindricalType(construction, matptr);
	return retval;
}

nDetImplant* nDetDetectorTypes::getImplantType(const G4String &geom, nDetConstruction *construction, nDetMaterials *matptr){
	nDetImplant *retval = NULL;
	if(geom == "hagrid")
		retval = new hagridType(construction, matptr);
	else if (geom == "dome") // "Dome" implant
		retval = new domeType(construction, matptr);
	else
		retval = new implantType(construction, matptr);
	return retval;
}

///////////////////////////////////////////////////////////////////////////////
// class nextModuleType
///////////////////////////////////////////////////////////////////////////////

void nextModuleType::prepareToBuild(){
	// Set the maximum detector body size
	maxBodySize = G4ThreeVector(fDetectorWidth, fDetectorHeight, fDetectorLength);
}

void nextModuleType::buildDetector(){
	const G4double cellWidth = GetSegmentWidth();
	const G4double cellHeight = GetSegmentHeight();
	
	// Get the number of rows and columns for this segmented detector
	int Ncol = fNumColumns;
	int Nrow = fNumRows;

	// Update the detector's copy numbers
	firstSegmentCopyNum = scintCopyNum; 
	lastSegmentCopyNum = firstSegmentCopyNum; //startCopyNum+col*row;

    // Construct the scintillator cell
    G4Box *cellScint = new G4Box("scintillator", cellWidth/2, cellHeight/2, fDetectorLength/2);
    G4LogicalVolume *cellScint_logV = new G4LogicalVolume(cellScint, scintMaterial, "scint_log");
    cellScint_logV->SetVisAttributes(scintVisAtt);

	G4Box *mylarVertLayer = NULL;
	G4Box *mylarHorizLayer = NULL;
	
	G4LogicalVolume *mylarVertLayer_logV = NULL;
	G4LogicalVolume *mylarHorizLayer_logV = NULL;

	// Build the wrapping.
	G4PVPlacement *wrapping_physV = NULL;
	if(WrappingEnabled()){
		// Construct the outer wrapping.
		G4Box *wrappingBox = new G4Box("wrappingBox", fDetectorWidth/2+fWrappingThickness, fDetectorHeight/2+fWrappingThickness, fDetectorLength/2);
		G4Box *scintBox = new G4Box("scintBox", fDetectorWidth/2, fDetectorHeight/2, fDetectorLength/2);
		
		G4SubtractionSolid *wrappingBody = new G4SubtractionSolid("wrapping", wrappingBox, scintBox);
		G4LogicalVolume *wrapping_logV = new G4LogicalVolume(wrappingBody, wrappingMaterial, "wrapping_logV");
		wrapping_logV->SetVisAttributes(wrappingVisAtt);
		
		// Place the outer wrapping into the assembly.
		wrapping_physV = addToDetectorBody(wrapping_logV, "Wrapping");
		
		// Construct vertical and horizontal reflector layers for later use.
		mylarVertLayer = new G4Box("mylarVertLayer", fWrappingThickness/2, fDetectorHeight/2, fDetectorLength/2);
		mylarHorizLayer = new G4Box("mylarHorizLayer", cellWidth/2, fWrappingThickness/2, fDetectorLength/2);

		mylarVertLayer_logV = new G4LogicalVolume(mylarVertLayer, wrappingMaterial, "mylarVertLayer_logV");
		mylarHorizLayer_logV = new G4LogicalVolume(mylarHorizLayer, wrappingMaterial, "mylarHorizLayer_logV");
		
		mylarVertLayer_logV->SetVisAttributes(wrappingVisAtt);
		mylarHorizLayer_logV->SetVisAttributes(wrappingVisAtt);
	}

	// Place the scintillator segments into the assembly.
	std::vector<G4PVPlacement*> mylarVertLayer_physV(Ncol, NULL);
	std::vector<std::vector<G4PVPlacement*> > mylarHorizLayer_physV(Ncol, std::vector<G4PVPlacement*>(Nrow, NULL));
	std::vector<std::vector<G4PVPlacement*> > cellScint_physV(Ncol, std::vector<G4PVPlacement*>(Nrow, NULL));	
	for(int col = 0; col < Ncol; col++){
		for(int row = 0; row < Nrow; row++){
			G4ThreeVector cellCenter(-fDetectorWidth/2 + col*fWrappingThickness + (col+0.5)*cellWidth, -fDetectorHeight/2 + row*fWrappingThickness + (row+0.5)*cellHeight, 0);

			// Copy numbers (segment IDs), indexed from 1
			std::stringstream stream; stream << "Scint-" << col << "," << row;
			cellScint_physV[col][row] = addSegmentToBody(cellScint_logV, stream.str(), cellCenter);
			scintBody_physV.push_back(cellScint_physV[col][row]);
		
			// Place vertical and horizontal reflectors.
			if(WrappingEnabled()){ 
				if(row == 0 && col != Ncol-1){ // New vertical reflector layer.
					std::stringstream stream2; stream2 << "Wrapping-" << col;
					mylarVertLayer_physV[col] = addToDetectorBody(mylarVertLayer_logV, stream2.str().c_str(), G4ThreeVector(cellCenter.getX()+cellWidth/2+fWrappingThickness/2, 0, 0));
				}
				if(row != Nrow-1){ // New horizontal reflector layer.
					std::stringstream stream2; stream2 << "Wrapping-" << col << "," << row;
					mylarHorizLayer_physV[col][row] = addToDetectorBody(mylarHorizLayer_logV, stream2.str().c_str(), G4ThreeVector(cellCenter.getX(), cellCenter.getY()+cellHeight/2+fWrappingThickness/2, 0));
				}
			}
		}
	}

	// Define logical reflector surfaces.
	if(WrappingEnabled()){ 
		for(int col = 0; col < Ncol; col++){
			for(int row = 0; row < Nrow; row++){
				G4PVPlacement *cellPhysical = cellScint_physV[col][row];
				
				int leftCol = col-1;
				int rightCol = col+1;
				int downRow = row-1;
				int upRow = row+1;
				
				// Border with the outer wrapping.
				if((col == 0 || row == 0) || (col == Ncol-1 || row == Nrow-1)) 
					new G4LogicalBorderSurface("Wrapping", cellPhysical, wrapping_physV, wrappingOpSurf);
				
				// Internal reflector layers.
				if(leftCol >= 0 && leftCol < Ncol) // Left side vertical layer.
					new G4LogicalBorderSurface("Wrapping", cellPhysical, mylarVertLayer_physV.at(col-1), wrappingOpSurf);
				if(rightCol >= 0 && rightCol < Ncol) // Right side vertical layer.
					new G4LogicalBorderSurface("Wrapping", cellPhysical, mylarVertLayer_physV.at(col), wrappingOpSurf);
				if(downRow >= 0 && downRow < Nrow) // Bottom side horizontal layer.
					new G4LogicalBorderSurface("Wrapping", cellPhysical, mylarHorizLayer_physV.at(col).at(row-1), wrappingOpSurf);
				if(upRow >= 0 && upRow < Nrow) // Top side vertical layer.
					new G4LogicalBorderSurface("Wrapping", cellPhysical, mylarHorizLayer_physV.at(col).at(row), wrappingOpSurf);
			}
		}
	}
	
	// Update the scintillator copy number.
	scintCopyNum += Nrow*Ncol;
	
	// Update the Z offset and layer width/height
	layerSizeX = fDetectorWidth;
	layerSizeY = fDetectorHeight;
	offsetZ = fDetectorLength/2;
}

///////////////////////////////////////////////////////////////////////////////
// class ellipticalType
///////////////////////////////////////////////////////////////////////////////

void ellipticalType::prepareToBuild(){
	if(fTrapezoidLength <= 0) // User failed to specify the trapezoid length, use a default
		fTrapezoidLength = 10*mm;

	// Width of the detector (defined by the trapezoid length, trapezoid angle, and SiPM dimensions).
	fDetectorWidth = 2*(pmtWidth/2+(fTrapezoidLength/std::tan(fTrapezoidAngle*deg)))*mm;
	fDetectorLength += 2*fTrapezoidLength;
	
	// Set the maximum detector body size
	maxBodySize = G4ThreeVector(fDetectorWidth, fDetectorHeight, fDetectorLength);
}

void ellipticalType::buildDetector(){
	const G4double deltaWidth = fWrappingThickness/std::sin(fTrapezoidAngle*deg);
	const G4double bodyLength = fDetectorLength - 2*fTrapezoidLength;

	// Create the geometry.
	G4Trd *innerTrapezoid = new G4Trd("innerTrapezoid", fDetectorWidth/2, pmtWidth/2, fDetectorHeight/2, pmtHeight/2, fTrapezoidLength/2);
	G4Trd *outerTrapezoid = new G4Trd("outerTrapezoid", fDetectorWidth/2+deltaWidth, pmtWidth/2+deltaWidth, fDetectorHeight/2+deltaWidth, pmtHeight/2+deltaWidth, fTrapezoidLength/2);
	G4Box *innerBody = new G4Box("innerBody", fDetectorWidth/2, fDetectorHeight/2, bodyLength/2);
	G4Box *outerBody = new G4Box("outerBody", fDetectorWidth/2+deltaWidth, fDetectorHeight/2+deltaWidth, bodyLength/2);

	// Build the detector body using unions.
	G4RotationMatrix *trapRot = new G4RotationMatrix;
	trapRot->rotateY(180*deg);
	G4UnionSolid *scintBody1 = new G4UnionSolid("", innerBody, innerTrapezoid, 0, G4ThreeVector(0, 0, +bodyLength/2+fTrapezoidLength/2));
	G4UnionSolid *scintBody2 = new G4UnionSolid("scint", scintBody1, innerTrapezoid, trapRot, G4ThreeVector(0, 0, -bodyLength/2-fTrapezoidLength/2));
	G4LogicalVolume *scint_logV = new G4LogicalVolume(scintBody2, scintMaterial, "scint_logV");	
    scint_logV->SetVisAttributes(scintVisAtt);

	// Place the scintillator inside the assembly.
	G4PVPlacement *ellipseBody_physV = addToDetectorBody(scint_logV, "Scint");

	scintBody_physV.push_back(ellipseBody_physV);

	// Build the wrapping.
	G4PVPlacement *ellipseWrapping_physV = NULL;
	if(WrappingEnabled()){
		G4UnionSolid *wrappingBody1 = new G4UnionSolid("", outerBody, outerTrapezoid, 0, G4ThreeVector(0, 0, +bodyLength/2+fTrapezoidLength/2));
		G4UnionSolid *wrappingBody2 = new G4UnionSolid("", wrappingBody1, outerTrapezoid, trapRot, G4ThreeVector(0, 0, -bodyLength/2-fTrapezoidLength/2));
		G4SubtractionSolid *wrappingBody3 = new G4SubtractionSolid("wrapping", wrappingBody2, scintBody2);
		G4LogicalVolume *wrapping_logV = new G4LogicalVolume(wrappingBody3, wrappingMaterial, "wrapping_logV");		
		wrapping_logV->SetVisAttributes(wrappingVisAtt);

		// Place the wrapping around the scintillator.
		ellipseWrapping_physV = addToDetectorBody(wrapping_logV, "Wrapping");
	}

	// Account for the trapezoids and update the layer width/height
	layerSizeX = fDetectorWidth;
	layerSizeY = fDetectorHeight;
	offsetZ = fDetectorLength/2;
	
	// Directly modify the size of the grease layer.
	layerSizeX = pmtWidth;

    // Reflective wrapping.
    if(WrappingEnabled())
	    new G4LogicalBorderSurface("Wrapping", ellipseBody_physV, ellipseWrapping_physV, wrappingOpSurf);
}

///////////////////////////////////////////////////////////////////////////////
// class rectangularType
///////////////////////////////////////////////////////////////////////////////

void rectangularType::prepareToBuild(){
	// Set the maximum detector body size
	maxBodySize = G4ThreeVector(fDetectorWidth, fDetectorHeight, fDetectorLength);
}

void rectangularType::buildDetector(){
	G4Box *plateBody = new G4Box("", fDetectorWidth/2, fDetectorHeight/2, fDetectorLength/2);
    G4LogicalVolume *plateBody_logV = new G4LogicalVolume(plateBody, scintMaterial, "plateBody_logV");
    plateBody_logV->SetVisAttributes(scintVisAtt);

	// Place the scintillator inside the assembly.
	G4PVPlacement *plateBody_physV = addToDetectorBody(plateBody_logV, "Scint");

	scintBody_physV.push_back(plateBody_physV);

	// Build the wrapping.
	if(WrappingEnabled()){
		G4Box *plateWrappingBox = new G4Box("", fDetectorWidth/2 + fWrappingThickness, fDetectorHeight/2 + fWrappingThickness, fDetectorLength/2);
		G4SubtractionSolid *plateWrapping = new G4SubtractionSolid("", plateWrappingBox, plateBody);
		G4LogicalVolume *plateWrapping_logV = new G4LogicalVolume(plateWrapping, wrappingMaterial, "plateWrapping_logV");
		plateWrapping_logV->SetVisAttributes(wrappingVisAtt);
		
		// Place the wrapping around the scintillator.
		G4PVPlacement *plateWrapping_physV = addToDetectorBody(plateWrapping_logV, "Wrapping");
		
		// Reflective wrapping.
		
		new G4LogicalBorderSurface("Wrapping", plateBody_physV, plateWrapping_physV, wrappingOpSurf);
	}
	
	// Update the Z offset and layer width/height
	layerSizeX = fDetectorWidth;
	layerSizeY = fDetectorHeight;
	offsetZ = fDetectorLength/2;
}

///////////////////////////////////////////////////////////////////////////////
// class cylindricalType
///////////////////////////////////////////////////////////////////////////////

void cylindricalType::prepareToBuild(){
	SetSquarePMTs(false); // Disable square PMTs
	
	// Set the maximum detector body size. Note here that X=Y=fDetectorWidth
	maxBodySize = G4ThreeVector(fDetectorWidth, fDetectorWidth, fDetectorLength);
}

void cylindricalType::buildDetector(){
	// Make sure the height and width match
	fDetectorHeight = fDetectorWidth;

    G4Tubs *cylinderBody = new G4Tubs("scintBody", 0, fDetectorWidth/2, fDetectorLength/2, 0, 2*CLHEP::pi);
    G4LogicalVolume *cylinderBody_logV = new G4LogicalVolume(cylinderBody, scintMaterial, "cylinderBody_logV");
    cylinderBody_logV->SetVisAttributes(scintVisAtt);

	// Place the scintillator inside the assembly.
	G4PVPlacement *cylinderBody_physV = addToDetectorBody(cylinderBody_logV, "Scint");

	scintBody_physV.push_back(cylinderBody_physV);

	// Build the wrapping.
	if(WrappingEnabled()){
		G4Tubs *cylinderWrappingBox = new G4Tubs("outerWrapping", 0, fDetectorWidth/2 + fWrappingThickness, fDetectorLength/2, 0, 2*CLHEP::pi);
		G4SubtractionSolid *cylinderWrapping = new G4SubtractionSolid("", cylinderWrappingBox, cylinderBody);
		G4LogicalVolume *cylinderWrapping_logV = new G4LogicalVolume(cylinderWrapping, wrappingMaterial, "cylinderWrapping_logV");
		cylinderWrapping_logV->SetVisAttributes(wrappingVisAtt);
		
		// Place the wrapping around the scintillator.
		G4PVPlacement *cylinderWrapping_physV = addToDetectorBody(cylinderWrapping_logV, "Wrapping");	
		
		// Reflective wrapping.
		new G4LogicalBorderSurface("Wrapping", cylinderBody_physV, cylinderWrapping_physV, wrappingOpSurf);
	}
	
	// Update the Z offset and layer width/height
	layerSizeX = fDetectorWidth;
	layerSizeY = fDetectorHeight;
	offsetZ = fDetectorLength/2;
}



///////////////////////////////////////////////////////////////////////////////
// class implantType
///////////////////////////////////////////////////////////////////////////////

void implantType::prepareToBuild(){
	// Set the maximum detector body size
	maxBodySize = G4ThreeVector(fDetectorWidth, fDetectorHeight, fDetectorLength);
}

void implantType::buildDetector(){
	const G4double cellWidth = GetSegmentWidth();
	const G4double cellHeight = GetSegmentHeight();
	
	// Get the number of rows and columns for this segmented detector
	int Ncol = fNumColumns;
	int Nrow = fNumRows;

	// Update the detector's copy numbers
	firstSegmentCopyNum = scintCopyNum; 
	lastSegmentCopyNum = firstSegmentCopyNum; //startCopyNum+col*row;

    // Construct the scintillator cell
    G4Box *cellScint = new G4Box("scintillator", cellWidth/2, cellHeight/2, fDetectorLength/2);
    G4LogicalVolume *cellScint_logV = new G4LogicalVolume(cellScint, scintMaterial, "scint_log");
    cellScint_logV->SetVisAttributes(scintVisAtt);

	G4Box *mylarVertLayer = NULL;
	G4Box *mylarHorizLayer = NULL;
	
	G4LogicalVolume *mylarVertLayer_logV = NULL;
	G4LogicalVolume *mylarHorizLayer_logV = NULL;

	// Build the wrapping.
	G4PVPlacement *wrapping_physV = NULL;
	G4PVPlacement* wrappingFacePhys = NULL;
	if(WrappingEnabled()){
		// Construct the outer wrapping.
		G4Box *wrappingBox = new G4Box("wrappingBox", fDetectorWidth/2+fWrappingThickness, fDetectorHeight/2+fWrappingThickness, fDetectorLength/2);
		G4Box *scintBox = new G4Box("scintBox", fDetectorWidth/2, fDetectorHeight/2, fDetectorLength/2);

		G4SubtractionSolid *wrappingBody = new G4SubtractionSolid("wrapping", wrappingBox, scintBox);
		G4LogicalVolume *wrapping_logV = new G4LogicalVolume(wrappingBody, outerMylar, "wrapping_logV");
		wrapping_logV->SetVisAttributes(wrappingVisAtt);

		// Place the outer wrapping into the assembly.
		wrapping_physV = addToDetectorBody(wrapping_logV, "Wrapping");

		if (UseFrontWrapping()) {
			// Place front wrapping unless specified otherwise
			auto wrappingFace = new G4Box("wrappingFace",fDetectorWidth/2+fWrappingThickness, fDetectorHeight/2+fWrappingThickness,fWrappingThickness/4.); 
			auto wrappingFace_log = new G4LogicalVolume(wrappingFace,outerMylar,"wrappingFace_log");
			wrappingFace_log->SetVisAttributes(windowVisAtt);	
			wrappingFacePhys = addFrontComponent(wrappingFace_log,-1*fDetectorLength/2-fWrappingThickness/4.,"wrappingFacePhys");
		}

		// Construct vertical and horizontal reflector layers for later use.
		mylarVertLayer = new G4Box("mylarVertLayer", fWrappingThickness/2, fDetectorHeight/2, fDetectorLength/2);
		mylarHorizLayer = new G4Box("mylarHorizLayer", cellWidth/2, fWrappingThickness/2, fDetectorLength/2);

		mylarVertLayer_logV = new G4LogicalVolume(mylarVertLayer, wrappingMaterial, "mylarVertLayer_logV");
		mylarHorizLayer_logV = new G4LogicalVolume(mylarHorizLayer, wrappingMaterial, "mylarHorizLayer_logV");
		
		mylarVertLayer_logV->SetVisAttributes(wrappingVisAtt);
		mylarHorizLayer_logV->SetVisAttributes(wrappingVisAtt);
	}

	// Place the scintillator segments into the assembly.
	std::vector<G4PVPlacement*> mylarVertLayer_physV(Ncol, NULL);
	std::vector<std::vector<G4PVPlacement*> > mylarHorizLayer_physV(Ncol, std::vector<G4PVPlacement*>(Nrow, NULL));
	std::vector<std::vector<G4PVPlacement*> > cellScint_physV(Ncol, std::vector<G4PVPlacement*>(Nrow, NULL));	
	for(int col = 0; col < Ncol; col++){
		for(int row = 0; row < Nrow; row++){
			G4ThreeVector cellCenter(-fDetectorWidth/2 + col*fWrappingThickness + (col+0.5)*cellWidth, -fDetectorHeight/2 + row*fWrappingThickness + (row+0.5)*cellHeight, 0);

			// Copy numbers (segment IDs), indexed from 1
			std::stringstream stream; stream << "Scint-" << col << "," << row;
			cellScint_physV[col][row] = addSegmentToBody(cellScint_logV, stream.str(), cellCenter);
			scintBody_physV.push_back(cellScint_physV[col][row]);
		
			// Place vertical and horizontal reflectors.
			if(WrappingEnabled()){ 
				if(row == 0 && col != Ncol-1){ // New vertical reflector layer.
					std::stringstream stream2; stream2 << "Wrapping-" << col;
					mylarVertLayer_physV[col] = addToDetectorBody(mylarVertLayer_logV, stream2.str().c_str(), G4ThreeVector(cellCenter.getX()+cellWidth/2+fWrappingThickness/2, 0, 0));
				}
				if(row != Nrow-1){ // New horizontal reflector layer.
					std::stringstream stream2; stream2 << "Wrapping-" << col << "," << row;
					mylarHorizLayer_physV[col][row] = addToDetectorBody(mylarHorizLayer_logV, stream2.str().c_str(), G4ThreeVector(cellCenter.getX(), cellCenter.getY()+cellHeight/2+fWrappingThickness/2, 0));
				}
			}
		}
		
	}

	// Define logical reflector surfaces.
	if(WrappingEnabled()){ 
		std::cout<<"Wrapping the scintillator"<<std::endl;
		for(int col = 0; col < Ncol; col++){
			for(int row = 0; row < Nrow; row++){
				G4PVPlacement *cellPhysical = cellScint_physV[col][row];
				
				int leftCol = col-1;
				int rightCol = col+1;
				int downRow = row-1;
				int upRow = row+1;
				
				// Border with the outer wrapping.
				if((col == 0 || row == 0) || (col == Ncol-1 || row == Nrow-1)) 
					new G4LogicalBorderSurface("Wrapping", cellPhysical, wrapping_physV, wrappingOuterOpSurf);
				// Place front reflector on every cell
				new G4LogicalBorderSurface("wrappingFaceLBS", cellPhysical, wrappingFacePhys, wrappingOuterOpSurf);

				// Internal reflector layers.
				if(leftCol >= 0 && leftCol < Ncol) // Left side vertical layer.
					new G4LogicalBorderSurface("Wrapping", cellPhysical, mylarVertLayer_physV.at(col-1), wrappingOpSurf);
				if(rightCol >= 0 && rightCol < Ncol) // Right side vertical layer.
					new G4LogicalBorderSurface("Wrapping", cellPhysical, mylarVertLayer_physV.at(col), wrappingOpSurf);
				if(downRow >= 0 && downRow < Nrow) // Bottom side horizontal layer.
					new G4LogicalBorderSurface("Wrapping", cellPhysical, mylarHorizLayer_physV.at(col).at(row-1), wrappingOpSurf);
				if(upRow >= 0 && upRow < Nrow) // Top side vertical layer.
					new G4LogicalBorderSurface("Wrapping", cellPhysical, mylarHorizLayer_physV.at(col).at(row), wrappingOpSurf);
			}
		}
	}
	
	// Update the scintillator copy number.
	scintCopyNum += Nrow*Ncol;
	
	// Update the Z offset and layer width/height
	layerSizeX = fDetectorWidth;
	layerSizeY = fDetectorHeight;
	offsetZ = fDetectorLength/2;
}



///////////////////////////////////////////////////////////////////////////////
// class domeType
///////////////////////////////////////////////////////////////////////////////

void domeType::prepareToBuild(){
	// Set the maximum detector body size
	maxBodySize = G4ThreeVector(fDetectorWidth, fDetectorHeight, fDetectorLength);
}

void domeType::buildDetector(){

	/*
	const G4double cellWidth = GetSegmentWidth();
	const G4double cellHeight = GetSegmentHeight();
	
	// Get the number of rows and columns for this segmented detector
	int Ncol = fNumColumns;
	int Nrow = fNumRows;

	// Update the detector's copy numbers
	firstSegmentCopyNum = scintCopyNum; 
	lastSegmentCopyNum = firstSegmentCopyNum; //startCopyNum+col*row;


    // Construct the scintillator cell
    G4Box *cellScint = new G4Box("scintillator", cellWidth/2, cellHeight/2, fDetectorLength/2);
    G4LogicalVolume *cellScint_logV = new G4LogicalVolume(cellScint, scintMaterial, "scint_log");
    cellScint_logV->SetVisAttributes(scintVisAtt);

	G4Box *mylarVertLayer = NULL;
	G4Box *mylarHorizLayer = NULL;
	
	G4LogicalVolume *mylarVertLayer_logV = NULL;
	G4LogicalVolume *mylarHorizLayer_logV = NULL;

	// Build the wrapping.
	G4PVPlacement *wrapping_physV = NULL;
	G4PVPlacement* wrappingFacePhys = NULL;
	*/



	const G4double cellWidth = GetSegmentWidth();
	const G4double cellHeight = GetSegmentHeight();
	
	// Get the number of rows and columns for this segmented detector
	int Ncol = fNumColumns;
	int Nrow = fNumRows;

	// Update the detector's copy numbers
	firstSegmentCopyNum = scintCopyNum; 
	lastSegmentCopyNum = firstSegmentCopyNum; //startCopyNum+col*row;


	G4Box *mylarVertLayer = NULL;
	G4Box *mylarHorizLayer = NULL;
	
	G4LogicalVolume *mylarVertLayer_logV = NULL;
	G4LogicalVolume *mylarHorizLayer_logV = NULL;

	// Build the wrapping.
	G4PVPlacement *wrapping_physV = NULL;
	G4PVPlacement* wrappingFacePhys = NULL;


	// NEED LOGIC TO DECIDE WHETHER IMPLANT IS PIXELATED,
	//  AND PRIORITIES (PIXELATED OR SHAPE MORE IMPORTANT?)



	// NEED TO CREATE A CALCULATOR/CONTROL PARAMETER FOR CURVATURE OF DOME

	// GIVEN RADIUS OF CURVATURE R (SUBJECT TO CONSTRAINT R >= Wdet/sqrt(2)),
	// WE HAVE PIXEL THICKNESS h = sqrt(R*R-2*x*x) = sqrt(R*R-2*y*y)


	// IF DISCRETE/SEGMENTED DOME, DO THE FOLLOWING....

	if (isPixelated && (domeShape == "spherical")) {

		/////////////////
		// NEED TO ALSO ADJUST WRAPPING THICKNESS/LENGTH
		/////////////////

		
		//if(WrappingEnabled()){
		//	// Construct the outer wrapping.
		//	G4Box *wrappingBox = new G4Box("wrappingBox", fDetectorWidth/2+fWrappingThickness, fDetectorHeight/2+fWrappingThickness, fDetectorLength/2);
		//	G4Box *scintBox = new G4Box("scintBox", fDetectorWidth/2, fDetectorHeight/2, fDetectorLength/2);

		//	G4SubtractionSolid *wrappingBody = new G4SubtractionSolid("wrapping", wrappingBox, scintBox);
		//	G4LogicalVolume *wrapping_logV = new G4LogicalVolume(wrappingBody, outerMylar, "wrapping_logV");
		//	wrapping_logV->SetVisAttributes(wrappingVisAtt);

		//	// Place the outer wrapping into the assembly.
		//	wrapping_physV = addToDetectorBody(wrapping_logV, "Wrapping");

		//	//if (UseFrontWrapping()) {
		//	//	// Place front wrapping unless specified otherwise
		//	//	auto wrappingFace = new G4Box("wrappingFace",fDetectorWidth/2+fWrappingThickness, fDetectorHeight/2+fWrappingThickness,fWrappingThickness/4.); 
		//	//	auto wrappingFace_log = new G4LogicalVolume(wrappingFace,outerMylar,"wrappingFace_log");
		//	//	wrappingFace_log->SetVisAttributes(windowVisAtt);	
		//	//	wrappingFacePhys = addFrontComponent(wrappingFace_log,-1*fDetectorLength/2-fWrappingThickness/4.,"wrappingFacePhys");
		//	//}

		//	// Construct vertical and horizontal reflector layers for later use.
		//	mylarVertLayer = new G4Box("mylarVertLayer", fWrappingThickness/2, fDetectorHeight/2, fDetectorLength/2);
		//	mylarHorizLayer = new G4Box("mylarHorizLayer", cellWidth/2, fWrappingThickness/2, fDetectorLength/2);

		//	mylarVertLayer_logV = new G4LogicalVolume(mylarVertLayer, wrappingMaterial, "mylarVertLayer_logV");
		//	mylarHorizLayer_logV = new G4LogicalVolume(mylarHorizLayer, wrappingMaterial, "mylarHorizLayer_logV");
			
		//	mylarVertLayer_logV->SetVisAttributes(wrappingVisAtt);
		//	mylarHorizLayer_logV->SetVisAttributes(wrappingVisAtt);
		//}
		

		// Place the scintillator segments into the assembly.
		std::vector<G4PVPlacement*> mylarVertLayer_physV(Ncol, NULL);
		std::vector<std::vector<G4PVPlacement*> > mylarHorizLayer_physV(Ncol, std::vector<G4PVPlacement*>(Nrow, NULL));
		std::vector<std::vector<G4PVPlacement*> > cellScint_physV(Ncol, std::vector<G4PVPlacement*>(Nrow, NULL));	
		for(int col = 0; col < Ncol; col++){
			for(int row = 0; row < Nrow; row++){
				double x = -fDetectorWidth/2 + col*fWrappingThickness + (col+0.5)*cellWidth;
				if (x < 0)
					x -= 0.5*cellWidth;
				else
					x += 0.5*cellWidth;
				double y = -fDetectorHeight/2 + row*fWrappingThickness + (row+0.5)*cellHeight;
				if (y < 0)
					y -= 0.5*cellHeight;
				else
					y += 0.5*cellHeight;
				double h = sqrt(domeDimension*domeDimension-x*x-y*y);
				double z_center = fDetectorLength/2-h/2; // where z_max is just the scintillator thickness as set by the user

				// Construct the scintillator cell
				G4Box *cellScint = new G4Box("scintillator", cellWidth/2, cellHeight/2, h/2);
				G4LogicalVolume *cellScint_logV = new G4LogicalVolume(cellScint, scintMaterial, "scint_log");
				cellScint_logV->SetVisAttributes(scintVisAtt);

				G4ThreeVector cellCenter(-fDetectorWidth/2 + col*fWrappingThickness + (col+0.5)*cellWidth, -fDetectorHeight/2 + row*fWrappingThickness + (row+0.5)*cellHeight, z_center);

				// Copy numbers (segment IDs), indexed from 1
				std::stringstream stream; stream << "Scint-" << col << "," << row;
				cellScint_physV[col][row] = addSegmentToBody(cellScint_logV, stream.str(), cellCenter);
				scintBody_physV.push_back(cellScint_physV[col][row]);
			
				// Place vertical and horizontal reflectors.
				//if(WrappingEnabled()){ 
				//	if(row == 0 && col != Ncol-1){ // New vertical reflector layer.
				//		std::stringstream stream2; stream2 << "Wrapping-" << col;
				//		mylarVertLayer_physV[col] = addToDetectorBody(mylarVertLayer_logV, stream2.str().c_str(), G4ThreeVector(cellCenter.getX()+cellWidth/2+fWrappingThickness/2, 0, (h-fDetectorLength/2)/2));
				//	}
				//	if(row != Nrow-1){ // New horizontal reflector layer.
				//		std::stringstream stream2; stream2 << "Wrapping-" << col << "," << row;
				//		mylarHorizLayer_physV[col][row] = addToDetectorBody(mylarHorizLayer_logV, stream2.str().c_str(), G4ThreeVector(cellCenter.getX(), cellCenter.getY()+cellHeight/2+fWrappingThickness/2, (h-fDetectorLength/2)/2));
				//	}
			}
		}
	
		/*
		// Define logical reflector surfaces.
		if(WrappingEnabled()){ 
			std::cout<<"Wrapping the scintillator"<<std::endl;
			for(int col = 0; col < Ncol; col++){
				for(int row = 0; row < Nrow; row++){
					G4PVPlacement *cellPhysical = cellScint_physV[col][row];
					
					int leftCol = col-1;
					int rightCol = col+1;
					int downRow = row-1;
					int upRow = row+1;
					
					// Border with the outer wrapping.
					if((col == 0 || row == 0) || (col == Ncol-1 || row == Nrow-1)) 
						new G4LogicalBorderSurface("Wrapping", cellPhysical, wrapping_physV, wrappingOuterOpSurf);
					// Place front reflector on every cell
					new G4LogicalBorderSurface("wrappingFaceLBS", cellPhysical, wrappingFacePhys, wrappingOuterOpSurf);

					// Internal reflector layers.
					if(leftCol >= 0 && leftCol < Ncol) // Left side vertical layer.
						new G4LogicalBorderSurface("Wrapping", cellPhysical, mylarVertLayer_physV.at(col-1), wrappingOpSurf);
					if(rightCol >= 0 && rightCol < Ncol) // Right side vertical layer.
						new G4LogicalBorderSurface("Wrapping", cellPhysical, mylarVertLayer_physV.at(col), wrappingOpSurf);
					if(downRow >= 0 && downRow < Nrow) // Bottom side horizontal layer.
						new G4LogicalBorderSurface("Wrapping", cellPhysical, mylarHorizLayer_physV.at(col).at(row-1), wrappingOpSurf);
					if(upRow >= 0 && upRow < Nrow) // Top side vertical layer.
						new G4LogicalBorderSurface("Wrapping", cellPhysical, mylarHorizLayer_physV.at(col).at(row), wrappingOpSurf);
				}
			}
		}
		*/
	}
	
	else if (isPixelated && (domeShape == "pyramidal")) {

		/////////////////
		// NEED TO ALSO ADJUST WRAPPING THICKNESS/LENGTH
		/////////////////

		
		//if(WrappingEnabled()){
		//	// Construct the outer wrapping.
		//	G4Box *wrappingBox = new G4Box("wrappingBox", fDetectorWidth/2+fWrappingThickness, fDetectorHeight/2+fWrappingThickness, fDetectorLength/2);
		//	G4Box *scintBox = new G4Box("scintBox", fDetectorWidth/2, fDetectorHeight/2, fDetectorLength/2);

		//	G4SubtractionSolid *wrappingBody = new G4SubtractionSolid("wrapping", wrappingBox, scintBox);
		//	G4LogicalVolume *wrapping_logV = new G4LogicalVolume(wrappingBody, outerMylar, "wrapping_logV");
		//	wrapping_logV->SetVisAttributes(wrappingVisAtt);

		//	// Place the outer wrapping into the assembly.
		//	wrapping_physV = addToDetectorBody(wrapping_logV, "Wrapping");

		//	//if (UseFrontWrapping()) {
		//	//	// Place front wrapping unless specified otherwise
		//	//	auto wrappingFace = new G4Box("wrappingFace",fDetectorWidth/2+fWrappingThickness, fDetectorHeight/2+fWrappingThickness,fWrappingThickness/4.); 
		//	//	auto wrappingFace_log = new G4LogicalVolume(wrappingFace,outerMylar,"wrappingFace_log");
		//	//	wrappingFace_log->SetVisAttributes(windowVisAtt);	
		//	//	wrappingFacePhys = addFrontComponent(wrappingFace_log,-1*fDetectorLength/2-fWrappingThickness/4.,"wrappingFacePhys");
		//	//}

		//	// Construct vertical and horizontal reflector layers for later use.
		//	mylarVertLayer = new G4Box("mylarVertLayer", fWrappingThickness/2, fDetectorHeight/2, fDetectorLength/2);
		//	mylarHorizLayer = new G4Box("mylarHorizLayer", cellWidth/2, fWrappingThickness/2, fDetectorLength/2);

		//	mylarVertLayer_logV = new G4LogicalVolume(mylarVertLayer, wrappingMaterial, "mylarVertLayer_logV");
		//	mylarHorizLayer_logV = new G4LogicalVolume(mylarHorizLayer, wrappingMaterial, "mylarHorizLayer_logV");
			
		//	mylarVertLayer_logV->SetVisAttributes(wrappingVisAtt);
		//	mylarHorizLayer_logV->SetVisAttributes(wrappingVisAtt);
		//}
		

		// Place the scintillator segments into the assembly.
		std::vector<G4PVPlacement*> mylarVertLayer_physV(Ncol, NULL);
		std::vector<std::vector<G4PVPlacement*> > mylarHorizLayer_physV(Ncol, std::vector<G4PVPlacement*>(Nrow, NULL));
		std::vector<std::vector<G4PVPlacement*> > cellScint_physV(Ncol, std::vector<G4PVPlacement*>(Nrow, NULL));	
		for(int col = 0; col < Ncol; col++){
			for(int row = 0; row < Nrow; row++){
				//double x = -fDetectorWidth/2 + col*fWrappingThickness + (col+0.5)*cellWidth;
				//if (x < 0)
				//	x -= 0.5*cellWidth;
				//else
				//	x += 0.5*cellWidth;
				//double y = -fDetectorHeight/2 + row*fWrappingThickness + (row+0.5)*cellHeight;
				//if (y < 0)
				//	y -= 0.5*cellHeight;
				//else
				//	y += 0.5*cellHeight;
				double midpoint = (Ncol-1)/2.0;
				int numSteps = std::max(std::floor(std::fabs(col-midpoint)),std::floor(std::fabs(row-midpoint)));
				double h = (fDetectorLength-numSteps*domeDimension > marginSize)? fDetectorLength-numSteps*domeDimension : marginSize;

				//double h = sqrt(domeDimension*domeDimension-x*x-y*y);
				double z_center = fDetectorLength/2-h/2; // where z_max is just the scintillator thickness as set by the user

				// Construct the scintillator cell
				G4Box *cellScint = new G4Box("scintillator", cellWidth/2, cellHeight/2, h/2);
				G4LogicalVolume *cellScint_logV = new G4LogicalVolume(cellScint, scintMaterial, "scint_log");
				cellScint_logV->SetVisAttributes(scintVisAtt);

				G4ThreeVector cellCenter(-fDetectorWidth/2 + col*fWrappingThickness + (col+0.5)*cellWidth, -fDetectorHeight/2 + row*fWrappingThickness + (row+0.5)*cellHeight, z_center);

				// Copy numbers (segment IDs), indexed from 1
				std::stringstream stream; stream << "Scint-" << col << "," << row;
				cellScint_physV[col][row] = addSegmentToBody(cellScint_logV, stream.str(), cellCenter);
				scintBody_physV.push_back(cellScint_physV[col][row]);
			
				// Place vertical and horizontal reflectors.
				//if(WrappingEnabled()){ 
				//	if(row == 0 && col != Ncol-1){ // New vertical reflector layer.
				//		std::stringstream stream2; stream2 << "Wrapping-" << col;
				//		mylarVertLayer_physV[col] = addToDetectorBody(mylarVertLayer_logV, stream2.str().c_str(), G4ThreeVector(cellCenter.getX()+cellWidth/2+fWrappingThickness/2, 0, (h-fDetectorLength/2)/2));
				//	}
				//	if(row != Nrow-1){ // New horizontal reflector layer.
				//		std::stringstream stream2; stream2 << "Wrapping-" << col << "," << row;
				//		mylarHorizLayer_physV[col][row] = addToDetectorBody(mylarHorizLayer_logV, stream2.str().c_str(), G4ThreeVector(cellCenter.getX(), cellCenter.getY()+cellHeight/2+fWrappingThickness/2, (h-fDetectorLength/2)/2));
				//	}
			}
		}
	
		/*
		// Define logical reflector surfaces.
		if(WrappingEnabled()){ 
			std::cout<<"Wrapping the scintillator"<<std::endl;
			for(int col = 0; col < Ncol; col++){
				for(int row = 0; row < Nrow; row++){
					G4PVPlacement *cellPhysical = cellScint_physV[col][row];
					
					int leftCol = col-1;
					int rightCol = col+1;
					int downRow = row-1;
					int upRow = row+1;
					
					// Border with the outer wrapping.
					if((col == 0 || row == 0) || (col == Ncol-1 || row == Nrow-1)) 
						new G4LogicalBorderSurface("Wrapping", cellPhysical, wrapping_physV, wrappingOuterOpSurf);
					// Place front reflector on every cell
					new G4LogicalBorderSurface("wrappingFaceLBS", cellPhysical, wrappingFacePhys, wrappingOuterOpSurf);

					// Internal reflector layers.
					if(leftCol >= 0 && leftCol < Ncol) // Left side vertical layer.
						new G4LogicalBorderSurface("Wrapping", cellPhysical, mylarVertLayer_physV.at(col-1), wrappingOpSurf);
					if(rightCol >= 0 && rightCol < Ncol) // Right side vertical layer.
						new G4LogicalBorderSurface("Wrapping", cellPhysical, mylarVertLayer_physV.at(col), wrappingOpSurf);
					if(downRow >= 0 && downRow < Nrow) // Bottom side horizontal layer.
						new G4LogicalBorderSurface("Wrapping", cellPhysical, mylarHorizLayer_physV.at(col).at(row-1), wrappingOpSurf);
					if(upRow >= 0 && upRow < Nrow) // Top side vertical layer.
						new G4LogicalBorderSurface("Wrapping", cellPhysical, mylarHorizLayer_physV.at(col).at(row), wrappingOpSurf);
				}
			}
		}
		*/
	}


	// ELSE, CONTINUOUS DOME, DO THE FOLLOWING:


	// CASE 1: SPHERICAL DOME (CONTROL PARAMETER IS RADIUS OF CURVATURE)

	// GIVEN RADIUS OF CURVATURE R (SUBJECT TO CONSTRAINT R >= Wdet/sqrt(2)),
	// WE HAVE PIXEL THICKNESS h = sqrt(R*R-2*x*x) = sqrt(R*R-2*y*y)



	// CASE 2: PYRAMIDAL DOME (NO CONTROL PARAMETER)



	// CASE 3: IMPORT GDML GEOMETRY (EXCESS THICKNESS?)


	// NOTE: WILL NEED TO FORCE THE NUMBER OF SEGMENTS TO 1 FOR CONTINUOUS DOME


	// Update the scintillator copy number.
	scintCopyNum += Nrow*Ncol;
	
	// Update the Z offset and layer width/height
	layerSizeX = fDetectorWidth;
	layerSizeY = fDetectorHeight;
	offsetZ = fDetectorLength/2;
}




///////////////////////////////////////////////////////////////////////////////
// class hagridType
///////////////////////////////////////////////////////////////////////////////

void hagridType::prepareToBuild(){
	SetSquarePMTs(false); // Disable square PMTs
	
	// Set the maximum detector body size. Note here that X=Y=fDetectorWidth
	maxBodySize = G4ThreeVector(fDetectorWidth, fDetectorWidth, fDetectorLength);
}

void hagridType::buildDetector(){
	// Make sure the height and width match
	fDetectorHeight = fDetectorWidth;

    G4Tubs *cylinderBody = new G4Tubs("scintBody", 0, fDetectorWidth/2, fDetectorLength/2, 0, 2*CLHEP::pi);
    G4LogicalVolume *cylinderBody_logV = new G4LogicalVolume(cylinderBody, scintMaterial, "cylinderBody_logV");
    cylinderBody_logV->SetVisAttributes(scintVisAtt);

	// Place the scintillator inside the assembly.
	G4PVPlacement *cylinderBody_physV = addToDetectorBody(cylinderBody_logV, "Scint");

	scintBody_physV.push_back(cylinderBody_physV);
	G4PVPlacement *wrapping_physV = NULL;
	G4PVPlacement* wrappingFacePhys = NULL;

	// Build the wrapping.
	if(WrappingEnabled()){
		std::cout<<"Wrapping Thickness: "<<fWrappingThickness/cm<<" detWidth "<<fDetectorWidth/cm<<std::endl;
		G4Tubs *cylinderWrappingBox = new G4Tubs("outerWrapping", 0, fDetectorWidth/2 + fWrappingThickness, fDetectorLength/2, 0, 2*CLHEP::pi);
		G4SubtractionSolid *cylinderWrapping = new G4SubtractionSolid("", cylinderWrappingBox, cylinderBody);
		G4LogicalVolume *cylinderWrapping_logV = new G4LogicalVolume(cylinderWrapping, wrappingMaterial, "cylinderWrapping_logV");
		cylinderWrapping_logV->SetVisAttributes(wrappingVisAtt);
		
		// Place the wrapping around the scintillator.
		G4PVPlacement *cylinderWrapping_physV = addToDetectorBody(cylinderWrapping_logV, "Wrapping");	

		auto wrappingFace = new G4Tubs("wrappingFace",0,fDetectorWidth/2+fWrappingThickness, fWrappingThickness/2,0,CLHEP::twopi); 
		auto wrappingFace_log = new G4LogicalVolume(wrappingFace,outerMylar,"wrappingFace_log");
		wrappingFace_log->SetVisAttributes(windowVisAtt);	
		wrappingFacePhys = addFrontComponent(wrappingFace_log,-1*fDetectorLength/2-fWrappingThickness/2.,"wrappingFacePhys");
		
		// Reflective wrapping.
		new G4LogicalBorderSurface("Wrapping", cylinderBody_physV, cylinderWrapping_physV, wrappingOpSurf);
	}
	
	// Update the Z offset and layer width/height
	layerSizeX = fDetectorWidth;
	layerSizeY = fDetectorHeight;
	offsetZ = fDetectorLength/2;
}



///////////////////////////////////////////////////////////////////////////////
// class cubeType
///////////////////////////////////////////////////////////////////////////////
/*

void cubeType::prepareToBuild(){
	// Set the maximum detector body size
	maxBodySize = G4ThreeVector(fDetectorWidth, fDetectorHeight, fDetectorLength);
}

void cylindricalType::buildDetector(){
	const G4double cellWidth = GetSegmentWidth();
	const G4double cellHeight = GetSegmentHeight();
	
	// Get the number of rows and columns for this segmented detector
	int Ncol = fNumColumns;
	int Nrow = fNumRows;

	// Update the detector's copy numbers
	firstSegmentCopyNum = scintCopyNum; 
	lastSegmentCopyNum = firstSegmentCopyNum; //startCopyNum+col*row;

    // Construct the scintillator cell
    G4Box *cellScint = new G4Box("scintillator", cellWidth/2, cellHeight/2, fDetectorLength/2);
    G4LogicalVolume *cellScint_logV = new G4LogicalVolume(cellScint, scintMaterial, "scint_log");
    cellScint_logV->SetVisAttributes(scintVisAtt);

	G4Box *mylarVertLayer = NULL;
	G4Box *mylarHorizLayer = NULL;
	
	G4LogicalVolume *mylarVertLayer_logV = NULL;
	G4LogicalVolume *mylarHorizLayer_logV = NULL;

	// Build the wrapping.
	G4PVPlacement *wrapping_physV = NULL;
	G4PVPlacement* wrappingFacePhys = NULL;
	if(WrappingEnabled()){
		// Construct the outer wrapping.
		G4Box *wrappingBox = new G4Box("wrappingBox", fDetectorWidth/2+fWrappingThickness, fDetectorHeight/2+fWrappingThickness, fDetectorLength/2);
		G4Box *scintBox = new G4Box("scintBox", fDetectorWidth/2, fDetectorHeight/2, fDetectorLength/2);
		
		G4SubtractionSolid *wrappingBody = new G4SubtractionSolid("wrapping", wrappingBox, scintBox);
		G4LogicalVolume *wrapping_logV = new G4LogicalVolume(wrappingBody, outerMylar, "wrapping_logV");
		wrapping_logV->SetVisAttributes(wrappingVisAtt);
		
		// Place the outer wrapping into the assembly.
		wrapping_physV = addToDetectorBody(wrapping_logV, "Wrapping");

		if (UseFrontWrapping()) {
			// Place front wrapping unless specified otherwise
			auto wrappingFace = new G4Box("wrappingFace",fDetectorWidth/2+fWrappingThickness, fDetectorHeight/2+fWrappingThickness,fWrappingThickness/4.); 
			auto wrappingFace_log = new G4LogicalVolume(wrappingFace,outerMylar,"wrappingFace_log");
			wrappingFace_log->SetVisAttributes(windowVisAtt);	
			wrappingFacePhys = addFrontComponent(wrappingFace_log,-1*fDetectorLength/2-fWrappingThickness/4.,"wrappingFacePhys");
		}

		// Construct vertical and horizontal reflector layers for later use.
		mylarVertLayer = new G4Box("mylarVertLayer", fWrappingThickness/2, fDetectorHeight/2, fDetectorLength/2);
		mylarHorizLayer = new G4Box("mylarHorizLayer", cellWidth/2, fWrappingThickness/2, fDetectorLength/2);

		mylarVertLayer_logV = new G4LogicalVolume(mylarVertLayer, wrappingMaterial, "mylarVertLayer_logV");
		mylarHorizLayer_logV = new G4LogicalVolume(mylarHorizLayer, wrappingMaterial, "mylarHorizLayer_logV");
		
		mylarVertLayer_logV->SetVisAttributes(wrappingVisAtt);
		mylarHorizLayer_logV->SetVisAttributes(wrappingVisAtt);
	}

	// Place the scintillator segments into the assembly.
	std::vector<G4PVPlacement*> mylarVertLayer_physV(Ncol, NULL);
	std::vector<std::vector<G4PVPlacement*> > mylarHorizLayer_physV(Ncol, std::vector<G4PVPlacement*>(Nrow, NULL));
	std::vector<std::vector<G4PVPlacement*> > cellScint_physV(Ncol, std::vector<G4PVPlacement*>(Nrow, NULL));	
	for(int col = 0; col < Ncol; col++){
		for(int row = 0; row < Nrow; row++){
			G4ThreeVector cellCenter(-fDetectorWidth/2 + col*fWrappingThickness + (col+0.5)*cellWidth, -fDetectorHeight/2 + row*fWrappingThickness + (row+0.5)*cellHeight, 0);

			// Copy numbers (segment IDs), indexed from 1
			std::stringstream stream; stream << "Scint-" << col << "," << row;
			cellScint_physV[col][row] = addSegmentToBody(cellScint_logV, stream.str(), cellCenter);
			scintBody_physV.push_back(cellScint_physV[col][row]);
		
			// Place vertical and horizontal reflectors.
			if(WrappingEnabled()){ 
				if(row == 0 && col != Ncol-1){ // New vertical reflector layer.
					std::stringstream stream2; stream2 << "Wrapping-" << col;
					mylarVertLayer_physV[col] = addToDetectorBody(mylarVertLayer_logV, stream2.str().c_str(), G4ThreeVector(cellCenter.getX()+cellWidth/2+fWrappingThickness/2, 0, 0));
				}
				if(row != Nrow-1){ // New horizontal reflector layer.
					std::stringstream stream2; stream2 << "Wrapping-" << col << "," << row;
					mylarHorizLayer_physV[col][row] = addToDetectorBody(mylarHorizLayer_logV, stream2.str().c_str(), G4ThreeVector(cellCenter.getX(), cellCenter.getY()+cellHeight/2+fWrappingThickness/2, 0));
				}
			}
		}
		
	}

	// Define logical reflector surfaces.
	if(WrappingEnabled()){ 
		std::cout<<"Wrapping the scintillator"<<std::endl;
		for(int col = 0; col < Ncol; col++){
			for(int row = 0; row < Nrow; row++){
				G4PVPlacement *cellPhysical = cellScint_physV[col][row];
				
				int leftCol = col-1;
				int rightCol = col+1;
				int downRow = row-1;
				int upRow = row+1;
				
				// Border with the outer wrapping.
				if((col == 0 || row == 0) || (col == Ncol-1 || row == Nrow-1)) 
					new G4LogicalBorderSurface("Wrapping", cellPhysical, wrapping_physV, wrappingOuterOpSurf);
				// Place front reflector on every cell
				new G4LogicalBorderSurface("wrappingFaceLBS", cellPhysical, wrappingFacePhys, wrappingOuterOpSurf);

				// Internal reflector layers.
				if(leftCol >= 0 && leftCol < Ncol) // Left side vertical layer.
					new G4LogicalBorderSurface("Wrapping", cellPhysical, mylarVertLayer_physV.at(col-1), wrappingOpSurf);
				if(rightCol >= 0 && rightCol < Ncol) // Right side vertical layer.
					new G4LogicalBorderSurface("Wrapping", cellPhysical, mylarVertLayer_physV.at(col), wrappingOpSurf);
				if(downRow >= 0 && downRow < Nrow) // Bottom side horizontal layer.
					new G4LogicalBorderSurface("Wrapping", cellPhysical, mylarHorizLayer_physV.at(col).at(row-1), wrappingOpSurf);
				if(upRow >= 0 && upRow < Nrow) // Top side vertical layer.
					new G4LogicalBorderSurface("Wrapping", cellPhysical, mylarHorizLayer_physV.at(col).at(row), wrappingOpSurf);
			}
		}
	}
	
	// Update the scintillator copy number.
	scintCopyNum += Nrow*Ncol;
	
	// Update the Z offset and layer width/height
	layerSizeX = fDetectorWidth;
	layerSizeY = fDetectorHeight;
	offsetZ = fDetectorLength/2;
}

*/