#include <fstream>
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

#include "termColors.hh"
#include "nDetMaterials.hh"
#include "nDetMaterialsMessenger.hh"
#include "nDetDynamicMaterial.hh"

#include "optionHandler.hh"

nDetMaterials::nDetMaterials() : messenger(), isInitialized(false), scintsAreDefined(false), lightYieldScale(1) { 
	messenger = new nDetMaterialsMessenger(this);
}

nDetMaterials::~nDetMaterials(){
	if(isInitialized){
		// Optical surfaces
		delete fTeflonOpSurf;
		delete fAl2O3OpSurf;
		delete fSiliconOpSurf;
		delete fPerfectOpSurf;
		delete fGreaseOpSurf;	
		delete fNOA61OpSurf;
		delete fNOA68OpSurf;
		delete fNOA136OpSurf;
		delete fNOA170OpSurf;
		delete fMylarOpSurf;
		delete fEsrOpSurf;
		delete fAluminum;
		
		// Materials
		delete fGrease;
		delete fNOA61;
		delete fNOA68;
		delete fNOA136;
		delete fNOA170;		
		delete fABS;
		delete fAl2O3;
		
		// Material properties tables
		delete fPerfectMPT;
		delete fGreaseMPT;
		delete fNOA61MPT;
		delete fNOA68MPT;
		delete fNOA136MPT;
		delete fNOA170MPT;
		delete fABSMPT;
		delete fAl2O3MPT;
		delete fEsrMPT;

		// Vis attributes
		delete visAssembly;
		delete visSensitive;
		delete visWindow;
		delete visGrease;
		delete visWrapping;
		delete visScint;
		delete visShadow;
	}
	if(scintsAreDefined){
		delete fEJ200;
		delete fEJ276;
		delete fOGS;
		delete fEJ200MPT;
		delete fEJ276MPT;
		delete fOGSMPT;
		delete fYSO;
		delete fLYSO;
		delete fYSOMPT;
		delete fLYSOMPT;
		delete fLaBr3;
		delete fLaBr3MPT;
		delete fYAP;
		delete fYAPMPT;
		delete fGAGG;
		delete fGAGGMPT;
		delete fCeBr3;
		delete fCeBr3MPT;
		//delete fCsI;
		//delete fCsIMPT
	}
	delete messenger;
}

void nDetMaterials::initialize(){
	if(isInitialized) return;
	
	defineMaterials();
	defineScintillators();

	elementList["H"] = fH;
	elementList["C"] = fC;
	elementList["O"] = fO;
	elementList["N"] = fN;
	elementList["S"] = fS;
	elementList["F"] = fF;
	elementList["Si"] = fSi;
	elementList["Al"] = fAl;
 	elementList["Y"] = fY;
	elementList["La"] = fLa;
	elementList["Br"] = fBr;
	//elementList["Cs"] = fCs;
	//elementList["I"] = fI;
	elementList["Gd"] = fGd;
	elementList["Ga"] = fGa;
	elementList["Ce"] = fCe;
	elementList["Zr"] = fZr;
	elementList["Lu"] = fLu;
	
	materialList["air"] = fAir;
	materialList["vacuum"] = fVacuum;
	materialList["teflon"] = fTeflon;
	materialList["al2o3"] = fAl2O3;
	materialList["ej200"] = fEJ200;
	materialList["ej276"] = fEJ276;
	materialList["ogs"] = fOGS;
	materialList["grease"] = fGrease;
	materialList["noa61"] = fNOA61;
	materialList["noa68"] = fNOA68;
	materialList["noa136"] = fNOA136;
	materialList["noa170"] = fNOA170;	
	materialList["abs"] = fABS;
 	materialList["yso"] = fYSO;
	materialList["lyso"] = fLYSO;
	materialList["labr3"] = fLaBr3;
	materialList["yap"] = fYAP;
	materialList["cebr3"] = fCeBr3;
	//materialList["csi"] = fCsI;
	materialList["quartz"] = fSiO2;
	materialList["silicon"] = fSilicon;
	materialList["mylar"] = fMylar;
	materialList["acrylic"] = fAcrylic;
	//materialList["mercapto"] = fMercapto;
	//materialList["triallyl"] = fTriallyl;
	//materialList["tetrahydrofurfuryl"] = fTetrahydroff;
	//materialList["urethaneAcrylate"] = fUrethaneAcrylate;
	//materialList["acrylate"] = fAcrylate;
	//materialList["zro2"] = fZrO2;
	materialList["aluminum"] = fAluminum;

	opticalSurfaceList["teflon"] = fTeflonOpSurf;
	opticalSurfaceList["al2o3"] = fAl2O3OpSurf;
	opticalSurfaceList["silicon"] = fSiliconOpSurf;
	opticalSurfaceList["mylar"] = fMylarOpSurf;
	opticalSurfaceList["esr"] = fEsrOpSurf;
	opticalSurfaceList["perfect"] = fPerfectOpSurf;
	opticalSurfaceList["grease"] = fGreaseOpSurf;
	opticalSurfaceList["noa61"] = fNOA61OpSurf;
	opticalSurfaceList["noa68"] = fNOA68OpSurf;
	opticalSurfaceList["noa136"] = fNOA136OpSurf;
	opticalSurfaceList["noa170"] = fNOA170OpSurf;

	opticalSurfaceList["air"] = fAirOpSurf;
	opticalSurfaceList["aluminum"] = fAluminumOpSurf;

	visAttributesList["air"] = visWrapping;
	//visAttributesList["vacuum"] = &G4VisAttributes::Invisible;
	visAttributesList["teflon"] = visWrapping;
	visAttributesList["al2o3"] = visWrapping;
	visAttributesList["ej200"] = visScint;
	visAttributesList["ej276"] = visScint;
	visAttributesList["ogs"] = visScint;
	visAttributesList["yso"] = visScint;
	visAttributesList["lyso"] = visScint;
	visAttributesList["labr3"] = visScint;
	visAttributesList["yap"] = visScint;
	visAttributesList["cebr3"] = visScint;
	//visAttributesList["csi"] = visScint;
	visAttributesList["grease"] = visGrease;
	visAttributesList["noa61"] = visGrease;
	visAttributesList["noa68"] = visGrease;
	visAttributesList["noa136"] = visGrease;
	visAttributesList["noa170"] = visGrease;	
	visAttributesList["quartz"] = visWindow;
	visAttributesList["silicon"] = visSensitive;
	visAttributesList["mylar"] = visWrapping;
	visAttributesList["acrylic"] = visWindow;
	//visAttributesList["mercapto"] = visGrease;
	//visAttributesList["triallyl"] = visGrease;
	//visAttributesList["urethaneAcrylate"] = visGrease;
	//visAttributesList["acrylate"] = visGrease;
	//visAttributesList["tetrahydrofurfuryl"] = visGrease;
	//visAttributesList["zro2"] = visGrease;
	visAttributesList["aluminum"] = visShadow;
}

G4Material* nDetMaterials::getUserDetectorMaterial(const G4String &name){
	G4Material *mat = getMaterial(name);
	
	if(!mat) // Material was not found
		std::cout << Display::ErrorStr("nDetDynamicMaterial") << "Detector material named \"" << name << "\" was not found in list!\n" << Display::ResetStr();
		
	return mat;
}

G4Material* nDetMaterials::getUserSurfaceMaterial(const G4String &name){
	G4Material *mat = NULL;
	if(name == "mylar" || name == "esr" || name == "perfect")
		mat = fMylar;
	else
		mat = getMaterial(name);

	if(!mat) // Material was not found
		std::cout << Display::ErrorStr("nDetDynamicMaterial") << "Wrapping material named \"" << name << "\" was not found in list!\n" << Display::ResetStr();

	return mat;
}

G4OpticalSurface* nDetMaterials::getUserOpticalSurface(const G4String &name){
	G4OpticalSurface *opt = getOpticalSurface(name);

	if(!opt) // Surface was not found
		std::cout << Display::ErrorStr("nDetDynamicMaterial") << "Optical surface named \"" << name << "\" was not found in list!\n" << Display::ResetStr();

  return opt; // default
}

G4VisAttributes* nDetMaterials::getUserVisAttributes(const G4String &name){
	G4VisAttributes *visatt = NULL;
	if(name == "mylar" || name == "teflon" || name == "al2o3" || name == "esr" || name == "perfect" || "air")
		visatt = visWrapping;
	else
		visatt = getVisualAttributes(name);

	if(!visatt){ // Visual attributes was not found (use a default)
		std::cout << Display::WarningStr("nDetDynamicMaterial") << "Visual attributes named \"" << name << "\" was not found in list!\n" << Display::ResetStr();
		visatt = visAssembly;
	}
		
	return visatt;
}

void nDetMaterials::setLightYield(const G4double &yield){ 
	lightYieldScale = yield; 
	
	if(scintsAreDefined) // Re-define the scintillators as the light yield has changed
		defineScintillators();
}

bool nDetMaterials::searchForElement(const G4String &name){
	size_t dummy;
	return nist.searchElementList(name, dummy);
}

G4Element* nDetMaterials::getElement(const G4String &name){
	std::map<G4String, G4Element*>::iterator iter = elementList.find(name);
	if(iter != elementList.end()) // Found element in the dictionary
		return iter->second;
		
	// No element found in the dictionary, search for it in the NIST database
	G4Element *ptr = nist.searchForElement(name);
	if(ptr){ // Found element in the database, add it to the dictionary
		elementList[name] = ptr;
		return ptr;
	}
	
	return NULL;
}

bool nDetMaterials::searchForMaterial(const G4String &name){
	size_t dummy;
	return nist.searchMaterialList(name, dummy);
}

G4Material* nDetMaterials::getMaterial(const G4String &name){
	std::map<G4String, G4Material*>::iterator iter = materialList.find(name);
	if(iter != materialList.end()) // Found material in the dictionary
		return iter->second;
	// No material found in the dictionary, search for it in the NIST database
	G4Material *ptr = nist.searchForMaterial(name);
	if(ptr){ // Found material in the database, add it to the dictionary
		materialList[name] = ptr;
		return ptr;
	}
	
	return NULL;
}

G4OpticalSurface* nDetMaterials::getOpticalSurface(const G4String &name){
	std::map<G4String, G4OpticalSurface*>::iterator iter = opticalSurfaceList.find(name);
	if(iter != opticalSurfaceList.end()) // Found optical surface in the dictionary
		return iter->second;
	return NULL;
}
	
G4VisAttributes* nDetMaterials::getVisualAttributes(const G4String &name){
	std::map<G4String, G4VisAttributes*>::iterator iter = visAttributesList.find(name);
	if(iter != visAttributesList.end()) // Found optical surface in the dictionary
		return iter->second;
	return NULL;
}

void nDetMaterials::listMaterials() const {
	std::cout << "/////////////////////////////\n";
	std::cout << "// Defined Materials\n";
	std::cout << "/////////////////////////////\n\n";
	for(auto material : materialList){
		std::cout << "  " << material.first << std::endl;
	}
	std::cout << std::endl;
}

void nDetMaterials::listVisAttributes() const {
	std::cout << "/////////////////////////////\n";
	std::cout << "// Defined Visual Attributes\n";
	std::cout << "/////////////////////////////\n\n";
	for(auto visatt : visAttributesList){
		std::cout << "  " << visatt.first << std::endl;
	}
	std::cout << std::endl;
}

void nDetMaterials::listOptSurfaces() const {
	std::cout << "/////////////////////////////\n";
	std::cout << "// Defined Optical Surfaces\n";
	std::cout << "/////////////////////////////\n\n";
	for(auto optsurf : opticalSurfaceList){
		std::cout << "  " << optsurf.first << std::endl;
	}
	std::cout << std::endl;
}

bool nDetMaterials::buildNewMaterial(const G4String &filename){
	nDetDynamicMaterial* dmat = new nDetDynamicMaterial();
	dmat->setScalingFactor(lightYieldScale);
	if(!dmat->read(filename, this)) // Read the material file
		return false;
	std::string name = dmat->getFilePrefix();
	int nameCounter = 2;
	while(materialList.find(name) != materialList.end()){ // Check for material existing in material list
		std::stringstream stream;
		stream << name << "-" << nameCounter++;
		
		name = stream.str();
	}
	if(nameCounter != 2){ // Print a warning about existing material name
		std::cout << Display::WarningStr("nDetDynamicMaterial") << "Material named \"" << dmat->getFilePrefix() << "\" is already in material list!\n" << Display::ResetStr();
		std::cout << Display::WarningStr("nDetDynamicMaterial") << " Renaming new material to \"" << name << "\"\n" << Display::ResetStr();
	}
	materialList[name] = dmat->getMaterial();
	visAttributesList[name] = dmat->getVisAttributes();
	opticalSurfaceList[name] = dmat->getOpticalSurface();
	return true;
}

void nDetMaterials::printMaterial(const G4String &name){
	G4Material *mat = getMaterial(name);
	if(mat){
		std::cout << "/////////////////////////////\n";
		std::cout << "// Material (" << name << ")\n";
		std::cout << "/////////////////////////////\n\n";
		std::cout << (*mat);
		// Re-writing G4MaterialPropertiesTable::DumpTable() because it prints to G4cout only
		// This is messy, but I have to do a lot of workarounds to deal with the strange setup of the G4MaterialPropertiesTable class
		std::cout << "\n/////////////////////////////\n";
		std::cout << "// Properties\n";
		std::cout << "/////////////////////////////\n\n";
		G4MaterialPropertiesTable* MPT = mat->GetMaterialPropertiesTable();
		if(MPT){
#ifndef GEANT_OLDER_VERSION
			std::vector<G4String> propNames = MPT->GetMaterialPropertyNames(); // Ugh...
			std::vector<G4String> cpropNames = MPT->GetMaterialConstPropertyNames(); // UGH... WHY, GEANT???
			for(auto prop : (*MPT->GetPropertyMap())){ // Iterate over variable properties
				G4PhysicsVector *vec = prop.second;
				G4String propertyName = propNames[prop.first];
				std::cout << std::string(propertyName.length(), '-') << std::endl;
				std::cout << propertyName << std::endl;
				std::cout << std::string(propertyName.length(), '-') << std::endl;
				if(vec != NULL){
					// Now I re-write G4PhysicsVector::DumpValues() for the same reason
					for (size_t i = 0; i < vec->GetVectorLength(); i++) // Iterate over all values in the vectors
						std::cout << vec->Energy(i) << "\t" << (*vec)[i] << std::endl;
				}
			}
			std::cout << "\n/////////////////////////////\n";
			std::cout << "// Constant Properties\n";
			std::cout << "/////////////////////////////\n\n";
			for(auto cprop : (*MPT->GetConstPropertyMap())) // Iterate over constant properties
				std::cout << cpropNames[cprop.first] << " = " << cprop.second << std::endl;
#else
			const std::map<G4String, G4MaterialPropertyVector*, std::less< G4String > >* pMap = MPT->GetPropertiesMap();
			const std::map<G4String, G4double, std::less< G4String > >* cpMap = MPT->GetPropertiesCMap();
			for(auto prop : (*pMap)){ // Iterate over variable properties
				G4PhysicsVector *vec = prop.second;
				std::cout << std::string(prop.first.length(), '-') << std::endl;
				std::cout << prop.first << std::endl;
				std::cout << std::string(prop.first.length(), '-') << std::endl;
				if(vec != NULL){
					// Now I re-write G4PhysicsVector::DumpValues() for the same reason
					for (size_t i = 0; i < vec->GetVectorLength(); i++) // Iterate over all values in the vectors
						std::cout << vec->Energy(i) << "\t" << (*vec)[i] << std::endl;
				}
			}
			std::cout << "\n/////////////////////////////\n";
			std::cout << "// Constant Properties\n";
			std::cout << "/////////////////////////////\n\n";
			for(auto cprop : (*cpMap)) // Iterate over constant properties
				std::cout << cprop.first << " = " << cprop.second << std::endl;
#endif
		}
		else{
			std::cout << " No properties table defined for this material\n";
		}
		std::cout << std::endl;
	}
}

void nDetMaterials::defineMaterials(){
	if(isInitialized) return;

	/////////////////////////////////////////////////////////////////
	// Visual Attributes
	/////////////////////////////////////////////////////////////////

	visAssembly = new G4VisAttributes();

	visSensitive = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75)); // grey
	visSensitive->SetForceSolid(true);

	visWindow = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0)); // cyan	

	visGrease = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // red

	visWrapping = new G4VisAttributes();
	visWrapping->SetColor(1, 0, 1, 0.5); // Alpha=50%
	visWrapping->SetForceSolid(true);

	visScint = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); // blue

	visShadow = new G4VisAttributes();
	visShadow->SetColor(0, 1, 0, 0.5); // green, Alpha=50%
	visShadow->SetForceSolid(true);	

	/////////////////////////////////////////////////////////////////
	// Elements & Molecular Components
	/////////////////////////////////////////////////////////////////

	// Elements
	fH = nist.searchForElement("H");
	fC = nist.searchForElement("C");
	fO = nist.searchForElement("O");
	fN = nist.searchForElement("N");
	fS = nist.searchForElement("S");
	fF = nist.searchForElement("F");
	fSi = nist.searchForElement("Si");
	fAl = nist.searchForElement("Al");
 	fY = nist.searchForElement("Y");
	fLa = nist.searchForElement("La");
	fBr = nist.searchForElement("Br");
	//fCs = nist.searchForElement("Cs");
	//fI = nist.searchForElement("I");
	fGd = nist.searchForElement("Gd");
	fGa = nist.searchForElement("Ga");
	fCe = nist.searchForElement("Ce");
	fZr = nist.searchForElement("Zr");
	fLu = nist.searchForElement("Lu");

	// Air
    fAir = nist.searchForMaterial("G4_AIR");
		G4double airPhotonEnergy[2] = {2.*eV,3.*eV};
		G4double airRefIndex[2] = {1.0003,1.0003};
		G4double airAbsorption[2] = {1*m,1*m};
		fAirMPT = new G4MaterialPropertiesTable();
		fAirMPT->AddProperty("RINDEX",airPhotonEnergy,airRefIndex,2);
		fAirMPT->AddProperty("ABSLENGTH",airPhotonEnergy,airAbsorption,2);
		std::cout<<"Getting air absorption length "<<fAirMPT->GetProperty("RINDEX")<<std::endl;
		fAir->SetMaterialPropertiesTable(fAirMPT);

	// Lab vacuum
	fVacuum = nist.searchForMaterial("G4_Galactic");

	/////////////////////////////////////////////////////////////////
	// Teflon (C2F4)n
	/////////////////////////////////////////////////////////////////

    fTeflon = nist.searchForMaterial("G4_TEFLON");

    G4double photonEnergy_teflon[9] = {1.607*eV, 1.743*eV, 1.908*eV, 2.108*eV, 2.354*eV, 2.664*eV, 3.070*eV, 3.621*eV, 4.413*eV};
    G4double reflectivity_teflon[9] = {0.514, 0.583, 0.656, 0.727, 0.789, 0.836, 0.868, 0.887, 0.892}; // https://www.osti.gov/servlets/purl/1184400 (1 layer)
    G4double efficiency_teflon[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    G4double absorption_teflon[9] =  {0.333*cm, 0.333*cm, 0.333*cm, 0.333*cm, 0.333*cm, 0.333*cm, 0.333*cm, 0.333*cm, 0.333*cm};
    G4double refIndex_teflon[9] = {1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315, 1.315};
    
    fTeflonMPT = new G4MaterialPropertiesTable();
    fTeflonMPT->AddProperty("REFLECTIVITY", photonEnergy_teflon, reflectivity_teflon, 9);
    fTeflonMPT->AddProperty("EFFICIENCY", photonEnergy_teflon, efficiency_teflon, 9);
    fTeflonMPT->AddProperty("RINDEX", photonEnergy_teflon, refIndex_teflon, 9);
    fTeflonMPT->AddProperty("ABSLENGTH", photonEnergy_teflon, absorption_teflon, 9);
    fTeflon->SetMaterialPropertiesTable(fTeflonMPT);

	/////////////////////////////////////////////////////////////////
	// Aluminum Oxide (Al2O3) (COME BACK AND CHANGE THESE NUMBERS LATER)
	/////////////////////////////////////////////////////////////////

    fAl2O3 = new G4Material("Al2O3", 3.95*g/cm3, 2);

	fAl2O3->AddElement(fAl, 2);
	fAl2O3->AddElement(fO, 3);

	// taken from Teflon above
    G4double photonEnergy_al2o3[9] = {1.607*eV, 1.743*eV, 1.908*eV, 2.108*eV, 2.354*eV, 2.664*eV, 3.070*eV, 3.621*eV, 4.413*eV};
    G4double reflectivity_al2o3[9] = {0.514, 0.583, 0.656, 0.727, 0.789, 0.836, 0.868, 0.887, 0.892};
    G4double efficiency_al2o3[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    G4double absorption_al2o3[9] =  {0.333*cm, 0.333*cm, 0.333*cm, 0.333*cm, 0.333*cm, 0.333*cm, 0.333*cm, 0.333*cm, 0.333*cm};

    G4double refIndex_al2o3[9] = {1.77, 1.77, 1.77, 1.77, 1.77, 1.77, 1.77, 1.77, 1.77};
    
    fAl2O3MPT = new G4MaterialPropertiesTable();
    fAl2O3MPT->AddProperty("REFLECTIVITY", photonEnergy_al2o3, reflectivity_al2o3, 9);
    fAl2O3MPT->AddProperty("EFFICIENCY", photonEnergy_al2o3, efficiency_al2o3, 9);
    fAl2O3MPT->AddProperty("RINDEX", photonEnergy_al2o3, refIndex_al2o3, 9);
    fAl2O3MPT->AddProperty("ABSLENGTH", photonEnergy_al2o3, absorption_al2o3, 9);
    fAl2O3->SetMaterialPropertiesTable(fAl2O3MPT);

	/////////////////////////////////////////////////////////////////
	// Silicone Optical Grease (C2H6OSi)n
	/////////////////////////////////////////////////////////////////

    fGrease = new G4Material("Grease", 1.06*g/cm3, 4);

    fGrease->AddElement(fC, 2);
    fGrease->AddElement(fH, 6);
    fGrease->AddElement(fO, 1);
    fGrease->AddElement(fSi, 1);

    /*const G4int nEntries_Grease = 5;
    G4double Photon_Energy[nEntries_Grease] = { 2.757*eV, 3.102*eV, 3.312*eV, 3.545*eV, 4.136*eV };
    G4double RefractiveIndex_grease[nEntries_Grease] = {1.468, 1.473, 1.477, 1.482, 1.496};
    G4double Absorption_grease[nEntries_Grease] = { 195*mm,  195*mm, 195*mm, 195*mm, 195*mm };
    fGreaseMPT=new G4MaterialPropertiesTable();

    fGreaseMPT->AddProperty("RINDEX", Photon_Energy, RefractiveIndex_grease, nEntries_Grease);
    fGreaseMPT->AddProperty("ABSLENGTH", Photon_Energy, Absorption_grease, nEntries_Grease);*/

    G4double greasePhotonEnergy[3] = { 2*eV, 3*eV, 4*eV};
    G4double refIndexGrease[3] = {1.465, 1.465, 1.465};
	G4double absorptionGrease[3] = {195*mm, 195*mm, 195*mm};

    fGreaseMPT = new G4MaterialPropertiesTable();
    fGreaseMPT->AddProperty("RINDEX", greasePhotonEnergy, refIndexGrease, 3);
    fGreaseMPT->AddProperty("ABSLENGTH", greasePhotonEnergy, absorptionGrease, 3);
    fGrease->SetMaterialPropertiesTable(fGreaseMPT);

	/*********************************************
	 NOTE: the NOA grease options are implemented
	  with the same material properties as ordinary
	  optical grease; only the index of refraction
	  is changed.
	 *********************************************/

	/////////////////////////////////////////////////////////////////
	// NOA-61 (Norland optical adhesive)
	/////////////////////////////////////////////////////////////////

    fNOA61 = new G4Material("NOA-61", 1.06*g/cm3, 4);

    fNOA61->AddElement(fC, 2);
    fNOA61->AddElement(fH, 6);
    fNOA61->AddElement(fO, 1);
    fNOA61->AddElement(fSi, 1);

    G4double NOA61PhotonEnergy[3] = { 2*eV, 3*eV, 4*eV};
    G4double refIndexNOA61[3] = {1.56, 1.56, 1.56};
    G4double absorptionNOA61[3] = {195*mm, 195*mm, 195*mm};

    fNOA61MPT = new G4MaterialPropertiesTable();
    fNOA61MPT->AddProperty("RINDEX", NOA61PhotonEnergy, refIndexNOA61, 3);
    fNOA61MPT->AddProperty("ABSLENGTH", NOA61PhotonEnergy, absorptionNOA61, 3);
    fNOA61->SetMaterialPropertiesTable(fNOA61MPT);


	/////////////////////////////////////////////////////////////////
	// NOA-68 (Norland optical adhesive)
	/////////////////////////////////////////////////////////////////

    fNOA68 = new G4Material("NOA-68", 1.06*g/cm3, 4);

    fNOA68->AddElement(fC, 2);
    fNOA68->AddElement(fH, 6);
    fNOA68->AddElement(fO, 1);
    fNOA68->AddElement(fSi, 1);

    G4double NOA68PhotonEnergy[3] = { 2*eV, 3*eV, 4*eV};
    G4double refIndexNOA68[3] = {1.54, 1.54, 1.54};
    G4double absorptionNOA68[3] = {195*mm, 195*mm, 195*mm};

    fNOA68MPT = new G4MaterialPropertiesTable();
    fNOA68MPT->AddProperty("RINDEX", NOA68PhotonEnergy, refIndexNOA68, 3);
    fNOA68MPT->AddProperty("ABSLENGTH", NOA68PhotonEnergy, absorptionNOA68, 3);
    fNOA68->SetMaterialPropertiesTable(fNOA68MPT);


	/////////////////////////////////////////////////////////////////
	// NOA-136 (Norland optical adhesive)
	/////////////////////////////////////////////////////////////////

    fNOA136 = new G4Material("NOA-136", 1.06*g/cm3, 4);

    fNOA136->AddElement(fC, 2);
    fNOA136->AddElement(fH, 6);
    fNOA136->AddElement(fO, 1);
    fNOA136->AddElement(fSi, 1);

    G4double NOA136PhotonEnergy[3] = { 2*eV, 3*eV, 4*eV};
    G4double refIndexNOA136[3] = {1.36, 1.36, 1.36};
    G4double absorptionNOA136[3] = {195*mm, 195*mm, 195*mm};

    fNOA136MPT = new G4MaterialPropertiesTable();
    fNOA136MPT->AddProperty("RINDEX", NOA136PhotonEnergy, refIndexNOA136, 3);
    fNOA136MPT->AddProperty("ABSLENGTH", NOA136PhotonEnergy, absorptionNOA136, 3);
    fNOA136->SetMaterialPropertiesTable(fNOA136MPT);


	/////////////////////////////////////////////////////////////////
	// NOA-170 (Norland optical adhesive)
	/////////////////////////////////////////////////////////////////

    fNOA170 = new G4Material("NOA-170", 1.06*g/cm3, 4);

    fNOA170->AddElement(fC, 2);
    fNOA170->AddElement(fH, 6);
    fNOA170->AddElement(fO, 1);
    fNOA170->AddElement(fSi, 1);

    G4double NOA170PhotonEnergy[3] = { 2*eV, 3*eV, 4*eV};
    G4double refIndexNOA170[3] = {1.70, 1.70, 1.70};
    G4double absorptionNOA170[3] = {195*mm, 195*mm, 195*mm};

    fNOA170MPT = new G4MaterialPropertiesTable();
    fNOA170MPT->AddProperty("RINDEX", NOA170PhotonEnergy, refIndexNOA170, 3);
    fNOA170MPT->AddProperty("ABSLENGTH", NOA170PhotonEnergy, absorptionNOA170, 3);
    fNOA170->SetMaterialPropertiesTable(fNOA170MPT);


	/////////////////////////////////////////////////////////////////
	// ABS plastic ((C8H8)(C4H6)(C3H3N))
	/////////////////////////////////////////////////////////////////

	fABS = new G4Material("ABS Plastic", 1.07*g/cm3, 3);

	fABS->AddElement(fC,15);
	fABS->AddElement(fH,17);
	fABS->AddElement(fN,1);

    G4double ABSPhotonEnergy[3] = { 2*eV, 3*eV, 4*eV};
    G4double refIndexABS[3] = {1.54, 1.54, 1.54};
    G4double absorptionABS[3] = {195*mm, 195*mm, 195*mm};

    fABSMPT = new G4MaterialPropertiesTable();
    fABSMPT->AddProperty("RINDEX", ABSPhotonEnergy, refIndexABS, 3);
    fABSMPT->AddProperty("ABSLENGTH", ABSPhotonEnergy, absorptionABS, 3);
    fABS->SetMaterialPropertiesTable(fABSMPT);

	/////////////////////////////////////////////////////////////////
	// Quartz (SiO2)
	/////////////////////////////////////////////////////////////////

    fSiO2 = nist.searchForMaterial("G4_SILICON_DIOXIDE");

    //optical properties of SiO2 - fused silica or fused quartz
    G4double PhotonEnergy[5] = { 2.484*eV, 2.615*eV, 2.760*eV, 2.922*eV, 3.105*eV };
    /*G4double RefractiveIndex_SiO2[5] = { 1.54, 1.54, 1.54, 1.54, 1.54 };
    G4double Absorption_SiO2[5] = {125.*cm, 123.5*cm, 122.*cm, 121.*cm, 120.*cm};

    fSiO2MPT = new G4MaterialPropertiesTable();
    fSiO2MPT->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex_SiO2, 5);
    fSiO2MPT->AddProperty("ABSLENGTH", PhotonEnergy, Absorption_SiO2, 5);*/
    
    G4double refIndexSiO2[3] = {1.458, 1.458, 1.458};
    G4double absorptionSiO2[3] = {125*cm, 125*cm, 125*cm};

    fSiO2MPT = new G4MaterialPropertiesTable();
    fSiO2MPT->AddProperty("RINDEX", greasePhotonEnergy, refIndexSiO2, 3);
    fSiO2MPT->AddProperty("ABSLENGTH", greasePhotonEnergy, absorptionSiO2, 3);   
    fSiO2->SetMaterialPropertiesTable(fSiO2MPT);


	/////////////////////////////////////////////////////////////////
	// Silicon (Si)
	/////////////////////////////////////////////////////////////////

    fSilicon = nist.searchForMaterial("G4_Si");

    // optical properties,    
    fSiliconMPT = new G4MaterialPropertiesTable();
    
    /*G4double RefractiveReal_Si[nEntries_Sil] = { 4.293, 4.453, 4.676, 5.008, 5.587 };
    G4double RefractiveImg_Si[nEntries_Sil] = { 0.045, 0.060, 0.091, 0.150, 0.303 };
    fSiliconMPT->AddProperty("REALRINDEX", PhotonEnergy, RefractiveReal_Si, nEntries_Sil);
    fSiliconMPT->AddProperty("IMAGINARYRINDEX", PhotonEnergy, RefractiveImg_Si, nEntries_Sil);*/
    
    // Ideal detector
    G4double EfficiencyIndex_Si[5] = { 1., 1., 1., 1., 1. };
    G4double Reflective_Si[5] = { 0., 0., 0., 0., 0.};
    
    // Non-ideal
    //G4double EfficiencyIndex_Si[nEntries_Sil] = { 0.37, 0.42, 0.39, 0.36, 0.32 };        
    //G4double Reflective_Si[nEntries_Sil] = { 0.49, 0.45, 0.42, 0.40, 0.39};

    fSiliconMPT->AddProperty("EFFICIENCY", PhotonEnergy, EfficiencyIndex_Si, 5);
    fSiliconMPT->AddProperty("REFLECTIVITY", PhotonEnergy, Reflective_Si, 5);

    fSilicon->SetMaterialPropertiesTable(fSiliconMPT);

	/////////////////////////////////////////////////////////////////
	// ACRYLIC (C5O2H8)n -- CRT
	/////////////////////////////////////////////////////////////////

    fAcrylic = nist.searchForMaterial("G4_GLASS_PLATE");

	// Photon energy (eV)
	G4double ENERGY_ACRYLIC[11] = {6.1992*eV, 4.95936*eV, 4.1328*eV, 3.5424*eV, 3.0996*eV, 2.7552*eV, 2.47968*eV, 2.25426*eV, 2.0664*eV, 1.90745*eV, 1.7712*eV};
	
	// Refractive index
	G4double RINDEX_ACRYLIC[11] = {1.57237, 1.54724, 1.5286, 1.51533, 1.50629, 1.50033, 1.49633, 1.49313, 1.4896, 1.48461, 1.47702};
	
	G4MaterialPropertiesTable *MPT_Acrylic = new G4MaterialPropertiesTable();
	MPT_Acrylic->AddProperty("RINDEX", ENERGY_ACRYLIC, RINDEX_ACRYLIC, 11);

	// Photon energy (eV)
	G4double energy2[25] = {12.3984*eV, 5.0292*eV, 4.75755*eV, 4.69897*eV, 4.66072*eV, 4.61377*eV, 4.56776*eV, 4.5316*eV, 4.48721*eV, 
	                        4.43507*eV, 4.13105*eV, 3.87258*eV, 3.64454*eV, 3.43671*eV, 3.25129*eV, 3.10158*eV, 2.94219*eV, 2.81212*eV, 
	                        2.69307*eV, 2.58078*eV, 2.47479*eV, 2.38212*eV, 2.29384*eV, 2.21614*eV, 1.7712*eV};
	                       
	// Acrylic absorption length (mm)
	G4double abslength[25] = {0, 0, 1.02102, 1.28345, 1.8604, 2.44271, 3.20801, 4.15751, 5.55214, 6.99127, 13.0334, 19.3961, 27.6547, 
	                          32.87, 34.1447, 35.5173, 34.1447, 32.87, 32.87, 34.1447, 34.1447, 35.5173, 35.5173, 34.1447, 33.7719};
	
	for(int i=0;i<(sizeof(abslength)/sizeof(abslength[0]));i++)
		abslength[i] *=10.;
	MPT_Acrylic->AddProperty("ABSLENGTH", energy2, abslength, 25);

	fAcrylic->SetMaterialPropertiesTable(MPT_Acrylic);


	/////////////////////////////////////////////////////////////////
	// Natural Aluminum
	/////////////////////////////////////////////////////////////////

	fAluminum = nist.searchForMaterial("G4_Al");

	double AlEnergies[3] = {2.0*eV, 3.0*eV, 4.0*eV};	                       
	double AlRefIndex[3] = {0.86, 0.50, 0.28};
	double AlAbsLength[3] = {1*mm, 1*mm, 1*mm};

	// H. EHRENREICH et al. Phys. Rev. 132, 5 (1963)
	double aluminumEnergy[8] = {0.000*eV, 0.697*eV, 1.172*eV, 1.350*eV, 1.504*eV, 2.305*eV, 3.174*eV, 4.000*eV};
	double aluminumReflectivity[8] = {1.000, 0.977, 0.950, 0.911, 0.869, 0.914, 0.921, 0.922};

	fAluminumMPT = new G4MaterialPropertiesTable();
	fAluminumMPT->AddProperty("RINDEX", AlEnergies, AlRefIndex, 3);
	fAluminumMPT->AddProperty("ABSLENGTH", AlEnergies, AlAbsLength, 3);
	fAluminumMPT->AddProperty("REFLECTIVITY", aluminumEnergy, aluminumReflectivity, 8);
	fAluminum->SetMaterialPropertiesTable(fAluminumMPT);

	/////////////////////////////////////////////////////////////////
	// Optical Surfaces
	/////////////////////////////////////////////////////////////////

	// Teflon	
    fTeflonOpSurf = new G4OpticalSurface("TeflonSurface", glisur, ground, dielectric_metal, 0.1); // polish level

    fTeflonOpSurf->SetFinish(ground);
    //fTeflonOpSurf->SetType(dielectric_metal);
    fTeflonOpSurf->SetType(dielectric_dielectric);
    fTeflonOpSurf->SetMaterialPropertiesTable(fTeflonMPT);

	// Aluminum Oxide (Al2O3) (COME BACK AND CHANGE THESE NUMBERS LATER)
	fAl2O3OpSurf = new G4OpticalSurface("Al2O3Surface", glisur, ground, dielectric_metal, 0.1);
	fAl2O3OpSurf->SetFinish(ground);
	fAl2O3OpSurf->SetType(dielectric_metal);
	fAl2O3OpSurf->SetMaterialPropertiesTable(fAl2O3MPT);

    // Silicon
    fSiliconOpSurf = new G4OpticalSurface("SiPMSurface");
    fSiliconOpSurf->SetType(dielectric_metal);
    fSiliconOpSurf->SetFinish(polished);
    fSiliconOpSurf->SetModel(glisur);
    fSiliconOpSurf->SetMaterialPropertiesTable(fSiliconMPT);

	// Aluminized mylar
    G4Material *Mylar = nist.searchForMaterial("G4_MYLAR");

    fMylar = new G4Material("AluminizedMylar", 1.39*g/cm3, 2);
    fMylar->AddMaterial(Mylar, 0.8);
    fMylar->AddMaterial(fAluminum, 0.2);

    //G4double RefractiveReal_Mylar[5] = {0.81257,0.72122,0.63324,0.55571,0.48787};
    //G4double RefractiveImg_Mylar[5] = {6.0481,5.7556,5.4544,5.1464,4.8355};

	//fMylarMPT = new G4MaterialPropertiesTable();
	//fMylarMPT->AddProperty("REALRINDEX", PhotonEnergy, RefractiveReal_Mylar, 5);
	//fMylarMPT->AddProperty("IMAGINARYRINDEX", PhotonEnergy, RefractiveImg_Mylar, 5);

    fMylarOpSurf = new G4OpticalSurface("MylarSurface");
	fMylarOpSurf->SetType(dielectric_metal);
	fMylarOpSurf->SetFinish(polished); // dielectric_metal only allows polished or ground. Polished dielectric_metal uses only reflectivity or absorption.
	fMylarOpSurf->SetModel(glisur);
	fMylarOpSurf->SetMaterialPropertiesTable(fAluminumMPT);

	// 100% reflectivity
	double perfectEfficiency[3] = {0.0, 0.0, 0.0};
	double perfectReflectivity[3] = {1.0, 1.0, 1.0};
	double perfectSpecularSpike[3] = {1.0, 1.0, 1.0};
	
	// ESR (98% reflectivity)
	double esrReflectivity[3] = {0.98, 0.98, 0.98};
    
    fEsrMPT = new G4MaterialPropertiesTable();
	fEsrMPT->AddProperty("EFFICIENCY", AlEnergies, perfectEfficiency, 3);
	fEsrMPT->AddProperty("REFLECTIVITY", AlEnergies, esrReflectivity, 3);    
    
    // ESR film (built in look-up-table)
    fEsrOpSurf = new G4OpticalSurface("EsrSurface");
    fEsrOpSurf->SetType(dielectric_LUT);
    fEsrOpSurf->SetModel(LUT);    
    //fEsrOpSurf->SetFinish(polishedvm2000air);
    fEsrOpSurf->SetFinish(polishedvm2000glue);
    fEsrOpSurf->SetMaterialPropertiesTable(fEsrMPT);

	fPerfectMPT = new G4MaterialPropertiesTable();
	fPerfectMPT->AddProperty("EFFICIENCY", AlEnergies, perfectEfficiency, 3);
	fPerfectMPT->AddProperty("REFLECTIVITY", AlEnergies, perfectReflectivity, 3);
	fPerfectMPT->AddProperty("SPECULARSPIKECONSTANT", AlEnergies, perfectSpecularSpike, 3);

	fPerfectOpSurf = new G4OpticalSurface("PerfectReflector");
	fPerfectOpSurf->SetType(dielectric_metal);
	fPerfectOpSurf->SetFinish(polished);
	fPerfectOpSurf->SetModel(glisur);
	fPerfectOpSurf->SetMaterialPropertiesTable(fPerfectMPT);

	fGreaseOpSurf = new G4OpticalSurface("GreaseSurface");
	fGreaseOpSurf->SetType(dielectric_dielectric);
	fGreaseOpSurf->SetFinish(ground);
	fGreaseOpSurf->SetModel(unified); // Defaults to Lambertian reflection (i.e. rough surface) --CRT
	fGreaseOpSurf->SetMaterialPropertiesTable(fGreaseMPT);	

	fNOA61OpSurf = new G4OpticalSurface("NOA61Surface");
	fNOA61OpSurf->SetType(dielectric_dielectric);
	fNOA61OpSurf->SetFinish(ground);
	fNOA61OpSurf->SetModel(unified);
	fNOA61OpSurf->SetMaterialPropertiesTable(fNOA61MPT);	

	fNOA68OpSurf = new G4OpticalSurface("NOA68Surface");
	fNOA68OpSurf->SetType(dielectric_dielectric);
	fNOA68OpSurf->SetFinish(ground);
	fNOA68OpSurf->SetModel(unified);
	fNOA68OpSurf->SetMaterialPropertiesTable(fNOA68MPT);

	fNOA136OpSurf = new G4OpticalSurface("NOA136Surface");
	fNOA136OpSurf->SetType(dielectric_dielectric);
	fNOA136OpSurf->SetFinish(ground);
	fNOA136OpSurf->SetModel(unified);
	fNOA136OpSurf->SetMaterialPropertiesTable(fNOA136MPT);

	fNOA170OpSurf = new G4OpticalSurface("NOA170Surface");
	fNOA170OpSurf->SetType(dielectric_dielectric);
	fNOA170OpSurf->SetFinish(ground);
	fNOA170OpSurf->SetModel(unified);
	fNOA170OpSurf->SetMaterialPropertiesTable(fNOA170MPT);

	fAirOpSurf = new G4OpticalSurface("AirSurface");
	//fAirOpSurf->SetType;
	//fAirOpSurf->SetFinish(polished);
	//fAirOpSurf->SetModel(glisur);
	fAirOpSurf->SetMaterialPropertiesTable(fAirMPT);

		fAluminumOpSurf = new G4OpticalSurface("PerfectReflector");
	fAluminumOpSurf->SetType(dielectric_metal);
	fAluminumOpSurf->SetFinish(polished);
	fAluminumOpSurf->SetModel(glisur);
	fAluminumOpSurf->SetMaterialPropertiesTable(fAluminumMPT);
	

	isInitialized = true;
}

void nDetMaterials::defineScintillators(){
	if(scintsAreDefined){
		delete fEJ200;
		delete fEJ276;
		delete fOGS;
		delete fEJ200MPT;
		delete fEJ276MPT;
		delete fOGSMPT;
		delete fYSO;
		delete fYSOMPT;
		delete fLYSO;
		delete fLYSOMPT;
		delete fLaBr3;
		delete fLaBr3MPT;
		delete fCeBr3;
		delete fCeBr3MPT;
		//delete fCsI;
		//delete fCsIMPT;
	}

	/////////////////////////////////////////////////////////////////
	// EJ200 N(H)=52.4%, N(C)=47.6%
	/////////////////////////////////////////////////////////////////

    fEJ200 = new G4Material("EJ200", 1.023*g/cm3, 2);
    fEJ200->AddElement(fH, 0.08457);
    fEJ200->AddElement(fC, 0.91543);

	G4double photonEnergy_Ej200[44] = {2.004*eV, 2.058*eV, 2.112*eV, 2.166*eV, 2.220*eV, 2.274*eV, 2.328*eV, 2.382*eV, 2.436*eV, 2.490*eV, 
		                               2.517*eV, 2.552*eV, 2.585*eV, 2.613*eV, 2.635*eV, 2.656*eV, 2.686*eV, 2.720*eV, 2.749*eV, 2.772*eV, 
		                               2.791*eV, 2.809*eV, 2.826*eV, 2.842*eV, 2.861*eV, 2.884*eV, 2.919*eV, 2.946*eV, 2.954*eV, 2.961*eV, 
		                               2.967*eV, 2.974*eV, 2.981*eV, 2.987*eV, 2.994*eV, 3.001*eV, 3.009*eV, 3.018*eV, 3.029*eV, 3.041*eV, 
		                               3.056*eV, 3.083*eV, 3.137*eV, 3.191*eV};

	G4double ScintilFast_EJ200[44] = {0.000, 0.001, 0.001, 0.002, 0.003, 0.006, 0.010, 0.018, 0.033, 0.060, 
		                              0.084, 0.122, 0.175, 0.234, 0.294, 0.356, 0.416, 0.473, 0.533, 0.594, 
		                              0.657, 0.720, 0.784, 0.846, 0.903, 0.962, 1.000, 0.917, 0.857, 0.798, 
		                              0.732, 0.669, 0.604, 0.542, 0.480, 0.422, 0.359, 0.297, 0.237, 0.170, 
		                              0.105, 0.028, 0.004, 0.000};
                                  
	G4double photonEnergy_Ej200_2[2] = {2.004*eV, 3.191*eV};
	G4double RefIndex_EJ200[2] = {1.580, 1.580};
	G4double Absorption_EJ200[2] = {400*cm, 400*cm};

    fEJ200MPT = new G4MaterialPropertiesTable();
    fEJ200MPT->AddProperty("RINDEX", photonEnergy_Ej200_2, RefIndex_EJ200, 2);
    fEJ200MPT->AddProperty("ABSLENGTH", photonEnergy_Ej200_2, Absorption_EJ200, 2);
    fEJ200MPT->AddProperty("FASTCOMPONENT", photonEnergy_Ej200, ScintilFast_EJ200, 44);

    //fEJ200MPT->AddConstProperty("SCINTILLATIONYIELD", 0.64*17400/MeV); // 64% of Anthracene
    fEJ200MPT->AddConstProperty("SCINTILLATIONYIELD", 10000/MeV); // Scintillation efficiency as per Eljen specs
    fEJ200MPT->AddConstProperty("RESOLUTIONSCALE", 1.0); // Intrinsic resolution

    //fEJ200MPT->AddConstProperty("RISETIMECONSTANT", 0.9*ns); Geant4 10.1 TODO
    fEJ200MPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 0.9*ns);
    fEJ200MPT->AddConstProperty("FASTTIMECONSTANT", 2.1*ns);
    fEJ200MPT->AddConstProperty("YIELDRATIO",1);// the strength of the fast component as a function of total scintillation yield

	std::cout << "nDetConstruction: Photon yield is set to " << lightYieldScale << "x scale\n";
    G4double pEF = lightYieldScale; 
    G4double pSF = pEF * 1.35;

    //light yield - data taken form V.V. Verbinski et al, Nucl. Instrum. & Meth. 65 (1968) 8-25
	G4double particleEnergy[36] = {1E-3, 0.10*MeV, 0.13*MeV, 0.17*MeV, 0.20*MeV, 0.24*MeV, 0.30*MeV, 0.34*MeV, 0.40*MeV, 
		                           0.48*MeV, 0.60*MeV, 0.72*MeV, 0.84*MeV, 1.00*MeV, 1.30*MeV, 1.70*MeV, 2.00*MeV, 2.40*MeV, 
		                           3.00*MeV, 3.40*MeV, 4.00*MeV, 4.80*MeV, 6.00*MeV, 7.20*MeV, 8.40*MeV, 10.00*MeV, 11.00*MeV, 
		                           12.00*MeV, 13.00*MeV, 14.00*MeV, 15.00*MeV, 16.00*MeV, 17.00*MeV, 18.00*MeV, 19.00*MeV, 20.00*MeV};
		                           
	G4double electronYield[36] = {0.0E+00*pEF, 1.0E+03*pEF, 1.3E+03*pEF, 1.7E+03*pEF, 2.0E+03*pEF, 2.4E+03*pEF, 3.0E+03*pEF, 3.4E+03*pEF, 4.0E+03*pEF, 
		                          4.8E+03*pEF, 6.0E+03*pEF, 7.2E+03*pEF, 8.4E+03*pEF, 1.0E+04*pEF, 1.3E+04*pEF, 1.7E+04*pEF, 2.0E+04*pEF, 2.4E+04*pEF, 
		                          3.0E+04*pEF, 3.4E+04*pEF, 4.0E+04*pEF, 4.8E+04*pEF, 6.0E+04*pEF, 7.2E+04*pEF, 8.4E+04*pEF, 1.0E+05*pEF, 1.1E+05*pEF, 
		                          1.2E+05*pEF, 1.3E+05*pEF, 1.4E+05*pEF, 1.5E+05*pEF, 1.6E+05*pEF, 1.7E+05*pEF, 1.8E+05*pEF, 1.9E+05*pEF, 2.0E+05*pEF};

    fEJ200MPT->AddProperty("ELECTRONSCINTILLATIONYIELD", particleEnergy, electronYield, 36)->SetSpline(true);
	G4double protonYield[36] = {0.6*pSF, 67.1*pSF, 88.6*pSF, 120.7*pSF, 146.5*pSF, 183.8*pSF, 246.0*pSF, 290.0*pSF, 365.0*pSF, 
		                        483.0*pSF, 678.0*pSF, 910.0*pSF, 1175.0*pSF, 1562.0*pSF, 2385.0*pSF, 3660.0*pSF, 4725.0*pSF, 6250.0*pSF, 
		                        8660.0*pSF, 10420.0*pSF, 13270.0*pSF, 17180.0*pSF, 23100.0*pSF, 29500.0*pSF, 36200.0*pSF, 45500.0*pSF, 51826.7*pSF, 
		                        58313.7*pSF, 65047.2*pSF, 72027.4*pSF, 79254.2*pSF, 86727.6*pSF, 94447.6*pSF, 102414.2*pSF, 110627.4*pSF, 119087.2*pSF};

    fEJ200MPT->AddProperty("PROTONSCINTILLATIONYIELD", particleEnergy, protonYield, 36)->SetSpline(true);
    fEJ200MPT->AddProperty("ALPHASCINTILLATIONYIELD", particleEnergy, protonYield, 36)->SetSpline(true);
	fEJ200MPT->AddProperty("TRITONSCINTILLATIONYIELD", particleEnergy, protonYield, 36)->SetSpline(true);
    fEJ200MPT->AddProperty("DEUTERONSCINTILLATIONYIELD", particleEnergy, protonYield, 36)->SetSpline(true);
	fEJ200MPT->AddProperty("TRITONSCINTILLATIONYIELD", particleEnergy, protonYield, 36)->SetSpline(true);

	G4double ionYield[36] = {0.2*pEF, 10.4*pEF, 12.7*pEF, 15.7*pEF, 17.9*pEF, 20.8*pEF, 25.1*pEF, 27.9*pEF, 31.9*pEF, 
		                     36.8*pEF, 43.6*pEF, 50.2*pEF, 56.9*pEF, 65.7*pEF, 81.3*pEF, 101.6*pEF, 116.5*pEF, 136.3*pEF, 
		                     166.2*pEF, 187.1*pEF, 218.6*pEF, 260.5*pEF, 323.5*pEF, 387.5*pEF, 451.5*pEF, 539.9*pEF, 595.5*pEF, 
		                     651.8*pEF, 708.7*pEF, 766.2*pEF, 824.2*pEF, 882.9*pEF, 942.2*pEF, 1002.1*pEF, 1062.6*pEF, 1123.7*pEF}; 
                         
    fEJ200MPT->AddProperty("IONSCINTILLATIONYIELD", particleEnergy, ionYield, 36)->SetSpline(true);

    fEJ200->SetMaterialPropertiesTable(fEJ200MPT);

	/////////////////////////////////////////////////////////////////
	// EJ276 N(H)=48.1%, N(C)=51.9%
	/////////////////////////////////////////////////////////////////

    fEJ276 = new G4Material("EJ276", 1.096*g/cm3, 2);
    fEJ276->AddElement(fH, 0.07216);
    fEJ276->AddElement(fC, 0.92784);

	G4double photonEnergy_Ej276[36] = {3.131*eV, 3.087*eV, 3.060*eV, 3.044*eV, 3.029*eV, 3.017*eV, 3.010*eV, 3.001*eV, 2.993*eV, 2.984*eV, 
		                               2.976*eV, 2.967*eV, 2.959*eV, 2.950*eV, 2.941*eV, 2.910*eV, 2.857*eV, 2.838*eV, 2.821*eV, 2.802*eV, 
		                               2.784*eV, 2.764*eV, 2.739*eV, 2.705*eV, 2.671*eV, 2.646*eV, 2.625*eV, 2.599*eV, 2.567*eV, 2.533*eV, 
		                               2.500*eV, 2.468*eV, 2.437*eV, 2.406*eV, 2.377*eV, 2.350*eV};

	G4double ScintilFast_EJ276[36] = {0.000, 0.010, 0.088, 0.157, 0.225, 0.293, 0.354, 0.415, 0.492, 0.570, 
	                                  0.649, 0.730, 0.807, 0.882, 0.934, 1.000, 0.890, 0.826, 0.761, 0.692, 
	                                  0.629, 0.569, 0.509, 0.445, 0.388, 0.326, 0.263, 0.200, 0.144, 0.100, 
	                                  0.068, 0.038, 0.024, 0.012, 0.002, 0.000};

    fEJ276MPT = new G4MaterialPropertiesTable();
    fEJ276MPT->AddProperty("RINDEX", photonEnergy_Ej200_2, RefIndex_EJ200, 2);
    fEJ276MPT->AddProperty("ABSLENGTH", photonEnergy_Ej200_2, Absorption_EJ200, 2);
    fEJ276MPT->AddProperty("FASTCOMPONENT", photonEnergy_Ej276, ScintilFast_EJ276, 36);

    fEJ276MPT->AddConstProperty("SCINTILLATIONYIELD", 8600/MeV); // Scintillation efficiency as per Eljen specs
    fEJ276MPT->AddConstProperty("RESOLUTIONSCALE", 1.0); // Intrinsic resolution

    fEJ276MPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 0.9*ns);
    fEJ276MPT->AddConstProperty("FASTTIMECONSTANT", 2.1*ns);
    fEJ276MPT->AddConstProperty("YIELDRATIO",1);// the strength of the fast component as a function of total scintillation yield

	G4double electronYield_EJ276[36];
	G4double protonYield_EJ276[36];
	G4double ionYield_EJ276[36];

	// Produce the scaled light-yield for EJ276 (scaled down from EJ200 by 14%).
	for(size_t i = 0; i < 36; i++){
		electronYield_EJ276[i] = 0.86 * electronYield[i];
		protonYield_EJ276[i] = 0.86 * protonYield[i];
		ionYield_EJ276[i] = 0.86 * ionYield[i];
	}

    fEJ276MPT->AddProperty("ELECTRONSCINTILLATIONYIELD", particleEnergy, electronYield_EJ276, 36)->SetSpline(true);
    fEJ276MPT->AddProperty("PROTONSCINTILLATIONYIELD", particleEnergy, protonYield_EJ276, 36)->SetSpline(true);
    fEJ276MPT->AddProperty("IONSCINTILLATIONYIELD", particleEnergy, ionYield_EJ276, 36)->SetSpline(true);

	fEJ276->SetMaterialPropertiesTable(fEJ276MPT);


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// OGS (Data from "Organic glass scintillator (OGS) property comparisons to Stilbene, EJ-276 and BC-404", https://www.sciencedirect.com/science/article/pii/S0168900221007634)
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    fOGS = new G4Material("OGS", 1.096*g/cm3, 2);
    fOGS->AddElement(fH, 0.07216);
    fOGS->AddElement(fC, 0.92784);


	// From OGS paper
	G4double photonEnergy_OGS[27] = {2.10*eV, 2.16719*eV, 2.21073*eV, 2.26614*eV, 2.33082*eV, 2.43387*eV, 2.53619*eV, 2.59857*eV, 2.64731*eV, 2.70073*eV, 2.74456*eV, 2.81126*eV, 2.84578*eV, 2.86505*eV, 2.87811*eV, 2.89972*eV, 2.91841*eV, 2.94258*eV, 2.97066*eV, 2.99193*eV, 3.02781*eV, 3.06821*eV, 3.11734*eV, 3.17601*eV, 3.23693*eV, 3.29593*eV, 3.35*eV};
	G4double ScintilFast_OGS[27] = {0.000, 0.00053091, 0.00074703, 0.00101013, 0.00896905, 0.0477394, 0.132478, 0.213213, 0.343732, 0.512598, 0.581748, 0.811983, 0.969273, 1, 0.984702, 0.873582, 0.697279, 0.448147, 0.187523, 0.0725652, 0.0228196, 0.00758775, 0.00387609, 0.00401704, 0.00415799, 0.00428954, 0.000};

    fOGSMPT = new G4MaterialPropertiesTable();
    fOGSMPT->AddProperty("RINDEX", photonEnergy_Ej200_2, RefIndex_EJ200, 2);
    fOGSMPT->AddProperty("ABSLENGTH", photonEnergy_Ej200_2, Absorption_EJ200, 2);
    fOGSMPT->AddProperty("FASTCOMPONENT", photonEnergy_OGS, ScintilFast_OGS, 27);

    fOGSMPT->AddConstProperty("SCINTILLATIONYIELD", 15800/MeV); // Scintillation efficiency as per Eljen specs
    fOGSMPT->AddConstProperty("RESOLUTIONSCALE", 1.0); // Intrinsic resolution

    fOGSMPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 0.82*ns);
    fOGSMPT->AddConstProperty("FASTTIMECONSTANT", 2.2*ns);
    fOGSMPT->AddConstProperty("YIELDRATIO",1);// the strength of the fast component as a function of total scintillation yield

	G4double electronYield_OGS[36];
	G4double protonYield_OGS[36];
	G4double ionYield_OGS[36];

	// Produce the scaled light-yield for OGS (scaled up from EJ276).
	for(size_t i = 0; i < 36; i++){
		electronYield_OGS[i] = 2.8 * electronYield_EJ276[i];
		protonYield_OGS[i] = 2.8 * protonYield_EJ276[i];
		ionYield_OGS[i] = 2.8 * ionYield_EJ276[i];
	}

    fOGSMPT->AddProperty("ELECTRONSCINTILLATIONYIELD", particleEnergy, electronYield_OGS, 36)->SetSpline(true);
    fOGSMPT->AddProperty("PROTONSCINTILLATIONYIELD", particleEnergy, protonYield_OGS, 36)->SetSpline(true);
    fOGSMPT->AddProperty("IONSCINTILLATIONYIELD", particleEnergy, ionYield_OGS, 36)->SetSpline(true);
	fOGSMPT->AddProperty("DEUTERONSCINTILLATIONYIELD", particleEnergy, protonYield_OGS, 36)->SetSpline(true);
	fOGSMPT->AddProperty("TRITONSCINTILLATIONYIELD", particleEnergy, protonYield_OGS, 36)->SetSpline(true);

	fOGS->SetMaterialPropertiesTable(fOGSMPT);


	/////////////////////////////////////////////////////////////////
	// YSO Y2SiO5
	/////////////////////////////////////////////////////////////////

    fYSO = new G4Material("YSO", 4.5*g/cm3, 3);
    fYSO->AddElement(fY, 0.25);
    fYSO->AddElement(fSi, 0.125);
	fYSO->AddElement(fO, 0.625);

	G4double photonEnergy_YSO[44] = {2.004*eV, 2.058*eV, 2.112*eV, 2.166*eV, 2.220*eV, 2.274*eV, 2.328*eV, 2.382*eV, 2.436*eV, 2.490*eV, 
		                               2.517*eV, 2.552*eV, 2.585*eV, 2.613*eV, 2.635*eV, 2.656*eV, 2.686*eV, 2.720*eV, 2.749*eV, 2.772*eV, 
		                               2.791*eV, 2.809*eV, 2.826*eV, 2.842*eV, 2.861*eV, 2.884*eV, 2.919*eV, 2.946*eV, 2.954*eV, 2.961*eV, 
		                               2.967*eV, 2.974*eV, 2.981*eV, 2.987*eV, 2.994*eV, 3.001*eV, 3.009*eV, 3.018*eV, 3.029*eV, 3.041*eV, 
		                               3.056*eV, 3.083*eV, 3.137*eV, 3.191*eV};

	G4double ScintilFast_YSO[44] = {0.000, 0.001, 0.001, 0.002, 0.003, 0.006, 0.010, 0.018, 0.033, 0.060, 
		                              0.084, 0.122, 0.175, 0.234, 0.294, 0.356, 0.416, 0.473, 0.533, 0.594, 
		                              0.657, 0.720, 0.784, 0.846, 0.903, 0.962, 1.000, 0.917, 0.857, 0.798, 
		                              0.732, 0.669, 0.604, 0.542, 0.480, 0.422, 0.359, 0.297, 0.237, 0.170, 
		                              0.105, 0.028, 0.004, 0.000};

	//G4double photonEnergy_YSO[13] = {16.8179*keV, 23.0707*keV, 30.8922*keV, 50.0014*keV, 59.1033*keV, 80.4609*keV, 122.4170*keV, 280.3698*keV, 504.5552*keV, 659.8267*keV, 830.5057*keV, 1188.7809*keV, 1316.9683*keV};

	//G4double ScintilFast_YSO[13] = {0.0482746, 0.0585062, 0.0648312, 0.0774812, 0.0748768, 0.0707841, 0.0776672, 0.0854804, 0.0878988, 0.0888289, 0.089387, 0.0880848, 0.0878988}; //This should be the efficiency curve..?

	//G4double electronYield_YSO[13] = {54.3455, 65.8639, 72.9843, 87.2251, 84.2932, 79.6859, 87.4346, 96.2304, 98.9529, 100, 100.628, 99.1623, 98.9529}; //I am not sure if this is the correct units

	G4double photonEnergy_YSO_2[2] = {16.8179*keV, 23.0707*keV};

	G4double RefIndex_YSO[2] = {1.80, 1.80}; //Data taken from https://www.advatech-uk.co.uk/yso_ce.html 
	G4double Absorption_YSO[2] = {2.57*cm, 2.57*cm}; // this is found in: Large size LSO:Ce and YSO:Ce scintillators for 50 MeV range /spl gamma/-ray detector

    fYSOMPT = new G4MaterialPropertiesTable();
    fYSOMPT->AddProperty("RINDEX", photonEnergy_YSO_2, RefIndex_YSO, 2);
    fYSOMPT->AddProperty("ABSLENGTH", photonEnergy_YSO_2, Absorption_YSO, 2);
    fYSOMPT->AddProperty("FASTCOMPONENT", photonEnergy_YSO, ScintilFast_YSO, 13);

    fYSOMPT->AddConstProperty("SCINTILLATIONYIELD", 24000/MeV); // Photon yield as found in paper above for 137Cs
    fYSOMPT->AddConstProperty("RESOLUTIONSCALE", 1.0); // Intrinsic resolution

    fYSOMPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 2.0*ns);
    fYSOMPT->AddConstProperty("FASTTIMECONSTANT", 50.0*ns);
    fYSOMPT->AddConstProperty("YIELDRATIO",1);// the strength of the fast component as a function of total scintillation yield

	G4double electronYield_YSO[36];
	G4double protonYield_YSO[36];
	G4double ionYield_YSO[36];

	// Produce the scaled light-yield for YAP (scaled down from EJ200 by 14%). This will need to be adjusted for YAP as it has a different light yield than EJ276.
	for(size_t i = 0; i < 36; i++){
		electronYield_YSO[i] = 0.86 * electronYield[i];
		protonYield_YSO[i] = 0.86 * protonYield[i];
		ionYield_YSO[i] = 0.86 * ionYield[i];
	}

	//TEMPORARY CHANGE
	/*
	G4double particleEnergy_YSO[130];
	G4double electronYield_YSO[130];
	G4double protonYield_YSO[130];
	G4double ionYield_YSO[130];
	std::ifstream fin[3];
	fin[0].open("lightYieldYSO_electron.dat");
	fin[1].open("lightYieldYSO_proton.dat");
	fin[2].open("lightYieldYSO_024.dat");
	if(!fin[0].good())
		std::cout<<"!!!!! ERROR LOADING LIGHT YIELD FILES !!!!!!"<<std::endl;
	else
		std::cout<<"Able to load light yield data for YSO"<<std::endl;
	double en, light;
	for(size_t i=0;i<130;i++){
		fin[0]>>en>>light;
		electronYield_YSO[i] = light*pEF;
		particleEnergy_YSO[i] = en*MeV;
		fin[1]>>en>>light;
		protonYield_YSO[i] = light*pEF;
		fin[2]>>en>>light;
		ionYield_YSO[i] = light*peF;
	}
	*/

	// G4double electronYield_YSO[36];
	// G4double protonYield_YSO[36];
	// G4double ionYield_YSO[36];

	// Produce the scaled light-yield for YSO (scaled down from EJ200 by 14%). This will need to be adjusted for YSO as it has a different light yield than EJ276.
	// for(size_t i = 0; i < 36; i++){
	// 	electronYield_YSO[i] = 0.86 * electronYield[i];
	// 	protonYield_YSO[i] = 0.86 * protonYield[i];
	// 	ionYield_YSO[i] = 0.86 * ionYield[i];
	// }

	// TEMPORARY CHANGE; CHANGE BACK TO YSO PARTICLE ENERGY LATER
    fYSOMPT->AddProperty("ELECTRONSCINTILLATIONYIELD", particleEnergy, electronYield_YSO, 130)->SetSpline(true);
    fYSOMPT->AddProperty("PROTONSCINTILLATIONYIELD", particleEnergy, protonYield_YSO, 130)->SetSpline(true);
    fYSOMPT->AddProperty("IONSCINTILLATIONYIELD", particleEnergy, ionYield_YSO, 130)->SetSpline(true);
	fYSOMPT->AddProperty("ALPHASCINTILLATIONYIELD", particleEnergy, protonYield_YSO, 130)->SetSpline(true);

	fYSO->SetMaterialPropertiesTable(fYSOMPT);
	

	/////////////////////////////////////////////////////////////////
	//  LYSO
	/////////////////////////////////////////////////////////////////

	fLYSO = new G4Material("LYSO", 7.25*g/cm3, 4);
  	fLYSO->AddElement(fY,0.0125);
  	fLYSO->AddElement(fLu,0.2375);
	fLYSO->AddElement(fSi,0.125);
	fLYSO->AddElement(fO,0.625);

   /* Old values, copy-pasted from somewhere
	G4double photonEnergy_LYSO[44] = {2.004*eV, 2.058*eV, 2.112*eV, 2.166*eV, 2.220*eV, 2.274*eV, 2.328*eV, 2.382*eV, 2.436*eV, 2.490*eV, 
		                               2.517*eV, 2.552*eV, 2.585*eV, 2.613*eV, 2.635*eV, 2.656*eV, 2.686*eV, 2.720*eV, 2.749*eV, 2.772*eV, 
		                               2.791*eV, 2.809*eV, 2.826*eV, 2.842*eV, 2.861*eV, 2.884*eV, 2.919*eV, 2.946*eV, 2.954*eV, 2.961*eV, 
		                               2.967*eV, 2.974*eV, 2.981*eV, 2.987*eV, 2.994*eV, 3.001*eV, 3.009*eV, 3.018*eV, 3.029*eV, 3.041*eV, 
		                               3.056*eV, 3.083*eV, 3.137*eV, 3.191*eV};
	G4double ScintilFast_LYSO[44] = {0.000, 0.001, 0.001, 0.002, 0.003, 0.006, 0.010, 0.018, 0.033, 0.060, 
		                              0.084, 0.122, 0.175, 0.234, 0.294, 0.356, 0.416, 0.473, 0.533, 0.594, 
		                              0.657, 0.720, 0.784, 0.846, 0.903, 0.962, 1.000, 0.917, 0.857, 0.798, 
		                              0.732, 0.669, 0.604, 0.542, 0.480, 0.422, 0.359, 0.297, 0.237, 0.170, 
		                              0.105, 0.028, 0.004, 0.000};
   */
									
	// Taken from LYSO paper
	G4double photonEnergy_LYSO[29] = {2.00*eV, 2.07314*eV, 2.09979*eV, 2.15677*eV, 2.23207*eV, 2.31647*eV, 2.41342*eV, 2.52101*eV, 2.59669*eV, 2.65531*eV, 2.70912*eV, 2.76647*eV, 2.80472*eV, 2.85371*eV, 2.90158*eV, 3.0023*eV, 3.05851*eV, 3.10698*eV, 3.15362*eV, 3.18091*eV, 3.22274*eV, 3.24945*eV, 3.28026*eV, 3.30607*eV, 3.33797*eV, 3.35894*eV, 3.41353*eV, 3.47401*eV, 3.55*eV};
	G4double ScintilFast_LYSO[29] = {0.000, 0.00054604, 0.00490535, 0.00935039, 0.0223922, 0.0429297, 0.0720368, 0.130025, 0.207184, 0.301415, 0.412744, 0.5583, 0.667456, 0.778768, 0.860128, 0.937286, 0.981218, 1.000, 0.97814, 0.920422, 0.789995, 0.650992, 0.505575, 0.370848, 0.258587, 0.154868, 0.0522575, 0.012753, 0.00};

	G4double photonEnergy_LYSO2[2] = {16.8179*keV, 23.0707*keV};
	G4double RefIndex_LYSO[2] = {1.82,1.82}; 
	G4double Absorption_LYSO[2] = {2.7*cm, 2.7*cm}; 

	fLYSOMPT = new G4MaterialPropertiesTable();
	fLYSOMPT->AddProperty("RINDEX", photonEnergy_LYSO2, RefIndex_LYSO, 2);
	fLYSOMPT->AddProperty("ABSLENGTH", photonEnergy_LYSO2, Absorption_LYSO, 2);
	fLYSOMPT->AddProperty("FASTCOMPONENT", photonEnergy_LYSO, ScintilFast_LYSO, 29);

	fLYSOMPT->AddConstProperty("SCINTILLATIONYIELD", 25000/MeV); 
	fLYSOMPT->AddConstProperty("RESOLUTIONSCALE", 1.0); // Intrinsic resolution

	fLYSOMPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 2.0*ns);
	fLYSOMPT->AddConstProperty("FASTTIMECONSTANT", 42.0*ns);
	fLYSOMPT->AddConstProperty("YIELDRATIO",1); // the strength of the fast component as a function of total scintillation yield

/*
	// Original version, not trustworthy because copy-pasted.
	G4double electronYield_LYSO[36];
	G4double protonYield_LYSO[36];
	G4double ionYield_LYSO[36];

	for(size_t i = 0; i < 36; i++){
		electronYield_LYSO[i] = 0.86 * electronYield[i];
		protonYield_LYSO[i] = 0.86 * protonYield[i];
		ionYield_LYSO[i] = 0.86 * ionYield[i];
	}

	fLYSOMPT->AddProperty("ELECTRONSCINTILLATIONYIELD", particleEnergy, electronYield_LYSO, 36)->SetSpline(true);
	fLYSOMPT->AddProperty("PROTONSCINTILLATIONYIELD", particleEnergy, protonYield_LYSO, 36)->SetSpline(true);
	fLYSOMPT->AddProperty("IONSCINTILLATIONYIELD", particleEnergy, ionYield_LYSO, 36)->SetSpline(true);
	fLYSOMPT->AddProperty("ALPHASCINTILLATIONYIELD", particleEnergy, protonYield_LYSO, 36)->SetSpline(true);
	fLYSOMPT->AddProperty("DEUTERONSCINTILLATIONYIELD", particleEnergy, protonYield_LYSO, 36)->SetSpline(true);
	fLYSOMPT->AddProperty("TRITONSCINTILLATIONYIELD", particleEnergy, protonYield_LYSO, 36)->SetSpline(true);
*/

	// Average Value Approximation (Birks formula) coefficients for LYSO
	//G4double a_LYSO = 0.596, b_LYSO = 0.00877; // OG
	G4double a_LYSO = 6.07798e-01, b_LYSO = 1.34628e-02, c_LYSO = 2.26845e-02; // Formula with offset

	// Use a different photon emission factor since AVA/Birks formula is calibrated to energy
	G4double pEF_LYSO = 10000;

	// Just assume here that: energy = light yield, i.e. no quenching occurs
	std::vector<G4double> particleEnergyElectron_LYSO = {0.001*MeV, 0.01*MeV, 0.011*MeV, 0.012*MeV, 0.013*MeV, 0.014*MeV, 0.015*MeV, 0.016*MeV, 0.017*MeV, 0.018*MeV, 0.02*MeV, 0.0225*MeV, 0.025*MeV, 0.0275*MeV, 0.03*MeV, 0.0325*MeV, 0.035*MeV, 0.0375*MeV, 0.04*MeV, 0.045*MeV, 0.05*MeV, 0.055*MeV, 0.06*MeV, 0.065*MeV, 0.07*MeV, 0.08*MeV, 0.09*MeV, 0.1*MeV, 0.11*MeV, 0.12*MeV, 0.13*MeV, 0.14*MeV, 0.15*MeV, 0.16*MeV, 0.17*MeV, 0.18*MeV, 0.2*MeV, 0.225*MeV, 0.25*MeV, 0.275*MeV, 0.3*MeV, 0.325*MeV, 0.35*MeV, 0.375*MeV, 0.4*MeV, 0.45*MeV, 0.5*MeV, 0.55*MeV, 0.6*MeV, 0.65*MeV, 0.7*MeV, 0.8*MeV, 0.9*MeV, 1*MeV, 1.1*MeV, 1.2*MeV, 1.3*MeV, 1.4*MeV, 1.5*MeV, 1.6*MeV, 1.7*MeV, 1.8*MeV, 2*MeV, 2.25*MeV, 2.5*MeV, 2.75*MeV, 3*MeV, 3.25*MeV, 3.5*MeV, 3.75*MeV, 4*MeV, 4.5*MeV, 5*MeV, 5.5*MeV, 6*MeV, 6.5*MeV, 7*MeV, 8*MeV, 9*MeV, 10*MeV, 11*MeV, 12*MeV, 13*MeV, 14*MeV, 15*MeV, 16*MeV, 17*MeV, 18*MeV, 20*MeV, 22.5*MeV, 25*MeV, 27.5*MeV, 30*MeV, 32.5*MeV, 35*MeV, 37.5*MeV, 40*MeV, 45*MeV, 50*MeV, 55*MeV, 60*MeV, 65*MeV, 70*MeV, 80*MeV, 90*MeV, 100*MeV, 110*MeV, 120*MeV, 130*MeV, 140*MeV, 150*MeV, 160*MeV, 170*MeV, 180*MeV, 200*MeV, 225*MeV, 250*MeV, 275*MeV, 300*MeV, 325*MeV, 350*MeV, 375*MeV, 400*MeV, 450*MeV, 500*MeV, 550*MeV, 600*MeV, 650*MeV, 700*MeV, 800*MeV, 900*MeV, 1000*MeV, 1100*MeV, 1200*MeV, 1300*MeV, 1400*MeV, 1500*MeV, 1600*MeV, 1700*MeV, 1800*MeV, 2000*MeV, 2250*MeV, 2500*MeV, 2750*MeV, 3000*MeV, 3250*MeV, 3500*MeV, 3750*MeV, 4000*MeV, 4500*MeV, 5000*MeV, 5500*MeV, 6000*MeV, 6500*MeV, 7000*MeV, 8000*MeV, 9000*MeV, 10000*MeV};
	std::vector<G4double> electronYield_LYSO = {0.00*MeV, 0.01*MeV, 0.011*MeV, 0.012*MeV, 0.013*MeV, 0.014*MeV, 0.015*MeV, 0.016*MeV, 0.017*MeV, 0.018*MeV, 0.02*MeV, 0.0225*MeV, 0.025*MeV, 0.0275*MeV, 0.03*MeV, 0.0325*MeV, 0.035*MeV, 0.0375*MeV, 0.04*MeV, 0.045*MeV, 0.05*MeV, 0.055*MeV, 0.06*MeV, 0.065*MeV, 0.07*MeV, 0.08*MeV, 0.09*MeV, 0.1*MeV, 0.11*MeV, 0.12*MeV, 0.13*MeV, 0.14*MeV, 0.15*MeV, 0.16*MeV, 0.17*MeV, 0.18*MeV, 0.2*MeV, 0.225*MeV, 0.25*MeV, 0.275*MeV, 0.3*MeV, 0.325*MeV, 0.35*MeV, 0.375*MeV, 0.4*MeV, 0.45*MeV, 0.5*MeV, 0.55*MeV, 0.6*MeV, 0.65*MeV, 0.7*MeV, 0.8*MeV, 0.9*MeV, 1*MeV, 1.1*MeV, 1.2*MeV, 1.3*MeV, 1.4*MeV, 1.5*MeV, 1.6*MeV, 1.7*MeV, 1.8*MeV, 2*MeV, 2.25*MeV, 2.5*MeV, 2.75*MeV, 3*MeV, 3.25*MeV, 3.5*MeV, 3.75*MeV, 4*MeV, 4.5*MeV, 5*MeV, 5.5*MeV, 6*MeV, 6.5*MeV, 7*MeV, 8*MeV, 9*MeV, 10*MeV, 11*MeV, 12*MeV, 13*MeV, 14*MeV, 15*MeV, 16*MeV, 17*MeV, 18*MeV, 20*MeV, 22.5*MeV, 25*MeV, 27.5*MeV, 30*MeV, 32.5*MeV, 35*MeV, 37.5*MeV, 40*MeV, 45*MeV, 50*MeV, 55*MeV, 60*MeV, 65*MeV, 70*MeV, 80*MeV, 90*MeV, 100*MeV, 110*MeV, 120*MeV, 130*MeV, 140*MeV, 150*MeV, 160*MeV, 170*MeV, 180*MeV, 200*MeV, 225*MeV, 250*MeV, 275*MeV, 300*MeV, 325*MeV, 350*MeV, 375*MeV, 400*MeV, 450*MeV, 500*MeV, 550*MeV, 600*MeV, 650*MeV, 700*MeV, 800*MeV, 900*MeV, 1000*MeV, 1100*MeV, 1200*MeV, 1300*MeV, 1400*MeV, 1500*MeV, 1600*MeV, 1700*MeV, 1800*MeV, 2000*MeV, 2250*MeV, 2500*MeV, 2750*MeV, 3000*MeV, 3250*MeV, 3500*MeV, 3750*MeV, 4000*MeV, 4500*MeV, 5000*MeV, 5500*MeV, 6000*MeV, 6500*MeV, 7000*MeV, 8000*MeV, 9000*MeV, 10000*MeV};
	for (int i=0; i < electronYield_LYSO.size(); i++)
		electronYield_LYSO[i] *= pEF_LYSO/MeV;

	// Calculate light output with AVA/Birks formula. NOTE that the stopping power dE/dx data is taken
	//  from SRIM!
	std::vector<G4double> particleEnergy_LYSO = {0.01*MeV, 0.011*MeV, 0.012*MeV, 0.013*MeV, 0.014*MeV, 0.015*MeV, 0.016*MeV, 0.017*MeV, 0.018*MeV, 0.02*MeV, 0.0225*MeV, 0.025*MeV, 0.0275*MeV, 0.03*MeV, 0.0325*MeV, 0.035*MeV, 0.0375*MeV, 0.04*MeV, 0.045*MeV, 0.05*MeV, 0.055*MeV, 0.06*MeV, 0.065*MeV, 0.07*MeV, 0.08*MeV, 0.09*MeV, 0.1*MeV, 0.11*MeV, 0.12*MeV, 0.13*MeV, 0.14*MeV, 0.15*MeV, 0.16*MeV, 0.17*MeV, 0.18*MeV, 0.2*MeV, 0.225*MeV, 0.25*MeV, 0.275*MeV, 0.3*MeV, 0.325*MeV, 0.35*MeV, 0.375*MeV, 0.4*MeV, 0.45*MeV, 0.5*MeV, 0.55*MeV, 0.6*MeV, 0.65*MeV, 0.7*MeV, 0.8*MeV, 0.9*MeV, 1*MeV, 1.1*MeV, 1.2*MeV, 1.3*MeV, 1.4*MeV, 1.5*MeV, 1.6*MeV, 1.7*MeV, 1.8*MeV, 2*MeV, 2.25*MeV, 2.5*MeV, 2.75*MeV, 3*MeV, 3.25*MeV, 3.5*MeV, 3.75*MeV, 4*MeV, 4.5*MeV, 5*MeV, 5.5*MeV, 6*MeV, 6.5*MeV, 7*MeV, 8*MeV, 9*MeV, 10*MeV, 11*MeV, 12*MeV, 13*MeV, 14*MeV, 15*MeV, 16*MeV, 17*MeV, 18*MeV, 20*MeV, 22.5*MeV, 25*MeV, 27.5*MeV, 30*MeV, 32.5*MeV, 35*MeV, 37.5*MeV, 40*MeV, 45*MeV, 50*MeV, 55*MeV, 60*MeV, 65*MeV, 70*MeV, 80*MeV, 90*MeV, 100*MeV, 110*MeV, 120*MeV, 130*MeV, 140*MeV, 150*MeV, 160*MeV, 170*MeV, 180*MeV, 200*MeV, 225*MeV, 250*MeV, 275*MeV, 300*MeV, 325*MeV, 350*MeV, 375*MeV, 400*MeV, 450*MeV, 500*MeV, 550*MeV, 600*MeV, 650*MeV, 700*MeV, 800*MeV, 900*MeV, 1000*MeV, 1100*MeV, 1200*MeV, 1300*MeV, 1400*MeV, 1500*MeV, 1600*MeV, 1700*MeV, 1800*MeV, 2000*MeV, 2250*MeV, 2500*MeV, 2750*MeV, 3000*MeV, 3250*MeV, 3500*MeV, 3750*MeV, 4000*MeV, 4500*MeV, 5000*MeV, 5500*MeV, 6000*MeV, 6500*MeV, 7000*MeV, 8000*MeV, 9000*MeV, 10000*MeV};
	std::vector<G4double> protonStopping_LYSO = {118.687, 123.895, 128.913, 133.64, 138.273, 142.713, 146.858, 151.007, 154.86, 162.277, 170.687, 178.41, 185.444, 191.885, 197.834, 203.287, 208.246, 212.708, 220.742, 227.487, 233.039, 237.697, 241.461, 244.528, 248.773, 250.928, 251.49, 250.858, 249.331, 247.106, 244.385, 241.166, 237.749, 234.134, 230.42, 222.896, 213.471, 204.35, 195.832, 187.917, 180.504, 173.693, 167.382, 161.573, 151.358, 142.645, 135.235, 128.826, 123.318, 118.611, 110.8, 104.691, 99.9231, 95.9168, 91.9215, 88.1169, 84.6929, 81.5694, 78.7063, 76.0535, 73.591, 69.1767, 64.4423, 60.4087, 56.9157, 53.8732, 51.191, 48.7991, 46.6674, 44.7359, 41.3835, 38.5714, 36.1597, 34.0783, 32.2571, 30.646, 27.9243, 25.6929, 23.8417, 22.2708, 20.92, 19.7393, 18.7087, 17.7882, 16.9678, 16.2274, 15.557, 14.3864, 13.1858, 12.1952, 11.3548, 10.6445, 10.0241, 9.48088, 9.00365, 8.57945, 7.8571, 7.26382, 6.76759, 6.3454, 5.98223, 5.66609, 5.14185, 4.72367, 4.38251, 4.09939, 3.85928, 3.65319, 3.47512, 3.31905, 3.18099, 3.05894, 2.94889, 2.76081, 2.57073, 2.41666, 2.29061, 2.18456, 2.09552, 2.01849, 1.95146, 1.89343, 1.79639, 1.72035, 1.65832, 1.6073, 1.56528, 1.52926, 1.47423, 1.4322, 1.40219, 1.37817, 1.36016, 1.34615, 1.33514, 1.32713, 1.32112, 1.31611, 1.31311, 1.3091, 1.30809, 1.31008, 1.31407, 1.31907, 1.32406, 1.33106, 1.33705, 1.34405, 1.35805, 1.37104, 1.38504, 1.39804, 1.41003, 1.42203, 1.44403, 1.46502, 1.48402};
	std::vector<G4double> alphaStopping_LYSO = {118.69, 122.91, 127.16, 131.26, 135.21, 139.2, 143.13, 146.99, 150.89, 158.45, 167.65, 176.54, 185.33, 193.78, 202.08, 210.14, 218.03, 225.76, 240.705, 255.042, 268.852, 282.221, 295.138, 307.595, 331.405, 353.914, 374.998, 395.04, 413.928, 431.953, 448.909, 465.091, 480.494, 495.115, 508.953, 534.766, 563.567, 589.117, 611.604, 631.421, 648.861, 664.12, 677.595, 689.283, 708.09, 721.83, 731.496, 737.68, 740.98, 741.992, 738.845, 730.626, 719.129, 705.647, 691.177, 676.116, 660.963, 646.016, 631.375, 617.138, 603.404, 577.447, 547.987, 521.438, 497.598, 476.163, 456.733, 439.207, 423.285, 408.664, 383.13, 361.402, 342.779, 326.559, 312.342, 299.827, 278.703, 261.083, 244.167, 229.954, 217.443, 206.534, 196.825, 188.118, 180.212, 173.106, 166.601, 155.192, 143.183, 133.176, 124.67, 117.264, 110.86, 105.156, 100.153, 95.62, 87.8351, 81.3711, 75.9078, 71.215, 67.1326, 63.5605, 57.5671, 52.7244, 48.7322, 45.3804, 42.5088, 40.0375, 37.8664, 35.9554, 34.2646, 32.7438, 31.3731, 28.9919, 26.5607, 24.5697, 22.9089, 21.4983, 20.2877, 19.2272, 18.2967, 17.4864, 16.0957, 14.9652, 14.0248, 13.2344, 12.5641, 11.9738, 11.0134, 10.253, 9.64276, 9.13853, 8.71733, 8.36017, 8.05403, 7.7879, 7.55579, 7.3527, 7.17161, 6.86646, 6.56631, 6.33119, 6.14309, 5.991, 5.86593, 5.76187, 5.67582, 5.60377, 5.49069, 5.40963, 5.35157, 5.30953, 5.27949, 5.25846, 5.23741, 5.23436, 5.24233};
	std::vector<G4double> C12Stopping_LYSO = {340.1, 344.4, 348.5, 352.3, 356, 359.5, 362.9, 366.3, 369.5, 375.7, 383.1, 390.8, 398.7, 406.2, 413.2, 419.8, 426.2, 432.5, 444.4, 456.2, 467.6, 478.9, 490, 501.1, 522.8, 544.42, 565.79, 586.91, 607.98, 628.93, 649.69, 670.41, 690.88, 711.14, 731.28, 771.02, 819.47, 866.73, 912.65, 957.15, 1000.46, 1042.42, 1083, 1123.26, 1198.22, 1270.66, 1338.46, 1402.55, 1462.87, 1521.38, 1629.86, 1727.8, 1819.07, 1901.61, 1978.34, 2049.24, 2114.27, 2174.41, 2230.64, 2282.95, 2330.32, 2416.22, 2506.08, 2579.14, 2639.35, 2689.67, 2729.09, 2761.57, 2788.12, 2809.71, 2839.02, 2854.46, 2861.98, 2861.57, 2855.22, 2845.92, 2816.41, 2779, 2736.67, 2692.39, 2646.16, 2600.95, 2555.78, 2510.62, 2467.48, 2425.36, 2385.25, 2308.06, 2219.86, 2140.7, 2059.57, 1975.46, 1899.36, 1829.28, 1764.2, 1704.14, 1596.03, 1501.94, 1419.87, 1346.8, 1282.75, 1224.7, 1124.63, 1042.56, 971.914, 911.373, 858.438, 811.709, 770.083, 732.56, 698.54, 667.622, 639.406, 589.479, 537.452, 494.329, 458.61, 429.095, 404.681, 384.97, 366.759, 349.251, 319.435, 294.923, 274.413, 257.005, 241.997, 228.891, 207.281, 189.973, 175.866, 164.161, 154.256, 145.652, 138.349, 131.846, 126.143, 121.041, 116.439, 108.635, 100.732, 94.3487, 89.0663, 84.6343, 80.8726, 77.6311, 74.8198, 72.3487, 68.2468, 64.9752, 62.314, 60.1129, 58.272, 56.7112, 54.2399, 52.3789, 50.9481};
	std::vector<G4double> O24Stopping_LYSO = {529.6, 534.1, 537.8, 541.1, 543.9, 546.4, 548.7, 550.7, 552.6, 555.8, 559.3, 562.4, 565.1, 567.6, 570.1, 572.4, 574.6, 576.9, 581.3, 584.1, 585.2, 587, 589.2, 591.8, 598.2, 605.5, 613.7, 622.5, 631.7, 641.2, 650.9, 660.7, 670.6, 680.7, 690.8, 710.9, 736.4, 761.8, 787.3, 812.9, 838.3, 863.8, 889.4, 914.7, 965.6, 1016.1, 1066.5, 1116.2, 1166.46, 1215.23, 1312.01, 1407.05, 1500.02, 1590.73, 1679, 1764.74, 1847.86, 1928.28, 2006.97, 2081.89, 2155.99, 2295.67, 2458.21, 2607.34, 2743.9, 2868.82, 2983, 3087.41, 3183, 3270.74, 3421.58, 3546.8, 3651.31, 3737.03, 3807.92, 3865.95, 3953.34, 4012.04, 4050.98, 4074.09, 4086.33, 4090.67, 4088.1, 4080.6, 4068.16, 4053.76, 4036.4, 3997.78, 3942.14, 3882.63, 3822.19, 3762.83, 3703.51, 3645.24, 3590, 3535.78, 3434.42, 3343.13, 3241.88, 3132.68, 3033.5, 2939.34, 2770.09, 2619.89, 2486.72, 2367.58, 2261.47, 2165.37, 2078.28, 2000.21, 1927.14, 1861.08, 1800.03, 1689.94, 1572.85, 1471.77, 1384.71, 1307.66, 1238.61, 1177.57, 1122.54, 1071.51, 982.958, 907.617, 842.983, 787.154, 738.53, 696.209, 629.074, 575.746, 531.924, 495.106, 463.89, 436.877, 413.465, 392.755, 374.447, 358.039, 343.332, 317.92, 292.008, 270.798, 253.29, 238.383, 225.677, 214.672, 205.068, 196.564, 182.257, 170.652, 161.148, 153.144, 146.241, 140.438, 130.834, 123.43, 117.528};
	std::vector<G4double> Ar40Stopping_LYSO = {1204, 1225.3, 1245.2, 1262.8, 1278.1, 1293.2, 1306, 1317.6, 1328, 1347.4, 1368, 1384.8, 1399.9, 1411.4, 1422.4, 1432, 1440.2, 1447, 1459.6, 1469.1, 1477.6, 1484.4, 1490.4, 1495.7, 1503.8, 1510.5, 1517, 1522, 1526, 1529.3, 1532.1, 1534.8, 1537.1, 1539.6, 1542.3, 1547.8, 1556, 1565.9, 1577.4, 1590.8, 1605.8, 1622.1, 1640.1, 1659.3, 1699, 1743.7, 1788.4, 1835.8, 1883.3, 1931.7, 2027.2, 2121.6, 2213.7, 2301.9, 2387.6, 2470.3, 2550.6, 2628.4, 2704.4, 2778.4, 2850.4, 2989.6, 3158.7, 3322.8, 3484.2, 3644.3, 3802.9, 3959.7, 4117.5, 4274.1, 4584.3, 4891.7, 5193.6, 5487.8, 5774, 6052, 6574.91, 7055.75, 7492.99, 7890.28, 8248.39, 8571.16, 8862.45, 9125.17, 9362.24, 9576.62, 9769.25, 10097.1, 10422.9, 10679.4, 10866.5, 11014, 11121.9, 11200, 11248.4, 11286.9, 11304.4, 11292.4, 11250.7, 11199.2, 11138, 11066.9, 10925.1, 10753.7, 10502.5, 10271.5, 10050.7, 9843.02, 9644.4, 9456.86, 9279.39, 9110.96, 8951.58, 8655.92, 8324.25, 8025.71, 7753.26, 7502.87, 7269.55, 7051.26, 6845.01, 6648.79, 6284.42, 5951.12, 5643.87, 5360.66, 5100.48, 4860.32, 4439.06, 4092.86, 3813.69, 3598.55, 3441.44, 3260.34, 3100.25, 2957.18, 2829.11, 2714.05, 2609, 2419.91, 2220.82, 2055.75, 1917.69, 1800.63, 1699.59, 1611.55, 1534.52, 1465.49, 1349.44, 1254.4, 1176.37, 1109.34, 1052.31, 1003.29, 922.461, 858.634, 807.113};
	std::vector<G4double> Ti56Stopping_LYSO = {1435.1, 1465.3, 1492.2, 1515.9, 1538.4, 1557.8, 1576, 1593.1, 1608.1, 1634.7, 1663.7, 1687.2, 1707.4, 1724.1, 1739.6, 1751.7, 1762.6, 1771.3, 1786, 1798.1, 1805.6, 1811.6, 1816.1, 1819.2, 1821.4, 1821.4, 1819.5, 1815.7, 1829.6, 1843.8, 1856.3, 1866.4, 1874.5, 1881.9, 1887.9, 1898, 1907.6, 1915.3, 1921.1, 1925.7, 1929.7, 1934.6, 1938.6, 1942.2, 1951.7, 1962.1, 1974.4, 1988.1, 2002.7, 2020, 2059, 2104, 2153.6, 2207.7, 2265.6, 2325.5, 2388.1, 2452.9, 2519.7, 2586.1, 2654.9, 2793.4, 2969.8, 3145, 3321.1, 3495.2, 3666.8, 3836.4, 4003.7, 4167.4, 4486.3, 4794, 5089.5, 5374.2, 5647.7, 5911.7, 6409.1, 6869.9, 7296.3, 7692.7, 8061.5, 8404.5, 8724.5, 9024.3, 9304.8, 9567.8, 9816.34, 10271.5, 10773.4, 11216.8, 11601.2, 11946.4, 12252.3, 12518.7, 12765.5, 12982.7, 13337.9, 13634, 13870.7, 14047.9, 14195.5, 14313.3, 14479.8, 14567, 14614.8, 14632.9, 14621.2, 14519.9, 14358.7, 14207.6, 14046.7, 13885.8, 13715.1, 13383.8, 12972.5, 12591.4, 12220.5, 11879.7, 11549.1, 11248.5, 10958, 10677.6, 10176.8, 9725.24, 9310.74, 8931.32, 8579.96, 8255.65, 7671.14, 7159.73, 6710.4, 6315.12, 5967.89, 5662.7, 5396.52, 5166.37, 4968.24, 4794.12, 4613.02, 4293.84, 3960.65, 3683.5, 3448.38, 3246.28, 3069.19, 2908.11, 2767.05, 2640.99, 2425.89, 2250.81, 2103.74, 1978.69, 1872.64, 1779.6, 1625.53, 1504.47, 1405.43};
	std::vector<G4double> Cu82Stopping_LYSO = {1806.8, 1853.72, 1895.51, 1934.18, 1969.75, 2002.24, 2031.64, 2059.97, 2085.23, 2131.57, 2181.7, 2225.55, 2262.16, 2295.6, 2323.8, 2348.9, 2371.8, 2391.6, 2424.9, 2452.7, 2474.2, 2492.5, 2506.5, 2517.2, 2534.1, 2543.4, 2547.1, 2548.3, 2546.2, 2541.6, 2535.8, 2528.6, 2521.2, 2526.2, 2534.1, 2537.3, 2525, 2503.6, 2479.7, 2454, 2429.8, 2406.8, 2386.4, 2367.7, 2335.5, 2312.1, 2293.9, 2282.4, 2275.1, 2271.6, 2274.7, 2286.9, 2306.2, 2331, 2359, 2389, 2422, 2457, 2494, 2531, 2570, 2649.8, 2752.5, 2858.6, 2967.5, 3077.8, 3190.7, 3305.2, 3420, 3536.4, 3772.9, 4011.5, 4251.7, 4493.2, 4734.3, 4975.2, 5453.1, 5923, 6382.1, 6828.2, 7261.5, 7680.3, 8084.1, 8472.6, 8845.5, 9204.7, 9549.8, 10199.6, 10940.7, 11615.9, 12233.5, 12792.8, 13303.6, 13785.5, 14218.3, 14632, 15361.1, 15992.1, 16534.6, 17018.2, 17442.7, 17807.9, 18419.8, 18913.4, 19298.2, 19603.8, 19850, 20046.8, 20214, 20341.5, 20449.3, 20547.4, 20585.7, 20472.6, 20279.5, 20027, 19754.9, 19473.1, 19181.6, 18880.3, 18589.1, 18308.1, 17756.3, 17234.9, 16753.7, 16292.7, 15851.9, 15441.1, 14679.9, 13989, 13348.2, 12767.5, 12227, 11716.5, 11256.1, 10825.7, 10425.4, 10055.1, 9709.87, 9101.43, 8474.99, 7970.64, 7496.34, 7083.09, 6718.88, 6396.7, 6108.53, 5849.39, 5402.15, 5029.96, 4715.8, 4445.66, 4211.55, 4001.45, 3645.28, 3361.15, 3128.05};
	std::vector<G4double> Sn103Stopping_LYSO = {2529.57, 2613.09, 2689.4, 2760.5, 2826.5, 2888.4, 2946.1, 2999.7, 3051.2, 3145.9, 3250.9, 3343.4, 3426.5, 3502.3, 3570.9, 3633.2, 3691.3, 3744.1, 3838.4, 3920, 3991.1, 4053.8, 4110, 4159.9, 4243.8, 4312.7, 4368.7, 4416, 4454.7, 4486.8, 4513.4, 4536.5, 4555.3, 4570.7, 4582.8, 4601, 4637.3, 4656.3, 4659.7, 4654.8, 4644.6, 4631.7, 4618, 4605, 4578.4, 4554.9, 4534.8, 4518.4, 4503.9, 4491.7, 4472, 4458, 4448, 4441, 4435, 4433, 4431, 4431, 4433, 4435, 4440, 4453, 4478, 4512, 4555, 4608, 4667, 4733, 4807, 4886, 5061, 5255, 5463, 5682, 5909, 6144, 6623, 7111, 7597, 8076.5, 8547.3, 9005.6, 9451.2, 9882.3, 10300.1, 10704.2, 11092, 11832.4, 12690.6, 13487.4, 14220.7, 14919.1, 15561.5, 16187.3, 16775.8, 17336.5, 18393.5, 19386.2, 20313.1, 21193.4, 22016.3, 22801.3, 24256.3, 25566.1, 26769.5, 27865.6, 28873.8, 29803.5, 30654.6, 31456.7, 32189.7, 32883.4, 33537.8, 34718, 36008, 36769.8, 37383, 37857.2, 38222.2, 38497.8, 38704, 38850.6, 38994.9, 39000.2, 38896.3, 38713, 38480.1, 38197.7, 37533.6, 36770.4, 35957.7, 35095.6, 34223.7, 33332.1, 32450.7, 31579.5, 30718.5, 29887.5, 29076.6, 27555.2, 25833.7, 24352.5, 23121.5, 22120.6, 21209.9, 20339.3, 19538.7, 18818.2, 17537.4, 16446.8, 15496.2, 14675.8, 13955.4, 13315, 12224.5, 11334, 10603.6};
	std::vector<G4double> Te107Stopping_LYSO = {2580.85, 2666.99, 2746.95, 2820.74, 2889.4, 2953.9, 3014.3, 3070.6, 3123.8, 3222, 3332.3, 3430.2, 3517.7, 3597, 3669, 3734.7, 3795.3, 3851.7, 3951, 4037.7, 4114, 4181.8, 4241.3, 4294.5, 4385, 4459.5, 4521.3, 4571.4, 4614.9, 4649.9, 4680.5, 4705.6, 4727.4, 4744.8, 4759, 4781.5, 4828.1, 4878.4, 4898.5, 4901.2, 4892.7, 4879.4, 4862.6, 4844.5, 4808.2, 4773.8, 4744.7, 4720.3, 4699.9, 4683.7, 4659.2, 4645, 4636, 4633, 4632, 4633, 4636, 4640, 4646, 4652, 4657, 4673, 4696, 4722, 4755, 4793, 4835, 4883, 4937, 4997, 5131, 5282, 5451, 5634, 5828, 6032, 6461, 6912, 7372, 7835, 8298, 8755.2, 9203, 9642.6, 10069.3, 10485.4, 10892.6, 11658.3, 12551.6, 13374.2, 14133.9, 14849.2, 15518.9, 16142.3, 16738.6, 17297.5, 18351.1, 19331, 20245.6, 21103.8, 21915, 22688.4, 24120.8, 25428.5, 26630.2, 27734.8, 28761.7, 29710.4, 30600.5, 31431.8, 32214, 32947.1, 33640.8, 34920, 36358.9, 37369.8, 38072.2, 38675.8, 39160.2, 39545.4, 39851.2, 40087.4, 40401.1, 40545.8, 40571.5, 40497.8, 40344.7, 40142, 39607.4, 38943.8, 38190.9, 37388.4, 36556.4, 35694.6, 34833.1, 33971.8, 33110.6, 32279.5, 31458.6, 29906.9, 28145.3, 26603.9, 25302.8, 24241.9, 23371, 22430.4, 21569.7, 20779.2, 19388.3, 18197.5, 17166.9, 16266.4, 15476, 14765.6, 13565, 12584.5, 11774.1};
	std::vector<G4double> Xe108Stopping_LYSO = {2701, 2793.6, 2878.8, 2958.8, 3032.5, 3102, 3167.3, 3228.3, 3286.2, 3393.6, 3513.9, 3620.6, 3717.6, 3805.1, 3886.2, 3959.8, 4028, 4090.9, 4203.8, 4303.7, 4391.8, 4470, 4539.6, 4603.6, 4713.1, 4804.9, 4882.2, 4948.3, 5005.3, 5054.3, 5096.6, 5134.1, 5166.9, 5196.1, 5221.7, 5263.5, 5296.5, 5309.9, 5318.5, 5324.8, 5326.1, 5326, 5322.5, 5317.5, 5304.4, 5288, 5270, 5252, 5235, 5218, 5191, 5171, 5157, 5151, 5150, 5155, 5164, 5177, 5194, 5214, 5235, 5285, 5355, 5430, 5510, 5592, 5674, 5758, 5843, 5927, 6098, 6271, 6446, 6625, 6807, 6993, 7379, 7781, 8199, 8632, 9078, 9531, 9992, 10456.9, 10927.1, 11391, 11858.1, 12780.2, 13899.6, 14979, 16005.9, 16978.8, 17916.4, 18787.8, 19622.5, 20419.8, 21880.9, 23208.6, 24411.3, 25497.9, 26497.7, 27419.9, 29040.2, 30446.3, 31656.6, 32740.1, 33695.9, 34563.7, 35343.1, 36063.7, 36725.4, 37337.9, 37911.2, 38949.5, 40127.6, 40927.8, 41399.6, 41802.7, 42126.7, 42371.5, 42556.9, 42702.8, 42866, 42920.4, 42875.7, 42771.7, 42618.3, 42415.3, 41900.5, 41286.6, 40593.4, 39830.8, 39028.5, 38196.6, 37355, 36493.5, 35642.2, 34791.1, 33950.1, 32348.3, 30486.5, 28845.1, 27433.9, 26272.8, 25342, 24331.2, 23400.5, 22560, 21049, 19758.2, 18647.5, 17667, 16806.5, 16036.1, 14735.4, 13674.8, 12784.4};

	// Eventually want to be able to interpolate this based on the ion specified as the source, 
	//  but don't worry for now. Currently set to: alpha particle.
	//G4double ionStopping_LYSO[157] = {118.69, 122.91, 127.16, 131.26, 135.21, 139.2, 143.13, 146.99, 150.89, 158.45, 167.65, 176.54, 185.33, 193.78, 202.08, 210.14, 218.03, 225.76, 240.705, 255.042, 268.852, 282.221, 295.138, 307.595, 331.405, 353.914, 374.998, 395.04, 413.928, 431.953, 448.909, 465.091, 480.494, 495.115, 508.953, 534.766, 563.567, 589.117, 611.604, 631.421, 648.861, 664.12, 677.595, 689.283, 708.09, 721.83, 731.496, 737.68, 740.98, 741.992, 738.845, 730.626, 719.129, 705.647, 691.177, 676.116, 660.963, 646.016, 631.375, 617.138, 603.404, 577.447, 547.987, 521.438, 497.598, 476.163, 456.733, 439.207, 423.285, 408.664, 383.13, 361.402, 342.779, 326.559, 312.342, 299.827, 278.703, 261.083, 244.167, 229.954, 217.443, 206.534, 196.825, 188.118, 180.212, 173.106, 166.601, 155.192, 143.183, 133.176, 124.67, 117.264, 110.86, 105.156, 100.153, 95.62, 87.8351, 81.3711, 75.9078, 71.215, 67.1326, 63.5605, 57.5671, 52.7244, 48.7322, 45.3804, 42.5088, 40.0375, 37.8664, 35.9554, 34.2646, 32.7438, 31.3731, 28.9919, 26.5607, 24.5697, 22.9089, 21.4983, 20.2877, 19.2272, 18.2967, 17.4864, 16.0957, 14.9652, 14.0248, 13.2344, 12.5641, 11.9738, 11.0134, 10.253, 9.64276, 9.13853, 8.71733, 8.36017, 8.05403, 7.7879, 7.55579, 7.3527, 7.17161, 6.86646, 6.56631, 6.33119, 6.14309, 5.991, 5.86593, 5.76187, 5.67582, 5.60377, 5.49069, 5.40963, 5.35157, 5.30953, 5.27949, 5.25846, 5.23741, 5.23436, 5.24233};

	std::vector<G4double> protonYield_LYSO(particleEnergy_LYSO.size()+1);
	protonYield_LYSO[0] = 1;
	std::vector<G4double> alphaYield_LYSO(particleEnergy_LYSO.size()+1);
	alphaYield_LYSO[0] = 4;
	std::vector<G4double> C12Yield_LYSO(particleEnergy_LYSO.size()+1);
	C12Yield_LYSO[0] = 12;
	std::vector<G4double> O24Yield_LYSO(particleEnergy_LYSO.size()+1);
	O24Yield_LYSO[0] = 24;
	std::vector<G4double> Ar40Yield_LYSO(particleEnergy_LYSO.size()+1);
	Ar40Yield_LYSO[0] = 40;
	std::vector<G4double> Ti56Yield_LYSO(particleEnergy_LYSO.size()+1);
	Ti56Yield_LYSO[0] = 56;
	std::vector<G4double> Cu82Yield_LYSO(particleEnergy_LYSO.size()+1);
	Cu82Yield_LYSO[0] = 82;
	std::vector<G4double> Sn103Yield_LYSO(particleEnergy_LYSO.size()+1);
	Sn103Yield_LYSO[0] = 103;
	std::vector<G4double> Te107Yield_LYSO(particleEnergy_LYSO.size()+1);
	Te107Yield_LYSO[0] = 107;
	std::vector<G4double> Xe108Yield_LYSO(particleEnergy_LYSO.size()+1);
	Xe108Yield_LYSO[0] = 108;


	//G4double ionYield_LYSO[157]; // Eventually want to be able to interpolate this based on the ion specified as the source, but don't worry for now
	for (int i=1; i < protonYield_LYSO.size(); i++) {
	/*	// OG
		protonYield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i-1]*a_LYSO/(1+b_LYSO*protonStopping_LYSO[i-1]);
		alphaYield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i-1]*a_LYSO/(1+b_LYSO*alphaStopping_LYSO[i-1]);
		C12Yield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i-1]*a_LYSO/(1+b_LYSO*C12Stopping_LYSO[i-1]);
		O24Yield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i-1]*a_LYSO/(1+b_LYSO*O24Stopping_LYSO[i-1]);
		Ar40Yield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i-1]*a_LYSO/(1+b_LYSO*Ar40Stopping_LYSO[i-1]);
		Ti56Yield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i-1]*a_LYSO/(1+b_LYSO*Ti56Stopping_LYSO[i-1]);
		Cu82Yield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i-1]*a_LYSO/(1+b_LYSO*Cu82Stopping_LYSO[i-1]);
		Sn103Yield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i-1]*a_LYSO/(1+b_LYSO*Sn103Stopping_LYSO[i-1]);
		Xe108Yield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i-1]*a_LYSO/(1+b_LYSO*Xe108Stopping_LYSO[i-1]);
		//ionYield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i]*a_LYSO/(1+b_LYSO*ionStopping_LYSO[i]);
	*/

		// New formula with offset
		protonYield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i-1]*dLdE(protonStopping_LYSO[i-1],a_LYSO,b_LYSO,c_LYSO);
		alphaYield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i-1]*dLdE(alphaStopping_LYSO[i-1],a_LYSO,b_LYSO,c_LYSO);
		C12Yield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i-1]*dLdE(C12Stopping_LYSO[i-1],a_LYSO,b_LYSO,c_LYSO);
		O24Yield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i-1]*dLdE(O24Stopping_LYSO[i-1],a_LYSO,b_LYSO,c_LYSO);
		Ar40Yield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i-1]*dLdE(Ar40Stopping_LYSO[i-1],a_LYSO,b_LYSO,c_LYSO);
		Ti56Yield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i-1]*dLdE(Ti56Stopping_LYSO[i-1],a_LYSO,b_LYSO,c_LYSO);
		Cu82Yield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i-1]*dLdE(Cu82Stopping_LYSO[i-1],a_LYSO,b_LYSO,c_LYSO);
		Sn103Yield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i-1]*dLdE(Sn103Stopping_LYSO[i-1],a_LYSO,b_LYSO,c_LYSO);
		Te107Yield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i-1]*dLdE(Te107Stopping_LYSO[i-1],a_LYSO,b_LYSO,c_LYSO);
		Xe108Yield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i-1]*dLdE(Xe108Stopping_LYSO[i-1],a_LYSO,b_LYSO,c_LYSO);
		//ionYield_LYSO[i] = pEF_LYSO*particleEnergy_LYSO[i-1]*dLdE(ionStopping_LYSO[i-1],a_LYSO,b_LYSO,c_LYSO);

	}

	// Create a matrix of light yields. NOTE that the first column gives the mass
	//  of the ion in a particular row
	std::vector< std::vector<G4double> > lightMatrix_LYSO;
	lightMatrix_LYSO.push_back(protonYield_LYSO);
	lightMatrix_LYSO.push_back(alphaYield_LYSO);
	lightMatrix_LYSO.push_back(C12Yield_LYSO);
	lightMatrix_LYSO.push_back(O24Yield_LYSO);
	lightMatrix_LYSO.push_back(Ar40Yield_LYSO);
	lightMatrix_LYSO.push_back(Ti56Yield_LYSO);
	lightMatrix_LYSO.push_back(Cu82Yield_LYSO);
	lightMatrix_LYSO.push_back(Sn103Yield_LYSO);
	lightMatrix_LYSO.push_back(Te107Yield_LYSO);
	lightMatrix_LYSO.push_back(Xe108Yield_LYSO);

	lightYieldTable["lyso"] = lightMatrix_LYSO;


	// Be careful about whether 0 or 1 is used as the index!!!
	fLYSOMPT->AddProperty("ELECTRONSCINTILLATIONYIELD", &particleEnergyElectron_LYSO[0], &electronYield_LYSO[0], 158)->SetSpline(true);
	fLYSOMPT->AddProperty("PROTONSCINTILLATIONYIELD", &particleEnergy_LYSO[0], &protonYield_LYSO[1], 157)->SetSpline(true);
	fLYSOMPT->AddProperty("IONSCINTILLATIONYIELD", &particleEnergy_LYSO[0], &alphaYield_LYSO[1], 157)->SetSpline(true);		// set ion light to alpha particle as default
	fLYSOMPT->AddProperty("ALPHASCINTILLATIONYIELD", &particleEnergy_LYSO[0], &alphaYield_LYSO[1], 157)->SetSpline(true);
	fLYSOMPT->AddProperty("DEUTERONSCINTILLATIONYIELD", &particleEnergy_LYSO[0], &protonYield_LYSO[1], 157)->SetSpline(true);
	fLYSOMPT->AddProperty("TRITONSCINTILLATIONYIELD", &particleEnergy_LYSO[0], &protonYield_LYSO[1], 157)->SetSpline(true);

	fLYSO->SetMaterialPropertiesTable(fLYSOMPT);


	/////////////////////////////////////////////////////////////////
	// LaBr3 
	/////////////////////////////////////////////////////////////////

    fLaBr3 = new G4Material("LaBr3", 5.08*g/cm3, 2);
    fLaBr3->AddElement(fLa, 0.25);
    fLaBr3->AddElement(fBr, 0.75);

	G4double photonEnergy_LaBr3[44] = {2.004*eV, 2.058*eV, 2.112*eV, 2.166*eV, 2.220*eV, 2.274*eV, 2.328*eV, 2.382*eV, 2.436*eV, 2.490*eV, 
		                               2.517*eV, 2.552*eV, 2.585*eV, 2.613*eV, 2.635*eV, 2.656*eV, 2.686*eV, 2.720*eV, 2.749*eV, 2.772*eV, 
		                               2.791*eV, 2.809*eV, 2.826*eV, 2.842*eV, 2.861*eV, 2.884*eV, 2.919*eV, 2.946*eV, 2.954*eV, 2.961*eV, 
		                               2.967*eV, 2.974*eV, 2.981*eV, 2.987*eV, 2.994*eV, 3.001*eV, 3.009*eV, 3.018*eV, 3.029*eV, 3.041*eV, 
		                               3.056*eV, 3.083*eV, 3.137*eV, 3.191*eV};

	G4double ScintilFast_LaBr3[44] = {0.000, 0.001, 0.001, 0.002, 0.003, 0.006, 0.010, 0.018, 0.033, 0.060, 
		                              0.084, 0.122, 0.175, 0.234, 0.294, 0.356, 0.416, 0.473, 0.533, 0.594, 
		                              0.657, 0.720, 0.784, 0.846, 0.903, 0.962, 1.000, 0.917, 0.857, 0.798, 
		                              0.732, 0.669, 0.604, 0.542, 0.480, 0.422, 0.359, 0.297, 0.237, 0.170, 
		                              0.105, 0.028, 0.004, 0.000};
                                  
	G4double photonEnergy_LaBr3_2[2] = {2.004*eV, 3.191*eV};
	G4double RefIndex_LaBr3[2] = {1.90, 1.90};
	G4double Absorption_LaBr3[2] = {1.881*cm, 1.881*cm};

	fLaBr3MPT = new G4MaterialPropertiesTable();
	fLaBr3MPT->AddProperty("RINDEX", photonEnergy_LaBr3_2, RefIndex_LaBr3, 2);
	fLaBr3MPT->AddProperty("ABSLENGTH", photonEnergy_LaBr3_2, Absorption_LaBr3, 2);
	fLaBr3MPT->AddProperty("FASTCOMPONENT", photonEnergy_LaBr3, ScintilFast_LaBr3, 44);

	fLaBr3MPT->AddConstProperty("SCINTILLATIONYIELD", 630000/MeV); 
	fLaBr3MPT->AddConstProperty("RESOLUTIONSCALE", 1.0); // Intrinsic resolution

	//fLaBr3MPT->AddConstProperty("RISETIMECONSTANT", 0.9*ns); Geant4 10.1 TODO
	fLaBr3MPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 0.9*ns);
	fLaBr3MPT->AddConstProperty("FASTTIMECONSTANT", 2.1*ns);
	fLaBr3MPT->AddConstProperty("YIELDRATIO",1);// the strength of the fast component as a function of total scintillation yield

	std::cout << "nDetConstruction: Photon yield is set to " << lightYieldScale << "x scale\n";

  //light yield - data taken form V.V. Verbinski et al, Nucl. Instrum. & Meth. 65 (1968) 8-25
	G4double particleEnergy_LaBr3[36] = {1E-3, 0.10*MeV, 0.13*MeV, 0.17*MeV, 0.20*MeV, 0.24*MeV, 0.30*MeV, 0.34*MeV, 0.40*MeV, 
		                           0.48*MeV, 0.60*MeV, 0.72*MeV, 0.84*MeV, 1.00*MeV, 1.30*MeV, 1.70*MeV, 2.00*MeV, 2.40*MeV, 
		                           3.00*MeV, 3.40*MeV, 4.00*MeV, 4.80*MeV, 6.00*MeV, 7.20*MeV, 8.40*MeV, 10.00*MeV, 11.00*MeV, 
		                           12.00*MeV, 13.00*MeV, 14.00*MeV, 15.00*MeV, 16.00*MeV, 17.00*MeV, 18.00*MeV, 19.00*MeV, 20.00*MeV};
		                           
	G4double electronYield_LaBr3[36] = {0.0E+00*pEF, 1.0E+03*pEF, 1.3E+03*pEF, 1.7E+03*pEF, 2.0E+03*pEF, 2.4E+03*pEF, 3.0E+03*pEF, 3.4E+03*pEF, 4.0E+03*pEF, 
		                          4.8E+03*pEF, 6.0E+03*pEF, 7.2E+03*pEF, 8.4E+03*pEF, 1.0E+04*pEF, 1.3E+04*pEF, 1.7E+04*pEF, 2.0E+04*pEF, 2.4E+04*pEF, 
		                          3.0E+04*pEF, 3.4E+04*pEF, 4.0E+04*pEF, 4.8E+04*pEF, 6.0E+04*pEF, 7.2E+04*pEF, 8.4E+04*pEF, 1.0E+05*pEF, 1.1E+05*pEF, 
		                          1.2E+05*pEF, 1.3E+05*pEF, 1.4E+05*pEF, 1.5E+05*pEF, 1.6E+05*pEF, 1.7E+05*pEF, 1.8E+05*pEF, 1.9E+05*pEF, 2.0E+05*pEF};

	fLaBr3MPT->AddProperty("ELECTRONSCINTILLATIONYIELD", particleEnergy_LaBr3, electronYield_LaBr3, 36)->SetSpline(true);
	G4double protonYield_LaBr3[36] = {0.6*pSF, 67.1*pSF, 88.6*pSF, 120.7*pSF, 146.5*pSF, 183.8*pSF, 246.0*pSF, 290.0*pSF, 365.0*pSF, 
		                        483.0*pSF, 678.0*pSF, 910.0*pSF, 1175.0*pSF, 1562.0*pSF, 2385.0*pSF, 3660.0*pSF, 4725.0*pSF, 6250.0*pSF, 
		                        8660.0*pSF, 10420.0*pSF, 13270.0*pSF, 17180.0*pSF, 23100.0*pSF, 29500.0*pSF, 36200.0*pSF, 45500.0*pSF, 51826.7*pSF, 
		                        58313.7*pSF, 65047.2*pSF, 72027.4*pSF, 79254.2*pSF, 86727.6*pSF, 94447.6*pSF, 102414.2*pSF, 110627.4*pSF, 119087.2*pSF};

	// fLaBr3MPT->AddProperty("PROTONSCINTILLATIONYIELD", particleEnergy_LaBr3, protonYield_LaBr3, 36)->SetSpline(true);

	G4double ionYield_LaBr3[36] = {0.2*pEF, 10.4*pEF, 12.7*pEF, 15.7*pEF, 17.9*pEF, 20.8*pEF, 25.1*pEF, 27.9*pEF, 31.9*pEF, 
		                     36.8*pEF, 43.6*pEF, 50.2*pEF, 56.9*pEF, 65.7*pEF, 81.3*pEF, 101.6*pEF, 116.5*pEF, 136.3*pEF, 
		                     166.2*pEF, 187.1*pEF, 218.6*pEF, 260.5*pEF, 323.5*pEF, 387.5*pEF, 451.5*pEF, 539.9*pEF, 595.5*pEF, 
		                     651.8*pEF, 708.7*pEF, 766.2*pEF, 824.2*pEF, 882.9*pEF, 942.2*pEF, 1002.1*pEF, 1062.6*pEF, 1123.7*pEF}; 
                         
	// fLaBr3MPT->AddProperty("IONSCINTILLATIONYIELD", particleEnergy_LaBr3, ionYield_LaBr3, 36)->SetSpline(true);

	fLaBr3->SetMaterialPropertiesTable(fLaBr3MPT);


	/////////////////////////////////////////////////////////////////
	// YAP YAlO3
	/////////////////////////////////////////////////////////////////

  fYAP = new G4Material("YAP", 5.35*g/cm3, 3);
  fYAP->AddElement(fY, 0.2);
  fYAP->AddElement(fAl, 0.2);
	fYAP->AddElement(fO, 0.6);

	G4double photonEnergy_YAP[44] = {2.004*eV, 2.058*eV, 2.112*eV, 2.166*eV, 2.220*eV, 2.274*eV, 2.328*eV, 2.382*eV, 2.436*eV, 2.490*eV, 
		                               2.517*eV, 2.552*eV, 2.585*eV, 2.613*eV, 2.635*eV, 2.656*eV, 2.686*eV, 2.720*eV, 2.749*eV, 2.772*eV, 
		                               2.791*eV, 2.809*eV, 2.826*eV, 2.842*eV, 2.861*eV, 2.884*eV, 2.919*eV, 2.946*eV, 2.954*eV, 2.961*eV, 
		                               2.967*eV, 2.974*eV, 2.981*eV, 2.987*eV, 2.994*eV, 3.001*eV, 3.009*eV, 3.018*eV, 3.029*eV, 3.041*eV, 
		                               3.056*eV, 3.083*eV, 3.137*eV, 3.191*eV};

	G4double ScintilFast_YAP[44] = {0.000, 0.001, 0.001, 0.002, 0.003, 0.006, 0.010, 0.018, 0.033, 0.060, 
		                              0.084, 0.122, 0.175, 0.234, 0.294, 0.356, 0.416, 0.473, 0.533, 0.594, 
		                              0.657, 0.720, 0.784, 0.846, 0.903, 0.962, 1.000, 0.917, 0.857, 0.798, 
		                              0.732, 0.669, 0.604, 0.542, 0.480, 0.422, 0.359, 0.297, 0.237, 0.170, 
		                              0.105, 0.028, 0.004, 0.000};

	//G4double photonEnergy_YAP[13] = {16.8179*keV, 23.0707*keV, 30.8922*keV, 50.0014*keV, 59.1033*keV, 80.4609*keV, 122.4170*keV, 280.3698*keV, 504.5552*keV, 659.8267*keV, 830.5057*keV, 1188.7809*keV, 1316.9683*keV};

	//G4double ScintilFast_YAP[13] = {0.0482746, 0.0585062, 0.0648312, 0.0774812, 0.0748768, 0.0707841, 0.0776672, 0.0854804, 0.0878988, 0.0888289, 0.089387, 0.0880848, 0.0878988}; //This should be the efficiency curve..?

	//G4double electronYield_YAP[13] = {54.3455, 65.8639, 72.9843, 87.2251, 84.2932, 79.6859, 87.4346, 96.2304, 98.9529, 100, 100.628, 99.1623, 98.9529}; //I am not sure if this is the correct units

	G4double photonEnergy_YAP2[2] = {16.8179*keV, 23.0707*keV};

	G4double RefIndex_YAP[2] = {1.94,1.94};
	G4double Absorption_YAP[2] = {2.7*cm, 2.7*cm}; 

	fYAPMPT = new G4MaterialPropertiesTable();
	fYAPMPT->AddProperty("RINDEX", photonEnergy_YAP2, RefIndex_YAP, 2);
	fYAPMPT->AddProperty("ABSLENGTH", photonEnergy_YAP2, Absorption_YAP, 2);
	fYAPMPT->AddProperty("FASTCOMPONENT", photonEnergy_YAP, ScintilFast_YAP, 44);

	fYAPMPT->AddConstProperty("SCINTILLATIONYIELD", 26000/MeV); 
	fYAPMPT->AddConstProperty("RESOLUTIONSCALE", 1.0); // Intrinsic resolution

	fYAPMPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 2.0*ns);
	fYAPMPT->AddConstProperty("FASTTIMECONSTANT", 50.0*ns);
	fYAPMPT->AddConstProperty("YIELDRATIO",1);// the strength of the fast component as a function of total scintillation yield

	G4double electronYield_YAP[36];
	G4double protonYield_YAP[36];
	G4double ionYield_YAP[36];

	// Produce the scaled light-yield for YAP (scaled down from EJ200 by 14%). This will need to be adjusted for YAP as it has a different light yield than EJ276.
	for(size_t i = 0; i < 36; i++){
		electronYield_YAP[i] = 0.86 * electronYield[i];
		protonYield_YAP[i] = 0.86 * protonYield[i];
		ionYield_YAP[i] = 0.86 * ionYield[i];
	}

	fYAPMPT->AddProperty("ELECTRONSCINTILLATIONYIELD", particleEnergy, electronYield_YAP, 36)->SetSpline(true);
	fYAPMPT->AddProperty("PROTONSCINTILLATIONYIELD", particleEnergy, protonYield_YAP, 36)->SetSpline(true);
	fYAPMPT->AddProperty("IONSCINTILLATIONYIELD", particleEnergy, ionYield_YAP, 36)->SetSpline(true);
	fYAPMPT->AddProperty("ALPHASCINTILLATIONYIELD", particleEnergy, protonYield_YAP, 36)->SetSpline(true);
	fYAPMPT->AddProperty("DEUTERONSCINTILLATIONYIELD", particleEnergy, protonYield_YAP, 36)->SetSpline(true);
    fYAPMPT->AddProperty("TRITONSCINTILLATIONYIELD", particleEnergy, protonYield_YAP, 36)->SetSpline(true);

	fYAP->SetMaterialPropertiesTable(fYAPMPT);

	/////////////////////////////////////////////////////////////////
	// GAGG Gd3Al2Ga3O12
	/////////////////////////////////////////////////////////////////

  fGAGG = new G4Material("GAGG", 6.63*g/cm3, 4);
  fGAGG->AddElement(fGd, 0.15);
  fGAGG->AddElement(fAl, 0.1);
  fGAGG->AddElement(fGa, 0.15);
	fGAGG->AddElement(fO, 0.6);

	G4double photonEnergy_GAGG[24] = {2.657150587*eV,2.602196113*eV,2.546601968*eV,2.515488356*eV,2.493333596*eV,2.468871385*eV,2.444884509*eV,2.418773247*eV,2.398282404*eV,2.370667847*eV,2.35098064*eV,2.32682673*eV,2.286884568*eV,2.237186055*eV,2.183269296*eV,2.135911364*eV,2.077142785*eV,2.058265389*eV,2.008969264*eV,1.953517553*eV,1.913896916*eV,1.868114944*eV,1.818611898*eV,1.775832214*eV};

	G4double ScintilFast_GAGG[24] = {0.004895105,0.097202797,0.298601399,0.55034965,0.713986014,0.858741259,0.951048951,0.990909091,0.999300699,0.999300699,0.988811189,0.96993007,0.902797203,0.781118881,0.648951049,0.520979021,0.361538462,0.327972028,0.258741259,0.204195804,0.174825175,0.137062937,0.101398601,0.084615385};

	//G4double photonEnergy_GAGG[13] = {16.8179*keV, 23.0707*keV, 30.8922*keV, 50.0014*keV, 59.1033*keV, 80.4609*keV, 122.4170*keV, 280.3698*keV, 504.5552*keV, 659.8267*keV, 830.5057*keV, 1188.7809*keV, 1316.9683*keV};

	//G4double ScintilFast_GAGG[13] = {0.0482746, 0.0585062, 0.0648312, 0.0774812, 0.0748768, 0.0707841, 0.0776672, 0.0854804, 0.0878988, 0.0888289, 0.089387, 0.0880848, 0.0878988}; //This should be the efficiency curve..?

	//G4double electronYield_GAGG[13] = {54.3455, 65.8639, 72.9843, 87.2251, 84.2932, 79.6859, 87.4346, 96.2304, 98.9529, 100, 100.628, 99.1623, 98.9529}; //I am not sure if this is the correct units

	G4double photonEnergy_GAGG2[2] = {16.8179*keV, 23.0707*keV};
	G4double RefIndex_GAGG[2] = {1.9,1.9}; 
	G4double Absorption_GAGG[2] = {2.7*cm, 2.7*cm}; 

	fGAGGMPT = new G4MaterialPropertiesTable();
	fGAGGMPT->AddProperty("RINDEX", photonEnergy_GAGG2, RefIndex_GAGG, 2);
	fGAGGMPT->AddProperty("ABSLENGTH", photonEnergy_GAGG2, Absorption_GAGG, 2);
	fGAGGMPT->AddProperty("FASTCOMPONENT", photonEnergy_GAGG, ScintilFast_GAGG, 24);

	fGAGGMPT->AddConstProperty("SCINTILLATIONYIELD", 30000/MeV); 
	fGAGGMPT->AddConstProperty("RESOLUTIONSCALE", 1.0); // Intrinsic resolution

	fGAGGMPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 2.0*ns);
	fGAGGMPT->AddConstProperty("FASTTIMECONSTANT", 50.0*ns);
	fGAGGMPT->AddConstProperty("YIELDRATIO",1);// the strength of the fast component as a function of total scintillation yield

	G4double electronYield_GAGG[36];
	G4double protonYield_GAGG[36];
	G4double ionYield_GAGG[36];

	// Produce the scaled light-yield for GAGG (scaled down from EJ200 by 14%). This will need to be adjusted for GAGG as it has a different light yield than EJ276.
	for(size_t i = 0; i < 36; i++){
		electronYield_GAGG[i] = 0.86 * electronYield[i];
		protonYield_GAGG[i] = 0.86 * protonYield[i];
		ionYield_GAGG[i] = 0.86 * ionYield[i];
	}

	fGAGGMPT->AddProperty("ELECTRONSCINTILLATIONYIELD", particleEnergy, electronYield_GAGG, 36)->SetSpline(true);
	fGAGGMPT->AddProperty("PROTONSCINTILLATIONYIELD", particleEnergy, protonYield_GAGG, 36)->SetSpline(true);
	fGAGGMPT->AddProperty("IONSCINTILLATIONYIELD", particleEnergy, ionYield_GAGG, 36)->SetSpline(true);

	fGAGG->SetMaterialPropertiesTable(fGAGGMPT);

	/////////////////////////////////////////////////////////////////
	// Cerium Bromide CeBr3
	/////////////////////////////////////////////////////////////////

    fCeBr3 = new G4Material("CeBr3", 5.2*g/cm3, 2);
    fCeBr3->AddElement(fCe, 0.25);
    fCeBr3->AddElement(fBr, 0.75);


	// taken from LaBr3
	G4double photonEnergy_CeBr3[44] = {2.004*eV, 2.058*eV, 2.112*eV, 2.166*eV, 2.220*eV, 2.274*eV, 2.328*eV, 2.382*eV, 2.436*eV, 2.490*eV, 
		                               2.517*eV, 2.552*eV, 2.585*eV, 2.613*eV, 2.635*eV, 2.656*eV, 2.686*eV, 2.720*eV, 2.749*eV, 2.772*eV, 
		                               2.791*eV, 2.809*eV, 2.826*eV, 2.842*eV, 2.861*eV, 2.884*eV, 2.919*eV, 2.946*eV, 2.954*eV, 2.961*eV, 
		                               2.967*eV, 2.974*eV, 2.981*eV, 2.987*eV, 2.994*eV, 3.001*eV, 3.009*eV, 3.018*eV, 3.029*eV, 3.041*eV, 
		                               3.056*eV, 3.083*eV, 3.137*eV, 3.191*eV};

	// taken from LaBr3
	G4double ScintilFast_CeBr3[44] = {0.000, 0.001, 0.001, 0.002, 0.003, 0.006, 0.010, 0.018, 0.033, 0.060, 
		                              0.084, 0.122, 0.175, 0.234, 0.294, 0.356, 0.416, 0.473, 0.533, 0.594, 
		                              0.657, 0.720, 0.784, 0.846, 0.903, 0.962, 1.000, 0.917, 0.857, 0.798, 
		                              0.732, 0.669, 0.604, 0.542, 0.480, 0.422, 0.359, 0.297, 0.237, 0.170, 
		                              0.105, 0.028, 0.004, 0.000};


	G4double photonEnergy_CeBr3_2[2] = {2.004*eV, 3.191*eV};
	G4double RefIndex_CeBr3[2] = {2.09, 2.09}; //Data taken from https://www.advatech-uk.co.uk/cebr3.html
	G4double Absorption_CeBr3[2] = {1.881*cm, 1.881*cm}; // taken from LaBr3

	fCeBr3MPT = new G4MaterialPropertiesTable();
	fCeBr3MPT->AddProperty("RINDEX", photonEnergy_CeBr3, RefIndex_CeBr3, 2);
	fCeBr3MPT->AddProperty("ABSLENGTH", photonEnergy_CeBr3, Absorption_CeBr3, 2);
	fCeBr3MPT->AddProperty("FASTCOMPONENT", photonEnergy_CeBr3, ScintilFast_CeBr3, 44);

	fCeBr3MPT->AddConstProperty("SCINTILLATIONYIELD", 60000/MeV); //Data taken from https://www.advatech-uk.co.uk/cebr3.html
	fCeBr3MPT->AddConstProperty("RESOLUTIONSCALE", 1.0); // Intrinsic resolution

	fCeBr3MPT->AddConstProperty("FASTTIMECONSTANT", 19.0*ns);
	fCeBr3MPT->AddConstProperty("YIELDRATIO",1);

	// electron, proton, and ion yields taken from LaBr3
	G4double particleEnergy_CeBr3[36] = {1E-3, 0.10*MeV, 0.13*MeV, 0.17*MeV, 0.20*MeV, 0.24*MeV, 0.30*MeV, 0.34*MeV, 0.40*MeV, 
		                           0.48*MeV, 0.60*MeV, 0.72*MeV, 0.84*MeV, 1.00*MeV, 1.30*MeV, 1.70*MeV, 2.00*MeV, 2.40*MeV, 
		                           3.00*MeV, 3.40*MeV, 4.00*MeV, 4.80*MeV, 6.00*MeV, 7.20*MeV, 8.40*MeV, 10.00*MeV, 11.00*MeV, 
		                           12.00*MeV, 13.00*MeV, 14.00*MeV, 15.00*MeV, 16.00*MeV, 17.00*MeV, 18.00*MeV, 19.00*MeV, 20.00*MeV};
		                           
	G4double electronYield_CeBr3[36] = {0.0E+00*pEF, 1.0E+03*pEF, 1.3E+03*pEF, 1.7E+03*pEF, 2.0E+03*pEF, 2.4E+03*pEF, 3.0E+03*pEF, 3.4E+03*pEF, 4.0E+03*pEF, 
		                          4.8E+03*pEF, 6.0E+03*pEF, 7.2E+03*pEF, 8.4E+03*pEF, 1.0E+04*pEF, 1.3E+04*pEF, 1.7E+04*pEF, 2.0E+04*pEF, 2.4E+04*pEF, 
		                          3.0E+04*pEF, 3.4E+04*pEF, 4.0E+04*pEF, 4.8E+04*pEF, 6.0E+04*pEF, 7.2E+04*pEF, 8.4E+04*pEF, 1.0E+05*pEF, 1.1E+05*pEF, 
		                          1.2E+05*pEF, 1.3E+05*pEF, 1.4E+05*pEF, 1.5E+05*pEF, 1.6E+05*pEF, 1.7E+05*pEF, 1.8E+05*pEF, 1.9E+05*pEF, 2.0E+05*pEF};
	fCeBr3MPT->AddProperty("ELECTRONSCINTILLATIONYIELD", particleEnergy_CeBr3, electronYield_CeBr3, 36)->SetSpline(true);

	G4double protonYield_CeBr3[36] = {0.6*pSF, 67.1*pSF, 88.6*pSF, 120.7*pSF, 146.5*pSF, 183.8*pSF, 246.0*pSF, 290.0*pSF, 365.0*pSF, 
		                        483.0*pSF, 678.0*pSF, 910.0*pSF, 1175.0*pSF, 1562.0*pSF, 2385.0*pSF, 3660.0*pSF, 4725.0*pSF, 6250.0*pSF, 
		                        8660.0*pSF, 10420.0*pSF, 13270.0*pSF, 17180.0*pSF, 23100.0*pSF, 29500.0*pSF, 36200.0*pSF, 45500.0*pSF, 51826.7*pSF, 
		                        58313.7*pSF, 65047.2*pSF, 72027.4*pSF, 79254.2*pSF, 86727.6*pSF, 94447.6*pSF, 102414.2*pSF, 110627.4*pSF, 119087.2*pSF};
	fCeBr3MPT->AddProperty("PROTONSCINTILLATIONYIELD", particleEnergy_CeBr3, protonYield_CeBr3, 36)->SetSpline(true);

	G4double ionYield_CeBr3[36] = {0.2*pEF, 10.4*pEF, 12.7*pEF, 15.7*pEF, 17.9*pEF, 20.8*pEF, 25.1*pEF, 27.9*pEF, 31.9*pEF, 
		                     36.8*pEF, 43.6*pEF, 50.2*pEF, 56.9*pEF, 65.7*pEF, 81.3*pEF, 101.6*pEF, 116.5*pEF, 136.3*pEF, 
		                     166.2*pEF, 187.1*pEF, 218.6*pEF, 260.5*pEF, 323.5*pEF, 387.5*pEF, 451.5*pEF, 539.9*pEF, 595.5*pEF, 
		                     651.8*pEF, 708.7*pEF, 766.2*pEF, 824.2*pEF, 882.9*pEF, 942.2*pEF, 1002.1*pEF, 1062.6*pEF, 1123.7*pEF};                   
	fCeBr3MPT->AddProperty("IONSCINTILLATIONYIELD", particleEnergy_CeBr3, ionYield_CeBr3, 36)->SetSpline(true);
	fCeBr3MPT->AddProperty("ALPHASCINTILLATIONYIELD", particleEnergy_CeBr3, protonYield_CeBr3, 36)->SetSpline(true);
	fCeBr3MPT->AddProperty("GAMMASCINTILLATIONYIELD", particleEnergy_CeBr3, electronYield_CeBr3, 36)->SetSpline(true);
	fCeBr3->SetMaterialPropertiesTable(fCeBr3MPT);


 /*
	/////////////////////////////////////////////////////////////////
	// Cesium Iodide CsI
	/////////////////////////////////////////////////////////////////

    fCsI = new G4Material("CsI", 4.51*g/cm3, 2);
    fCsI->AddElement(fCs, 0.5);
    fCsI->AddElement(fI, 0.5);


	// taken from LaBr3 NEED TO UPDATE
	G4double photonEnergy_CsI[44] = {2.004*eV, 2.058*eV, 2.112*eV, 2.166*eV, 2.220*eV, 2.274*eV, 2.328*eV, 2.382*eV, 2.436*eV, 2.490*eV, 
		                               2.517*eV, 2.552*eV, 2.585*eV, 2.613*eV, 2.635*eV, 2.656*eV, 2.686*eV, 2.720*eV, 2.749*eV, 2.772*eV, 
		                               2.791*eV, 2.809*eV, 2.826*eV, 2.842*eV, 2.861*eV, 2.884*eV, 2.919*eV, 2.946*eV, 2.954*eV, 2.961*eV, 
		                               2.967*eV, 2.974*eV, 2.981*eV, 2.987*eV, 2.994*eV, 3.001*eV, 3.009*eV, 3.018*eV, 3.029*eV, 3.041*eV, 
		                               3.056*eV, 3.083*eV, 3.137*eV, 3.191*eV};

	// taken from LaBr3 NEED TO UPDATE
	G4double ScintilFast_CsI[44] = {0.000, 0.001, 0.001, 0.002, 0.003, 0.006, 0.010, 0.018, 0.033, 0.060, 
		                              0.084, 0.122, 0.175, 0.234, 0.294, 0.356, 0.416, 0.473, 0.533, 0.594, 
		                              0.657, 0.720, 0.784, 0.846, 0.903, 0.962, 1.000, 0.917, 0.857, 0.798, 
		                              0.732, 0.669, 0.604, 0.542, 0.480, 0.422, 0.359, 0.297, 0.237, 0.170, 
		                              0.105, 0.028, 0.004, 0.000};


	// NEED TO UPDATE THESE THREE
	G4double photonEnergy_CsI_2[2] = {2.004*eV, 3.191*eV};
	G4double RefIndex_CsI[2] = {2.09, 2.09}; //Data taken from https://www.advatech-uk.co.uk/cebr3.html
	G4double Absorption_CsI[2] = {1.881*cm, 1.881*cm}; // taken from LaBr3

	fCsIMPT = new G4MaterialPropertiesTable();
	fCsIMPT->AddProperty("RINDEX", photonEnergy_CsI, RefIndex_CsI, 2);
	fCsIMPT->AddProperty("ABSLENGTH", photonEnergy_CsI, Absorption_CsI, 2);
	fCsIMPT->AddProperty("FASTCOMPONENT", photonEnergy_CsI, ScintilFast_CsI, 44);

	// NEED TO UPDATE SCINTILLATION YIELD
	fCsIMPT->AddConstProperty("SCINTILLATIONYIELD", 60000/MeV); //Data taken from https://www.advatech-uk.co.uk/cebr3.html
	fCsIMPT->AddConstProperty("RESOLUTIONSCALE", 1.0); // Intrinsic resolution

	// NEED TO UPDATE FAST TIME CONSTANT
	fCsIMPT->AddConstProperty("FASTTIMECONSTANT", 19.0*ns);
	fCsIMPT->AddConstProperty("YIELDRATIO",1);

	// electron, proton, and ion yields taken from LaBr3 NEED TO UPDATE!!
	G4double particleEnergy_CsI[36] = {1E-3, 0.10*MeV, 0.13*MeV, 0.17*MeV, 0.20*MeV, 0.24*MeV, 0.30*MeV, 0.34*MeV, 0.40*MeV, 
		                           0.48*MeV, 0.60*MeV, 0.72*MeV, 0.84*MeV, 1.00*MeV, 1.30*MeV, 1.70*MeV, 2.00*MeV, 2.40*MeV, 
		                           3.00*MeV, 3.40*MeV, 4.00*MeV, 4.80*MeV, 6.00*MeV, 7.20*MeV, 8.40*MeV, 10.00*MeV, 11.00*MeV, 
		                           12.00*MeV, 13.00*MeV, 14.00*MeV, 15.00*MeV, 16.00*MeV, 17.00*MeV, 18.00*MeV, 19.00*MeV, 20.00*MeV};
		                           
	G4double electronYield_CsI[36] = {0.0E+00*pEF, 1.0E+03*pEF, 1.3E+03*pEF, 1.7E+03*pEF, 2.0E+03*pEF, 2.4E+03*pEF, 3.0E+03*pEF, 3.4E+03*pEF, 4.0E+03*pEF, 
		                          4.8E+03*pEF, 6.0E+03*pEF, 7.2E+03*pEF, 8.4E+03*pEF, 1.0E+04*pEF, 1.3E+04*pEF, 1.7E+04*pEF, 2.0E+04*pEF, 2.4E+04*pEF, 
		                          3.0E+04*pEF, 3.4E+04*pEF, 4.0E+04*pEF, 4.8E+04*pEF, 6.0E+04*pEF, 7.2E+04*pEF, 8.4E+04*pEF, 1.0E+05*pEF, 1.1E+05*pEF, 
		                          1.2E+05*pEF, 1.3E+05*pEF, 1.4E+05*pEF, 1.5E+05*pEF, 1.6E+05*pEF, 1.7E+05*pEF, 1.8E+05*pEF, 1.9E+05*pEF, 2.0E+05*pEF};
	fCsIMPT->AddProperty("ELECTRONSCINTILLATIONYIELD", particleEnergy_CsI, electronYield_CsI, 36)->SetSpline(true);

	G4double protonYield_CsI[36] = {0.6*pSF, 67.1*pSF, 88.6*pSF, 120.7*pSF, 146.5*pSF, 183.8*pSF, 246.0*pSF, 290.0*pSF, 365.0*pSF, 
		                        483.0*pSF, 678.0*pSF, 910.0*pSF, 1175.0*pSF, 1562.0*pSF, 2385.0*pSF, 3660.0*pSF, 4725.0*pSF, 6250.0*pSF, 
		                        8660.0*pSF, 10420.0*pSF, 13270.0*pSF, 17180.0*pSF, 23100.0*pSF, 29500.0*pSF, 36200.0*pSF, 45500.0*pSF, 51826.7*pSF, 
		                        58313.7*pSF, 65047.2*pSF, 72027.4*pSF, 79254.2*pSF, 86727.6*pSF, 94447.6*pSF, 102414.2*pSF, 110627.4*pSF, 119087.2*pSF};
	fCsIMPT->AddProperty("PROTONSCINTILLATIONYIELD", particleEnergy_CsI, protonYield_CsI, 36)->SetSpline(true);

	G4double ionYield_CsI[36] = {0.2*pEF, 10.4*pEF, 12.7*pEF, 15.7*pEF, 17.9*pEF, 20.8*pEF, 25.1*pEF, 27.9*pEF, 31.9*pEF, 
		                     36.8*pEF, 43.6*pEF, 50.2*pEF, 56.9*pEF, 65.7*pEF, 81.3*pEF, 101.6*pEF, 116.5*pEF, 136.3*pEF, 
		                     166.2*pEF, 187.1*pEF, 218.6*pEF, 260.5*pEF, 323.5*pEF, 387.5*pEF, 451.5*pEF, 539.9*pEF, 595.5*pEF, 
		                     651.8*pEF, 708.7*pEF, 766.2*pEF, 824.2*pEF, 882.9*pEF, 942.2*pEF, 1002.1*pEF, 1062.6*pEF, 1123.7*pEF};                   
	fCsIMPT->AddProperty("IONSCINTILLATIONYIELD", particleEnergy_CsI, ionYield_CsI, 36)->SetSpline(true);
	fCsIMPT->AddProperty("ALPHASCINTILLATIONYIELD", particleEnergy_CsI, protonYield_CsI, 36)->SetSpline(true);
	fCsIMPT->AddProperty("GAMMASCINTILLATIONYIELD", particleEnergy_CsI, electronYield_CsI, 36)->SetSpline(true);
	fCsI->SetMaterialPropertiesTable(fCsIMPT);
 */


	// Update the material dictionary
	materialList["ej200"] = fEJ200;
	materialList["ej276"] = fEJ276;
	materialList["ogs"] = fOGS;
 	materialList["yso"] = fYSO;
	materialList["lyso"] = fLYSO;
	materialList["labr3"] = fLaBr3;
	materialList["yap"] = fYAP;
	materialList["gagg"] = fGAGG;
	materialList["cebr3"] = fCeBr3;
	//materialList["csi"] = fCsI;
	
	scintsAreDefined = true;
}


// Formula with offset
double nDetMaterials::dLdE(G4double dEdx, G4double a, G4double b, G4double c) {
	return a/(1+b*dEdx) + c;
}


void nDetMaterials::setIonSourceType(G4int mass) {

	// need to make sure the scintillator actually has a light yield matrix first
	if(lightYieldTable.find(detectorMaterialName) == lightYieldTable.end()){
		std::cout << "No light yield matrix has been defined for " << detectorMaterialName << ". Source will default to alpha particle light yield." << std::endl;
		lightYieldTable.clear();
		return;
	}

	std::vector<G4double> particleEnergy = {0.01*MeV, 0.011*MeV, 0.012*MeV, 0.013*MeV, 0.014*MeV, 0.015*MeV, 0.016*MeV, 0.017*MeV, 0.018*MeV, 0.02*MeV, 0.0225*MeV, 0.025*MeV, 0.0275*MeV, 0.03*MeV, 0.0325*MeV, 0.035*MeV, 0.0375*MeV, 0.04*MeV, 0.045*MeV, 0.05*MeV, 0.055*MeV, 0.06*MeV, 0.065*MeV, 0.07*MeV, 0.08*MeV, 0.09*MeV, 0.1*MeV, 0.11*MeV, 0.12*MeV, 0.13*MeV, 0.14*MeV, 0.15*MeV, 0.16*MeV, 0.17*MeV, 0.18*MeV, 0.2*MeV, 0.225*MeV, 0.25*MeV, 0.275*MeV, 0.3*MeV, 0.325*MeV, 0.35*MeV, 0.375*MeV, 0.4*MeV, 0.45*MeV, 0.5*MeV, 0.55*MeV, 0.6*MeV, 0.65*MeV, 0.7*MeV, 0.8*MeV, 0.9*MeV, 1*MeV, 1.1*MeV, 1.2*MeV, 1.3*MeV, 1.4*MeV, 1.5*MeV, 1.6*MeV, 1.7*MeV, 1.8*MeV, 2*MeV, 2.25*MeV, 2.5*MeV, 2.75*MeV, 3*MeV, 3.25*MeV, 3.5*MeV, 3.75*MeV, 4*MeV, 4.5*MeV, 5*MeV, 5.5*MeV, 6*MeV, 6.5*MeV, 7*MeV, 8*MeV, 9*MeV, 10*MeV, 11*MeV, 12*MeV, 13*MeV, 14*MeV, 15*MeV, 16*MeV, 17*MeV, 18*MeV, 20*MeV, 22.5*MeV, 25*MeV, 27.5*MeV, 30*MeV, 32.5*MeV, 35*MeV, 37.5*MeV, 40*MeV, 45*MeV, 50*MeV, 55*MeV, 60*MeV, 65*MeV, 70*MeV, 80*MeV, 90*MeV, 100*MeV, 110*MeV, 120*MeV, 130*MeV, 140*MeV, 150*MeV, 160*MeV, 170*MeV, 180*MeV, 200*MeV, 225*MeV, 250*MeV, 275*MeV, 300*MeV, 325*MeV, 350*MeV, 375*MeV, 400*MeV, 450*MeV, 500*MeV, 550*MeV, 600*MeV, 650*MeV, 700*MeV, 800*MeV, 900*MeV, 1000*MeV, 1100*MeV, 1200*MeV, 1300*MeV, 1400*MeV, 1500*MeV, 1600*MeV, 1700*MeV, 1800*MeV, 2000*MeV, 2250*MeV, 2500*MeV, 2750*MeV, 3000*MeV, 3250*MeV, 3500*MeV, 3750*MeV, 4000*MeV, 4500*MeV, 5000*MeV, 5500*MeV, 6000*MeV, 6500*MeV, 7000*MeV, 8000*MeV, 9000*MeV, 10000*MeV};


	if (mass < 1) {
		std::cout << "Mass must be at least 1. Source will default to alpha particle light yield." << std::endl;
	}

	G4int upper = 999, upper_index = 0, lower = -1, lower_index = 0;
	for (int i=0; i < lightYieldTable[detectorMaterialName].size(); i++) {
		if (lightYieldTable[detectorMaterialName][i][0] == mass) {

			// if the mass is in the light yield matrix, then use the data from the table
			if (detectorMaterialName == "lyso") {
				fLYSOMPT->AddProperty("IONSCINTILLATIONYIELD", &particleEnergy[0], &lightYieldTable[detectorMaterialName][i][1], 157)->SetSpline(true);
				fLYSO->SetMaterialPropertiesTable(fLYSOMPT);

			}

			// clear light yield table since it isn't needed anymore
			lightYieldTable.clear();
			return;
		}
		else if (lightYieldTable[detectorMaterialName][i][0] < mass && lightYieldTable[detectorMaterialName][i][0] > lower) {
			lower = lightYieldTable[detectorMaterialName][i][0];
			lower_index = i;
		}
		else if (lightYieldTable[detectorMaterialName][i][0] > mass && lightYieldTable[detectorMaterialName][i][0] < upper) {
			upper = lightYieldTable[detectorMaterialName][i][0];
			upper_index = i;
		}
	}

	// if mass not in matrix, perform linear interpolation
	std::vector<G4double> ionYield(particleEnergy.size());

	if (upper == 999) { // Mass larger than largest mass for which light yield is defined
		std::cout << mass << " is larger mass than all ions with defined light yield. Linear interpolation cannot be performed. Source will default to alpha particle light yield." << std::endl;
		lightYieldTable.clear();
		return;
	}

	std::cout << "Light yield for A = " << mass << " has not been defined. A linear interpolation will be performed." << std::endl;
	for (int i=0; i < particleEnergy.size(); i++) {
		ionYield[i] = lightYieldTable[detectorMaterialName][upper_index][i+1] + ((double)upper-mass)/(upper-lower) * (lightYieldTable[detectorMaterialName][lower_index][i+1]-lightYieldTable[detectorMaterialName][upper_index][i+1]);
	}


	if (detectorMaterialName == "lyso") {
		fLYSOMPT->AddProperty("IONSCINTILLATIONYIELD", &particleEnergy[0], &ionYield[0], 157)->SetSpline(true);
		fLYSO->SetMaterialPropertiesTable(fLYSOMPT);
	}

	// clear light yield table since it isn't needed anymore
	lightYieldTable.clear();
}