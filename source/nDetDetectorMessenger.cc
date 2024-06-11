#include "nDetDetectorMessenger.hh"

#include "nDetDetector.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4SystemOfUnits.hh"
#include <regex>
using namespace std;

void nDetDetectorMessenger::addAllCommands(){
///////////////////////////////////////////////////////////////////////////////
// Detector commands
///////////////////////////////////////////////////////////////////////////////

	addDirectory("/nDet/detector/", "Detector geometry control");

	addCommand(new G4UIcmdWithAString("/nDet/detector/setPmtDimensions", this));
	addGuidance("Defines the size of the SiPMs in mm. SYNTAX: setPmtDimensions <sizeX> [sizeY]");

	addCommand(new G4UIcmdWithADouble("/nDet/detector/setDetectorLength", this));
	addGuidance("Defines the length (Z) of the detector in cm");

	addCommand(new G4UIcmdWithADouble("/nDet/detector/setDetectorWidth", this));
	addGuidance("Defines the width (X) of the detector in cm");

	addCommand(new G4UIcmdWithADouble("/nDet/detector/setDetectorThickness", this));
	addGuidance("Defines the thickness of the plastic in mm");

	addCommand(new G4UIcmdWithADouble("/nDet/detector/setMylarThickness", this));
	addGuidance("Defines the thickness of the plastic the mylar in mm (0 for no mylar)");

	addCommand(new G4UIcmdWithADouble("/nDet/detector/setGreaseThickness", this));
	addGuidance("Defines the thickness of the optical grease layer in mm (0 for no grease)");

	addCommand(new G4UIcmdWithADouble("/nDet/detector/setTrapezoidLength", this));
	addGuidance("Defines the length of the trapezoidal part of ellipse in cm");

	addCommand(new G4UIcmdWithAnInteger("/nDet/detector/setNumColumns", this));
	addGuidance("Set the number of columns in a segmented scintillator.");
	
	addCommand(new G4UIcmdWithAnInteger("/nDet/detector/setNumRows", this));
	addGuidance("Set the number of rows in a segmented scintillator.");
	
	addCommand(new G4UIcmdWithAnInteger("/nDet/detector/setPmtColumns", this));
	addGuidance("Set the number of anode columns in a segmented PSPMT.");
	
	addCommand(new G4UIcmdWithAnInteger("/nDet/detector/setPmtRows", this));
	addGuidance("Set the number of anode rows in a segmented PSPMT.");
	
	addCommand(new G4UIcmdWithADouble("/nDet/detector/setDiffuserLength", this));
	addGuidance("Defines the length of the straight light diffuser in cm");
	
	addCommand(new G4UIcmdWithAString("/nDet/detector/setWrapping", this));
	addGuidance("Set the material to use for reflective wrapping");	

	addCommand(new G4UIcmdWithAString("/nDet/detector/setMaterial", this));
	addGuidance("Set the material to use for detector construction.");
		
	addCommand(new G4UIcmdWith3VectorAndUnit("/nDet/detector/setPosition", this)); // position of source in cartesian coordinates (x, y, z).
	addGuidance("Set the position of the detector in cartesian coordinates (x, y, z)");

	addCommand(new G4UIcmdWith3Vector("/nDet/detector/setCylindrical", this)); // position of source in cylindrical coordinates (r, theta, y).
	addGuidance("Set the position of the detector in cylindrical coordinates (r, theta, y) where r and y are in cm and theta is in degrees");

	addCommand(new G4UIcmdWith3Vector("/nDet/detector/setSpherical", this)); // position of source in spherical coordinates (r, theta, phi).
	addGuidance("Set the position of the detector in spherical coordinates (r, theta, phi) where r is in cm and theta and phi are in degrees");
	
	addCommand(new G4UIcmdWith3Vector("/nDet/detector/setRotation", this)); // rotation of detector.
	addGuidance("Set the rotation of the detector by specifying angles about the x, y, and z axes (in deg)");

	addCommand(new G4UIcmdWithADouble("/nDet/detector/setWindowThickness", this));
	addGuidance("Defines the thickness of the pspmt quartz window in mm (0 for no window)");

	addCommand(new G4UIcmdWithAString("/nDet/detector/setPolished", this));
	addGuidance("Enable or disable polished optical grease faces (disabled by default)");
	addCandidates("true false");

	addCommand(new G4UIcmdWithoutParameter("/nDet/detector/setStart", this));
	addGuidance("Mark the current detector as a start detector");

	addCommand(new G4UIcmdWithAString("/nDet/detector/setSquarePmt", this));
	addGuidance("Enable or disable square PMTs (enabled by default). If disabled, circular PMTs will be used");
	addCandidates("true false");

	addCommand(new G4UIcmdWithADouble("/nDet/detector/setSegmentWidth", this));
	addGuidance("Set the width of scintillator segments for segmented detectors (in mm)");

	addCommand(new G4UIcmdWithADouble("/nDet/detector/setSegmentHeight", this));
	addGuidance("Set the height of scintillator segments for segmented detectors (in mm)");

	addCommand(new G4UIcmdWithADouble("/nDet/detector/setTrapezoidAngle", this));
	addGuidance("Defines the angle of the trapezoidal part of ellipse wrt the edge of the rectangular body (in degrees)");

	addCommand(new G4UIcmdWithADouble("/nDet/detector/setDetectorHeight", this));
	addGuidance("Defines the height (Y) of the detector in cm");	

///////////////////////////////////////////////////////////////////////////////
// Implant commands
///////////////////////////////////////////////////////////////////////////////

	addDirectory("/nDet/implant/", "Detector geometry control");

	addCommand(new G4UIcmdWithAString("/nDet/implant/setPmtDimensions", this));
	addGuidance("Defines the size of the SiPMs in mm. SYNTAX: setPmtDimensions <sizeX> [sizeY]");

	addCommand(new G4UIcmdWithADouble("/nDet/implant/setDetectorLength", this));
	addGuidance("Defines the length (Z) of the implant in cm");

	addCommand(new G4UIcmdWithADouble("/nDet/implant/setDetectorWidth", this));
	addGuidance("Defines the width (X) of the implant in cm");

	addCommand(new G4UIcmdWithADouble("/nDet/implant/setDetectorThickness", this));
	addGuidance("Defines the thickness of the plastic in mm");

	addCommand(new G4UIcmdWithADouble("/nDet/implant/setMylarThickness", this));
	addGuidance("Defines the thickness of the plastic the mylar in mm (0 for no mylar)");

	addCommand(new G4UIcmdWithADouble("/nDet/implant/setGreaseThickness", this));
	addGuidance("Defines the thickness of the optical grease layer in mm (0 for no grease)");

	addCommand(new G4UIcmdWithADouble("/nDet/implant/setTrapezoidLength", this));
	addGuidance("Defines the length of the trapezoidal part of ellipse in cm");

	addCommand(new G4UIcmdWithAnInteger("/nDet/implant/setNumColumns", this));
	addGuidance("Set the number of columns in a segmented scintillator.");
	
	addCommand(new G4UIcmdWithAnInteger("/nDet/implant/setNumRows", this));
	addGuidance("Set the number of rows in a segmented scintillator.");
	
	addCommand(new G4UIcmdWithAnInteger("/nDet/implant/setPmtColumns", this));
	addGuidance("Set the number of anode columns in a segmented PSPMT.");
	
	addCommand(new G4UIcmdWithAnInteger("/nDet/implant/setPmtRows", this));
	addGuidance("Set the number of anode rows in a segmented PSPMT.");

	addCommand(new G4UIcmdWithADouble("/nDet/implant/setDiffuserLength", this));
	addGuidance("Defines the length of the straight light diffuser in cm");
	
	addCommand(new G4UIcmdWithAString("/nDet/implant/setWrapping", this));
	addGuidance("Set the material to use for reflective wrapping");	

	addCommand(new G4UIcmdWithAString("/nDet/implant/setMaterial", this));
	addGuidance("Set the material to use for implant construction.");
		
	addCommand(new G4UIcmdWith3VectorAndUnit("/nDet/implant/setPosition", this)); // position of source in cartesian coordinates (x, y, z).
	addGuidance("Set the position of the implant in cartesian coordinates (x, y, z)");

	addCommand(new G4UIcmdWith3Vector("/nDet/implant/setCylindrical", this)); // position of source in cylindrical coordinates (r, theta, y).
	addGuidance("Set the position of the implant in cylindrical coordinates (r, theta, y) where r and y are in cm and theta is in degrees");

	addCommand(new G4UIcmdWith3Vector("/nDet/implant/setSpherical", this)); // position of source in spherical coordinates (r, theta, phi).
	addGuidance("Set the position of the implant in spherical coordinates (r, theta, phi) where r is in cm and theta and phi are in degrees");
	
	addCommand(new G4UIcmdWith3Vector("/nDet/implant/setRotation", this)); // rotation of implant.
	addGuidance("Set the rotation of the implant by specifying angles about the x, y, and z axes (in deg)");

	addCommand(new G4UIcmdWithADouble("/nDet/implant/setWindowThickness", this));
	addGuidance("Defines the thickness of the pspmt quartz window in mm (0 for no window)");

	addCommand(new G4UIcmdWithAString("/nDet/implant/setPolished", this));
	addGuidance("Enable or disable polished optical grease faces (disabled by default)");
	addCandidates("true false");

	addCommand(new G4UIcmdWithoutParameter("/nDet/implant/setStart", this));
	addGuidance("Mark the current implant as a start implant");

	addCommand(new G4UIcmdWithAString("/nDet/implant/setSquarePmt", this));
	addGuidance("Enable or disable square PMTs (enabled by default). If disabled, circular PMTs will be used");
	addCandidates("true false");

	addCommand(new G4UIcmdWithADouble("/nDet/implant/setSegmentWidth", this));
	addGuidance("Set the width of scintillator segments for segmented implants (in mm)");

	addCommand(new G4UIcmdWithADouble("/nDet/implant/setSegmentHeight", this));
	addGuidance("Set the height of scintillator segments for segmented implants (in mm)");
	
	addCommand(new G4UIcmdWithADouble("/nDet/implant/setTrapezoidAngle", this));
	addGuidance("Defines the angle of the trapezoidal part of ellipse wrt the edge of the rectangular body (in degrees)");

	addCommand(new G4UIcmdWithADouble("/nDet/implant/setDetectorHeight", this));
	addGuidance("Defines the height (Y) of the implant in cm");	

	addCommand(new G4UIcmdWithABool("/nDet/implant/setReconLogic",this));
	addGuidance("Sets the logic used for calculating the reconstructed position, 0 = anger logic, 1 = sipm logic");

	addCommand(new G4UIcmdWithADouble("/nDet/implant/setPmtGapThickness", this));
	addGuidance("Set the gap thickness in mm of a segmented PMT");

	addCommand(new G4UIcmdWithABool("/nDet/implant/useFrontWrapping",this));
	addGuidance("Set whether to apply reflective wrapping to the front of the detector");

	/////////////////////////////////////////////
	// Segmented Light Guide commands
	/////////////////////////////////////////////

	addDirectory("/nDet/segLightGuide/","Segmented Light Guide Geometry");

	addCommand(new G4UIcmdWithADouble("/nDet/segLightGuide/thickness",this));
	addGuidance("Defines the z thickness of the segmented light guide");

	addCommand(new G4UIcmdWithADouble("/nDet/segLightGuide/frontLength",this));
	addGuidance("Defines the length (x dimension) of the front of the segmented light guide");

	addCommand(new G4UIcmdWithADouble("/nDet/segLightGuide/frontWidth",this));
	addGuidance("Defines the width (y dimension) of the front of the segmented light guide");

	addCommand(new G4UIcmdWithADouble("/nDet/segLightGuide/backLength",this));
	addGuidance("Defines the length (x dimension) of the back of the segmented light guide");

	addCommand(new G4UIcmdWithADouble("/nDet/segLightGuide/backWidth",this));
	addGuidance("Defines the width (y dimension) of the back of the segmented light guide");

	addCommand(new G4UIcmdWithAnInteger("/nDet/segLightGuide/segx",this));
	addGuidance("Defines the number of segments along the x direction");

	addCommand(new G4UIcmdWithAnInteger("/nDet/segLightGuide/segy",this));
	addGuidance("Defines the number of segments along the y direction");

	addCommand(new G4UIcmdWithAString("/nDet/segLightGuide/setWrapping", this));
	addGuidance("Set the material to use for reflective wrapping on the light guide");	

}

void nDetDetectorMessenger::SetNewChildValue(G4UIcommand* command, G4String newValue){
	size_t index;
	if(!findCommand(command, newValue, index)) return;

	else if(index == 0) {
		fDetector->SetPmtDimension(newValue);
	}
	else if(index == 1) {
		G4double length = command->ConvertToDouble(newValue);
		fDetector->SetDetectorLength(length*cm);
	}
	else if(index == 2) {
		G4double length = command->ConvertToDouble(newValue);
		fDetector->SetDetectorWidth(length*cm);
	}
	else if(index == 3) {
		G4double thickness = command->ConvertToDouble(newValue);
		fDetector->SetDetectorHeight(thickness*mm);
	}
	else if(index == 4) {
		G4double val = command->ConvertToDouble(newValue);
		fDetector->SetMylarThickness(val*mm);
	}	
	else if(index == 5) {
		G4double val = command->ConvertToDouble(newValue);
		fDetector->SetGreaseThickness(val*mm);
	}
	else if(index == 6) {
		G4double val = command->ConvertToDouble(newValue);
		fDetector->SetTrapezoidLength(val*cm);
	}
	else if(index == 7){
		G4int val = command->ConvertToInt(newValue);
		fDetector->SetNumColumns(val);
	}
	else if(index == 8){
		G4int val = command->ConvertToInt(newValue);
		fDetector->SetNumRows(val);
	}
	else if(index == 9){
		G4int val = command->ConvertToInt(newValue);
		fDetector->SetNumPmtColumns(val);
	}
	else if(index == 10){
		G4int val = command->ConvertToInt(newValue);
		fDetector->SetNumPmtRows(val);
	}
	else if(index == 11){
		G4double val = command->ConvertToDouble(newValue);
		fDetector->SetDiffuserLength(val*cm);
	}
	else if(index == 12){
		fDetector->SetWrappingMaterial(newValue);
	}
	else if(index == 13){
		fDetector->SetDetectorMaterial(newValue);
	}	
	else if(index == 14){
		G4ThreeVector val = command->ConvertToDimensioned3Vector(newValue);
		fDetector->SetPosition(val);
	}
	else if(index == 15){
		G4ThreeVector val = command->ConvertTo3Vector(newValue);
		fDetector->SetPositionCylindrical(val);
	}
	else if(index == 16){
		G4ThreeVector val = command->ConvertTo3Vector(newValue);
		fDetector->SetPositionSpherical(val);
	}
	else if(index == 17){
		G4ThreeVector val = command->ConvertTo3Vector(newValue);
		fDetector->SetRotation(val);
	}
	else if(index == 18) {
		G4double val = command->ConvertToDouble(newValue);
		fDetector->SetWindowThickness(val*mm);
	}
	else if(index == 19){
		fDetector->SetPolishedInterface((newValue == "true") ? true : false);
	}
	else if(index == 20){
		fDetector->SetAsStart(true);
	}
	else if(index == 21){
		fDetector->SetSquarePMTs((newValue == "true") ? true : false);
	}
	else if(index == 22){
		fDetector->SetSegmentWidth(command->ConvertToDouble(newValue));
	}
	else if(index == 23){
		fDetector->SetSegmentHeight(command->ConvertToDouble(newValue));
	}
	else if(index == 24){
		fDetector->SetTrapezoidAngle(command->ConvertToDouble(newValue));
	}
	else if(index == 25) {
		G4double thickness = command->ConvertToDouble(newValue);
		fDetector->SetDetectorHeight(thickness*cm);
	}
	
	//For implant values

	else if(index == 26) {
		fDetector->SetPmtDimension(newValue);
	}
	else if(index == 27) {
		G4double length = command->ConvertToDouble(newValue);
		fDetector->SetDetectorLength(length*cm);
	}
	else if(index == 28) {
		G4double length = command->ConvertToDouble(newValue);
		fDetector->SetDetectorWidth(length*cm);
	}
	else if(index == 29) {
		G4double thickness = command->ConvertToDouble(newValue);
		fDetector->SetDetectorHeight(thickness*mm);
	}
	else if(index == 30) {
		G4double val = command->ConvertToDouble(newValue);
		fDetector->SetMylarThickness(val*mm);
	}	
	else if(index == 31) {
		G4double val = command->ConvertToDouble(newValue);
		fDetector->SetGreaseThickness(val*mm);
	}
	else if(index == 32) {
		G4double val = command->ConvertToDouble(newValue);
		fDetector->SetTrapezoidLength(val*cm);
	}
	else if(index == 33){
		G4int val = command->ConvertToInt(newValue);
		fDetector->SetNumColumns(val);
	}
	else if(index == 34){
		G4int val = command->ConvertToInt(newValue);
		fDetector->SetNumRows(val);
	}
	else if(index == 35){
		G4int val = command->ConvertToInt(newValue);
		fDetector->SetNumPmtColumns(val);
	}
	else if(index == 36){
		G4int val = command->ConvertToInt(newValue);
		fDetector->SetNumPmtRows(val);
	}
	else if(index == 37){
		G4double val = command->ConvertToDouble(newValue);
		fDetector->SetDiffuserLength(val*cm);
	}
	else if(index == 38){
		fDetector->SetWrappingMaterial(newValue);
	}
	else if(index == 39){
		fDetector->SetDetectorMaterial(newValue);
	}	
	else if(index == 40){
		G4ThreeVector val = command->ConvertToDimensioned3Vector(newValue);
		fDetector->SetPosition(val);
	}
	else if(index == 41){
		G4ThreeVector val = command->ConvertTo3Vector(newValue);
		fDetector->SetPositionCylindrical(val);
	}
	else if(index == 42){
		G4ThreeVector val = command->ConvertTo3Vector(newValue);
		fDetector->SetPositionSpherical(val);
	}
	else if(index == 43){
		G4ThreeVector val = command->ConvertTo3Vector(newValue);
		fDetector->SetRotation(val);
	}
	else if(index == 44) {
		G4double val = command->ConvertToDouble(newValue);
		fDetector->SetWindowThickness(val*mm);
	}
	else if(index == 45){
		fDetector->SetPolishedInterface((newValue == "true") ? true : false);
	}
	else if(index == 46){
		fDetector->SetAsStart(true);
	}
	else if(index == 47){
		fDetector->SetSquarePMTs((newValue == "true") ? true : false);
	}
	else if(index == 48){
		fDetector->SetSegmentWidth(command->ConvertToDouble(newValue));
	}
	else if(index == 49){
		fDetector->SetSegmentHeight(command->ConvertToDouble(newValue));
	}
	else if(index == 50){
		fDetector->SetTrapezoidAngle(command->ConvertToDouble(newValue));
	}
	else if(index == 51) {
		G4double thickness = command->ConvertToDouble(newValue);
		fDetector->SetDetectorHeight(thickness*cm);
	}
	else if(index == 52){
		bool val = command->ConvertToBool(newValue);
		std::cout<<"Setting Recon Logic to "<<val<<std::endl;
		fDetector->setReconLogic(val);
	}
	else if(index == 53) {
		G4double val = command->ConvertToDouble(newValue);
		std::cout<<"Setting PMT gap thickness to "<<val<<" mm"<<std::endl;
		fDetector->SetPmtGapThickness(val*mm);
	}
	else if (index == 54) {
		bool val = command->ConvertToBool(newValue);
		std::cout<<"Setting front wrapping to "<<val<<std::endl;
		fDetector->SetFrontWrapping(val);
	}

	else if(index == 55){
		G4double val = command->ConvertToDouble(newValue);
		fDetector->setSegZThick(val);
		std::cout<<"Setting Z thickness to "<<val<<std::endl;
	}
	else if(index == 56){
		G4double val = command->ConvertToDouble(newValue);
		fDetector->setSegTopThick(val);
	}
	else if(index == 57){
		G4double val = command->ConvertToDouble(newValue);
		fDetector->setSegTopWidth(val);
	}
	else if(index == 58){
		G4double val = command->ConvertToDouble(newValue);
		fDetector->setSegBotThick(val);
	}
	else if(index == 59){
		G4double val = command->ConvertToDouble(newValue);
		fDetector->setSegBotWidth(val);
	}
	else if(index == 60){
		G4int val = command->ConvertToInt(newValue);
		fDetector->setSegX(val);
	}
	else if(index == 61){
		G4int val = command->ConvertToInt(newValue);
		fDetector->setSegY(val);
	}
	else if(index == 62){

	}
}
