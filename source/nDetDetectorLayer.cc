#include "G4SystemOfUnits.hh"

#include "nDetDetector.hh"
#include "nDetDetectorLayer.hh"
#include "optionHandler.hh" // split_str


///////////////////////////////////////////////////////////////////////////////
// class greaseLayer
///////////////////////////////////////////////////////////////////////////////

bool greaseLayer::decodeArgs(){
	// Expects a space-delimited string of the form:
	// "addGreaseLayer width(mm) height(mm) thickness(mm) material"
	x = strtod(args.at(0).c_str(), NULL);
	y = strtod(args.at(1).c_str(), NULL);
	if (nUserArgs >= 4) {
		thickness = strtod(args.at(2).c_str(), NULL);
		grMaterial = args.at(3).c_str();
	}
	if(nUserArgs >= 3) {
		if (strtod(args.at(2).c_str(), NULL) > 0) {
			thickness = strtod(args.at(2).c_str(), NULL);
		}
		else
			grMaterial = args.at(2).c_str();
	}
	size = G4ThreeVector(x, y, thickness);
	return true;
}

void greaseLayer::construct(nDetDetector *obj){
	obj->applyGreaseLayer(x, y, thickness);
}

void greaseLayer::construct(nDetImplant *obj){
	obj->applyGreaseLayer(x, y, thickness, grMaterial);
}

std::string greaseLayer::syntaxStr() const {
	return std::string("addGreaseLayer <width> <height> [thickness] [material]");
}

///////////////////////////////////////////////////////////////////////////////
// class diffuserLayer
///////////////////////////////////////////////////////////////////////////////

diffuserLayer::diffuserLayer(const G4String &arg_) : nDetWorldObject(arg_, 3), x(0), y(0), thickness(0) { 
	material = "G4_SILICON_DIOXIDE"; 
}

bool diffuserLayer::decodeArgs(){
	// Expects a space-delimited string of the form:
	//  "addDiffuserLayer width(mm) height(mm) thickness(mm) material"
	x = strtod(args.at(0).c_str(), NULL);
	y = strtod(args.at(1).c_str(), NULL);
	thickness = strtod(args.at(2).c_str(), NULL);
	size = G4ThreeVector(x, y, thickness);
	return true;
}

void diffuserLayer::construct(nDetDetector *obj){
	obj->applyDiffuserLayer(x, y, thickness);
}

void diffuserLayer::construct(nDetImplant *obj){
	obj->applyDiffuserLayer(x, y, thickness);
}

std::string diffuserLayer::syntaxStr() const {
	return std::string("addDiffuserLayer <width> <height> <thickness> [material=G4_SILICON_DIOXIDE]");
}

///////////////////////////////////////////////////////////////////////////////
// class lightGuideLayer
///////////////////////////////////////////////////////////////////////////////

lightGuideLayer::lightGuideLayer(const G4String &arg_) : nDetWorldObject(arg_, 5), x1(0), x2(0), y1(0), y2(0), thickness(0) { 
	material = "G4_SILICON_DIOXIDE"; 
}

bool lightGuideLayer::decodeArgs(){
	// Expects a space-delimited string of the form:
	//  "addLightGuide width1(mm) width2(mm) height1(mm) height2(mm) thickness(mm) material"
	x1 = strtod(args.at(0).c_str(), NULL);
	x2 = strtod(args.at(1).c_str(), NULL);
	y1 = strtod(args.at(2).c_str(), NULL);
	y2 = strtod(args.at(3).c_str(), NULL);
	thickness = strtod(args.at(4).c_str(), NULL);
	size = G4ThreeVector(std::max(x1, x2), std::max(y1, y2), thickness);
	return true;
}

void lightGuideLayer::construct(nDetDetector *obj){
	obj->applyLightGuide(x1, x2, y1, y2, thickness);
}

void lightGuideLayer::construct(nDetImplant *obj){
	obj->applyLightGuide(x1, x2, y1, y2, thickness);
}

std::string lightGuideLayer::syntaxStr() const {
	return std::string("addLightGuide <width1> <width2> <height1> <height2> <thickness> [material=G4_SILICON_DIOXIDE]");
}

///////////////////////////////////////////////////////////////////////////////
// class segLightGuideLayer
///////////////////////////////////////////////////////////////////////////////

segLightGuideLayer::segLightGuideLayer(const 
G4String &arg_) : nDetWorldObject(arg_,8), fXseg(0),fYseg(0),fspacing(0.),ftopWidth(0.),ftopThick(0.),fbotWidth(0.),fbotThick(0.),fzThick(0.),fSegMaterial("G4_Pyrex_Glass"){
	material = "G4_Pyrex_Glass";
}

segLightGuideLayer::~segLightGuideLayer(){
	G4cout<<"Deleting segLightGuideLayer"<<G4endl;
}

bool segLightGuideLayer::decodeArgs(){
	fXseg = int(strtod(args.at(0).c_str(), NULL));
	fYseg = int(strtod(args.at(1).c_str(), NULL));
	fzThick = strtod(args.at(2).c_str(), NULL);
	ftopWidth = strtod(args.at(3).c_str(), NULL);
	ftopThick = strtod(args.at(4).c_str(), NULL);
	fbotWidth = strtod(args.at(5).c_str(), NULL);
	fbotThick = strtod(args.at(6).c_str(), NULL);
	fspacing = strtod(args.at(7).c_str(), NULL);
	return true;
}

void segLightGuideLayer::construct(nDetDetector *obj){
	G4cout<<"Placing segmented light guide on nDetDetector which is not defined"<<G4endl;
}

void segLightGuideLayer::construct(nDetImplant *obj){
	obj->applyGreaseLayer(ftopWidth,ftopThick);
	obj->applySegmentedLightGuide(fXseg, fYseg, fspacing, ftopWidth, ftopThick, fbotWidth, fbotThick, fzThick);
	obj->setSegZThick(fzThick);
}

std::string segLightGuideLayer::syntaxStr() const{
	return std::string("<# x seg> <# y seg> <z thickness> <front x width> <front y width> <back x width> <back y width> <spacing> <material>");
}


///////////////////////////////////////////////////////////////////////////////
// class phoswichLayer
///////////////////////////////////////////////////////////////////////////////

phoswichLayer::phoswichLayer(const G4String &arg_) : nDetWorldObject(arg_,8), fXseg(0),fYseg(0),fThick(0.),fWidth(0.),fHeight(0.),fWrapping(0.), fSegMaterial("") 
{}

phoswichLayer::~phoswichLayer(){
	G4cout<<"Deleting phoswichLayer"<<G4endl;
}

bool phoswichLayer::decodeArgs(){
	fXseg = int(strtod(args.at(0).c_str(), NULL));
	fYseg = int(strtod(args.at(1).c_str(), NULL));
	fThick = strtod(args.at(2).c_str(), NULL);
	fWidth = strtod(args.at(3).c_str(), NULL);
	fHeight = strtod(args.at(4).c_str(), NULL);
	fWrapping = strtod(args.at(5).c_str(), NULL);
	fSegMaterial = args.at(6);
	fWrapMat = args.at(7);
	return true;
}


void phoswichLayer::construct(nDetDetector *obj){
	G4cout<<"Placing phoswich scintillator on nDetDetector which is not defined"<<G4endl;
}

void phoswichLayer::construct(nDetImplant *obj){
	obj->applyPhoswich(fXseg,fYseg,fThick,fWidth,fHeight,fWrapping,fSegMaterial,fWrapMat);
}

std::string phoswichLayer::syntaxStr() const{
	return std::string("<# x seg> <# y seg> <z thickness> <x width> <y height> <wrapping thickness> <material>");
}

///////////////////////////////////////////////////////////////////////////////
// class gdmlLightGuideLayer
///////////////////////////////////////////////////////////////////////////////

bool gdmlLightGuideLayer::decodeArgs(){
	// Expects a space-delimited string of the form:
	//  "filename rotX(deg) rotY(deg) rotZ(deg) material"
	filename = args.at(0);
	rotVector = G4ThreeVector(strtod(args.at(1).c_str(), NULL)*deg, strtod(args.at(2).c_str(), NULL)*deg, strtod(args.at(3).c_str(), NULL)*deg);
	material = args.at(4);

	// Load the model
	solid.read(filename, material, false);
	solid.setRotation(rotVector);
	
	// Get the size of the model
	size = G4ThreeVector(solid.getWidth(), solid.getThickness(), solid.getLength());
	
	return true;
}

void gdmlLightGuideLayer::construct(nDetDetector *obj){
	obj->loadLightGuide(&solid, rotVector);
}

void gdmlLightGuideLayer::construct(nDetImplant *obj){
	obj->loadLightGuide(&solid, rotVector);
}

std::string gdmlLightGuideLayer::syntaxStr() const {
	return std::string("loadLightGuide <filename> <rotX> <rotY> <rotZ> <matString>");
}
