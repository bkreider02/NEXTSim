////
// $Id: nDetConstruction.cc,v1.0 Sept., 2015 $
// Written by Dr. Xiaodong Zhang
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "nDetConstruction.hh"
#include "nDetConstructionMessenger.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4VPhysicalVolume.hh"
#include "G4AssemblyVolume.hh"
#include "G4VSolid.hh"

#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"

#include "G4GeometryManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4SubtractionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4UnionSolid.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"
#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"
#include "G4PhysicalConstants.hh"

static const G4double inch = 2.54*cm;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


nDetConstruction::nDetConstruction(const G4double &scale/*=1*/){
    expHall_logV = NULL;
    expHall_physV = NULL;
    
    assembly_logV = NULL;
    assembly_physV = NULL;
    
    assemblyBoxX = 15.11*mm; // reducing size to eliminate non-physical PTFE.  Noticed excessive scattering of flourine and carbon.
    assemblyBoxX = 16.66667*mm; //10x30 bars KS 5/23/16
    assemblyBoxX = 50*mm; // 30X30 scintillator bars, each has a dimension of 3X100X3 mm3
    assemblyBoxY = 51.2*mm; // length of scintillator: 100 mm + 0.1 mm thick grease + 1 mm thick SiPM window + 0.1 mm thick sensitive part of SiPM on both ends of scintillator
    assemblyBoxZ = 45.31*mm; // closing gaps.  Teflon tape measures ~.02 mm thick.
    assemblyBoxZ = 50*mm; // the gaps of scintillators are filled with reflectors, such as Teflon
    
    ej200_logV = NULL;
    ej200X = 1.5*mm;
    ej200Y = 50*mm;
    ej200Z = 1.5*mm;
    
    // 100 micron thick grease on both ends of EJ200 bar
    grease_logV = NULL;
    greaseX = 1.5*mm;
    greaseY = 0.05*mm;
    greaseZ = 1.5*mm; 
    
    // 1 mm thick SiPM window
    qwSiPM_logV = NULL;
    qwSiPMx = 1.5*mm;
    qwSiPMy = 0.5*mm;
    qwSiPMz = 1.5*mm; 
    
    // note that two sensing parts on both ends of scintillator are the same
    psSiPM_logV = NULL;
    psSiPMx = 1.5*mm;
    psSiPMy = 0.05*mm;
    psSiPMz = 1.5*mm; 

	fNumColumns = 1;
	fNumRows = 1;

	fNumColumnsPmt = -1;
	fNumRowsPmt = -1;

    // the world volume is 10 mm bigger than assembly volume in three dimensions

    G4cout<<"nDetConstruction::nDetConstruction()->"<<this<<G4endl;
    G4double margin = 500*mm;

    fDetectorMessenger=new nDetConstructionMessenger(this);

    fGeometry="rectangle";

    fCheckOverlaps = false;
    fTeflonThickness = 0.11*mm;
    fMylarThickness = 0.*mm;
    fDetectorLength = 3.94*inch;
    fDetectorWidth = 6*mm;
    fTrapezoidLength = 1*inch;
    fHexagonRadius = 5*cm;
    fDetectorThickness = 0.24*inch;
    SiPM_dimension=3*mm;    
    
    fLightYieldScale = scale;

    expHallX = assemblyBoxX + margin;
    expHallY = assemblyBoxY + margin;
    expHallZ = assemblyBoxZ + margin;

    //Build the materials
    DefineMaterials();

    //Make Sensitive Detectors

    //G4SDManager* SDMan = G4SDManager::GetSDMpointer();
  return;
} // end of construction function.

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

nDetConstruction::~nDetConstruction()
{} // end of deconstruction function.

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* nDetConstruction::Construct() {

    G4cout << "nDetConstruction::Construct()-->" << G4endl;

    if (expHall_physV) {
        G4GeometryManager::GetInstance()->OpenGeometry();
        G4PhysicalVolumeStore::GetInstance()->Clean();
        G4LogicalVolumeStore::GetInstance()->Clean();
        G4SolidStore::GetInstance()->Clean();
        //G4LogicalSkinSurface::CleanSurfaceTable();
        //G4LogicalBorderSurface::CleanSurfaceTable();
    }


    return ConstructDetector();
}

G4VPhysicalVolume* nDetConstruction::ConstructDetector(){

    // build experiment hall
  buildExpHall();

  // build assembly volume
  //buildAssembly();

  //build Disk volume
    if(fGeometry == "disk")
        buildDisk();
    if(fGeometry == "hexagon")
        buildDisk2();
    //builds the Ellipse
    if(fGeometry == "ellipse")
        buildEllipse2();
    if(fGeometry == "rectangle")
        buildRectangle();
    if(fGeometry == "array")
    buildArray();
  return expHall_physV;
}

void nDetConstruction::setSegmentedPmt(const short &col_, const short &row_, const double &width_, const double &height_){
  center[0].setSegmentedPmt(col_, row_, width_, height_);
  center[1].setSegmentedPmt(col_, row_, width_, height_);
  std::cout << " nDetConstruction: Setting segmented PMTs with WxH=(" << width_ << " x " << height_ << ") and " << col_ << " columns and " << row_ << " rows.\n";
}

bool nDetConstruction::setPmtSpectralResponse(const char *fname){
  bool retval = (center[0].loadSpectralResponse(fname) && center[1].loadSpectralResponse(fname));
  if(retval)
    std::cout << " nDetConstruction: Successfully loaded PMT spectral response function\n";
  else
    std::cout << " nDetConstruction: ERROR! Failed to load PMT spectral response from \"" << fname << "\"!\n";
  return retval;
}

bool nDetConstruction::setPmtGainMatrix(const char *fname){
  bool retval = (center[0].loadGainMatrix(fname) && center[1].loadGainMatrix(fname));
  if(retval)
    std::cout << " nDetConstruction: Successfully loaded PMT anode gain matrix\n";
  else
    std::cout << " nDetConstruction: ERROR! Failed to load PMT anode gain matrix from \"" << fname << "\"!\n";
  return retval;
}

void nDetConstruction::buildExpHall()
{
  // **************** expHalll **************
  //      ------------- Volumes --------------
  
  G4Box* expHall_solidV = new G4Box("expHall_solidV",expHallX,expHallY,expHallZ);

  expHall_logV  = new G4LogicalVolume(expHall_solidV, fAir, "expHall_logV",0,0,0);
 // expHall_logV  = new G4LogicalVolume(expHall_solidV, Vacuum, "expHall_logV",0,0,0);
    expHall_logV->SetVisAttributes(G4VisAttributes::Invisible);

  expHall_physV = new G4PVPlacement(0,G4ThreeVector(0., 0., 0.),expHall_logV,
                                   "expHall_physV",0,false,0);
 
  return;
} // end of buildExpHall function


// build scintillator array that consists of 10X30 EJ200 scintillator bars
// Each bar is read out by two SiPMs at the two ends
// Teflon tape is used to wrap the bars and take care diffusing light back into the scintillator
void nDetConstruction::buildAssembly()
{ 
  //************* Wrapping materials (such as Teflon) and EJ200  *************

  //-----------------------------------------------------------
  //     The assembly volume is defined as a wrapping material
  //-----------------------------------------------------------

   G4Box* assemblyBox = new G4Box("assembly", assemblyBoxX, assemblyBoxY, assemblyBoxZ);
  assembly_logV = new G4LogicalVolume(assemblyBox, fTeflon, "assembly");
  assembly_physV = new G4PVPlacement(0,
				     G4ThreeVector(assemblyPx, assemblyPy, assemblyPz),
 				     assembly_logV,
				     "assembly_physV",
				     expHall_logV,
				     false, 0);

  fWrapSkinSurface = new G4LogicalSkinSurface("teflonSurface", assembly_logV, fTeflonOpticalSurface);


  //----------------------------------------------------------------
  //    scintillator bars
  //----------------------------------------------------------------

  G4Box* ej200_solidV = new G4Box("ej200_solidV", ej200X, ej200Y, ej200Z);

  ej200_logV = new G4LogicalVolume(ej200_solidV, fEJ200, "ej200_logV", 0, 0, 0);

  G4VisAttributes* ej200_VisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0));//blue
  ej200_logV->SetVisAttributes(ej200_VisAtt);

  // create assembly of scintillator bars
  //
  G4AssemblyVolume* ej200Assembly = new G4AssemblyVolume();
  // Rotation and translation of the Fiber logical volume
  G4RotationMatrix rotationMatrix_ej200;
  G4ThreeVector position_ej200(0., 0., 0.);
  G4Transform3D transform3D_ej200;

//  for(G4int x_index = 0; x_index<30; x_index++){
  for(G4int x_index = 0; x_index<10; x_index++){ //reduced array size to 10x30 KS 5/23/16
    for(G4int z_index = 0; z_index<30; z_index++){
      //position_ej200.setX(assemblyBoxX - assemblyBoxX/15*(1+x_index*2)/2);
      position_ej200.setX(assemblyBoxX - assemblyBoxX/5*(1+x_index*2)/2); //reduced KS 5/23/16
      position_ej200.setZ(assemblyBoxZ - assemblyBoxZ/15*(1+z_index*2)/2);
      position_ej200.setY(0.0);

      transform3D_ej200 = G4Transform3D(rotationMatrix_ej200, position_ej200);
      ej200Assembly->AddPlacedVolume(ej200_logV, transform3D_ej200);
    }
  }


   /*   position_ej200.setX(0.0); //reduced KS 5/23/16
      position_ej200.setZ(0.0);
      position_ej200.setY(0.0);
   */
      transform3D_ej200 = G4Transform3D(rotationMatrix_ej200, position_ej200);
      ej200Assembly->AddPlacedVolume(ej200_logV, transform3D_ej200);

  //  place the assembly
  //
  position_ej200.setX(0.0); position_ej200.setY(0.0); position_ej200.setZ(0.0);
  transform3D_ej200 = G4Transform3D(rotationMatrix_ej200, position_ej200);
  ej200Assembly->MakeImprint(assembly_logV, transform3D_ej200);

  // grease Sanit-Gobain BC-630 silicone grease
  // optical grease

  // grease box ...
  G4Box* grease_solidV = new G4Box("grease", greaseX, greaseY, greaseZ);
  grease_logV = new G4LogicalVolume(grease_solidV, fGrease, "grease_logV", 0, 0, 0);
  G4VisAttributes* grease_VisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0));//red
  grease_logV->SetVisAttributes(grease_VisAtt);

  // create assembly of grease
  //
  G4AssemblyVolume* greaseAssembly = new G4AssemblyVolume();
  // Rotation and translation of the grease logical volume
  G4RotationMatrix rotationMatrix_grease;
  G4ThreeVector position_grease(0., 0., 0.);
  G4Transform3D transform3D_grease;

  //for(G4int x_index = 0; x_index<30; x_index++){
  for(G4int x_index = 0; x_index<10; x_index++){
    for(G4int z_index = 0; z_index<30; z_index++){
      //position_grease.setX(assemblyBoxX - assemblyBoxX/15*(1+x_index*2)/2);
      position_grease.setX(assemblyBoxX - assemblyBoxX/5*(1+x_index*2)/2);
      position_grease.setZ(assemblyBoxZ - assemblyBoxZ/15*(1+z_index*2)/2);
      position_grease.setY(ej200Y+greaseY);

      transform3D_grease = G4Transform3D(rotationMatrix_grease, position_grease);
      greaseAssembly->AddPlacedVolume(grease_logV, transform3D_grease);

      position_grease.setY(-1*(ej200Y+greaseY));

      transform3D_grease = G4Transform3D(rotationMatrix_grease, position_grease);
      greaseAssembly->AddPlacedVolume(grease_logV, transform3D_grease);
    }
  }

  //  place the assembly
  //
  position_grease.setX(0.0); position_grease.setY(0.0); position_grease.setZ(0.0);
  transform3D_grease = G4Transform3D(rotationMatrix_grease, position_grease);
  greaseAssembly->MakeImprint(assembly_logV, transform3D_grease);


  G4Box* qwSiPM_solidV = new G4Box("qwSiPM_solidV", qwSiPMx, qwSiPMy, qwSiPMz);

  qwSiPM_logV = new G4LogicalVolume(qwSiPM_solidV, fSiO2, "qwSiPM_logV", 0, 0, 0);
  G4VisAttributes* qwSiPM_VisAtt= new G4VisAttributes(G4Colour(0.0,1.0,0.0));//green
  qwSiPM_logV->SetVisAttributes(qwSiPM_VisAtt);

  // create assembly of SiPM window
  //
  G4AssemblyVolume* SiPMwindowAssembly = new G4AssemblyVolume();
  // Rotation and translation of the grease logical volume
  G4RotationMatrix rotationMatrix_SiPMwindow;
  G4ThreeVector position_SiPMwindow(0., 0., 0.);
  G4Transform3D transform3D_SiPMwindow;

//  for(G4int x_index = 0; x_index<30; x_index++){
  for(G4int x_index = 0; x_index<10; x_index++){
    for(G4int z_index = 0; z_index<30; z_index++){
      //position_SiPMwindow.setX(assemblyBoxX - assemblyBoxX/15*(1+x_index*2)/2);
      position_SiPMwindow.setX(assemblyBoxX - assemblyBoxX/5*(1+x_index*2)/2);
      position_SiPMwindow.setZ(assemblyBoxZ - assemblyBoxZ/15*(1+z_index*2)/2);
      position_SiPMwindow.setY(ej200Y+2*greaseY+qwSiPMy);

      transform3D_SiPMwindow = G4Transform3D(rotationMatrix_SiPMwindow, position_SiPMwindow);
      SiPMwindowAssembly->AddPlacedVolume(qwSiPM_logV, transform3D_SiPMwindow);

      position_SiPMwindow.setY(-1*(ej200Y+2*greaseY+qwSiPMy));

      transform3D_SiPMwindow = G4Transform3D(rotationMatrix_SiPMwindow, position_SiPMwindow);
      SiPMwindowAssembly->AddPlacedVolume(qwSiPM_logV, transform3D_SiPMwindow);
    }
  }

  //  place the assembly
  //
  position_SiPMwindow.setX(0.0); position_SiPMwindow.setY(0.0); position_SiPMwindow.setZ(0.0);
  transform3D_SiPMwindow = G4Transform3D(rotationMatrix_SiPMwindow, position_SiPMwindow);
  SiPMwindowAssembly->MakeImprint(assembly_logV, transform3D_SiPMwindow);

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // the portion of SiPM just following the front quartz window is used for sensing the photons
  // the material of this portion or photon sensing part is made of silicon.
  //


  G4Box* psSiPM_solidV = new G4Box("psSiPM_solidV", psSiPMx, psSiPMy, psSiPMz);


  psSiPM_logV = new G4LogicalVolume(psSiPM_solidV, fSil, "psSiPM_logV", 0, 0, 0);
  G4VisAttributes* psSiPM_VisAtt= new G4VisAttributes(G4Colour(0.0,1.0,1.0));//
  psSiPM_logV->SetVisAttributes(psSiPM_VisAtt);


    fSiPMSkinSurface = new G4LogicalSkinSurface("SiPMSurface", psSiPM_logV,fSiliconPMOpticalSurface);

  // create assembly of SiPM sensitive part
  //
  G4AssemblyVolume* psSiPMAssembly = new G4AssemblyVolume();
  // Rotation and translation of the grease logical volume
  G4RotationMatrix rotationMatrix_psSiPM;
  G4ThreeVector position_psSiPM(0., 0., 0.);
  G4Transform3D transform3D_psSiPM;

//  for(G4int x_index = 0; x_index<30; x_index++){
  for(G4int x_index = 0; x_index<10; x_index++){
    for(G4int z_index = 0; z_index<30; z_index++){
//      position_psSiPM.setX(assemblyBoxX - assemblyBoxX/15*(1+x_index*2)/2);
      position_psSiPM.setX(assemblyBoxX - assemblyBoxX/5*(1+x_index*2)/2);
      position_psSiPM.setZ(assemblyBoxZ - assemblyBoxZ/15*(1+z_index*2)/2);
      position_psSiPM.setY(ej200Y+2*greaseY+2*qwSiPMy+psSiPMy);

      transform3D_psSiPM = G4Transform3D(rotationMatrix_psSiPM, position_psSiPM);
      psSiPMAssembly->AddPlacedVolume(psSiPM_logV, transform3D_psSiPM);

      position_psSiPM.setY(-1*(ej200Y+2*greaseY+2*qwSiPMy+psSiPMy));

      transform3D_psSiPM = G4Transform3D(rotationMatrix_psSiPM, position_psSiPM);
      psSiPMAssembly->AddPlacedVolume(psSiPM_logV, transform3D_psSiPM);
    }
  }

  //  place the assembly
  //
  position_psSiPM.setX(0.0); position_psSiPM.setY(0.0); position_psSiPM.setZ(0.0);
  transform3D_psSiPM = G4Transform3D(rotationMatrix_psSiPM, position_psSiPM);


  psSiPMAssembly->MakeImprint(assembly_logV, transform3D_psSiPM);

  return;
} // end of function. //





void nDetConstruction::buildDisk()
{ 
  //************* Wrapping materials (such as Teflon) and EJ200  *************

  //-----------------------------------------------------------
  //     The assembly volume is defined as a wrapping material
  //-----------------------------------------------------------





    greaseX = SiPM_dimension;
    greaseZ = SiPM_dimension;
    qwSiPMx = SiPM_dimension;
    qwSiPMz = SiPM_dimension;
    psSiPMx = SiPM_dimension;
    psSiPMz = SiPM_dimension;

    G4double hexagonSize=50*mm;

    G4double teflonThickness=0.11*mm;

    G4double minDiskRadius=0*mm;
    G4double maxDiskRadius=hexagonSize+teflonThickness;

    G4double diskLength=SiPM_dimension*mm+teflonThickness;

    G4double startAngle=0*deg;
    G4double endAngle=360*deg;

    G4int nSides=6;
    G4int nPlanes=2;

    const G4double zplanes[]={-diskLength,diskLength};
    const G4double maxRadius[]={maxDiskRadius,maxDiskRadius};
    const G4double minRadius[]={minDiskRadius,minDiskRadius};

    G4VSolid *assemblyDisk=NULL;

    if(fGeometry=="disk")
    assemblyDisk = (G4Tubs*)new G4Tubs("Disk",minDiskRadius,maxDiskRadius,diskLength,startAngle,endAngle);

    else
    assemblyDisk= (G4Polyhedra*) new G4Polyhedra("Hexagon",startAngle,endAngle,nSides,nPlanes,zplanes,minRadius,maxRadius);

    /*
    G4double BoxX=1.52*mm;
    G4double BoxY=1.52*mm;
    G4double BoxZ=1.65*mm;
    */
    G4double BoxX=SiPM_dimension+2*teflonThickness;
    G4double BoxY=SiPM_dimension+2*teflonThickness;
    G4double BoxZ=SiPM_dimension+2*(teflonThickness+greaseY+qwSiPMy+psSiPMy);

    G4Box *Hole=new G4Box("Hole",BoxX,BoxY,BoxZ);

    G4SubtractionSolid *theSolid=new G4SubtractionSolid("Wrapping",assemblyDisk,Hole);


    assembly_logV = new G4LogicalVolume(theSolid, fTeflon, "Wrap");


    G4RotationMatrix *rotationMatrix_Disk=new G4RotationMatrix;
    rotationMatrix_Disk->rotateY(90*deg);

    assembly_physV = new G4PVPlacement(rotationMatrix_Disk,
                                       G4ThreeVector(assemblyPx, assemblyPy, assemblyPz),
                                       assembly_logV,
                                       "assembly_physV",
                                       expHall_logV,
                                       false, 0);

    //----------------------------------------------------------------
    //    scintillator bars
    //----------------------------------------------------------------


  G4LogicalSkinSurface* assemblySkin = new G4LogicalSkinSurface("teflonSUrface", assembly_logV, fTeflonOpticalSurface);


     minDiskRadius=0*mm;
     maxDiskRadius=hexagonSize;
    
     diskLength=SiPM_dimension;
    
     startAngle=0*deg;
     endAngle=360*deg;

    const G4double zplanes2[]={-diskLength,diskLength};
    const G4double maxRadius2[]={maxDiskRadius,maxDiskRadius};
    const G4double minRadius2[]={minDiskRadius,minDiskRadius};

    G4VSolid *scintillatorDisk=NULL;
    if(fGeometry=="disk")
        scintillatorDisk = (G4Tubs*)new G4Tubs("SciDisk",minDiskRadius,maxDiskRadius,diskLength,startAngle,endAngle);
    else
        scintillatorDisk= (G4Polyhedra*) new G4Polyhedra("SciHexagon",startAngle,endAngle,nSides,nPlanes,zplanes2,minRadius2,maxRadius2);


/*
    BoxX=2.70*mm;
    BoxY=2.70*mm;
    BoxZ=1.62*mm;
*/
    BoxX=BoxZ;
    BoxY=BoxZ;
    BoxZ=SiPM_dimension+2*teflonThickness;

    G4Box *Hole2=new G4Box("Hole",BoxX,BoxY,BoxZ);
    
    G4SubtractionSolid *theScint=new G4SubtractionSolid("Detector",scintillatorDisk,Hole2);
  
    
  ej200_logV = new G4LogicalVolume(theScint, fEJ200, "ej200_logV", 0, 0, 0);

  G4VisAttributes* ej200_VisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0));//blue
  //ej200_VisAtt->SetForceSolid(true);
  ej200_logV->SetVisAttributes(ej200_VisAtt);

  //ej200_logV->SetSensitiveDetector(fScintSD);

  // create assembly of scintillator
  //
  
  G4RotationMatrix *rotationMatrix_Scint=new G4RotationMatrix;

    
  assembly_physV = new G4PVPlacement(rotationMatrix_Scint,
                                       G4ThreeVector(assemblyPx, assemblyPy, assemblyPz),
                                       ej200_logV,
                                       "Scint_physV",
                                       assembly_logV,
                                       false, 0);
  
  // grease Sanit-Gobain BC-630 silicone grease
  // optical grease


  // grease box ...
  G4Box* grease_solidV = new G4Box("grease", greaseX, greaseY, greaseZ);
  grease_logV = new G4LogicalVolume(grease_solidV, fGrease, "grease_logV", 0, 0, 0);
  G4VisAttributes* grease_VisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0));//red
  //grease_VisAtt->SetForceSolid(true);
  grease_logV->SetVisAttributes(grease_VisAtt);

  // create assembly of grease
  //
  G4AssemblyVolume* greaseAssembly = new G4AssemblyVolume();
  // Rotation and translation of the grease logical volume
  G4RotationMatrix rotationMatrix_grease;
  G4ThreeVector position_grease(0., 0., 0.);
  G4Transform3D transform3D_grease;

  for(G4int i = 0; i<2; i++){
      //position_grease.setX(assemblyBoxX - assemblyBoxX/15*(1+x_index*2)/2);
      if(i==0){
          position_grease.setX(0);
          position_grease.setZ(0);
          position_grease.setY(BoxY-greaseY);

          transform3D_grease = G4Transform3D(rotationMatrix_grease, position_grease);
          greaseAssembly->AddPlacedVolume(grease_logV, transform3D_grease);

          position_grease.setY(-1*(BoxY-greaseY));

          transform3D_grease = G4Transform3D(rotationMatrix_grease, position_grease);
          greaseAssembly->AddPlacedVolume(grease_logV, transform3D_grease);
      }
      else{
          rotationMatrix_grease.rotateZ(90*deg);
          position_grease.setY(0);
          position_grease.setZ(0);
          position_grease.setX(BoxY-greaseY);
          
          transform3D_grease = G4Transform3D(rotationMatrix_grease, position_grease);
          greaseAssembly->AddPlacedVolume(grease_logV, transform3D_grease);
          
          position_grease.setX(-1*(BoxY-greaseY));
          
          transform3D_grease = G4Transform3D(rotationMatrix_grease, position_grease);
          greaseAssembly->AddPlacedVolume(grease_logV, transform3D_grease);
          
      }
  }

  //  place the assembly
  //
  position_grease.setX(0.0); position_grease.setY(0.0); position_grease.setZ(0.0);
  transform3D_grease = G4Transform3D(rotationMatrix_grease, position_grease);
  greaseAssembly->MakeImprint(assembly_logV, transform3D_grease);


  G4Box* qwSiPM_solidV = new G4Box("qwSiPM_solidV", qwSiPMx, qwSiPMy, qwSiPMz);

  qwSiPM_logV = new G4LogicalVolume(qwSiPM_solidV, fSiO2, "qwSiPM_logV", 0, 0, 0);
  G4VisAttributes* qwSiPM_VisAtt= new G4VisAttributes(G4Colour(0.0,1.0,0.0));//green
  //qwSiPM_VisAtt->SetForceSolid(true);
  qwSiPM_logV->SetVisAttributes(qwSiPM_VisAtt);

  // create assembly of SiPM window
  //
  G4AssemblyVolume* SiPMwindowAssembly = new G4AssemblyVolume();
  // Rotation and translation of the grease logical volume
  G4RotationMatrix rotationMatrix_SiPMwindow;
  G4ThreeVector position_SiPMwindow(0., 0., 0.);
  G4Transform3D transform3D_SiPMwindow;

  for(G4int i = 0; i<2; i++){
      if(i==0){
      position_SiPMwindow.setX(0);
      position_SiPMwindow.setZ(0);
      position_SiPMwindow.setY(BoxY-2*greaseY-qwSiPMy);

      transform3D_SiPMwindow = G4Transform3D(rotationMatrix_SiPMwindow, position_SiPMwindow);
      SiPMwindowAssembly->AddPlacedVolume(qwSiPM_logV, transform3D_SiPMwindow);

      position_SiPMwindow.setY(-1*(BoxY-2*greaseY-qwSiPMy));

      transform3D_SiPMwindow = G4Transform3D(rotationMatrix_SiPMwindow, position_SiPMwindow);
      SiPMwindowAssembly->AddPlacedVolume(qwSiPM_logV, transform3D_SiPMwindow);
      }
      else{
          rotationMatrix_SiPMwindow.rotateZ(90*deg);
          position_SiPMwindow.setY(0);
          position_SiPMwindow.setZ(0);
          position_SiPMwindow.setX(BoxY-2*greaseY-qwSiPMy);
          
          transform3D_SiPMwindow = G4Transform3D(rotationMatrix_SiPMwindow, position_SiPMwindow);
          SiPMwindowAssembly->AddPlacedVolume(qwSiPM_logV, transform3D_SiPMwindow);
          
          position_SiPMwindow.setX(-1*(BoxY-2*greaseY-qwSiPMy));
          
          transform3D_SiPMwindow = G4Transform3D(rotationMatrix_SiPMwindow, position_SiPMwindow);
          SiPMwindowAssembly->AddPlacedVolume(qwSiPM_logV, transform3D_SiPMwindow);
      
      }
  }
  

  //  place the assembly
  //
  position_SiPMwindow.setX(0.0); position_SiPMwindow.setY(0.0); position_SiPMwindow.setZ(0.0);
  transform3D_SiPMwindow = G4Transform3D(rotationMatrix_SiPMwindow, position_SiPMwindow);
  SiPMwindowAssembly->MakeImprint(assembly_logV, transform3D_SiPMwindow);

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // the portion of SiPM just following the front quartz window is used for sensing the photons
  // the material of this portion or photon sensing part is made of silicon.
  // 

  // define the optical surface for SiPM
  G4OpticalSurface* SiPMSurface = new G4OpticalSurface("SiPMSurface");
  SiPMSurface->SetType(dielectric_metal);
  SiPMSurface->SetFinish(polished);
  SiPMSurface->SetModel(glisur);
  SiPMSurface->SetMaterialPropertiesTable(fSilMPT);

  G4Box* psSiPM_solidV = new G4Box("psSiPM_solidV", psSiPMx, psSiPMy, psSiPMz);



    psSiPM_logV = new G4LogicalVolume(psSiPM_solidV,fSil, "psSiPM_logV", 0, 0, 0);
  G4VisAttributes* psSiPM_VisAtt= new G4VisAttributes(G4Colour(0.0,1.0,1.0));//
  psSiPM_VisAtt->SetForceSolid(true);
  psSiPM_logV->SetVisAttributes(psSiPM_VisAtt);

  //Add Sensitive Detector to SiPM
  G4LogicalSkinSurface* psSiPMSurface = new G4LogicalSkinSurface("SiPMSurface", psSiPM_logV, SiPMSurface);

  // create assembly of SiPM sensitive part
  //
  G4AssemblyVolume* psSiPMAssembly = new G4AssemblyVolume();
  // Rotation and translation of the grease logical volume
  G4RotationMatrix rotationMatrix_psSiPM;
  G4ThreeVector position_psSiPM(0., 0., 0.);
  G4Transform3D transform3D_psSiPM;

  for(G4int i = 0; i<2; i++){
      if(i==0){
      position_psSiPM.setX(0);
      position_psSiPM.setZ(0);
      position_psSiPM.setY(BoxY-2*greaseY-2*qwSiPMy-psSiPMy);

      transform3D_psSiPM = G4Transform3D(rotationMatrix_psSiPM, position_psSiPM);
      psSiPMAssembly->AddPlacedVolume(psSiPM_logV, transform3D_psSiPM);

      position_psSiPM.setY(-1*(BoxY-2*greaseY-2*qwSiPMy-psSiPMy));

      transform3D_psSiPM = G4Transform3D(rotationMatrix_psSiPM, position_psSiPM);
      psSiPMAssembly->AddPlacedVolume(psSiPM_logV, transform3D_psSiPM);
      }
      else{
          rotationMatrix_psSiPM.rotateZ(90*deg);
          position_psSiPM.setY(0);
          position_psSiPM.setZ(0);
          position_psSiPM.setX(BoxY-2*greaseY-2*qwSiPMy-psSiPMy);
          
          transform3D_psSiPM = G4Transform3D(rotationMatrix_psSiPM, position_psSiPM);
          psSiPMAssembly->AddPlacedVolume(psSiPM_logV, transform3D_psSiPM);
          
          position_psSiPM.setX(-1*(BoxY-2*greaseY-2*qwSiPMy-psSiPMy));
          
          transform3D_psSiPM = G4Transform3D(rotationMatrix_psSiPM, position_psSiPM);
          psSiPMAssembly->AddPlacedVolume(psSiPM_logV, transform3D_psSiPM);
      
      }
  }

  //  place the assembly
  //
  position_psSiPM.setX(0.0); position_psSiPM.setY(0.0); position_psSiPM.setZ(0.0);
  transform3D_psSiPM = G4Transform3D(rotationMatrix_psSiPM, position_psSiPM);
    psSiPMAssembly->MakeImprint(assembly_logV, transform3D_psSiPM);

  return;
} // end of function. //


void nDetConstruction::buildEllipse() {


    G4double teflonThickness=0.11*mm;

    G4double mylarThickness=0.0125*mm; //25 um mylar

    //First we build the box
    G4double xdimension=1.18/2*inch+teflonThickness;
    G4double ydimension=0.24/2*inch+teflonThickness;
    G4double zdimension=1.94/2*inch+teflonThickness;

    G4Box *theBox=new G4Box("theBox",xdimension,ydimension,zdimension);

    //Now the trapezoid

    G4double dy1=ydimension;
    G4double dy2=ydimension;
    G4double dx1=xdimension;
    G4double dx2=ydimension;
    G4double dz=0.5*inch;


    //Ellipse shaped
    G4Trd *theTrapezoid1=new G4Trd("theTrapezoid1", dx1, dx2, dy1, dy2, dz);

    G4Trd *theTrapezoid2=new G4Trd("theTrapezoid2", dx1, dx2, dy1, dy2, dz);

    G4Box *wrappinBox=new G4Box("theBox1",ydimension,ydimension,(greaseY+qwSiPMy+psSiPMy));

    G4Box *wrappinBox2=new G4Box("theBox2",ydimension,ydimension,(greaseY+qwSiPMy+psSiPMy));

    G4ThreeVector translation1(0,0,(1.94/2+1./2)*inch);

    G4UnionSolid *theWrapping0=new G4UnionSolid("theWrapping0",theBox,theTrapezoid1,0,translation1);

    G4ThreeVector translation11(0,0,(1.94/2+1.)*inch+1*(greaseY+qwSiPMy+psSiPMy));

    G4UnionSolid *theWrapping1=new G4UnionSolid("theWrapping1",theWrapping0,wrappinBox,0,translation11);

    G4ThreeVector translation2(0,0,-(1.94/2+1./2)*inch);
    G4RotationMatrix *rot=new G4RotationMatrix();

    rot->rotateX(180*deg);

    G4UnionSolid *theWrapping2=new G4UnionSolid("theWrapping",theWrapping1,theTrapezoid2,rot,translation2);

    G4ThreeVector translation22(0,0,-1*((1.94/2+1.)*inch+1*(greaseY+qwSiPMy+psSiPMy)));



    //theBox
    xdimension=1.18/2*inch+teflonThickness;
    ydimension=0.24/2*inch+teflonThickness;
    zdimension=(1.94/2+1)*inch+teflonThickness+2*(greaseY+qwSiPMy+psSiPMy);

    G4VSolid *theWrapping=NULL;

    if(fGeometry == "ellipse")
        theWrapping=new G4UnionSolid("theWrapping",theWrapping2,wrappinBox2,rot,translation22);
    else
        theWrapping=new G4Box("theWrapping",xdimension,ydimension,zdimension);



    assembly_logV=new G4LogicalVolume(theWrapping,fTeflon,"Wrapping_log");

    G4VPhysicalVolume *Wrapping_physVol=new G4PVPlacement(0,G4ThreeVector(0,0,0),assembly_logV,"Wrap",expHall_logV,0,0,true);

    G4LogicalSkinSurface *thePlasticSkin=new G4LogicalSkinSurface("Wrapping",assembly_logV,fTeflonOpticalSurface);




    //Building the Scintillator

    //First we build the box
    xdimension=1.18/2*inch;
    ydimension=0.24/2*inch;
    zdimension=1.94/2*inch;

    G4Box *theBoxScint=new G4Box("theBoxScint",xdimension,ydimension,zdimension);

    //Now the trapezoid

    dy1=ydimension;
    dy2=ydimension;
    dx1=xdimension;
    dx2=ydimension;
    dz=0.5*inch;

    G4Trd *theTrapezoidScint1=new G4Trd("theTrapezoidScint1", dx1, dx2, dy1, dy2, dz);

    G4Trd *theTrapezoidScint2=new G4Trd("theTrapezoidScint2", dx1, dx2, dy1, dy2, dz);

    G4UnionSolid *thePlastic0=new G4UnionSolid("thePlastic0",theBoxScint,theTrapezoidScint1,0,translation1);

    //translation2(0,0,-(1.94/2+1./2)*inches);

    //G4RotationMatrix *rot2=new G4RotationMatrix();

    //rot->rotateX(180*deg);



    //theBox
    xdimension=1.18/2*inch;
    ydimension=0.24/2*inch;
    zdimension=(1.94/2+1)*inch;

    G4VSolid *thePlastic=NULL;
    if(fGeometry=="ellipse")
       thePlastic=new G4UnionSolid("thePlastic",thePlastic0,theTrapezoidScint2,rot,translation2);
    else
       thePlastic=new G4Box("thePlastic",xdimension,ydimension,zdimension);


    ej200_logV=new G4LogicalVolume(thePlastic,fEJ200,"plastic_log");
    //thePlastic_log->SetSensitiveDetector(fScintSD);
    G4VisAttributes* ej200_VisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0));//blue
    ej200_logV->SetVisAttributes(ej200_VisAtt);
    G4VPhysicalVolume *physVol=new G4PVPlacement(0,G4ThreeVector(0,0,0),ej200_logV,"Scint",assembly_logV,0,0,true);

    buildSiPMs();
    return;

}

//****************************************** Material Definitions  *********************************************//

void nDetConstruction::DefineMaterials() {
    G4cout<<"nDetConstruction::DefineMaterials()"<<G4endl;
    //Elements
    G4double z;
    G4double a;

    fH=new G4Element("Hydrogen", "H", z=1., a=1.00794*g/mole);
    fC=new G4Element("Carbon", "C", z=6., a=12.01*g/mole);
    fO=new G4Element("Oxygen", "O", z=8., a=16.00*g/mole);
    fF=new G4Element("Fluorine", "F", z=9., a=18.9984*g/mole);
    fSi=new G4Element("Aluminium","Si",z=14.,a=28.09*g/mole);
    fAl=new G4Element("Aluminium","Al",z=13.,a=26.9815*g/mole);

    //Materials & Properties
    G4NistManager* manNist = G4NistManager::Instance();

    fAir=manNist->FindOrBuildMaterial("G4_AIR");

    const G4int nEntries_Air = 2;
    G4double photonEnergy_Air[nEntries_Air] = { 1.55*eV,  4.13*eV }; // from 800 nm to 300 nm
    G4double refractiveIndex_Air[nEntries_Air] = { 1.00, 1.00 };
    G4double Absorption_Air[nEntries_Air] = { 1700 * m ,1700 * m  };
    // note that the aborption data are not correct. If they are employed, the purpose is to kill photons quickly
    // to speed up the program...
    //G4double absorptionLength_Air[nEntries] = {90.2*mm, 90.2*mm};
    fAirMPT = new G4MaterialPropertiesTable();
    fAirMPT->AddProperty("RINDEX", photonEnergy_Air, refractiveIndex_Air, nEntries_Air);
    fAirMPT->AddProperty("ABSLENGTH", photonEnergy_Air, Absorption_Air, nEntries_Air);
    fAir->SetMaterialPropertiesTable(fAirMPT);

    G4double density;
    int natoms;
    int ncomponents;

	/////////////////////////////////////////////////////////////////
	// Teflon (C2F4)n
	/////////////////////////////////////////////////////////////////

    fTeflon= new G4Material("Teflon", density=2.2*g/cm3,2);
    fTeflon->AddElement(fC,natoms=2);
    fTeflon->AddElement(fF,natoms=4);

    fTeflonMPT = new G4MaterialPropertiesTable();

    const G4int nEntries_Teflon = 2;

    G4double photonEnergy_teflon[nEntries_Teflon] = {2.*eV, 3.47*eV};
    G4double refl_teflon[nEntries_Teflon] = {0.99, 0.99};
    //G4double refl_teflon[nEntries_Teflon] = {0.5, 0.5};
    G4double effi_teflon[nEntries_Teflon] = {0, 0};
    G4double Absorption_Teflon[nEntries_Teflon] =  { 0.333 * cm ,0.333 * cm  };
    G4double refractiveIndex_teflon[nEntries_Teflon] = {1.315, 1.315};
    fTeflonMPT->AddProperty("REFLECTIVITY", photonEnergy_teflon, refl_teflon, nEntries_Teflon);
    fTeflonMPT->AddProperty("EFFICIENCY", photonEnergy_teflon, effi_teflon, nEntries_Teflon);
    fTeflonMPT->AddProperty("RINDEX", photonEnergy_teflon,refractiveIndex_teflon,nEntries_Teflon);
    fTeflonMPT->AddProperty("ABSLENGTH", photonEnergy_teflon,Absorption_Teflon,nEntries_Teflon);

    fTeflonMPT->DumpTable();

    fTeflon->SetMaterialPropertiesTable(fTeflonMPT);

	/////////////////////////////////////////////////////////////////
	// EJ200 (C9H10)n
	/////////////////////////////////////////////////////////////////

    fEJ200 = new G4Material("EJ200", 1.023*g/cm3, 2);
    fEJ200->AddElement(fH, 10);
    fEJ200->AddElement(fC, 9);

    /*const G4int nEntries_EJ200 = 5;
    G4double PhotonEnergy[nEntries_EJ200] = { 2.484*eV, 2.615*eV, 2.760*eV, 2.922*eV, 3.105*eV };
    G4double RefractiveIndex_EJ200[nEntries_EJ200] = { 1.58, 1.58, 1.58, 1.58, 1.58 };
    G4double Absorption_EJ200[nEntries_EJ200] =  { 400*cm, 400*cm, 400*cm, 400*cm, 400*cm };
    G4double ScintilFast_EJ200[nEntries_EJ200] = { 0.05, 0.25, 0.55, 1.00, 0.0};*/

    G4double PhotonEnergy[33] = {2.490*eV, 2.517*eV, 2.552*eV, 2.585*eV, 2.613*eV, 2.635*eV, 2.656*eV, 2.686*eV, 2.720*eV, 2.749*eV,
                                 2.772*eV, 2.791*eV, 2.809*eV, 2.826*eV, 2.842*eV, 2.861*eV, 2.884*eV, 2.919*eV, 2.946*eV, 2.954*eV, 
                                 2.961*eV, 2.967*eV, 2.974*eV, 2.981*eV, 2.987*eV, 2.994*eV, 3.001*eV, 3.009*eV, 3.018*eV, 3.029*eV, 
                                 3.041*eV, 3.056*eV, 3.083*eV};

    G4double ScintilFast_EJ200[33] = {0.060, 0.084, 0.122, 0.175, 0.234, 0.294, 0.356, 0.416, 0.473, 0.533,
                                      0.594, 0.657, 0.720, 0.784, 0.846, 0.903, 0.962, 1.000, 0.917, 0.857, 
                                      0.798, 0.732, 0.669, 0.604, 0.542, 0.480, 0.422, 0.359, 0.297, 0.237, 
                                      0.170, 0.105, 0.028};

	G4double PhotonEnergy2[2] = {0.060, 0.028};
	G4double RefIndex_EJ200[2] = {1.58, 1.58};
	G4double Absorption_EJ200[2] = {400*cm, 400*cm};

    fEJ200MPT = new G4MaterialPropertiesTable();
    fEJ200MPT->AddProperty("RINDEX", PhotonEnergy2, RefIndex_EJ200, 2);
    fEJ200MPT->AddProperty("ABSLENGTH", PhotonEnergy2, Absorption_EJ200, 2);
    fEJ200MPT->AddProperty("FASTCOMPONENT", PhotonEnergy, ScintilFast_EJ200, 33);

    fEJ200MPT->AddConstProperty("SCINTILLATIONYIELD", 0.64*17400/MeV); // 64% of Anthracene
    fEJ200MPT->AddConstProperty("RESOLUTIONSCALE", 1.0); // Intrinsic resolution

    //fEJ200MPT->AddConstProperty("RISETIMECONSTANT", 0.9*ns); Geant4 10.1 TODO
    fEJ200MPT->AddConstProperty("FASTSCINTILLATIONRISETIME", 0.9*ns);
    fEJ200MPT->AddConstProperty("FASTTIMECONSTANT", 2.1*ns);
    fEJ200MPT->AddConstProperty("YIELDRATIO",1);// the strength of the fast component as a function of total scintillation yield

    G4double pEF = fLightYieldScale; 
    G4double pSF = pEF * 1.35;

    //light yield - data taken form V.V. Verbinski et al, Nucl. Instrum. & Meth. 65 (1968) 8-25
    G4double particleEnergy[26] = {0.001*MeV, 0.1*MeV, 0.13*MeV, 0.17*MeV, 0.2*MeV,
                                   0.24*MeV, 0.3*MeV, 0.34*MeV, 0.4*MeV, 0.48*MeV, 
                                   0.6*MeV, 0.72*MeV, 0.84*MeV, 1.*MeV, 1.3*MeV, 
                                   1.7*MeV, 2.*MeV, 2.4*MeV, 3.*MeV, 3.4*MeV, 4.*MeV, 
                                   4.8*MeV, 6.*MeV,  7.2*MeV, 8.4*MeV, 10.*MeV };

    G4double electronYield[26] = {0*pEF, 1000*pEF, 1300*pEF, 1700*pEF,
                                2000*pEF, 2400*pEF, 3000*pEF, 3400*pEF,
                                4000*pEF, 4800*pEF, 6000*pEF, 7200*pEF,
                                8400*pEF,10000*pEF, 13000*pEF, 17000*pEF,
                                20000*pEF, 24000*pEF, 30000*pEF, 34000*pEF,
                                40000*pEF, 48000*pEF, 60000*pEF, 72000*pEF,
                                84000*pEF, 100000*pEF };

    fEJ200MPT->AddProperty("ELECTRONSCINTILLATIONYIELD", particleEnergy, electronYield, 26)->SetSpline(true);
    G4double protonYield[26] = {0.6*pSF, 67.1*pSF, 88.6*pSF, 120.7*pSF,
                                146.5*pSF, 183.8*pSF, 246*pSF, 290*pSF,
                                365*pSF, 483*pSF, 678*pSF, 910*pSF,
                                1175*pSF, 562*pSF, 2385*pSF, 3660*pSF,
                                4725*pSF,6250*pSF, 8660*pSF, 10420*pSF,
                                13270*pSF,17180*pSF, 23100*pSF, 29500*pSF, 
                                36200*pSF, 45500*pSF};

    fEJ200MPT->AddProperty("PROTONSCINTILLATIONYIELD", particleEnergy, protonYield, 26)->SetSpline(true);

    G4double ionYield[26] = {0.2*pEF, 10.4*pEF, 12.7*pEF, 15.7*pEF,
                             17.9*pEF, 20.8*pEF, 25.1*pEF, 27.9*pEF,
                             31.9*pEF, 36.8*pEF, 43.6*pEF, 50.2*pEF,
                             56.9*pEF, 65.7*pEF, 81.3*pEF, 101.6*pEF,
                             116.5*pEF, 136.3*pEF, 166.15*pEF, 187.1*pEF,
                             218.6*pEF, 260.54*pEF, 323.5*pEF, 387.5*pEF,
                             451.54*pEF, 539.9*pEF };

    fEJ200MPT->AddProperty("IONSCINTILLATIONYIELD", particleEnergy, ionYield, 26)->SetSpline(true);

    fEJ200->SetMaterialPropertiesTable(fEJ200MPT);

	/////////////////////////////////////////////////////////////////
	// Silicone Optical Grease (C2H6OSi)n
	/////////////////////////////////////////////////////////////////
	
    density = 1.06*g/cm3;
    fGrease = new G4Material("Grease",density,ncomponents=4);

    fGrease->AddElement(fC, natoms=2);
    fGrease->AddElement(fH, natoms=6);
    fGrease->AddElement(fO, natoms=1);
    fGrease->AddElement(fSi, natoms=1);

    const G4int nEntries_Grease = 5;
    G4double Photon_Energy[nEntries_Grease] = { 2.757*eV, 3.102*eV, 3.312*eV, 3.545*eV, 4.136*eV };
    G4double RefractiveIndex_grease[nEntries_Grease] = {1.468, 1.473, 1.477, 1.482, 1.496};
    G4double Absorption_grease[nEntries_Grease] = { 195*mm,  195*mm, 195*mm, 195*mm, 195*mm };
    fGreaseMPT=new G4MaterialPropertiesTable();

    fGreaseMPT->AddProperty("RINDEX", Photon_Energy, RefractiveIndex_grease, nEntries_Grease);
    fGreaseMPT->AddProperty("ABSLENGTH", Photon_Energy, Absorption_grease, nEntries_Grease);

    fGrease->SetMaterialPropertiesTable(fGreaseMPT);

    fSiO2 = manNist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

    const G4int nEntries_SiO2 = 5;

    //optical properties of SiO2 - fused silica or fused quartz
    G4double RefractiveIndex_SiO2[nEntries_SiO2] = { 1.54, 1.54, 1.54, 1.54, 1.54 };
    G4double Absorption_SiO2[nEntries_SiO2] =  {125.*cm, 123.5*cm, 122.*cm, 121.*cm, 120.*cm};

    fSiO2MPT = new G4MaterialPropertiesTable();
    fSiO2MPT->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex_SiO2, nEntries_SiO2);
    fSiO2MPT->AddProperty("ABSLENGTH", PhotonEnergy, Absorption_SiO2,nEntries_SiO2);
    fSiO2->SetMaterialPropertiesTable(fSiO2MPT);

    fSil=manNist->FindOrBuildMaterial("G4_Si");

    // optical properties,
    const G4int nEntries_Sil = 5;
    G4double RefractiveReal_Si[nEntries_Sil] = { 4.293, 4.453, 4.676, 5.008, 5.587 };
    G4double RefractiveImg_Si[nEntries_Sil] = { 0.045, 0.060, 0.091, 0.150, 0.303 };
    //G4double EfficiencyIndex_Si[nEntries_Sil] = { 0.37, 0.42, 0.39, 0.36, 0.32 };
    //IDEAL DETECTOR
    G4double EfficiencyIndex_Si[nEntries_Sil] = { 1., 1., 1., 1., 1. };
    //G4double Reflective_Si[nEntries_Sil] = { 0.49, 0.45, 0.42, 0.40, 0.39};
    G4double Reflective_Si[nEntries_Sil] = { 0., 0., 0., 0., 0.};

    fSilMPT = new G4MaterialPropertiesTable();
    //fSilMPT->AddProperty("REALRINDEX", PhotonEnergy, RefractiveReal_Si, nEntries_Sil);
    //fSilMPT->AddProperty("IMAGINARYRINDEX", PhotonEnergy, RefractiveImg_Si, nEntries_Sil);
    fSilMPT->AddProperty("EFFICIENCY",   PhotonEnergy, EfficiencyIndex_Si, nEntries_Sil);
    fSilMPT->AddProperty("REFLECTIVITY",   PhotonEnergy, Reflective_Si, nEntries_Sil);

    fSil->SetMaterialPropertiesTable(fSilMPT);

	/////////////////////////////////////////////////////////////////
	// ACRYLIC (C5O2H8)n -- CRT
	/////////////////////////////////////////////////////////////////

	density = 1.19*g/cm3;
    fAcrylic = new G4Material("Acrylic", density, ncomponents=3);
    
    fAcrylic->AddElement(fC, natoms=5);
    fAcrylic->AddElement(fH, natoms=8);
    fAcrylic->AddElement(fO, natoms=2);

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
	
	MPT_Acrylic->AddProperty("ABSLENGTH", energy2, abslength, 25);

	fAcrylic->SetMaterialPropertiesTable(MPT_Acrylic);

	/////////////////////////////////////////////////////////////////
	// Optical Surfaces
	/////////////////////////////////////////////////////////////////
	
    fTeflonOpticalSurface = new G4OpticalSurface("TeflonSurface", glisur, ground, dielectric_metal, 0.1); // polish level
    //DPL Testing the other finishes for teflon

    fTeflonOpticalSurface->SetFinish(ground);
    //fTeflonOpticalSurface->SetType(dielectric_metal);
    fTeflonOpticalSurface->SetType(dielectric_dielectric);
    fTeflonOpticalSurface->SetMaterialPropertiesTable(fTeflonMPT);

    //SiliconOpticalSurface
    fSiliconPMOpticalSurface = new G4OpticalSurface("SiPMSurface");
    fSiliconPMOpticalSurface->SetType(dielectric_metal);
    fSiliconPMOpticalSurface->SetFinish(polished);
    fSiliconPMOpticalSurface->SetModel(glisur);
    fSiliconPMOpticalSurface->SetMaterialPropertiesTable(fSilMPT);

    G4Material *Al=manNist->FindOrBuildMaterial("G4_Al");

    G4cout<<"manNist->FindOrBuildMaterial(\"G4_Al\")"<<Al<<G4endl;

    G4Material *Mylar=manNist->FindOrBuildMaterial("G4_MYLAR");

    fMylar=new G4Material("AluninizedMylar",density=1.39*g/cm3,ncomponents=2);
    fMylar->AddMaterial(Mylar,0.8);
    fMylar->AddMaterial(Al,0.2);

    const G4int nEntries_Mylar = 5;
    G4double RefractiveReal_Mylar[nEntries_Mylar]={0.81257,0.72122,0.63324,0.55571,0.48787};
    G4double RefractiveImg_Mylar[nEntries_Mylar]={6.0481,5.7556,5.4544,5.1464,4.8355};

    fMylarMPT=new G4MaterialPropertiesTable();
    fMylarMPT->AddProperty("REALRINDEX", PhotonEnergy,RefractiveReal_Mylar,nEntries_Mylar);
    fMylarMPT->AddProperty("IMAGINARYRINDEX", PhotonEnergy,RefractiveImg_Mylar,nEntries_Mylar);

    fMylarOpticalSurface=new G4OpticalSurface("MylarSurface",glisur,polished,dielectric_metal,1.0);

    return;
}

void nDetConstruction::buildSiPMs(const G4double &zOffset/*=0*/) {


    if(fGeometry == "ellipse" || fGeometry == "rectangle") {
        greaseX = SiPM_dimension;
        greaseZ = SiPM_dimension;
        qwSiPMx = SiPM_dimension;
        qwSiPMz = SiPM_dimension;
        psSiPMx = SiPM_dimension;
        psSiPMz = SiPM_dimension;
    }

    fNdetectors=2;
    G4int factor=1;

    G4ThreeVector mother_position(0, 0, 0);


    //The Grease
    G4Box* grease_solidV = new G4Box("grease", greaseX, greaseY, greaseZ);
    grease_logV = new G4LogicalVolume(grease_solidV, fGrease, "grease_logV", 0, 0, 0);
    G4VisAttributes* grease_VisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0));//red
    //grease_VisAtt->SetForceSolid(true);
    grease_logV->SetVisAttributes(grease_VisAtt);

    G4AssemblyVolume *theGrease=new G4AssemblyVolume();

    for(G4int m=0;m<fNdetectors;m++) {

        if(m==1)
            factor=-1;
        G4ThreeVector position(0, 0, factor*(fDetectorLength/ 2+greaseY+zOffset));
        G4RotationMatrix *rot = new G4RotationMatrix();
        rot->rotateX(90 * deg);

        theGrease->AddPlacedVolume(grease_logV, position, rot);
    }
    theGrease->MakeImprint(assembly_logV,mother_position,NULL);


    //The SiPM window
    G4Box* qwSiPM_solidV = new G4Box("qwSiPM_solidV", qwSiPMx, qwSiPMy, qwSiPMz);

    qwSiPM_logV = new G4LogicalVolume(qwSiPM_solidV,fSiO2,"SiPMwindow_log");

    G4VisAttributes* window_VisAtt= new G4VisAttributes(G4Colour(0.0,1.0,1.0));//cyan
    //window_VisAtt->SetForceSolid(true);
    qwSiPM_logV->SetVisAttributes(window_VisAtt);

    factor=1;

    G4AssemblyVolume* theSiPMWindows=new G4AssemblyVolume();
    for(G4int m=0;m<fNdetectors;m++) {

        if(m==1)
            factor=-1;
        G4ThreeVector position(0, 0, factor*(fDetectorLength/2 + qwSiPMy+2*greaseY+zOffset));
        G4RotationMatrix *rot = new G4RotationMatrix();
        rot->rotateX(90 * deg);

        theSiPMWindows->AddPlacedVolume(qwSiPM_logV, position, rot);
    }
    theSiPMWindows->MakeImprint(assembly_logV,mother_position,NULL);

    //the Sensitive surface

    G4Box* psSiPM_solidV = new G4Box("psSiPM_solidV", psSiPMx, psSiPMy, psSiPMz);


    psSiPM_logV = new G4LogicalVolume(psSiPM_solidV,fSil, "psSiPM_logV", 0, 0, 0);

    G4LogicalSkinSurface *SiPMSurface=new G4LogicalSkinSurface("theSiPM",psSiPM_logV,fSiliconPMOpticalSurface);

    G4VisAttributes* psSiPM_VisAtt= new G4VisAttributes(G4Colour(0.1,1.0,0.0));//magenta
    psSiPM_VisAtt->SetForceSolid(true);
    psSiPM_logV->SetVisAttributes(psSiPM_VisAtt);
    factor=1;
    G4AssemblyVolume* psSiPMAssembly = new G4AssemblyVolume();
    for(G4int m=0;m<fNdetectors;m++) {

        if(m==1)
            factor=-1;
        G4ThreeVector position(0, 0, factor*(fDetectorLength / 2  + 2*qwSiPMy+2*greaseY+psSiPMy+zOffset));
        G4RotationMatrix *rot = new G4RotationMatrix();
        rot->rotateX(90 * deg);

        psSiPMAssembly->AddPlacedVolume(psSiPM_logV, position, rot);
    }
    psSiPMAssembly->MakeImprint(assembly_logV,mother_position,NULL);

    return;
}

void nDetConstruction::ConstructSDandField(){
}


G4VSolid* nDetConstruction::ConstructEllipse(G4String name,G4ThreeVector dimensions,G4double thickness) {

    //First we build the box
    G4double xdimensionBox=dimensions.x()/2.;
    G4double ydimensionBox=thickness/2.;
    G4double zdimensionBox=dimensions.z()/2.;

    G4Box *theBox = new G4Box("theBox"+name,xdimensionBox,ydimensionBox,zdimensionBox);

    //Now the trapezoid

    G4double trapezoidlength=fTrapezoidLength/2.;

    G4double dy1=thickness/2.;
    G4double dy2=thickness/2.;
    G4double dx1=dimensions.x()/2.;
    G4double dx2=dimensions.y()/2.;
    G4double dz=trapezoidlength;

    G4Trd *theTrapezoid1=new G4Trd("theTrapezoid1"+name, dx1, dx2, dy1, dy2, dz);

    G4Trd *theTrapezoid2=new G4Trd("theTrapezoid2"+name, dx1, dx2, dy1, dy2, dz);

    G4ThreeVector translation(0,0,zdimensionBox+dz);

    G4UnionSolid *part1=new G4UnionSolid("Box+Trapezoid"+name,theBox,theTrapezoid1,0,translation);

    translation=-1*translation;

    G4RotationMatrix *rot=new G4RotationMatrix();

    rot->rotateX(180*deg);


    G4UnionSolid *theSolid=new G4UnionSolid("Trapezoid+Box+Trapezoid"+name,part1,theTrapezoid2,rot,translation);

    return theSolid;

}



G4VSolid* nDetConstruction::ConstructHexagon(G4String name,G4double radius,G4double thickness,G4ThreeVector Boxdimensions){



    G4double minDiskRadius=0*mm;
    G4double maxDiskRadius=radius;


    G4double startAngle=0*deg;
    G4double endAngle=360*deg;

    G4int nSides=6;
    G4int nPlanes=2;

    const G4double zplanes[]={-thickness,thickness};
    const G4double maxRadius[]={maxDiskRadius,maxDiskRadius};
    const G4double minRadius[]={minDiskRadius,minDiskRadius};


    G4Polyhedra* theHexagon = new G4Polyhedra("Hexagon",startAngle,endAngle,nSides,nPlanes,zplanes,minRadius,maxRadius);

    G4double BoxX = Boxdimensions.x();
    G4double BoxY = Boxdimensions.y();
    G4double BoxZ = Boxdimensions.z();

    G4Box *theHole = new G4Box("Hole", BoxX, BoxY, BoxZ);

    G4SubtractionSolid *theSolid = new G4SubtractionSolid(name, theHexagon, theHole);

    //Define the Hole

    return theSolid;

}


G4VSolid* nDetConstruction::ConstructHexagon(G4String name,G4double radius,G4double thickness){


G4double minDiskRadius=0*mm;
G4double maxDiskRadius=radius;


G4double startAngle=0*deg;
G4double endAngle=360*deg;

G4int nSides=6;
G4int nPlanes=2;

const G4double zplanes[]={-thickness,thickness};
const G4double maxRadius[]={maxDiskRadius,maxDiskRadius};
const G4double minRadius[]={minDiskRadius,minDiskRadius};


G4Polyhedra* theHexagon = new G4Polyhedra(name,startAngle,endAngle,nSides,nPlanes,zplanes,minRadius,maxRadius);


return theHexagon;

}

void nDetConstruction::GetSegmentFromCopyNum(const G4int &copyNum, G4int &col, G4int &row) const {
	if(copyNum > 0){ // Scintillator
		col = ((copyNum-1) / fNumRows) + 1;
		row = ((copyNum-1) % fNumRows) + 1;
	}
	else{ // Outer wrappings (mylar, teflon, etc) (id==0)
		col = 0;
		row = 0;
	}
}

bool nDetConstruction::AddDetectedPhoton(const G4Step *step, const double &mass/*=1*/){
    if(step->GetPostStepPoint()->GetPosition().z() > 0){
        if(center[0].addPoint(step, mass))
            return true;
    }
    if(step->GetPostStepPoint()->GetPosition().z() < 0){
        if(center[1].addPoint(step, mass))
            return true;
    }
    return false;
}

void nDetConstruction::GetDetectedPhotons(size_t &numLeft, size_t &numRight){
	numLeft = center[0].getNumDetected();
	numRight = center[1].getNumDetected();
}

void nDetConstruction::Clear(){
    center[0].clear();
    center[1].clear();
}

void nDetConstruction::UpdateGeometry(){
    // clean-up previous geometry
    G4SolidStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    
    //define new one
    G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();

    if(PmtIsSegmented())
      setSegmentedPmt(fNumColumnsPmt, fNumRowsPmt, SiPM_dimension*2, SiPM_dimension*2);
}

void nDetConstruction::buildEllipse2() {


    G4double teflonThickness=0.22*mm;


    G4double plasticLength = fDetectorLength-2*fTrapezoidLength ;

    //fMylarThickness=0.025*mm; //25 um mylar



    //First we build the box
    G4double xdimension=fDetectorWidth+teflonThickness+fMylarThickness;
    G4double ydimension=0.24*inch+teflonThickness+2*fMylarThickness;
    G4double zdimension=plasticLength+teflonThickness;

    G4ThreeVector dimensions(xdimension,ydimension,zdimension);

    G4VSolid* theEllipse=ConstructEllipse("wrap",dimensions,ydimension);

    G4Box *wrappinBox=new G4Box("theBox1",ydimension/2,ydimension/2,(greaseY+qwSiPMy+psSiPMy));

    //G4Box *wrappinBox2=new G4Box("theBox2",ydimension,ydimension,(greaseY+qwSiPMy+psSiPMy));

    G4ThreeVector translation11(0,0,(plasticLength/2+fTrapezoidLength)+1*(greaseY+qwSiPMy+psSiPMy));


    G4UnionSolid *theWrapping0=new G4UnionSolid("Wrap",theEllipse,wrappinBox,0,translation11);


    translation11=-1*translation11;

    G4UnionSolid *theWrapping=new G4UnionSolid("Wrap",theWrapping0,wrappinBox,0,translation11);


    assembly_logV=new G4LogicalVolume(theWrapping,fTeflon,"Wrapping_log");

    //In Ellipse
    fWrapSkinSurface=new G4LogicalSkinSurface("Wrapping",assembly_logV,fTeflonOpticalSurface);
    //fWrapSkinSurface=new G4LogicalSkinSurface("Wrapping",assembly_logV,0); // TODO DPL  Equivalent to black tape

    G4VPhysicalVolume *Wrapping_physVol=new G4PVPlacement(0,G4ThreeVector(0,0,0),assembly_logV,"Wrap",expHall_logV,0,0,true);


    //Building the Scintillator

    xdimension=fDetectorWidth;
    ydimension=0.24*inch;
    zdimension=plasticLength;

    G4ThreeVector dimensions2(xdimension,ydimension,zdimension);

    G4VSolid *theScint=ConstructEllipse("Scint",dimensions2,ydimension);

    ej200_logV=new G4LogicalVolume(theScint,fEJ200,"Scint_log");

    G4VisAttributes* ej200_VisAtt= new G4VisAttributes(G4Colour(0.0,0.0,1.0));//blue
    ej200_logV->SetVisAttributes(ej200_VisAtt);

    G4VPhysicalVolume *scint_phys=new G4PVPlacement(0,G4ThreeVector(0,0,0),ej200_logV,"Scint",assembly_logV,0,0,true);

    //Building the Mylar covers

    if(fMylarThickness>0) {

        xdimension = fDetectorWidth;
        ydimension = 0.24 * inch;
        zdimension = plasticLength;

        G4ThreeVector dimensions3(xdimension, ydimension, zdimension);

        G4VSolid *theMylar = ConstructEllipse("mylar", dimensions3, fMylarThickness);

        mylar_logV = new G4LogicalVolume(theMylar, fMylar, "mylar");

        G4VisAttributes *mylar_VisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0)); //magenta
        mylar_VisAtt->SetForceSolid(true);
        mylar_logV->SetVisAttributes(mylar_VisAtt);
        //assembly_logV->SetVisAttributes(mylar_VisAtt);

        G4ThreeVector position(0, dimensions3.y() / 2 + fMylarThickness / 2, 0);
        G4PVPlacement *mylar_phys = new G4PVPlacement(0, position, mylar_logV, "Mylar1", assembly_logV, true, 0, true);
        G4ThreeVector position2(0, -dimensions3.y() / 2 - fMylarThickness / 2, 0);
        G4PVPlacement *mylar_phys2 = new G4PVPlacement(0, position2, mylar_logV, "Mylar2", assembly_logV, true, 0,
                                                       true);


        fMylarSurface = new G4LogicalBorderSurface("Mylar", scint_phys, mylar_phys, fMylarOpticalSurface);
        G4LogicalBorderSurface *log2 = new G4LogicalBorderSurface("Mylar2", scint_phys, mylar_phys2,
                                                                  fMylarOpticalSurface);

    }

    buildSiPMs();
    return;

}

void nDetConstruction::buildRectangle(){
	const G4double cellWidth = (fDetectorThickness-2*fNumColumns*fMylarThickness)/fNumColumns;
	const G4double cellHeight = (fDetectorThickness-2*fNumRows*fMylarThickness)/fNumRows;

    const G4double greaseThickness = 1*mm;
    const G4double windowThickness = 1*mm;

	G4double assemblyLength = fDetectorLength + 2*greaseThickness + 4*windowThickness;
	if(2*SiPM_dimension < fDetectorWidth) // Account for light guides
	    assemblyLength += 2*fTrapezoidLength + 2*greaseThickness;

    G4Box *theRectangle = new G4Box("rectangle", fDetectorWidth/2, fDetectorThickness/2, fDetectorLength/2);

	G4Box *assembly = new G4Box("assembly", fDetectorWidth/2, fDetectorThickness/2, assemblyLength/2);
    assembly_logV = new G4LogicalVolume(assembly, fMylar, "assembly_logV");
    G4VisAttributes* assembly_VisAtt = new G4VisAttributes();
    assembly_VisAtt->SetVisibility(false);
    assembly_logV->SetVisAttributes(assembly_VisAtt);

    fWrapSkinSurface = new G4LogicalSkinSurface("wrapping", assembly_logV, fMylarOpticalSurface); //Outside

    G4VPhysicalVolume *totalDetector_physVol = new G4PVPlacement(0, G4ThreeVector(0,0,0), assembly_logV, "Wrapping", expHall_logV, 0, 0, true);//fCheckOverlaps);

    //Building the Scintillator
    G4Box *theScint = new G4Box("scintillator", cellWidth/2, cellHeight/2, fDetectorLength/2);

    ej200_logV = new G4LogicalVolume(theScint, fEJ200, "scint_log");

    G4VisAttributes* ej200_VisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); //blue
    G4VisAttributes *mylar_VisAtt = new G4VisAttributes();
    mylar_VisAtt->SetColor(1, 0, 1, 0.5); // Alpha=50%
    mylar_VisAtt->SetForceSolid(true);
    
    ej200_logV->SetVisAttributes(ej200_VisAtt);

    G4Box *mylarCellBoxTop, *mylarCellBoxSide;
    G4LogicalVolume *mylarCellLogTop=NULL;
    G4LogicalVolume *mylarCellLogSide=NULL;
    G4PVPlacement *cellMylarWrap[4];

    //Building the Mylar covers
    if(fMylarThickness > 0){
		mylarCellBoxTop = new G4Box("mylar", cellWidth/2, fMylarThickness/2, fDetectorLength/2);
		mylarCellLogTop = new G4LogicalVolume(mylarCellBoxTop, fMylar, "mylar_log");

		mylarCellBoxSide = new G4Box("mylar", fMylarThickness/2, cellHeight/2, fDetectorLength/2);
		mylarCellLogSide = new G4LogicalVolume(mylarCellBoxSide, fMylar, "mylar_log");

		mylarCellLogTop->SetVisAttributes(mylar_VisAtt);
		mylarCellLogSide->SetVisAttributes(mylar_VisAtt);
	}
	
	for(int col = 0; col < fNumColumns; col++){
		for(int row = 0; row < fNumRows; row++){
			G4ThreeVector center(-fDetectorWidth/2 + (2*col+1)*fMylarThickness + (col+0.5)*cellWidth, -fDetectorThickness/2 + (2*row+1)*fMylarThickness + (row+0.5)*cellHeight, 0);
		
			// The physical scintillator bar
			G4PVPlacement *cellPhysical = new G4PVPlacement(0, center, ej200_logV, "Scint", assembly_logV, 0, col*fNumRows+row+1, fCheckOverlaps); // Copy numbers (segment IDs), indexed from 1
		
			if(fMylarThickness > 0){ // Wrap the bar in mylar
				cellMylarWrap[0] = new G4PVPlacement(0, G4ThreeVector(center.getX(), center.getY() + cellHeight/2 + fMylarThickness/2, 0), mylarCellLogTop, "Mylar", assembly_logV, true, 0, fCheckOverlaps); // Top layer (+Y)
				cellMylarWrap[1] = new G4PVPlacement(0, G4ThreeVector(center.getX(), center.getY() - cellHeight/2 - fMylarThickness/2, 0), mylarCellLogTop, "Mylar", assembly_logV, true, 0, fCheckOverlaps); // Bottom layer (-Y)
				cellMylarWrap[2] = new G4PVPlacement(0, G4ThreeVector(center.getX() - cellWidth/2 - fMylarThickness/2, center.getY(), 0), mylarCellLogSide, "Mylar", assembly_logV, true, 0, fCheckOverlaps); // -X layer
				cellMylarWrap[3] = new G4PVPlacement(0, G4ThreeVector(center.getX() + cellWidth/2 + fMylarThickness/2, center.getY(), 0), mylarCellLogSide, "Mylar", assembly_logV, true, 0, fCheckOverlaps); // +X layer
				for(size_t i = 0; i < 4; i++)
					new G4LogicalBorderSurface("Mylar", cellPhysical, cellMylarWrap[i], fMylarOpticalSurface);
			}
		}
	}

    G4VisAttributes* grease_VisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // red

    G4double greaseZ = fDetectorLength/2 + greaseThickness/2;
    G4double windowZ = fDetectorLength/2 + greaseThickness + windowThickness/2;
    G4double sensitiveZ = fDetectorLength/2 + greaseThickness + 1.5*windowThickness;

    if(2*SiPM_dimension < fDetectorWidth){ // Build the light guides (if needed)
        G4double trapezoidW1 = fDetectorWidth - 2*fMylarThickness;
        G4double trapezoidW2 = fDetectorThickness - 2*fMylarThickness;
    
        //The grease
        G4Box* grease_solidV = new G4Box("grease", trapezoidW1/2, trapezoidW2/2, greaseThickness/2);
        grease_logV = new G4LogicalVolume(grease_solidV, fGrease, "grease_logV");
        grease_logV->SetVisAttributes(grease_VisAtt);
    
        new G4PVPlacement(0, G4ThreeVector(0, 0, greaseZ), grease_logV, "Grease", assembly_logV, true, 0, fCheckOverlaps);
        new G4PVPlacement(0, G4ThreeVector(0, 0, -greaseZ), grease_logV, "Grease", assembly_logV, true, 0, fCheckOverlaps);
    
        // BE CAREFUL, for some reason SiPM_dimension is set to the user defined SiPM dimension divided by 2 in nDetConstructionMessenger!!!
        G4Trd *lightGuide = new G4Trd("lightGuide", trapezoidW1/2, SiPM_dimension, trapezoidW2/2, SiPM_dimension, fTrapezoidLength/2);
        G4LogicalVolume *lightGuideLog = new G4LogicalVolume(lightGuide, fAcrylic, "lightGuide_logV");
        
        G4double trapezoidZ = fDetectorLength/2 + fTrapezoidLength/2 + greaseThickness;
        
        // "Left" side (+z)
    	G4PVPlacement *trapPhysicalL = new G4PVPlacement(0, G4ThreeVector(0, 0, trapezoidZ), lightGuideLog, "Acrylic", assembly_logV, true, 0, fCheckOverlaps);
    	
    	// "Right" side (-z)
    	G4RotationMatrix *rot = new G4RotationMatrix;
        rot->rotateY(180*deg);
    	G4PVPlacement *trapPhysicalR = new G4PVPlacement(rot, G4ThreeVector(0, 0, -trapezoidZ), lightGuideLog, "Acrylic", assembly_logV, true, 0, fCheckOverlaps);
    	
    	// Offset the PSPMT to account for the light-guide and another layer of grease.
    	greaseZ += greaseThickness + fTrapezoidLength;
    	windowZ += greaseThickness + fTrapezoidLength;
    	sensitiveZ += greaseThickness + fTrapezoidLength;
    	
    	if(fMylarThickness > 0){ // Construct the mylar cover flaps
    	    G4double topFlapOffset = fDetectorThickness/2-fMylarThickness-SiPM_dimension;
    	    G4double sideFlapOffset = fDetectorWidth/2-fMylarThickness-SiPM_dimension;
    	
    	    // Top and bottom
    	    G4double zprime = std::sqrt(std::pow(topFlapOffset, 2.0) + std::pow(fTrapezoidLength, 2.0));
    	    std::vector<G4TwoVector> trapCoverPoints;
    	    trapCoverPoints.push_back(G4TwoVector(-fDetectorWidth/2+fMylarThickness, zprime/2));
    	    trapCoverPoints.push_back(G4TwoVector(fDetectorWidth/2-fMylarThickness, zprime/2));
    	    trapCoverPoints.push_back(G4TwoVector(SiPM_dimension, -zprime/2));
    	    trapCoverPoints.push_back(G4TwoVector(-SiPM_dimension, -zprime/2));

    	    G4ExtrudedSolid *trapMylarCover = new G4ExtrudedSolid("", trapCoverPoints, fMylarThickness/2, G4TwoVector(0, 0), 1, G4TwoVector(0, 0), 1);
    	    G4LogicalVolume *trapMylarCover_logV = new G4LogicalVolume(trapMylarCover, fMylar, "trapMylarCover_logV");
    	    trapMylarCover_logV->SetVisAttributes(mylar_VisAtt);

    	    // Sides
	        zprime = std::sqrt(std::pow(sideFlapOffset, 2.0) + std::pow(fTrapezoidLength, 2.0));
    	    std::vector<G4TwoVector> sideTrapCoverPoints;
    	    sideTrapCoverPoints.push_back(G4TwoVector(-fDetectorThickness/2+fMylarThickness, zprime/2));
    	    sideTrapCoverPoints.push_back(G4TwoVector(fDetectorThickness/2-fMylarThickness, zprime/2));
    	    sideTrapCoverPoints.push_back(G4TwoVector(SiPM_dimension, -zprime/2));
    	    sideTrapCoverPoints.push_back(G4TwoVector(-SiPM_dimension, -zprime/2));
    	    
    	    G4ExtrudedSolid *sideTrapMylarCover = new G4ExtrudedSolid("", sideTrapCoverPoints, fMylarThickness/2, G4TwoVector(0, 0), 1, G4TwoVector(0, 0), 1);
    	    G4LogicalVolume *sideTrapMylarCover_logV = new G4LogicalVolume(sideTrapMylarCover, fMylar, "sideTrapMylarCover_logV");
    	    sideTrapMylarCover_logV->SetVisAttributes(mylar_VisAtt);
    	    
    	    // Handle polygon rotation
    	    G4RotationMatrix *trapCoverRotL[4] = {new G4RotationMatrix(), new G4RotationMatrix(), new G4RotationMatrix(), new G4RotationMatrix()};
            G4RotationMatrix *trapCoverRotR[4] = {new G4RotationMatrix(), new G4RotationMatrix(), new G4RotationMatrix(), new G4RotationMatrix()};	       	    
    	    
    	    G4double trapAngle = std::atan2(fTrapezoidLength, topFlapOffset);

            trapCoverRotL[0]->rotateX(trapAngle);
            
            trapCoverRotR[0]->rotateY(pi);
            trapCoverRotR[0]->rotateX(trapAngle);
            
            trapCoverRotL[1]->rotateX(pi-trapAngle);
            
            trapCoverRotR[1]->rotateY(pi);
            trapCoverRotR[1]->rotateX(pi-trapAngle);
            
            trapAngle = std::atan2(fTrapezoidLength, sideFlapOffset);
            
            trapCoverRotL[2]->rotateY(trapAngle);
            trapCoverRotL[2]->rotateZ(270*deg);
            
            trapCoverRotR[2]->rotateX(pi);
            trapCoverRotR[2]->rotateY(trapAngle);
            trapCoverRotR[2]->rotateZ(270*deg);
            
            trapCoverRotL[3]->rotateY(pi-trapAngle);
            trapCoverRotL[3]->rotateZ(270*deg);

            trapCoverRotR[3]->rotateX(pi);
            trapCoverRotR[3]->rotateY(pi-trapAngle);
            trapCoverRotR[3]->rotateZ(270*deg);
    	    
    	    // Place the mylar covers
    	    zprime = 0.5*fMylarThickness/std::sqrt(1+std::pow(std::tan(trapAngle), 2.0));
    	    G4double hprime = zprime*std::tan(trapAngle);
    	    
    	    G4double verticalOffset = 0.5*(topFlapOffset) + SiPM_dimension + hprime;
    	    G4double horizontalOffset = 0.5*(sideFlapOffset) + SiPM_dimension + hprime;
    	    trapezoidZ += zprime;
    	    
    	    cellMylarWrap[0] = new G4PVPlacement(trapCoverRotL[0], G4ThreeVector(0, verticalOffset, trapezoidZ), trapMylarCover_logV, "Mylar", assembly_logV, true, 0, fCheckOverlaps);
    	    cellMylarWrap[1] = new G4PVPlacement(trapCoverRotL[1], G4ThreeVector(0, -verticalOffset, trapezoidZ), trapMylarCover_logV, "Mylar", assembly_logV, true, 0, fCheckOverlaps);
    	    cellMylarWrap[2] = new G4PVPlacement(trapCoverRotL[2], G4ThreeVector(-horizontalOffset, 0, trapezoidZ), sideTrapMylarCover_logV, "Mylar", assembly_logV, true, 0, fCheckOverlaps);
    	    cellMylarWrap[3] = new G4PVPlacement(trapCoverRotL[3], G4ThreeVector(horizontalOffset, 0, trapezoidZ), sideTrapMylarCover_logV, "Mylar", assembly_logV, true, 0, fCheckOverlaps);
    	    for(size_t i = 0; i < 4; i++)
    	   	    new G4LogicalBorderSurface("Mylar", trapPhysicalL, cellMylarWrap[i], fMylarOpticalSurface);

            // Now for the other side
    	    cellMylarWrap[0] = new G4PVPlacement(trapCoverRotR[0], G4ThreeVector(0, verticalOffset, -trapezoidZ), trapMylarCover_logV, "Mylar", assembly_logV, true, 0, fCheckOverlaps);
    	    cellMylarWrap[1] = new G4PVPlacement(trapCoverRotR[1], G4ThreeVector(0, -verticalOffset, -trapezoidZ), trapMylarCover_logV, "Mylar", assembly_logV, true, 0, fCheckOverlaps);
    	    cellMylarWrap[2] = new G4PVPlacement(trapCoverRotR[2], G4ThreeVector(-horizontalOffset, 0, -trapezoidZ), sideTrapMylarCover_logV, "Mylar", assembly_logV, true, 0, fCheckOverlaps);
    	    cellMylarWrap[3] = new G4PVPlacement(trapCoverRotR[3], G4ThreeVector(horizontalOffset, 0, -trapezoidZ), sideTrapMylarCover_logV, "Mylar", assembly_logV, true, 0, fCheckOverlaps);
    	    for(size_t i = 0; i < 4; i++)
    	        new G4LogicalBorderSurface("Mylar", trapPhysicalR, cellMylarWrap[i], fMylarOpticalSurface);
    	}
	}

	// Build the sensitive PMT surfaces.
	const G4String name = "psSiPM";

	// BE CAREFUL, for some reason SiPM_dimension is set to the user defined SiPM dimension divided by 2 in nDetConstructionMessenger!!!
	G4Box* grease_solidV2 = new G4Box("grease", SiPM_dimension, SiPM_dimension, greaseThickness/2);
	grease_logV = new G4LogicalVolume(grease_solidV2, fGrease, "grease_logV");
	grease_logV->SetVisAttributes(grease_VisAtt);
	
	new G4PVPlacement(0, G4ThreeVector(0, 0, greaseZ), grease_logV, "Grease", assembly_logV, true, 0, fCheckOverlaps);
    new G4PVPlacement(0, G4ThreeVector(0, 0, -greaseZ), grease_logV, "Grease", assembly_logV, true, 0, fCheckOverlaps);
    
    // The quartz window
    G4Box* window_solidV = new G4Box("window_solidV", SiPM_dimension, SiPM_dimension, windowThickness/2);
    G4LogicalVolume *window_logV = new G4LogicalVolume(window_solidV, fSiO2, "window_logV");
    G4VisAttributes* window_VisAtt= new G4VisAttributes(G4Colour(0.0, 1.0, 1.0)); // cyan
    window_logV->SetVisAttributes(window_VisAtt);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, windowZ), window_logV, "Quartz", assembly_logV, true, 0, fCheckOverlaps);
    new G4PVPlacement(0, G4ThreeVector(0, 0, -windowZ), window_logV, "Quartz", assembly_logV, true, 0, fCheckOverlaps);
    
    // The photon sensitive surface
    G4Box *sensitive_solidV = new G4Box(name+"_solidV", SiPM_dimension, SiPM_dimension, windowThickness/2);
    G4LogicalVolume *sensitive_logV = new G4LogicalVolume(sensitive_solidV, fSil, name+"_logV");
    G4LogicalSkinSurface *sensitiveSurface = new G4LogicalSkinSurface(name, sensitive_logV, fSiliconPMOpticalSurface);
    G4VisAttributes *sensitive_VisAtt = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75)); // grey
    sensitive_VisAtt->SetForceSolid(true);
    sensitive_logV->SetVisAttributes(sensitive_VisAtt);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, sensitiveZ), sensitive_logV, name, assembly_logV, true, 0, fCheckOverlaps);
    new G4PVPlacement(0,  G4ThreeVector(0, 0, -sensitiveZ), sensitive_logV, name, assembly_logV, true, 0, fCheckOverlaps);

    return;
}

void nDetConstruction::buildDisk2() {
    //************* Wrapping materials (such as Teflon) and EJ200  *************

    //-----------------------------------------------------------
    //     The assembly volume is defined as a wrapping material
    //-----------------------------------------------------------

    greaseX = SiPM_dimension;
    greaseZ = SiPM_dimension;
    qwSiPMx = SiPM_dimension;
    qwSiPMz = SiPM_dimension;
    psSiPMx = SiPM_dimension;
    psSiPMz = SiPM_dimension;

    fTeflonThickness = 0.11*mm;

    G4double DiskRadius = fHexagonRadius + fTeflonThickness;

    G4double thickness = SiPM_dimension + 2*fTeflonThickness + 2*fMylarThickness;


    G4double BoxX = SiPM_dimension + 2 * fTeflonThickness;
    G4double BoxY = SiPM_dimension + 2 * fTeflonThickness;
    G4double BoxZ = SiPM_dimension + 2 * (fTeflonThickness + greaseY + qwSiPMy + psSiPMy);

    G4ThreeVector BoxDimensions(BoxX,BoxY,BoxZ);

    G4VSolid *theWrapping = ConstructHexagon("theWrapping",DiskRadius,thickness,BoxDimensions);


    assembly_logV = new G4LogicalVolume(theWrapping,fTeflon,"Wrap");

    fWrapSkinSurface = new G4LogicalSkinSurface("WrapSkin",assembly_logV,fTeflonOpticalSurface);


    G4RotationMatrix *rotationMatrix_Disk = new G4RotationMatrix;
    rotationMatrix_Disk->rotateY(90 * deg);

    assembly_physV = new G4PVPlacement(rotationMatrix_Disk,
                                       G4ThreeVector(assemblyPx, assemblyPy, assemblyPz),
                                       assembly_logV,
                                       "assembly_physV",
                                       expHall_logV,
                                       false, 0);

    //----------------------------------------------------------------
    //    scintillator bars
    //----------------------------------------------------------------


    DiskRadius = fHexagonRadius;

    thickness = SiPM_dimension;

    BoxX = BoxZ;
    BoxY = BoxZ;
    BoxZ = SiPM_dimension + 2 * fTeflonThickness;

    G4ThreeVector BoxDimensions2(BoxX,BoxY,BoxZ);

    G4VSolid *theScint = ConstructHexagon("theScint",DiskRadius,thickness,BoxDimensions2);

    ej200_logV = new G4LogicalVolume(theScint, fEJ200, "ej200_logV", 0, 0, 0);

    G4VisAttributes *ej200_VisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));//blue
    //ej200_VisAtt->SetForceSolid(true);
    ej200_logV->SetVisAttributes(ej200_VisAtt);

    //ej200_logV->SetSensitiveDetector(fScintSD);

    // create assembly of scintillator
    //

    G4RotationMatrix *rotationMatrix_Scint = new G4RotationMatrix;


    G4PVPlacement *ej200_physV = new G4PVPlacement(rotationMatrix_Scint,
                                       G4ThreeVector(assemblyPx, assemblyPy, assemblyPz),
                                       ej200_logV,
                                       "Scint_physV",
                                       assembly_logV,
                                       false, 0);



    //the Mylar sheets

    if(fMylarThickness>0) {
        G4VSolid *theMylar = ConstructHexagon("mylar", DiskRadius, fMylarThickness, BoxDimensions2);

        mylar_logV = new G4LogicalVolume(theMylar, fMylar, "mylar");

        G4VisAttributes *mylar_VisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0)); //magenta
        mylar_VisAtt->SetForceSolid(true);
        mylar_logV->SetVisAttributes(mylar_VisAtt);
        //assembly_logV->SetVisAttributes(mylar_VisAtt);

        G4ThreeVector position(0, 0, SiPM_dimension + fMylarThickness);
        G4PVPlacement *mylar_phys = new G4PVPlacement(0, position, mylar_logV, "Mylar1", assembly_logV, true, 0, true);
        G4ThreeVector position2(0, 0, -SiPM_dimension - fMylarThickness);
        G4PVPlacement *mylar_phys2 = new G4PVPlacement(0, position2, mylar_logV, "Mylar2", assembly_logV, true, 0,
                                                       true);


        fMylarSurface = new G4LogicalBorderSurface("Mylar", assembly_physV, mylar_phys, fMylarOpticalSurface);
        G4LogicalBorderSurface *log2 = new G4LogicalBorderSurface("Mylar2", ej200_physV, mylar_phys2,
                                                                  fMylarOpticalSurface);

    }


    // grease Sanit-Gobain BC-630 silicone grease
    // optical grease

    // grease box ...
    G4Box *grease_solidV = new G4Box("grease", greaseX, greaseY, greaseZ);
    grease_logV = new G4LogicalVolume(grease_solidV, fGrease, "grease_logV", 0, 0, 0);
    G4VisAttributes *grease_VisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));//red
    //grease_VisAtt->SetForceSolid(true);
    grease_logV->SetVisAttributes(grease_VisAtt);

    // create assembly of grease
    //
    G4AssemblyVolume *greaseAssembly = new G4AssemblyVolume();
    // Rotation and translation of the grease logical volume
    G4RotationMatrix rotationMatrix_grease;
    G4ThreeVector position_grease(0., 0., 0.);
    G4Transform3D transform3D_grease;

    for (G4int i = 0; i < 2; i++) {
        //position_grease.setX(assemblyBoxX - assemblyBoxX/15*(1+x_index*2)/2);
        if (i == 0) {
            position_grease.setX(0);
            position_grease.setZ(0);
            position_grease.setY(BoxY - greaseY);

            transform3D_grease = G4Transform3D(rotationMatrix_grease, position_grease);
            greaseAssembly->AddPlacedVolume(grease_logV, transform3D_grease);

            position_grease.setY(-1 * (BoxY - greaseY));

            transform3D_grease = G4Transform3D(rotationMatrix_grease, position_grease);
            greaseAssembly->AddPlacedVolume(grease_logV, transform3D_grease);
        } else {
            rotationMatrix_grease.rotateZ(90 * deg);
            position_grease.setY(0);
            position_grease.setZ(0);
            position_grease.setX(BoxY - greaseY);

            transform3D_grease = G4Transform3D(rotationMatrix_grease, position_grease);
            greaseAssembly->AddPlacedVolume(grease_logV, transform3D_grease);

            position_grease.setX(-1 * (BoxY - greaseY));

            transform3D_grease = G4Transform3D(rotationMatrix_grease, position_grease);
            greaseAssembly->AddPlacedVolume(grease_logV, transform3D_grease);

        }
    }

    //  place the assembly
    //
    position_grease.setX(0.0);
    position_grease.setY(0.0);
    position_grease.setZ(0.0);
    transform3D_grease = G4Transform3D(rotationMatrix_grease, position_grease);
    greaseAssembly->MakeImprint(assembly_logV, transform3D_grease);


    G4Box *qwSiPM_solidV = new G4Box("qwSiPM_solidV", qwSiPMx, qwSiPMy, qwSiPMz);

    qwSiPM_logV = new G4LogicalVolume(qwSiPM_solidV, fSiO2, "qwSiPM_logV", 0, 0, 0);
    G4VisAttributes *qwSiPM_VisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));//green
    //qwSiPM_VisAtt->SetForceSolid(true);
    qwSiPM_logV->SetVisAttributes(qwSiPM_VisAtt);

    // create assembly of SiPM window
    //
    G4AssemblyVolume *SiPMwindowAssembly = new G4AssemblyVolume();
    // Rotation and translation of the grease logical volume
    G4RotationMatrix rotationMatrix_SiPMwindow;
    G4ThreeVector position_SiPMwindow(0., 0., 0.);
    G4Transform3D transform3D_SiPMwindow;

    for (G4int i = 0; i < 2; i++) {
        if (i == 0) {
            position_SiPMwindow.setX(0);
            position_SiPMwindow.setZ(0);
            position_SiPMwindow.setY(BoxY - 2 * greaseY - qwSiPMy);

            transform3D_SiPMwindow = G4Transform3D(rotationMatrix_SiPMwindow, position_SiPMwindow);
            SiPMwindowAssembly->AddPlacedVolume(qwSiPM_logV, transform3D_SiPMwindow);

            position_SiPMwindow.setY(-1 * (BoxY - 2 * greaseY - qwSiPMy));

            transform3D_SiPMwindow = G4Transform3D(rotationMatrix_SiPMwindow, position_SiPMwindow);
            SiPMwindowAssembly->AddPlacedVolume(qwSiPM_logV, transform3D_SiPMwindow);
        } else {
            rotationMatrix_SiPMwindow.rotateZ(90 * deg);
            position_SiPMwindow.setY(0);
            position_SiPMwindow.setZ(0);
            position_SiPMwindow.setX(BoxY - 2 * greaseY - qwSiPMy);

            transform3D_SiPMwindow = G4Transform3D(rotationMatrix_SiPMwindow, position_SiPMwindow);
            SiPMwindowAssembly->AddPlacedVolume(qwSiPM_logV, transform3D_SiPMwindow);

            position_SiPMwindow.setX(-1 * (BoxY - 2 * greaseY - qwSiPMy));

            transform3D_SiPMwindow = G4Transform3D(rotationMatrix_SiPMwindow, position_SiPMwindow);
            SiPMwindowAssembly->AddPlacedVolume(qwSiPM_logV, transform3D_SiPMwindow);

        }
    }


    //  place the assembly
    //
    position_SiPMwindow.setX(0.0);
    position_SiPMwindow.setY(0.0);
    position_SiPMwindow.setZ(0.0);
    transform3D_SiPMwindow = G4Transform3D(rotationMatrix_SiPMwindow, position_SiPMwindow);
    SiPMwindowAssembly->MakeImprint(assembly_logV, transform3D_SiPMwindow);

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // the portion of SiPM just following the front quartz window is used for sensing the photons
    // the material of this portion or photon sensing part is made of silicon.
    //

    // define the optical surface for SiPM
    G4OpticalSurface *SiPMSurface = new G4OpticalSurface("SiPMSurface");
    SiPMSurface->SetType(dielectric_metal);
    SiPMSurface->SetFinish(polished);
    SiPMSurface->SetModel(glisur);
    SiPMSurface->SetMaterialPropertiesTable(fSilMPT);

    G4Box *psSiPM_solidV = new G4Box("psSiPM_solidV", psSiPMx, psSiPMy, psSiPMz);


    psSiPM_logV = new G4LogicalVolume(psSiPM_solidV, fSil, "psSiPM_logV", 0, 0, 0);
    G4VisAttributes *psSiPM_VisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));//
    psSiPM_VisAtt->SetForceSolid(true);
    psSiPM_logV->SetVisAttributes(psSiPM_VisAtt);

    //Add Sensitive Detector to SiPM
    G4LogicalSkinSurface *psSiPMSurface = new G4LogicalSkinSurface("SiPMSurface", psSiPM_logV, SiPMSurface);

    // create assembly of SiPM sensitive part
    //
    G4AssemblyVolume *psSiPMAssembly = new G4AssemblyVolume();
    // Rotation and translation of the grease logical volume
    G4RotationMatrix rotationMatrix_psSiPM;
    G4ThreeVector position_psSiPM(0., 0., 0.);
    G4Transform3D transform3D_psSiPM;

    for (G4int i = 0; i < 2; i++) {
        if (i == 0) {
            position_psSiPM.setX(0);
            position_psSiPM.setZ(0);
            position_psSiPM.setY(BoxY - 2 * greaseY - 2 * qwSiPMy - psSiPMy);

            transform3D_psSiPM = G4Transform3D(rotationMatrix_psSiPM, position_psSiPM);
            psSiPMAssembly->AddPlacedVolume(psSiPM_logV, transform3D_psSiPM);

            position_psSiPM.setY(-1 * (BoxY - 2 * greaseY - 2 * qwSiPMy - psSiPMy));

            transform3D_psSiPM = G4Transform3D(rotationMatrix_psSiPM, position_psSiPM);
            psSiPMAssembly->AddPlacedVolume(psSiPM_logV, transform3D_psSiPM);
        } else {
            rotationMatrix_psSiPM.rotateZ(90 * deg);
            position_psSiPM.setY(0);
            position_psSiPM.setZ(0);
            position_psSiPM.setX(BoxY - 2 * greaseY - 2 * qwSiPMy - psSiPMy);

            transform3D_psSiPM = G4Transform3D(rotationMatrix_psSiPM, position_psSiPM);
            psSiPMAssembly->AddPlacedVolume(psSiPM_logV, transform3D_psSiPM);

            position_psSiPM.setX(-1 * (BoxY - 2 * greaseY - 2 * qwSiPMy - psSiPMy));

            transform3D_psSiPM = G4Transform3D(rotationMatrix_psSiPM, position_psSiPM);
            psSiPMAssembly->AddPlacedVolume(psSiPM_logV, transform3D_psSiPM);

        }
    }

    //  place the assembly
    //
    position_psSiPM.setX(0.0);
    position_psSiPM.setY(0.0);
    position_psSiPM.setZ(0.0);
    transform3D_psSiPM = G4Transform3D(rotationMatrix_psSiPM, position_psSiPM);
    psSiPMAssembly->MakeImprint(assembly_logV, transform3D_psSiPM);

    return;
} // end of function. //

void nDetConstruction::buildArray() {

    greaseX = SiPM_dimension;
    greaseZ = SiPM_dimension;
    qwSiPMx = SiPM_dimension;
    qwSiPMz = SiPM_dimension;
    psSiPMx = SiPM_dimension;
    psSiPMz = SiPM_dimension;

    qwSiPMy = 0.37/2*mm;
    psSiPMy = 0.09/2*mm;

    fNdetectors=8;

    fCheckOverlaps =false;

    //fDetectorLength = 20*cm;
    fDetectorWidth = 10*cm;
    G4double gap = 0.1*mm;
    G4double array_length = fNdetectors*2*SiPM_dimension+2*(fNdetectors-1)*gap;

    G4double offset = 2*(fTeflonThickness)+4*(greaseY+qwSiPMy+psSiPMy);


    G4double wrapping_length = fDetectorLength+offset;
    G4double wrapping_width = fDetectorWidth+2*fTeflonThickness;
    G4double wrapping_width2 = array_length+2*fTeflonThickness;
    G4double wrapping_thickness = 2*psSiPMx +2*fTeflonThickness;

    //G4VSolid *theWrapping = ConstructNextModule("Wrapping",wrapping_length,wrapping_width,wrapping_width2,wrapping_thickness);

    G4Box * theWrapping = new G4Box("Wrapping",wrapping_width2/2,wrapping_thickness/2,wrapping_length/2);

    assembly_logV = new G4LogicalVolume (theWrapping,fTeflon,"theWrapping_log");

    fWrapSkinSurface = new G4LogicalSkinSurface("WrapSkin",assembly_logV,fTeflonOpticalSurface);

    G4RotationMatrix *rot=new G4RotationMatrix();
    //rot->rotateY(90*deg);
    G4ThreeVector wrapping_position(0,0,0);
    //For NEXT ARRAYS
    //rot->rotateZ(90*deg);
    //G4ThreeVector wrapping_position(0,0,wrapping_length/4);

    assembly_physV= new G4PVPlacement(rot,wrapping_position,assembly_logV,"theWrapping_phys",expHall_logV,false,0,fCheckOverlaps);

    //G4VSolid *TheScint = ConstructNextModule("Scint",fDetectorLength,fDetectorWidth,array_length,2*psSiPMx);

    G4Box * TheScint = new G4Box("Wrapping",array_length/2,fDetectorThickness/2,fDetectorLength/2);


    ej200_logV=new G4LogicalVolume(TheScint,fEJ200,"ej200_logV");

    G4VisAttributes* Array_VisAtt= new G4VisAttributes(G4Colour(0.0,1.0,1.0));//green
    ej200_logV->SetVisAttributes(Array_VisAtt);
    //Array_VisAtt->SetForceSolid(true);

    //G4ThreeVector scint_position(0,0,-offset/4);
    G4ThreeVector scint_position(0,0,0);

    G4VPhysicalVolume *scint_phys= new G4PVPlacement(0,scint_position,ej200_logV,"theScint_phys",assembly_logV,false,0,fCheckOverlaps);

    G4AssemblyVolume *theArray_log=ConstructArray("Array",8);


    //G4ThreeVector position(0,0,fDetectorLength/4+greaseY-offset/4);
    G4ThreeVector position(0,0,fDetectorLength/2+greaseY);

    G4RotationMatrix *rotation=new G4RotationMatrix();
    rotation->rotateX(90*deg);



    theArray_log->MakeImprint(assembly_logV,position,rotation,0,fCheckOverlaps);

    //position.setZ(-3*fDetectorLength/4-offset/4-greaseY);
    position.setZ(-fDetectorLength/2-greaseY);
    rotation->rotateX(-180*deg);

    theArray_log->MakeImprint(assembly_logV,position,rotation,0,fCheckOverlaps);

    //G4VPhysicalVolume *theArray_phys=new G4PVPlacement(rotation,position,theArray_log,"Array_ph",assembly_logV,1,0,1);
    //G4VPhysicalVolume *theArray_phys=new G4PVPlacement(rotation,position,theArray_log,"Array_ph",expHall_logV,1,0,1);


    return;
}


//G4LogicalVolume *nDetConstruction::ConstructArray(G4String name, G4int NDetectors) {
    G4AssemblyVolume *nDetConstruction::ConstructArray(G4String name, G4int NDetectors) {

    greaseX = SiPM_dimension;
    greaseZ = SiPM_dimension;
    qwSiPMx = SiPM_dimension;
    qwSiPMz = SiPM_dimension;
    psSiPMx = SiPM_dimension;
    psSiPMz = SiPM_dimension;

    qwSiPMy = 0.37/2*mm;
    psSiPMy = 0.09/2*mm;
    greaseY = 0.05*mm;

    fNdetectors=NDetectors;

    G4double gap = 0.1*mm;
    G4double BoxX = NDetectors*SiPM_dimension+2*(NDetectors-1)*gap;
    G4double BoxY = (greaseY+qwSiPMy+psSiPMy);
    G4double BoxZ = SiPM_dimension+fTeflonThickness;
    G4Box *TheArrayBox=new G4Box("Array",BoxX,BoxY,BoxZ);

    G4LogicalVolume *theArray_log = new G4LogicalVolume(TheArrayBox,fTeflon,"Array",0);


    G4Box *GreaseBox = new G4Box("Grease",greaseX,greaseY,greaseZ);
    G4Box *SiPMGlass = new G4Box("Glass",qwSiPMx,qwSiPMy,qwSiPMz);
    G4Box *SiPMChip = new G4Box("Chip",psSiPMx,psSiPMy,psSiPMz);
    grease_logV = new G4LogicalVolume(GreaseBox,fGrease,"Grease_log");
    G4VisAttributes *grease_VisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));//red
    grease_VisAtt->SetForceSolid(true);
    grease_logV->SetVisAttributes(grease_VisAtt);

    qwSiPM_logV = new G4LogicalVolume(SiPMGlass,fSiO2,"Window_Log");
    G4VisAttributes* qwSiPM_VisAtt= new G4VisAttributes(G4Colour(0.0,1.0,0.0));//green
    qwSiPM_VisAtt->SetForceSolid(true);
    qwSiPM_logV->SetVisAttributes(qwSiPM_VisAtt);

    // define the optical surface for SiPM
    G4OpticalSurface *SiPMSurface = new G4OpticalSurface("SiPMSurface");
    SiPMSurface->SetType(dielectric_metal);
    SiPMSurface->SetFinish(polished);
    SiPMSurface->SetModel(glisur);
    SiPMSurface->SetMaterialPropertiesTable(fSilMPT);



    psSiPM_logV = new G4LogicalVolume(SiPMChip,fSil,"Chip_log");

    G4LogicalSkinSurface *psSiPMSurface = new G4LogicalSkinSurface("SiPMSurface", psSiPM_logV, SiPMSurface);

    G4VisAttributes *psSiPM_VisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));//
    psSiPM_VisAtt->SetForceSolid(true);
    psSiPM_logV->SetVisAttributes(psSiPM_VisAtt);


    G4AssemblyVolume *theArray = new G4AssemblyVolume();

    for(G4int j=0;j<fNdetectors;j++) {
        G4ThreeVector position((-4+j)*2*(greaseX + gap)+greaseX+gap , 0, 0);
        G4RotationMatrix *rot = new G4RotationMatrix();
        //rot->rotateZ(60 * deg);
        theArray->AddPlacedVolume(grease_logV,position,rot);
        position.setY(greaseY+qwSiPMy);
        theArray->AddPlacedVolume(qwSiPM_logV,position,rot);
        position.setY(greaseY+2*qwSiPMy+psSiPMy);
        theArray->AddPlacedVolume(psSiPM_logV,position,rot);
    }

    G4RotationMatrix *rotation=new G4RotationMatrix();
    G4ThreeVector pos(0,0,0);
    //theArray->MakeImprint(theArray_log,pos,rotation,0,1);


    //return theArray_log;
    return theArray;
}



G4VSolid* nDetConstruction::ConstructNextModule(G4String name, G4double length, G4double width1,G4double width2, G4double thickness) {


    //We build the first trapezoid
    fTrapezoidLength=length/2.;

    G4double trapezoidlength=length/2.;

    G4double dx1=width1/2.;
    G4double dx2=width2/2.;
    G4double dy1=thickness/2.;
    G4double dy2=thickness/2.;
    G4double dz=trapezoidlength/2;

    G4Trd *theTrapezoid1=new G4Trd("theTrapezoid1"+name, dx1, dx2, dy1, dy2, dz);

    G4Trd *theTrapezoid2=new G4Trd("theTrapezoid2"+name, dx1, dx2, dy1, dy2, dz);

    G4ThreeVector translation(0,0,-2*dz);

    G4RotationMatrix *rot=new G4RotationMatrix();

    rot->rotateX(180*deg);

    G4UnionSolid *themodule=new G4UnionSolid("NEXT"+name,theTrapezoid1,theTrapezoid2,rot,translation);



    return themodule;

}
