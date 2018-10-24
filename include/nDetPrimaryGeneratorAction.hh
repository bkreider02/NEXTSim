//
// $Id: nDetPrimaryGeneratorAction.hh,v 1.0 Sept., 2015 $
// Written by Dr. Xiaodong Zhang
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef nDetPrimaryGeneratorAction_h
#define nDetPrimaryGeneratorAction_h 1

#include "nDetRunAction.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

class G4ParticleGun;
class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class nDetPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    nDetPrimaryGeneratorAction() : particleGun(NULL), runAct(NULL) { }
  
    nDetPrimaryGeneratorAction(nDetRunAction*);
    
   ~nDetPrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);

	G4ParticleGun *GetParticleGun(){ return particleGun; }
	
	nDetRunAction *GetRunAction(){ return runAct; }
	
	void SetParticleGun(G4ParticleGun *gun){ particleGun = gun; }
	
	void SetRunAction(nDetRunAction *action){ runAct = action; }

  protected:
    G4ParticleGun *particleGun;
    
    nDetRunAction* runAct;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*nDetPrimaryGeneratorAction_h*/

// Class description of the G4VUserPrimaryGeneratorAction:
// 
// This is the abstract base class of the user's mandatory action class
// for primary vertex/particle generation. This class has only one pure 
// virtual method GeneratPrimaries() which is invoked from G4RunManager
// during the event loop.
// Note that this class is not intended for genrating primary vertex/particle
// by itself. This class should
//	- have one or more G4VUserPrimaryGenerator concrete classes such as 
//	  G4PaticleGun
//	- set/change properties of generator(s)
//	- pass G4Event object so that the generator(s) can generate primaries.
//
