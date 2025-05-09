#ifndef NDET_RUN_ACTION_HH
#define NDET_RUN_ACTION_HH

#include "globals.hh"
#include "G4UserRunAction.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"

#include "nDetConstruction.hh"
#include "nDetDataPack.hh"
#include "centerOfMass.hh"

#include <vector>

class G4Timer;
class G4Run;

class userActionManager;
class nDetEventAction;
class nDetStackingAction;
class nDetTrackingAction;
class nDetSteppingAction;
class nDetConstruction;
class pmtResponse;
class nDetDetectorParams;

class photonCounter;
class nDetParticleSource;

/** @class primaryTrackInfo
  * @brief Container for primary particle scatter information
  * @author Cory R. Thornsberry (cthornsb@vols.utk.edu)
  * @date July 5, 2019
  */

class primaryTrackInfo{
  public:
	G4ThreeVector pos; ///< 3d position of scatter (in mm)
	G4ThreeVector dir; ///< 3d momentum direction of primary particle prior to scattering
	G4double kE; ///< Kinetic energy of primary particle prior to scattering (in MeV)
	G4double dkE; ///< Change in kinetic energy of primary particle during scatter event (in MeV)
	G4double gtime; ///< Global time of scatter event (in ns)
	G4double plength; ///< Path length of primary particle between successive scatter events (in mm)
	G4double angle; ///< Angle between incoming and outgoing primary particle (in degrees)

	G4int copyNum; ///< Copy number of object inside which the scatter event occured
	G4int trackID; ///< Geant track ID of the primary particle track
	G4double atomicMass; ///< Atomic mass of the recoiling particle after the scatter event (NOT WORKING)
	
	G4bool inScint; ///< Flag indicating that the scatter event occured in a scintillator material

	const G4ParticleDefinition *part; ///< Pointer to the Geant particle definition of the primary particle

	/** G4Step constructor
	  */
	primaryTrackInfo(const G4Step *step);

	/** G4Track constructor
	  */
	primaryTrackInfo(const G4Track *track);

	/** Return true if this track and the track pointed to by @a rhs originate from the same position and return false otherwise
	  */
	bool compare(const G4Track *rhs);

	/** Return the angle between this track and the track pointed to by @a rhs and return -1 if the tracks do not originate at the same position
	  */
	double getAngle(const G4Track *rhs);

	/** Return the angle between the momentum direction of this track and vector @a rhs
	  */
	double getAngle(const G4ThreeVector &rhs);

	/** Get the path length between the position of this track and position @a rhs
	  */
	double getPathLength(const G4ThreeVector &rhs);

	/** Print all relevent information about this track to stdout
	  */
	void print();
	
  private:
	/** Copy relevant information from a track
	  */
	void setValues(const G4Track *track);
};

/** @class nDetRunAction
  * @brief Handles all processing of output variables and tracks primary particle scattering
  * @author Cory R. Thornsberry (cthornsb@vols.utk.edu)
  * @date July 5, 2019
  */

class nDetRunAction : public G4UserRunAction
{
  public:
	/** Default constructor
	  */
	nDetRunAction();

	/** Destructor
	  */
	virtual ~nDetRunAction();

	/** Action to perform at the beginning of each new run
	  */
	void BeginOfRunAction(const G4Run* aRun);

	/** Action to perform at the end of each run
	  */
	void EndOfRunAction(const G4Run* aRun);

	/** Initialize the pmtResponse objects in the center of mass (called in BeginOfRunAction())
	  */
	void initializePmtResponses();

	/** Get a pointer to the particle source object
	  */
	nDetParticleSource *getParticleSource(){ return source; }

	/** Get a pointer to the low level scatter info structure 
	  */
	nDetDebugStructure *getDebugData(){ return &debugData; }

	/** Set all user action pointers
	  * @param action Pointer to thread-local user action manager
	  */
	void setActions(userActionManager *action);

	/** Set the number of the current event
	  */
	void setEventNumber(const int &eventID){ evtData.eventID = eventID; }

	/** Enable or disable copying of light response traces to the output data structure
	  */
	void setOutputTraces(const bool &enabled){ outputTraces = enabled; }
	
	/** Enable or disable copying of debug variables to output debug data structure
	  */
	void setOutputDebug(const bool &enabled){ outputDebug = enabled; }

	/** Enable or disable copying of certain debug variables to multi detector data structure
	  */
	void setOutputMultiDebug(const bool &enabled){ outputMultiDebug = enabled; }


	/** Toggle the verbosity flag and return its state
	  */
	bool toggleVerboseMode(){ return (verbose = !verbose); }

	/** Segment the left and right PMTs of all defined detectors
	  * @param col_ The number of anode columns (horizontal) on the surface of the PMT
	  * @param row_ The number of anode rows (vertical) on the surface of the PMT
	  * @param width_ The width (horizontal) of the active area of the surface of the PMT
	  * @param height_ The height (vertical) of the active area of the surface of the PMT
	  */
    void setSegmentedPmt(const short &col_, const short &row_, const double &width_, const double &height_);

	/** Copy left and right PMT spectral response curves to all defined detectors
	  * @param left Pointer to the left PMT center-of-mass calculator from which to copy the PMT spectral response
	  * @param right Pointer to the right PMT center-of-mass calculator from which to copy the PMT spectral response
	  */
    void setPmtSpectralResponse(centerOfMass *left, centerOfMass *right);

	/** Copy left and right PMT gain matrices to all defined detectors
	  * @param left Pointer to the left PMT center-of-mass calculator from which to copy the PMT gain matrix
	  * @param right Pointer to the right PMT center-of-mass calculator from which to copy the PMT gain matrix
	  */
	void setPmtGainMatrix(centerOfMass *left, centerOfMass *right);

	/** Get the number of defined detectors
	  */
	size_t getNumDetectors() const { return userDetectors.size(); }

	/** Get a pointer to the left PMT dynode response object of one of the defined detectors
	  * @param index The index of the detector in the list of defined detectors
	  * @return The pointer to the dynode response if @a index corresponds to a defined detector and return NULL otherwise
	  */
	pmtResponse *getPmtResponseLeft(const size_t &index=0){ return (index < userDetectors.size() ? userDetectors.at(index).getCenterOfMassL()->getPmtResponse() : NULL); }

	/** Get a pointer to the right PMT dynode response object of one of the defined detectors
	  * @param index The index of the detector in the list of defined detectors
	  * @return The pointer to the dynode response if @a index corresponds to a defined detector and return NULL otherwise
	  */
	pmtResponse *getPmtResponseRight(const size_t &index=0){ return (index < userDetectors.size() ? userDetectors.at(index).getCenterOfMassR()->getPmtResponse() : NULL); }

		/** Get a pointer to the left PMT dynode response object of one of the defined detectors
	  * @param index The index of the detector in the list of defined detectors
	  * @return The pointer to the dynode response if @a index corresponds to a defined detector and return NULL otherwise
	  */
	pmtResponse *getPmtResponseImplant(const size_t &index=0){ return (index < userImplants.size() ? userImplants.at(index).getCenterOfMass()->getPmtResponse() : NULL); }

	/** Get vector of left PMT anode response objects of one of the defined detectors
	  * @param index The index of the detector in the list of defined detectors
	  * @return The vector of anode responses if @a index corresponds to a defined detector and return NULL otherwise
	  */
	std::vector<pmtResponse> *getAnodeResponseLeft(int eventID, const size_t &index=0){ 
		//return (index < userDetectors.size() ? userDetectors.at(index).getCenterOfMassL()->getAnodeResponse() : NULL); 
		return userDetectors.at(index).getCenterOfMassL()->getAnodeResponse();
	}

	/** Get vector of right PMT anode response objects of one of the defined detectors
	  * @param index The index of the detector in the list of defined detectors
	  * @return The vector of anode responses if @a index corresponds to a defined detector and return NULL otherwise
	  */
	std::vector<pmtResponse> *getAnodeResponseRight(int eventID, const size_t &index=0){
		//return (index < userDetectors.size() ? userDetectors.at(index).getCenterOfMassR()->getAnodeResponse() : NULL);
		return userDetectors.at(index).getCenterOfMassR()->getAnodeResponse();
	}	

	/** Get vector of left PMT anode response objects of one of the defined detectors
	  * @param index The index of the detector in the list of defined detectors
	  * @return The vector of anode responses if @a index corresponds to a defined detector and return NULL otherwise
	  */
	std::vector<pmtResponse> *getAnodeResponseImplant(int eventID, const size_t &index=0){
		//return (index < userImplants.size() ? userImplants.at(index).getCenterOfMass()->getAnodeResponse() : NULL);
		return userImplants.at(index).getCenterOfMass()->getAnodeResponse();
	}

	/** Get the total number of optical photons which have been produced during this run
	  */
    unsigned long long getNumPhotons() const { return numPhotonsTotal; }
 
	/** Get the total number of detected optical photons for this run
	  */   
    unsigned long long getNumPhotonsDet() const { return numPhotonsDetTotal; }

	/** Get the vector of detectors added by the user
	  */
	std::vector<nDetImplant> getUserImplants() {
		return userImplants;
	}

	/** Get the starting implant
	  */
	nDetImplant *getStartImplant() {
		return startImplant;
	}

	/** Get the detector parameters
	  */
	nDetDetectorParams getDetectorParameters() {
		return params;
	}

	/** Get a pointer to the singleton detector construction object
	  */
	nDetConstruction* getDetectorConstruction() {
		return detector;
	}

	/** Copy the list of defined detectors and set the start detector (if one exists)
	  * @param construction Pointer to the singleton detector constructor object
	  */
	void updateDetector(nDetConstruction *construction);

	/** Check if a scintillator copy number is part of a detector in the list of defined detectors
	  * @param num The copy number of a scintillator segment
	  * @return The parent copy number of the detector of which the segment is a member or return -1 if no match was found
	  */
	G4int checkCopyNumber(const G4int &num) const ;

	/** Get the column and row corresponding to a segment copy number in a segmented detector
	  * @param copyNum The copy number of a scintillator segment
	  * @param col The segmented detector column corresponding to the copy number
	  * @param row The segmented detector row corresponding to the copy number
	  * @return True if the copy number is found to be part of a detector and return false otherwise
	  */
	bool getSegmentFromCopyNum(const G4int &copyNum, G4int &col, G4int &row) const ;

	/** Process all detection events in all defined detectors, compute all output variables, fill output tree, and clear all data structures
	  */
	void process();

	/** Add a photon to the center-of-mass calculator corresponding to the detector inside of which the optical photon struck a sensitive surface
	  * @param step Pointer to the G4Step corresponding to the optical photon
	  * @param mass The mass of the detection event for the center-of-mass weighted average
	  * @return True if a matching detector was found and the optical photon was added successfully and return false otherwise
	  */
	bool AddDetectedPhoton(const G4Step *step, const double &mass=1);


	/** Add a detected secondary electron to the running total
	  * @param step Pointer to the G4Step corresponding to the optical photon
	  * @param energy The energy of the secondary electron (MeV)
	  * @return True if a matching detector was found and the secondary electron was added successfully and return false otherwise
	  */
	bool AddDetectedSecondary(const G4Step *step, const double &energy);


	/** Set initial primary particle scatter information with parameters from a G4Step
	  */
	void initializeNeutron(const G4Step *step);

	/** Set primary particle scatter inforamtion with parameters from a G4Step
	  * @return True if the change in primary particle kinetic energy is less than or equal to zero and return false otherwise
	  */
	bool scatterNeutron(const G4Step *step);

	/** Set final primary particle scatter information with parameters from a G4Step
	  */
	void finalizeNeutron(const G4Step *step);

	/** Copy the pmtResponse parameters specified in the macro from the responses contained in nDetConstruction
		to the ones contained in nDetImplant (need to do it for nDetDetector in the future as well).
	
	NOTE: this function was added because the centerOfMass and pmtResponse objects are cloned (using clone())
		before the response parameters are set. This is a very messy way of solving the problem, and it may 
		need to be changed in the future.
	
	* @param dest pointer to pmtResponse object whose parameters will be copied
	* @param src  pointer to pmtResponse object to which the parameters will be copied
	*/
	void copyResponseParameters(pmtResponse *dest, pmtResponse *src);


  private:
	G4Timer* timer; ///< Timer used to measure the total time for a run (in seconds)
	
	bool outputTraces; ///< Flag indicating that traces will be written to the output tree
	bool outputDebug; ///< Flag indicating that the user has requested low-level debug to be written to the output file
	bool outputMultiDebug; ///< Flag indicating reduced debug readout when using multiple detectors
	bool verbose; ///< Verbosity flag

	nDetEventAction *eventAction; ///< Pointer to the thread-local user event action
	nDetStackingAction *stacking; ///< Pointer to the thread-local user stacking action
	nDetTrackingAction *tracking; ///< Pointer to the thread-local user tracking action
	nDetSteppingAction *stepping; ///< Pointer to the thread-local user stepping action
	
	nDetConstruction *detector; ///< Pointer to the global detector singleton
	nDetParticleSource *source; ///< Pointer to the thread-local particle source
	nDetDetectorParams params; ///< not a Pointer to the parameters of the detector
	
	std::vector<primaryTrackInfo> primaryTracks; ///< Vector of primary particle scatter events
	
	G4ThreeVector prevDirection; ///< Momentum direction of the previous primary particle scatter
	G4ThreeVector prevPosition; ///< Position of the previous primary particle scatter

	nDetDataPack data; ///< Container object for all output data
	nDetEventStructure evtData; ///< Container object for output event information
	nDetOutputStructure outData; ///< Container object for single-detector output data
	nDetImplantOutputStructure outImplantData; ///< Container object for Implant output data
	nDetMultiOutputStructure multData; ///< Container object for multi-detector output data
	nDetDebugStructure debugData; ///< Container object for debug output data
	nDetTraceStructure traceData; ///< Container object for output traces

	photonCounter *counter; ///< Counter used to record the total number of optical photons produced by scattering

    unsigned long long numPhotonsTotal; ///< Total number of simulated optical photons (thread-local)
    unsigned long long numPhotonsDetTotal; ///< Total number of detected optical photons (thread-local)

	nDetDetector *startDetector; ///< Pointer to the detector used as a start signal for timing
	nDetImplant *startImplant; ///< Pointer to the detector used as a start signal for timing
	bool endImplant=0;
	bool endDetector=0;

	long numberOfEvents; ///< Number of events to be processed

	std::vector<nDetDetector> userDetectors; ///< Vector of detectors added by the user
	std::vector<nDetImplant> userImplants; ///< Vector of detectors added by the user

	/** Pop a primary scatter off the stack. Set all initial event conditions if this is the first scatter
	  * @return True if the stack of primary scatters is not empty after popping off a scatter and return false otherwise
	  */
	bool scatterEvent();
	
	/** Process a single event for a single detector
	  * @param det Pointer to the nDetDetector to process
	  * @return True if the detector has detected optical photons and return false otherwise
	  */
	bool processDetector(nDetDetector* det);

	/** Process a single event for a single detector
	  * @param det Pointer to the nDetDetector to process
	  * @return True if the detector has detected optical photons and return false otherwise
	  */
	bool processImplant(nDetImplant* imp);
	
	/** Process a single event for a single start detector
	  * @param det Pointer to the nDetDetector to process
	  * @param startTime Time-of-flight of the start detector event
	  * @return True if the detector has detected optical photons and return false otherwise
	  */
	bool processStartDetector(nDetDetector* det, double &startTime);

		/** Process a single event for a single start detector
	  * @param imp Pointer to the nDetImplant to process
	  * @param startTime Time-of-flight of the start detector event
	  * @return True if the detector has detected optical photons and return false otherwise
	  */
	bool processStartImplant(nDetImplant* imp, double &startTime);

};


#endif
