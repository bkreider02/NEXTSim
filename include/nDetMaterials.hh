#ifndef NDET_MATERIALS_HH
#define NDET_MATERIALS_HH

#include "globals.hh"

#include "nistDatabase.hh"

class G4Element;
class G4Material;
class G4MaterialPropertiesTable;
class G4OpticalSurface;
class G4VisAttributes;

class nDetMaterialsMessenger;

/** @class nDetMaterials
  * @brief Geant materials and elements used for detector construction
  * @author Cory R. Thornsberry (cthornsb@vols.utk.edu)
  * @date July 3, 2019
  *
  * This class is designed as a handler for common Geant materials which are used by
  * other NEXTSim classes. The materials and elements defined in this class are expandable
  * by using the nistDatabase class.
  */

class nDetMaterials{
  public:
    // Materials and elements
    G4Element* fH; ///< Hydrogen
    G4Element* fC; ///< Carbon
    G4Element* fO; ///< Oxygen
	G4Element* fN; ///< Nitrogen
	G4Element* fS; ///< Sulfur
    G4Element* fF; ///< Flourine
    G4Element* fSi; ///< Silicon
    G4Element* fAl; ///< Aluminum
 		G4Element* fY; ///< Yttrium
		G4Element* fLa; ///< Lanthinum
		G4Element* fBr; ///< Bromine
		//G4Element* fCs; ///< Cesium
		//G4Element* fI;  ///< Iodine
		G4Element* fGd; ///< Gadolinium
		G4Element* fGa; ///< Gallium
		G4Element* fCe; ///< Cerium
		G4Element* fZr; ///< Zirconium
		G4Element* fLu; ///< Lutetium

    G4Material* fAir; ///< Material corresponding to air
    G4Material* fVacuum; ///< Material corresponding to natural vacuum
    G4Material* fTeflon; ///< Material corresponding to teflon
	G4Material* fAl2O3; ///< Material corresponding to aluminum oxide
    G4Material* fEJ200; ///< Material corresponding to EJ-200 scintillator
    G4Material* fEJ276; ///< Material corresponding to EJ-276 scintillator
	G4Material* fOGS; ///< Material corresponding to organic glass scintillator (OGS)
    G4Material* fGrease; ///< Material corresponding to optical grease
	G4Material* fNOA61; ///< Material corresponding to NOA-61 (Norland optical adhesive)
	G4Material* fNOA68; ///< Material corresponding to NOA-68 (Norland optical adhesive)
	G4Material* fNOA136; ///< Material corresponding to NOA-136 (Norland optical adhesive)
	G4Material* fNOA170; ///< Material corresponding to NOA-170 (Norland optical adhesive)	
	G4Material* fABS; ///< Material corresponding to ABS plastic
    G4Material* fSiO2; ///< Material corresponding to quartz
    G4Material* fSilicon; ///< Material corresponding to silicon
    G4Material* fMylar; ///< Material corresponding to aluminized mylar
    G4Material* fAcrylic; ///< Material corresponding to acrylic
		G4Material* fAluminum; ///< Material corresponding to aluminum
 		G4Material* fYSO; ///< Material corresponding to yttrium orthosilicate
		G4Material* fLYSO; ///< Material corresponding to lutetium-yttrium oxyorthosilicate
		G4Material* fYAP; ///< yttrium aluminum perovskite
		G4Material* fGAGG; ///< Gadolinium Aluminum Gallium Garnet (Gd3Al2Ga3O12)
		G4Material* fLaBr3; ///< Lanthinum Bromide 
		G4Material* fCeBr3; ///< Cerium Bromide
		//G4Material* fCsI; ///< Cesium Iodide


    // Material table properties
    G4MaterialPropertiesTable* fAirMPT; ///< Material properties table for air
    G4MaterialPropertiesTable* fTeflonMPT; ///< Material properties table for teflon
	G4MaterialPropertiesTable* fAl2O3MPT; ///< Material properties table for aluminum oxide
    G4MaterialPropertiesTable* fEJ200MPT; ///< Material properties table for EJ-200 scintillator
    G4MaterialPropertiesTable* fEJ276MPT; ///< Material properties table for EJ-276 scintillator
	G4MaterialPropertiesTable* fOGSMPT; ///< Material properties table for OGS scintillator
		G4MaterialPropertiesTable* fYSOMPT; ///< Material properties table for YSO scintillator
		G4MaterialPropertiesTable* fLYSOMPT; ///< Material properties table for LYSO scintillator
		G4MaterialPropertiesTable* fYAPMPT; ///< MPT for YAP scintillator
		G4MaterialPropertiesTable* fGAGGMPT; ///< MPT for GAGG scintillator
		G4MaterialPropertiesTable* fLaBr3MPT; ///< MPT for LaBr3 scintillator
		G4MaterialPropertiesTable* fCeBr3MPT; ///< MPT for CeBr3 scintillator
		//G4MaterialPropertiesTable* fCsIMPT; ///< MPT for CsI scintillator
    G4MaterialPropertiesTable* fGreaseMPT; ///< Material properties table for optical grease
	//G4MaterialPropertiesTable* fMercaptoMPT; ///< Material properties table for Mercapto-ester
	//G4MaterialPropertiesTable* fTriallylMPT; ///< Material properties table for Triallyl Isocyanurate
	//G4MaterialPropertiesTable* fTetrahydroffMPT; ///< Material properties table for Tetrahydrofurfuryl
	//G4MaterialPropertiesTable* fUrethaneAcrylateMPT; ///< Material properties table for aliphatic urethane acrylate
	//G4MaterialPropertiesTable* fZrO2MPT; ///< Material properties table for zirconium dioxide
	//G4MaterialPropertiesTable* fAcrylateMPT; ///< Material properties table for acrylate	
    G4MaterialPropertiesTable* fNOA61MPT; ///< Material properties table for NOA-61
    G4MaterialPropertiesTable* fNOA68MPT; ///< Material properties table for NOA-68
    G4MaterialPropertiesTable* fNOA136MPT; ///< Material properties table for NOA-136
    G4MaterialPropertiesTable* fNOA170MPT; ///< Material properties table for NOA-170
	G4MaterialPropertiesTable* fABSMPT; /// Material properties table for ABS plastic
    G4MaterialPropertiesTable* fSiO2MPT; ///< Material properties table for quartz
    G4MaterialPropertiesTable* fSiliconMPT; ///< Material properties table for silicon
		G4MaterialPropertiesTable* fMylarMPT; ///< Material properties table for aluminized mylar
		G4MaterialPropertiesTable* fPerfectMPT; ///< Material properties table for a perfect reflector
		G4MaterialPropertiesTable* fAluminumMPT; ///< Material properties table for aluminum
		G4MaterialPropertiesTable* fEsrMPT; ///< Material properties table for 3M Enhanced Specular Reflector

    // Optical Surfaces
    G4OpticalSurface* fTeflonOpSurf; ///< Optical surface for teflon
	G4OpticalSurface* fAl2O3OpSurf; ///< Optical surface for aluminum oxide
    G4OpticalSurface* fSiliconOpSurf; ///< Optical surface for silicon
    G4OpticalSurface* fMylarOpSurf; ///< Optical surface for aluminized mylar
    G4OpticalSurface* fEsrOpSurf; ///< Optical surface for 3M Enhanced Specular Reflector
		G4OpticalSurface* fPerfectOpSurf; ///< Optical surface for a perfect reflector
		G4OpticalSurface* fGreaseOpSurf; ///< Optical surface for optical grease
		G4OpticalSurface* fNOA61OpSurf; ///< Optical surface for NOA-61
		G4OpticalSurface* fNOA68OpSurf; ///< Optical surface for NOA-68
		G4OpticalSurface* fNOA136OpSurf; ///< Optical surface for NOA-136
		G4OpticalSurface* fNOA170OpSurf; ///< Optical surface for NOA-170
		G4OpticalSurface* fAirOpSurf; ///< Optical surface for air
		G4OpticalSurface* fAluminumOpSurf; ///< Optical surface for aluminum

		// Visual attributes
		G4VisAttributes *visAssembly; ///< Visual attributes for the mother assembly
		G4VisAttributes *visSensitive; ///< Visual attributes for the photo-sensitive surface
		G4VisAttributes *visWindow; ///< Visual attributes for the optical window
		G4VisAttributes *visGrease; ///< Visual attributes for the optical grease	
		G4VisAttributes *visWrapping; ///< Visual attributes for the inner/outer detector wrappint
		G4VisAttributes *visScint; ///< Visual attributes for the scintillator material
		G4VisAttributes *visShadow; ///< Visual attributes for the shadow object
	
	// name of detector material
	std::string detectorMaterialName = "";

	/** Default constructor
	  */
	nDetMaterials();
	
	/** Destructor
	  */
	~nDetMaterials();

	/** Define all materials, scintillators, and optical surfaces
	  */
	void initialize();

	/** Return true if all materials and surfaces are defined and return false otherwise
	  */
	bool materialsAreDefined() const { return isInitialized; }

	/** Get a pointer to the scintillator material corresponding to @a name
	  * 
	  * Valid material names are shown in the table below
	  * | Name  | Description |
	  * |-------|-------------|
	  * | ej200 | <a href="https://eljentechnology.com/products/plastic-scintillators/ej-200-ej-204-ej-208-ej-212">Eljen EJ-200</a> plastic scintillator
	  * | ej276 | <a href="https://eljentechnology.com/products/plastic-scintillators/ej-276">Eljen EJ-276</a> plastic scintillator
	  *
	  * @return A pointer to the scintillator material if the name is valid or a pointer to EJ200 if it is invalid
	  */
	G4Material* getUserDetectorMaterial(const G4String &name);

	/** Get a pointer to the surface material corresponding to @a name
	  *
	  * Valid material names are shown in the table below
	  * | Name    | Description |
	  * |---------|-------------|
	  * | mylar   | Aluminized mylar (80% mylar, 20% aluminum by mass)
	  * | teflon  | Polytetrafluoroethylene (C2F4)n (PTFE)
	  * | esr     | <a href="https://www.3m.com/3M/en_US/company-us/all-3m-products/~/3M-Enhanced-Specular-Reflector-3M-ESR-/?N=5002385+3293061534&rt=rud">3M Enhanced Specular Reflector</a>
	  * | silicon | Natural silicon (from the NIST database)
	  * | perfect | Perfect optical reflector (reflectivity = 1, specular spike = 1)
	  * 
	  * @return A pointer to the surface material if the name is valid or a pointer to Mylar if it is invalid
	  */
	G4Material* getUserSurfaceMaterial(const G4String &name);

	/** Get a pointer to the optical surface corresponding to @a name
	  * @note See getUserSurfaceMaterial() for a description of valid optical surface names
	  * @return A pointer to the optical surface if the name is valid or a pointer to Mylar optical surface if it is invalid	
	  */
	G4OpticalSurface* getUserOpticalSurface(const G4String &name);

	/** Get a pointer to the visual attributes corresponding to @a name
	  */
	G4VisAttributes* getUserVisAttributes(const G4String &name);

	/** Get the current light yield multiplier used for scintillator photon production
	  */
	G4double getLightYield() const { return lightYieldScale; }

	/** Get a pointer to the database for pre-defined NIST element and materials
	  */
	nistDatabase* getNistDatabase(){ return &nist; }

	/** Get a pointer to the messenger for this class
	  */
	nDetMaterialsMessenger* getMessenger(){ return messenger; }

	/** Set the light yield multiplier for scintillator photon production
	  * @param yield The fraction of the light yield to use for optical photon production in scintillators (default is 1)
	  */	
	void setLightYield(const G4double &yield);

	/** Search for an element in the NIST database 
	  * @return True if the element is found in the NIST database and return false otherwise
	  */
	bool searchForElement(const G4String &name);

	/** Search for an element name in the pre-defined dictionary
	  * @note If the specified name corresponds to a pre-defined element, a pointer to that element is returned.
	  *       Otherwise the name is searched for within the Geant NIST element database.
	  * @return A pointer to the Geant element if it is found in the pre-defined dictionary OR the NIST database and return NULL otherwise
	  */
	G4Element* getElement(const G4String &name);

	/** Search for a material in the NIST database 
	  * @return True if the material is found in the NIST database and return false otherwise
	  */
	bool searchForMaterial(const G4String &name);
	
	/** Search for a material name in the pre-defined dictionary
	  * @note If the specified name corresponds to a pre-defined material, a pointer to that material is returned.
	  *       Otherwise the name is searched for within the Geant NIST material database.
	  * @return A pointer to the Geant material if it is found in the pre-defined dictionary OR the NIST database and return NULL otherwise
	  */
	G4Material* getMaterial(const G4String &name);

	/** Search for an optical surface name in the pre-defined dictionary
	  * @return A pointer to the Geant optical surface if it is found in the pre-defined dictionary and return NULL otherwise
	  */
	G4OpticalSurface* getOpticalSurface(const G4String &name);

	/** Search for an object visual attributes name in the pre-defined dictionary
	  * @return A pointer to the Geant visual attributes if it is found in the pre-defined dictionary and return NULL otherwise
	  */	
	G4VisAttributes* getVisualAttributes(const G4String &name);

	/** List all materials in the materials dictionary
	  */
	void listMaterials() const ;

	/** List all visual attributes in the dictionary
	  */
	void listVisAttributes() const ;

	/** List all optical surfaces in the dictionary
	  */
	void listOptSurfaces() const ;

	/** Calculate Average Scintillation Efficiency with AVA/Birks formula with offset
	 */
	G4double dLdE(G4double dEdx, G4double a, G4double b, G4double c);

	/** Define light yield for ion source if not already defined. Interpolation is performed if necessary
	  */
	void setIonSourceType(G4int mass);

	/** Read an input material file and build a new material to add to the list
	  * @param filename The path to the material file to read
	  * @return True if the file is read successfully and return false otherwise
	  */
	bool buildNewMaterial(const G4String &filename);

	/**
	  */
	void printMaterial(const G4String &name);

  private:
	nDetMaterialsMessenger* messenger; ///< Pointer to the messenger for this class

	bool isInitialized; ///< Flag indicating that all materials have been defined
	bool scintsAreDefined; ///< Flag indicating that the scintillator materials have been defined
	
    G4double lightYieldScale; ///< Scaling parameter used for scintillation light yield (default=1)	

	nistDatabase nist; ///< Database for pre-defined NIST element and material lookups

	std::map<std::string,std::vector< std::vector<G4double> > > lightYieldTable; ///< table of light yields; used to interpolate light for undefined ions

	std::map<G4String, G4Element*> elementList; ///< List of elements which are available for use
	
	std::map<G4String, G4Material*> materialList; ///< List of materials which are available for use
	
	std::map<G4String, G4OpticalSurface*> opticalSurfaceList; ///< List of optical surfaces which are available for use
	
	std::map<G4String, G4VisAttributes*> visAttributesList; ///< List of visual attributes which are available for use

	/** Define materials and optical surfaces
	  */
	void defineMaterials();

	/** Define scintillator materials EJ200 and EJ276
	  */
	void defineScintillators();
};

#endif
