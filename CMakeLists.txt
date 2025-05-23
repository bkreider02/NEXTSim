cmake_minimum_required (VERSION 2.8)

project (nextSim)

#Compile with c++11 or c++0x depending on the compiler version
if (CMAKE_COMPILER_IS_GNUCXX)
	if (${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "4.3" AND
            ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "4.7")
            #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x") # CHANGE IS HERE
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
    elseif(${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "4.6")
        #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") # CHANGE IS HERE
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    endif(${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "4.3" AND
            ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "4.7")

    #Make compiler messages nice with colored tags.
	if (${CMAKE_CXX_COMPILER_VERSION} VERSION_GREATER "4.9")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=auto")
	endif()

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -fPIC")
endif()

#if user does not specify prefix we assign it to the exec directory
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
   message(STATUS "Install Prefix not specified.")
   file(MAKE_DIRECTORY install)
   get_filename_component(INSTALL_DIR ${CMAKE_SOURCE_DIR}/install REALPATH)
   set(CMAKE_INSTALL_PREFIX ${INSTALL_DIR} CACHE PATH "Install Prefix" FORCE)
endif()
message(STATUS "Installing to ${CMAKE_INSTALL_PREFIX}")

#Define the default build type to be Release
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Build type, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
   message(STATUS "Build type not defined, using default: ${CMAKE_BUILD_TYPE}")
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

#Add additional CXX flags if we use the Debug option
if(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif(CMAKE_BUILD_TYPE MATCHES "Debug")

#------------------------------------------------------------------------------

#Install options
#option(BUILD_SHARED "Build and install shared libraries." OFF)
option(BUILD_TOOLS "Build and install tool programs." OFF)
option(GEANT4_MT "Enable multi-threading support (if available)." OFF)
option(GEANT4_GDML "Enable support for G4GDMLParser." OFF)
option(GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)

#------------------------------------------------------------------------------

#Set the name of the dictionary
set(DICTIONARY_NAME "NextSimEvent" CACHE STRING "Prefix to root dictionary.")

#Find required packages.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

#Find ROOT install.
find_package (ROOT REQUIRED)
mark_as_advanced(FORCE GENREFLEX_EXECUTABLE ROOTCINT_EXECUTABLE ROOT_CONFIG_EXECUTABLE)
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

#Find CADMesh
find_package (cadmesh)
message(${CADMESH_INCLUDE_DIRS})
include_directories(${CADMESH_INCLUDE_DIRS})
#link_directories(${cadmesh_LIBRARIES})

#Inform the user about the Geant4 install version
message(STATUS "Found ROOT install: v${ROOT_VERSION}")

#Find Geant install.
if(GEANT4_UIVIS)
	find_package(Geant4 REQUIRED ui_all vis_all)
else()
	find_package(Geant4 REQUIRED)
endif()

#Inform the user about the Geant4 install version
message(STATUS "Found Geant4 install: v${Geant4_VERSION}")

set(TOP_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

#Add the include directories.
include(${Geant4_USE_FILE})
include_directories(include)
include_directories(olaInclude)
include_directories(dict/include)
include_directories(/opt/CADmesh/geant4.10.3.3/include)

#Add the root dictionary directory.
add_subdirectory(dict)

#Add the source directories.
add_subdirectory(source)

if(BUILD_TOOLS)
	option(BUILD_TOOLS_NISTLIST "Build NIST library list program." OFF)
	option(BUILD_TOOLS_CONVERTER "Build nextSim to simpleScan output converter." OFF)
	option(BUILD_TOOLS_CMDSEARCH "Build nextSim macro command search program." OFF)
	option(BUILD_TOOLS_MACROREADER "Build nextSim macro file generator program." OFF)
	add_subdirectory(tools)
endif()

if(BUILD_TOOLS_CONVERTER)
	#Find simpleScan install.
	find_package (SimpleScan REQUIRED)
	include_directories(${SimpleScan_INCLUDE_DIR})
	mark_as_advanced(SimpleScan_SCAN_LIB SimpleScan_OPT_LIB SimpleScan_DICTIONARY)
endif(BUILD_TOOLS_CONVERTER)

#Build/install the miscellaneous stuff
add_subdirectory(share)