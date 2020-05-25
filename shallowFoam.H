#ifndef PRECICEADAPTER_H
#define PRECICEADAPTER_H
#include "INterfaceSF.H"

// the respective libraries in the Make/options.

// OpenFOAM header files
#include "Time.H"
#include "fvMesh.H"

// preCICE Solver Interface
#include "precice/SolverInterface.hpp"

namespace preciceAdapter
{

class Adapter
{

private:

    //- Structure of the configuration of each coupling interface.
    //  Every interface needs to know the coupling mesh, the OpenFOAM
    //  patches that consist the coupling surface and the kinds
    //  of data that are exchanged.
    struct InterfaceConfig
    {
        std::string meshName;
        std::string locationsType;
        bool meshConnectivity;
        std::vector<std::string> patchNames;
        std::vector<std::string> writeData;
        std::vector<std::string> readData;
    };

    //- Configuration interfaces
    std::vector<struct InterfaceConfig> interfacesConfig_;

    //- OpenFOAM runTime object
    const Foam::Time& runTime_;

    //- OpenFOAM fvMesh object
    const Foam::fvMesh& mesh_;

    // Configuration parameters used in the Adapter

        //- Remember if there were errors in the read() method
        bool errorsInConfigure = false;

        //- preCICE configuration file name
        Foam::word preciceConfigFilename_;

        //- preCICE participant name
        Foam::word participantName_;

        //- Read the adapter's configuration file
        bool configFileRead();

};

}

#endif

