/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    shallowFoam

Description
    Transient solver for depth averaged shallow water equations.

Authors
    KM-Turbulenz GmbH, 2009 (www.km-turbulenz.de)
    Florian Mintgen, 2012

\*---------------------------------------------------------------------------*/

#include <iostream>
#include <sstream>

#include "precice/SolverInterface.hpp"

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "createTimeControls.H"
#   include "CourantNo.H"
#   include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// This should be changed accordingly
char solverName[] = "rightSF";
char configFileName[] =	"../precice-config.xml";
int rank = 0;
int size = 1;
double meshHeight = 4;
double meshWidth = 0.1;

char meshName[] = "right-Mesh";
int nvertex1 = 1;

char meshName2[] = "left-Mesh";
int nvertex2 = 40;

precice::SolverInterface precice(solverName,configFileName,rank,size);
int dim = precice.getDimensions();

// Make data structures shallowFoam (1)
int meshID = precice.getMeshID(meshName);
int vertexSize = nvertex1;				// Set number of vertices at wet surface
double* coords = new double[vertexSize*dim]; 
int* vertexIDs = new int[vertexSize];
precice.setMeshVertices(meshID, vertexSize, coords, vertexIDs); 
delete[] coords;

int flowdID = precice.getDataID("FlowDepth", meshID); 	
int dischID = precice.getDataID("Discharge", meshID); 
int alphaID = precice.getDataID("Alpha", meshID); 
int velocID = precice.getDataID("Velocity", meshID); 
int prghID = precice.getDataID("Prgh", meshID);

double* flowdepth = new double[vertexSize*dim];
double* discharge = new double[vertexSize*dim];
double* alphaw = new double[vertexSize*dim];
double* velocity = new double[vertexSize*dim];
double* prgh = new double[vertexSize*dim];
// End make data structures shallowFoam

// Make data structures interFoam (2)
int meshID2 = precice.getMeshID(meshName2);
int vertexSize2 = nvertex2;				// Set number of vertices at wet surface
double* coords2 = new double[vertexSize2*dim]; 
int* vertexIDs2 = new int[vertexSize2];
delete[] coords2;

int flowdID2 = precice.getDataID("FlowDepth", meshID2); 	
int dischID2 = precice.getDataID("Discharge", meshID2); 
int alphaID2 = precice.getDataID("Alpha", meshID2); 
int velocID2 = precice.getDataID("Velocity", meshID2); 
int prghID2 = precice.getDataID("Prgh", meshID2);

double* flowdepth2 = new double[vertexSize2*dim];
double* discharge2 = new double[vertexSize2*dim];
double* alphaw2 = new double[vertexSize2*dim];
double* velocity2 = new double[vertexSize2*dim];
double* prgh2 = new double[vertexSize2*dim];
// End make data structures interFoam


// Defining Time step sizes
double dt; 		// solver timestep size	runTime.setDeltaT  or  runTime.deltaTValue()
double precice_dt; 	// maximum precice timestep size

// Starting Time Loop
    Info<< "\nStarting time loop\n" << endl;
    precice_dt = precice.initialize(); 

    precice.initializeData();

    while (precice.isCouplingOngoing())
    {
	// Data handling and conversion
	precice.readBlockScalarData(alphaID, vertexSize, vertexIDs, alphaw);
	precice.readBlockVectorData(velocID, vertexSize, vertexIDs, velocity);
	precice.readBlockScalarData(prghID, vertexSize, vertexIDs, prgh);

	// Alpha to H conversion
	int i;
	double H_left = 0;
	for ( i = 0; i<vertexSize ; i++ )
	{
	    H_left += alphaw[i] * meshHeight / vertexSize;
	}
		
        Info<< "H_left = " << H_left << nl << endl;

	// U to HU conversion
	double HUx_left = 0;
	double HUy_left = 0;
	for ( i = 0; i<vertexSize ; i++ )
	{
	    HUx_left += alpha[i] * velocity[i] / (meshWidth * meshHeight / vertexSize);
	    HUy_left += alpha[i] * velocity[i] / (meshWidth * meshHeight / vertexSize);
	}
		
        Info<< "H_left = " << H_left << nl << endl;


        #include "setDeltaT.H"
	dt = runTime.deltaTValue();
	dt = min(precice_dt, dt);
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

///////////////////////////////////////////////////////////////////////////////
///     Transport equations for H und HU
	fvScalarMatrix HEqn
        (
	    fvm::ddt(H)
	    + fvm::div(phi,H)
        );

        HEqn.solve();

        // Hclip = max (H, Hdry)
        Hclip = (H-Hdry)*pos(H-Hdry) + Hdry;

	// Bottom friction based on Strickler-equation
        alpha = (g*dim_s*dim_s/dim_m)*mag(U)/pow(kst,2.0)/pow(Hclip/dim_m,1.0/3.0)/Hclip;

        // nu_t = C_nu * u_tau * H
        nut = Cnu*sqrt(alpha*mag(HU))*H;
        nut = nut*pos(nutmax-nut) + nutmax*pos(nut-nutmax);
        nut.correctBoundaryConditions();

	// Zero-Gradient for faces at wet-dry interface
	#include "zeroGradientWetDry.H"

         fvVectorMatrix HUEqn
         (
             fvm::ddt(HU)          // d(HU_i)/dt                    Local derivative
	     + fvm::div(phi, HU)   // + d/dx_j ( U_j * HU_i )       Convection
	     + H*gradgSpHf          // + H * d/dx_i ( gSpH )        Downhill-slope force and acceleration due to change in flowdepth
	     + fvm::Sp(alpha,HU)   // + alpha * HU                  Bottom friction
	     - fvm::laplacian(nut,HU) // - d^2/dx_j^2 ( nut * HU )  Turbulent stresses
         );  

	HUEqn.solve();

///////////////////////////////////////////////////////////////////////////////
///     Calculation of U and phi

        // if Hclip < Hdry2, then U = 0
        HU = HU*pos(Hclip - Hdry2); 
        U = HU/Hclip;
       
        phi = (fvc::interpolate(U) & mesh.Sf());

///////////////////////////////////////////////////////////////////////////////

	//precice.writeBlockScalarData(flowdID, vertexSize, vertexIDs, flowdepth);
	//precice.writeBlockVectorData(dischID, vertexSize, vertexIDs, discharge);

	precice_dt = precice.advance(dt);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    precice.finalize(); // frees data structures and closes communication channels
// Ending Time Loop


    return(0);
}


// ************************************************************************* //
