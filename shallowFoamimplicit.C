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
#include <cmath>
#include <math.h> 

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
char meshName[] = "right-Mesh";
int nFaces = 40;
int nVertex = 2 * nFaces + 2;
double faceHeight = (meshHeight / nFaces);
double rhow3D = 1000;	// Same as fluid density in other solver
double rhoa3D = 1;	// Same as air density in other solver
double g3D = 9.81;	// Same as gravitational acceleration in other solver

precice::SolverInterface precice(solverName,configFileName,rank,size);

const std::string& coric = precice::constants::actionReadIterationCheckpoint(); 
const std::string& cowic = precice::constants::actionWriteIterationCheckpoint();

int dim = precice.getDimensions();


// Make data structures for shallowFoam
int meshID = precice.getMeshID(meshName);
int vertexSize = nVertex;				// Set number of vertices at wet surface
double* coords = new double[vertexSize*dim]; 

for ( int i = 0; i<vertexSize ; i++ )			// Hardcoding coords of coupling interface
{
    if ( i < vertexSize / 2 )
    {
	coords[0+3*i] = 4.0; 		// x
	coords[1+3*i] = 0.0+0.1*i; 	// y
	coords[2+3*i] = 0.0; 		// z
    }
    else
    {
	coords[0+3*i] = 4.0; 		// x
	coords[1+3*i] = 0.0+0.1*i; 	// y
	coords[2+3*i] = 1.0; 		// z
    }
}

int* vertexIDs = new int[vertexSize];
precice.setMeshVertices(meshID, vertexSize, coords, vertexIDs); 
 
int alphaID = precice.getDataID("Alpha", meshID); 
int velocID = precice.getDataID("Velocity", meshID); 
int prghID = precice.getDataID("Prgh", meshID);

double* alphaw = new double[vertexSize];
double* velocity = new double[vertexSize*dim];
double* prgh = new double[vertexSize];
// End make data structures for shallowFoam


// Make data structure for checkpoints
double* alphawcheck = alphaw;
double Hcheck;
// End make data structure for checkpoints


// Defining Time step sizes
double dt; 		// solver timestep size	runTime.setDeltaT  or  runTime.deltaTValue()
double precice_dt; 	// maximum precice timestep size

// Starting Time Loop
    Info<< "\nStarting time loop\n" << endl;
    precice_dt = precice.initialize(); 

    precice.initializeData();

    while (precice.isCouplingOngoing())
    {
	if(precice.isActionRequired(cowic))
	{
	    // save initial data
	    alphawcheck = alphaw;
	    Hcheck = H.boundaryField()[0][0];
	    precice.markActionFulfilled(cowic);
	}

	// Data reading
	precice.readBlockScalarData(alphaID, vertexSize, vertexIDs, alphaw);
	precice.readBlockVectorData(velocID, vertexSize, vertexIDs, velocity);
	precice.readBlockScalarData(prghID, vertexSize, vertexIDs, prgh);
	Info<< "alpha = " << alphaw[0] << endl;
	Info<< "velocity = " << velocity[0] << endl;
	// End Data reading

	// Alpha to H conversion and updating BC H
	double H_BC = 0;
	for ( int i = 1; i < nFaces + 1 ; i++ )				// For vertices above groud level
	{
	    H_BC += alphaw[i] * faceHeight;				// H = sum(alpha * faceheight)
	}

	H.boundaryFieldRef()[0][0] = H_BC;				// Set boundary condition H
	// End Alpha to H conversion and updating BC H

	Info<< "H_bound = " << H.boundaryField()[0][0] << endl;

	// U to HU conversion and updating BC HU
	double HUx_BC = 0;
	double HUy_BC = 0;
	double HUz_BC = 0;
	for ( int i = 0; i < nFaces + 1 ; i++ )				// For all vertices 
	{
	    if ( alphaw[i] > 0 )
	    {
		HUx_BC += alphaw[i] * velocity[0+i*3] / (nFaces+1);	// HUx = mean(alpha * Ux)
		HUy_BC += alphaw[i] * velocity[1+i*3] / (nFaces+1);	// HUy = mean(alpha * Uy)
		HUz_BC += alphaw[i] * velocity[2+i*3] / (nFaces+1);	// HUz = mean(alpha * Uz)
	    }
	}

        HU.boundaryFieldRef()[0][0].component(0) = HUx_BC;		// Set boundary condition HUx
        HU.boundaryFieldRef()[0][0].component(1) = HUy_BC;		// Set boundary condition HUy
        HU.boundaryFieldRef()[0][0].component(2) = HUz_BC;		// Set boundary condition HUz
	// End U to HU conversion and updating BC HU

	Info<< "HUx_BC = " << HUx_BC << endl;
	Info<< "HU = " << HU.boundaryField()[0][0] << endl;


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
             fvm::ddt(HU)          	// d(HU_i)/dt			Local derivative
	     + fvm::div(phi, HU)   	// + d/dx_j ( U_j * HU_i )	Convection
	     + H*gradgSpHf          	// + H * d/dx_i ( gSpH )	Downhill-slope force and acceleration due to change in flowdepth
	     + fvm::Sp(alpha,HU)   	// + alpha * HU			Bottom friction
	     - fvm::laplacian(nut,HU) 	// - d^2/dx_j^2 ( nut * HU )	Turbulent stresses
         );  

	HUEqn.solve();

///////////////////////////////////////////////////////////////////////////////
///     Calculation of U and phi

        // if Hclip < Hdry2, then U = 0
        HU = HU*pos(Hclip - Hdry2); 
        U = HU/Hclip;
       
        phi = (fvc::interpolate(U) & mesh.Sf());

///////////////////////////////////////////////////////////////////////////////

	
	double h = H.boundaryField()[0][0];		// Flow depth on 2D boundary
	double zb = S.boundaryField()[0][0];		// Bottom elevation on 2D boundary
	double zw = zb + h;				// Water surface height on 2D boundary

	// H to alpha conversion and updating BC alphaw
	for ( int i = 0; i<vertexSize ; i++ )			
	{
	    if ( zw < coords[1+3*i] - 0.5 * (meshHeight / nFaces))		// For H <
	    {
		alphaw[i] = 0;
	    }
	    else if ( zw > coords[1+3*(i)] + 0.5 * (meshHeight / nFaces))	// For H >
	    {
		alphaw[i] = 1;
	    }
	    else								// For < H <
	    {
		alphaw[i] = 0.5 + (zw - coords[1+3*i]) / faceHeight;		// Interpolate
	    }
	}
	
	//precice.writeBlockScalarData(flowdID, vertexSize, vertexIDs, alphaw);
	// End H to alpha conversion and updating BC alphaw
	
	Info<< "H_bound = " << H.boundaryField()[0][0] << endl;
	Info<< "new aplhaw = " << alphaw[1] << endl;


	// H to Prgh conversion and updating BC Prgh
	for ( int i = 0; i < nFaces + 1 ; i++ )	
	{
	    prgh[i] = (alphaw[i] * rhow3D + (1 - alphaw[i]) * rhoa3D) * g3D * zw;
	}
	// End H to Prgh conversion and updating BC Prgh

	Info<< "new prgh = " << prgh[1] << endl;


	// H and HU to U conversion and updating BC U
	for ( int i = 0; i < nFaces + 1 ; i++ )	
	{
	    //velocity[0+i*3] = U[0] + 1 * ( 1 + log ((coords[1+3*i]-zb)/h) );
	    //velocity[1+i*3] = U[0] + 1 * ( 1 + log ((coords[1+3*i]-zb)/h) );
	    //velocity[2+i*3] = U[0] + 1 * ( 1 + log ((coords[1+3*i]-zb)/h) );
	} 
	// End H and HU to U conversion and updating BC U

	Info<< "new U = " << velocity[1] << nl << endl;


	precice_dt = precice.advance(dt);
	
	if(precice.isActionRequired(coric))	// timestep not converged
	{ 
	    // load initial data
	    alphaw = alphawcheck;
	    H.boundaryFieldRef()[0][0] = Hcheck;
	    precice.markActionFulfilled(coric);
	}
	else{ 
	    runTime.setDeltaT(dt);
            runTime.write();
	}


        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    precice.finalize(); 			// Frees data structures and closes communication channels



    return(0);
}


// ************************************************************************* //
