/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     | Website:  https://openfoam.org
  \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
  \\/     M anipulation  |
  -------------------------------------------------------------------------------
  Authors: Eirik Fyhn and Hursanay Fyhn

  Inspired by turbulentDFSEMInlet.
  -------------------------------------------------------------------------------
  License
  This file is part of OpenFOAM.

  OpenFOAM is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

  Class
  Foam::syntheticTurbulenceInletFvPatchVectorField

  Description
  Boundary condition for creating synthetic turbulence at the inlet 

  Usage
  Example of the boundary condition specification:
  \verbatim
  <patchName>
  {
    type            syntheticTurbulenceInlet;
    value           uniform (0 0 0);    // Initial value
    geometryMode    annulus;
    Umean           1;
  }
  \endverbatim

  See also
  Foam::fixedValueFvPatchVectorField

  SourceFiles
  syntheticTurbulenceInletFvPatchVectorField.C

  \*---------------------------------------------------------------------------*/

#ifndef syntheticTurbulenceInletFvPatchVectorField_H
#define syntheticTurbulenceInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "eddy.H"
#include <random>
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  /*---------------------------------------------------------------------------*\
    Class syntheticTurbulenceInletFvPatchVectorField Declaration
    \*---------------------------------------------------------------------------*/

  class syntheticTurbulenceInletFvPatchVectorField
    :
      public fixedValueFvPatchVectorField
  {

    // Private Data
    
        //friction velocity
        scalar Utau;

        bool UtauIsInDict;

        //scaling factor to go from y to yplus (Utau/nu)
        scalar Sr;

        //average kinematic viscosity (nu=mu/rho)
        scalar avNu;

        //the geometry can be one of these: an annulus, a disk, a rectangle, and a channel (rectangle with PBC in one direction)
        word geometryMode;
        
        //For channel: Are the walls on the long or the short edge? The other direction is assumed to have PBC.
        word wallLongOrShort;

        //- Mean velocity field
        vectorField U_;

        // Two lengths where length2>length1; radii for annulus and disk, length and height for rectangle and channel;
        scalar length1;
        scalar length2;

        //the center coordinate
        vector centerCoord;
        
        //the normal vector of the patch
        vector patchNormal_;

        //the inplane vectors of the patch
        vector inplaneVec1_;
        vector inplaneVec2_;

        //- Normalisation factor for mean velocity
        autoPtr<Function1<scalar>> Umean;

        //Velocity scale to determine eddy size and intensity, defaults to Umean(time=startTime)
        scalar Uscale;

        //- Average velocity used for scaling
        scalar Uavg;

        //- Swirl number 
        scalar swirlNr;

        // Eddy information
        List<eddy> eddies;

        //- Number of eddies that is used simultaniously
        label nEddy_;

        //- Characteristic length scale
        scalar delta_;

        //- Ratio of sum of eddy volumes to eddy box volume, i.e. eddy density
        const scalar d_;

        //- von Karman constant
        const scalar kappa_;

        //- Normalisation factor for integral-length scale
        const scalar Lref;

        //- Heuristic scaling factor being applied on the normalisation factor
        const scalar scale_;

        //- Minimum eddy length in units of number of cells
        const label nCellPerEddy_;

        //- Minimum integral-length scale (across all processors)
        scalar minSigma_;

        //- Eddy box volume
        scalar v0_;

        //- Maximum integral-length scale (across all processors)
        scalar maxSigmaX_;

        //- have the initialization been done 
        bool notInitialized_;

        //the seed for the random number generator
        const unsigned int seed_;

        //random number generator
        std::mt19937 randGen;

    public:

      //- Runtime type information
      TypeName("syntheticTurbulenceInlet");


      // Constructors

      //- Construct from patch, internal field and dictionary
      syntheticTurbulenceInletFvPatchVectorField
        (
         const fvPatch&,
         const DimensionedField<vector, volMesh>&,
         const dictionary&
        );

      //- Construct by mapping given syntheticTurbulenceInletFvPatchVectorField
      //  onto a new patch
      syntheticTurbulenceInletFvPatchVectorField
        (
         const syntheticTurbulenceInletFvPatchVectorField&,
         const fvPatch&,
         const DimensionedField<vector, volMesh>&,
         const fieldMapper&
        );

      //- Disallow copy without setting internal field reference
      syntheticTurbulenceInletFvPatchVectorField
        (
         const syntheticTurbulenceInletFvPatchVectorField&
        ) = delete;

      //- Copy constructor setting internal field reference
      syntheticTurbulenceInletFvPatchVectorField
        (
         const syntheticTurbulenceInletFvPatchVectorField&,
         const DimensionedField<vector, volMesh>&
        );

      //- Construct and return a clone setting internal field reference
      virtual tmp<fvPatchVectorField> clone
        (
         const DimensionedField<vector, volMesh>& iF
        ) const
        {
          return tmp<fvPatchVectorField>
            (
             new syntheticTurbulenceInletFvPatchVectorField(*this, iF)
            );
        }


      // Member Functions
      
      //initialize geometry
      void initializeGeometry();

      // Get random numbers and intensities and increase random number counter by 1
      void getRandomPositionAndIntensities(vector& pos, vector& intensities);

      void getRandSigmaFromPos(vector pos, symmTensor& R, scalar& sigmax);

      //makes a new eddy
      eddy oneEddy(scalar depth); 

      // Initiliazing profiles
      void initializeUProfile();

      //- Initialize eddies
      void initializeEddies();

      //- Convect the eddies with the bulk velocity (Umean)
      void convectEddies();

      // contribution to velocity (U') at a point from all eddies
      vector uPrimeEddy(vector point);

      //- Update the coefficients associated with the patch field
      virtual void updateCoeffs();

      //- Write
      virtual void write(Ostream&) const;

  };
  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
