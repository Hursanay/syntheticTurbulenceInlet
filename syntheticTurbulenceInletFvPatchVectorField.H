/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     | Website:  https://openfoam.org
  \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
  \\/     M anipulation  |
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
  Velocity inlet boundary condition for OpenFOAM-dev that generates
  synthetic turbulent fluctuations consistent with a target
  Reynolds-stress tensor and integral-length-scale field, using a
  divergence-free synthetic eddy method (DFSEM).

  The implementation derives from ESI's turbulentDFSEMInlet with three
  classes of improvement:
    - equation-level corrections that bring the BC into line with the
      original DFSEM paper (Poletto, Craft & Revell, 2013);
    - re-derived normalisation constants that improve the match to DNS
      reference data (especially near-wall);
    - a substantially reduced configuration burden through
      geometry-aware sampling, automatic profile mapping, swirl, ghost
      eddies, and reproducible RNG.

  Supported geometryMode values: annulus, disk, rectangle, channel
  (rectangle with PBC in one direction), arbitrary. Bundled DNS
  profiles (U, R, L) live under URL/, and ready-to-run tutorials for
  each geometryMode are under tutorials/.

  See https://github.com/Hursanay/syntheticTurbulenceInlet/tree/main/documentation/doc_syntheticTurbulence.pdf for the full method
  description and the complete parameter table.

  If you use this boundary condition in academic work, please cite:
    1. H. Fyhn, A. Gruber, H. T. Nygard, J. Dawson, and K.-J.
       Nogenmyr, "Numerical and Experimental Investigation of a
       Non-Premixed Flame of Partially Decomposed Ammonia at
       Atmospheric Pressure," Journal of Engineering for Gas Turbines
       and Power (accepted, 2026).
    2. H. Fyhn, O. H. H. Meyer, Q. Wang, N. Worth, J. Koomen, T.
       Dammers, and A. Gruber, "Numerical and Experimental
       Investigation of a Laboratory-Scale Hydrogen-Fired FlameSheet
       Burner Operated at Atmospheric Pressure With Full Optical
       Access," Journal of Engineering for Gas Turbines and Power
       (accepted, 2026).

  Usage
  Load the library in system/controlDict:
  \verbatim
  libs ("libsyntheticTurbulenceInlet.so");
  \endverbatim

  Minimal patch entry in 0/U:
  \verbatim
  <patchName>
  {
    type            syntheticTurbulenceInlet;
    geometryMode    annulus;
    massFlowMode    true;
    massFlow        2.53e-4;
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
        scalar Utau_;

        bool UtauIsInDict_;

        //scaling factor to go from y to yplus (Utau_/nu)
        scalar Sr_;

        //average kinematic viscosity (nu=mu/rho)
        scalar avNu_;

        //the geometry can be one of these: annulus, disk, rectangle, channel (rectangle with PBC in one direction), arbitrary (the purpose is to generate turbulence that do not need to follow DNS profiles)
        word geometryMode_;
        
        //For channel: Are the walls on the long or the short edge? The other direction is assumed to have PBC.
        word wallLongOrShort_;

        //when geometryMode_ is arbitrary, its radius and its center are
        scalar arbitraryModeRadius_; //default largest distance from the patch center
        vector arbitraryModeCenter_; //default centerCoord_
        bool arbitraryModeCenterIsInDict_;

        //- Mean velocity field
        vectorField U_;

        // Two lengths where length2_>length1_; radii for annulus and disk, length and height for rectangle and channel;
        scalar length1_;
        scalar length2_;

        //the center coordinate
        vector centerCoord_;
        
        //the normal vector of the patch
        vector patchNormal_;

        //the inplane vectors of the patch
        vector inplaneVec1_;
        vector inplaneVec2_;

        //If user specifies "massFlow" or not (if not, user specifies "Umean")
        bool massFlowMode_;

        //- Normalisation factor for mean velocity or mass flow
        autoPtr<Function1<scalar>> inFlow_;

        //Mean velocity (if inFlow_ is mass flow, Umean_ is inFlow_/rho/area) 
        scalar Umean_;

        //- Average velocity from DNS data used for scaling
        scalar Uavg_;

        //- Swirl number 
        scalar swirlNr_;

        // Eddy information
        List<eddy> eddies_;

        //- Number of eddies_ that is used simultaniously
        label nEddy_;

        //- Characteristic length scale
        scalar delta_;

        //- Ratio of sum of eddy volumes to eddy box volume, i.e. eddy density
        const scalar d_;

        //- von Karman constant
        const scalar kappa_;

        //- Normalisation factor for integral-length scale
        const scalar Lref_;

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
        const label seed_;

        //random number generator
        std::mt19937 randGen_;

        //Density when incompressible solver is used (default 1)
        scalar rhoInf_;

    public:

      //- Runtime type information
      TypeName("syntheticTurbulenceInlet");


      // Constructors

      //- Construct from patch, internal field and dictionary
      syntheticTurbulenceInletFvPatchVectorField
        (
         const fvPatch&,
         const DimensionedField<vector, fvMesh>&,
         const dictionary&
        );

      //- Construct by mapping given syntheticTurbulenceInletFvPatchVectorField
      //  onto a new patch
      syntheticTurbulenceInletFvPatchVectorField
        (
         const syntheticTurbulenceInletFvPatchVectorField&,
         const fvPatch&,
         const DimensionedField<vector, fvMesh>&,
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
         const DimensionedField<vector, fvMesh>&
        );

      //- Construct and return a clone setting internal field reference
      virtual tmp<fvPatchVectorField> clone
        (
         const DimensionedField<vector, fvMesh>& iF
        ) const
        {
          return tmp<fvPatchVectorField>
            (
             new syntheticTurbulenceInletFvPatchVectorField(*this, iF)
            );
        }


      // Member Functions
      
      // Get the present mean U
      void updateUmean();

      //- Density on the patch faces: returns the looked-up rho patch field
      //  if phi is mass flow rate, otherwise a constant field of rhoInf_.
      tmp<scalarField> rhoOnPatch() const;

      //initialize geometry
      void initializeGeometry();

      // Get random numbers and intensities and increase random number counter by 1
      void getRandomPositionAndIntensities(vector& pos, vector& intensities);

      void getRandSigmaFromPos(const vector& pos, symmTensor& R, scalar& sigmax) const;

      //makes a new eddy
      eddy oneEddy(scalar depth); 

      // Initiliazing profiles
      void initializeUProfile();

      //- Initialize eddies_
      void initializeEddies();

      //- Convect the eddies_ with the bulk velocity (Umean_)
      void convectEddies();

      // contribution to velocity (U') at a point from all eddies
      vector uPrimeEddy(const vector& point) const;

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
