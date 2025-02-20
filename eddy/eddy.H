/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Authors: Eirik Fyhn and Hursanay Fyhn

Inspired by turbulentDFSEMInlet.
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Foam::eddy

Description
    Class to describe eddies for the syntheticTurbulenceInletInletFvPatchVectorField
    boundary condition.

SourceFiles
    eddy.C
    eddyI.H

\*---------------------------------------------------------------------------*/

#ifndef syntheticTurbulenceInletFvPatchVectorField_eddy_H
#define syntheticTurbulenceInletFvPatchVectorField_eddy_H

#include "vector.H"
#include "point.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class eddy;

bool operator==(const eddy& a, const eddy& b);
bool operator!=(const eddy& a, const eddy& b);

/*---------------------------------------------------------------------------*\
                            Class eddy Declaration
\*---------------------------------------------------------------------------*/

class eddy
{
    // Private Data

        static scalar GammaValues[8];
        static UList<scalar> Gamma;

        //- Reference position
        point position0_;

        //- Distance from reference position in normal direction
        scalar x_;

        //- Integral-length scales in 3-D space
        vector sigma_;

        //- Time-averaged intensity
        vector alpha_;

        //- Coordinate system transformation from local to global axes
        //  X-direction aligned with max stress eigenvalue
        tensor Rpg_;

        //- Model coefficient c1
        scalar c1_;

    // Private Member Functions

        //- Set the eddy scales: length, intensity
        bool setScales
        (
            const scalar sigmaX,
            const scalar gamma,
            const vector& e,
            const vector& lambda,
            vector& sigma,
            vector& alpha
        ) const;

public:

    // Constructors

        //- Construct null
        eddy();

        //- Construct from components
        eddy
        (
            const point& position0,     // reference position
            const scalar x,             // distance from reference position
            const scalar sigmaX,        // integral-length scale
            const symmTensor& R,        // Reynolds stress tensor
            const vector& e             // Random vector of 1's and -1's
        );

        //- Copy construct
        eddy(const eddy& e);


    // Public Data

        //- Flag to activate debug statements
        static int debug;


    // Public Member Functions

        // Access

            //- Return the reference position
            inline const point& position0() const noexcept;

            //- Return the distance from the reference position
            inline scalar x() const noexcept;

            //- Return the length scales in 3-D space
            inline const vector& sigma() const noexcept;

            //- Return the time-averaged intensity
            inline const vector& alpha() const noexcept;

            //- Return the coordinate system transformation from local
            //- principal to global axes
            inline const tensor& Rpg() const noexcept;

            //- Return the model coefficient c1
            inline scalar c1() const noexcept;

            //- Return the eddy position
            inline point position(const vector& n) const;

        // Helper functions

            //- Volume
            inline scalar volume() const;

            //- Move the eddy
            inline void move(const scalar dx);

        // Evaluation

            //- Return the fluctuating velocity contribution at local point xp
            vector uPrime(const point& xp, const vector& n) const;

    // Member Operators

        void operator=(const eddy& e){
            position0_ =e.position0_;
            x_ =e.x_;
            sigma_ =e.sigma_;
            alpha_ =e.alpha_;
            Rpg_ =e.Rpg_;
            c1_ =e.c1_;
        }


    // Friend Operators

        friend bool operator==(const eddy& a, const eddy& b)
        {
            return
             a.position0_ == b.position0_
             && a.x_ == b.x_
             && a.sigma_ == b.sigma_
             && a.alpha_ == b.alpha_
             && a.Rpg_ == b.Rpg_
             && a.c1_ == b.c1_;
        }

        friend bool operator!=(const eddy& a, const eddy& b)
        {
            return !(a == b);
        }
};
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "eddyI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
