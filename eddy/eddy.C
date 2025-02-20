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

\*---------------------------------------------------------------------------*/

#include "eddy.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::eddy::GammaValues[] = {1.0, Foam::sqrt(2.0), Foam::sqrt(3.0), Foam::sqrt(4.0), Foam::sqrt(5.0), Foam::sqrt(6.0), Foam::sqrt(7.0), Foam::sqrt(8.0)};
Foam::UList<Foam::scalar> Foam::eddy::Gamma(&GammaValues[0], 8);
int Foam::eddy::debug = 0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::eddy::setScales
(
    const scalar sigmaX,
    const scalar gamma,
    const vector& e,
    const vector& lambda,
    vector& sigma,
    vector& alpha
) const
{
    // Length scale in the largest eigenvalue direction
    const label d1 = 2; // Set the eddy orientation to position of max eigenvalue; (direction of eddy major axis, sigma_x in reference)
    const label d2 = (d1 + 1) % 3;
    const label d3 = (d1 + 2) % 3;

    sigma[d1] = sigmaX;

    // Note: sigma_average = 1/3*(sigma_x + sigma_y + sigma_z)
    // Substituting for sigma_y = sigma_x/gamma and sigma_z = sigma_y
    // Other length scales equal, as function of major axis length and gamma
    sigma[d2] = sigma[d1]/gamma;
    sigma[d3] = sigma[d2];

    // (PCR:Eq. 13)
    const vector sigma2(cmptMultiply(sigma, sigma));

    bool ok = true;

    const scalar x0 = lambda[2]/sigma2[2]
                     -lambda[0]/sigma2[0] 
                     +lambda[1]/sigma2[1] ;

    const scalar x1 = lambda[0]/sigma2[0]
                     -lambda[1]/sigma2[1] 
                     +lambda[2]/sigma2[2] ;

    const scalar x2 = lambda[1]/sigma2[1]
                     -lambda[2]/sigma2[2] 
                     +lambda[0]/sigma2[0] ;

    if (x0 <0 || x1 <0 || x2 <0) {
      alpha = vector (0,0,0);
      ok = false; //means will try another gamme in constructor
    } else {
      alpha[0] = e[0] * sqrt(x0);
      alpha[1] = e[1] * sqrt(x1);
      alpha[2] = e[2] * sqrt(x2);
    }

    if (debug > 1)
    {
        Pout<< ", gamma:" << gamma
            << ", lambda:" << lambda
            << ", sigma2: " << sigma2
            << ", sigmaX:" << sigmaX
            << ", sigma:" << sigma
            << ", alpha:" << alpha
            << endl;
    }

    return ok;
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::eddy::eddy()
:
    position0_(Zero),
    x_(0),
    sigma_(Zero),
    alpha_(Zero),
    Rpg_(tensor::I),
    c1_(-1)
{}


Foam::eddy::eddy
(
    const point& position0,
    const scalar x,
    const scalar sigmaX,
    const symmTensor& R,
    const vector& e 
)
:
    position0_(position0),
    x_(x),
    sigma_(Zero),
    alpha_(Zero),
    Rpg_(tensor::I),
    c1_(-1)
{
    // Principal stresses - eigenvalues returned in ascending order
    const vector lambda(eigenValues(R));

    // Eddy rotation from principal-to-global axes
    // - given by the 3 eigenvectors of the Reynold stress tensor as rows in
    //   the result tensor (transposed transformation tensor)
    // - returned in ascending eigenvalue order
    Rpg_ = eigenVectors(R, lambda).T();

    if (debug)
    {
        // Global->Principal transform = Rgp = Rpg.T()
        // Rgp & R & Rgp.T() should have eigenvalues on its diagonal and
        // zeros for all other components
        Pout<< "Rpg.T() & R & Rpg: " << (Rpg_.T() & R & Rpg_) << endl;
    }

    // Set intensities and length scales
    bool found = false;
    forAll(Gamma, i)
    {
        // Random length scale ratio, gamma = sigmax/sigmay = sigmax/sigmaz
        const scalar gamma = Gamma[i];

        if (setScales(sigmaX, gamma, e, lambda, sigma_, alpha_))
        {
            found = true;
            break;
        }
    }

    c1_ = 1.0/ Foam::sqrt(volume());

    if (found)
    {
        // Shuffle the gamma^2 values
        // Applying this shuffle: [2, 8, 6, 7, 4, 5, 1, 3]
        scalar G1 = Gamma[0];
        scalar G2 = Gamma[1];
        scalar G3 = Gamma[2];
        scalar G4 = Gamma[3];
        scalar G5 = Gamma[4];
        scalar G6 = Gamma[5];
        scalar G7 = Gamma[6];
        scalar G8 = Gamma[7];

        Gamma[0] = G2;
        Gamma[1] = G8;
        Gamma[2] = G6;
        Gamma[3] = G7;
        Gamma[4] = G4;
        Gamma[5] = G5;
        Gamma[6] = G1;
        Gamma[7] = G3;
    }
}


Foam::eddy::eddy(const eddy& e)
:
    position0_(e.position0_),
    x_(e.x_),
    sigma_(e.sigma_),
    alpha_(e.alpha_),
    Rpg_(e.Rpg_),
    c1_(e.c1_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::eddy::uPrime(const point& xp, const vector& n) const
{
    // This was different in DFSEM. To be consistent with the order of the components in sigma, it seems one should rather transform to eddy principal system before diving by sigma_.
    
    // Relative position inside eddy (eddy principal system)
    const vector rp(cmptDivide(Rpg_.T() & (xp - position(n)), sigma_));

    if (mag(rp) >= scalar(1))
    {
        return vector::zero;
    }

    // Shape function (eddy principal system)
    /* const vector q(cmptMultiply(sigma_, vector::one - cmptMultiply(rp, rp))); */
    const vector q(sigma_*(1 - magSqr(rp)));

    // Fluctuating velocity (eddy principal system) (PCR:Eq. 8)
    const vector uPrimep(cmptMultiply(q, rp^alpha_));

    // Convert into global system (PCR:Eq. 10)
    return c1_*(Rpg_ & uPrimep);
}

// ************************************************************************* //
