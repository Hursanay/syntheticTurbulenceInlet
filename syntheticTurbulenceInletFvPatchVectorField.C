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

\*---------------------------------------------------------------------------*/

#include "syntheticTurbulenceInletFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include <algorithm>

// DNS profile tables (Re_tau ~ 590 channel; constant for one dataset).
// Pulled in once at file scope rather than re-defined in each member
// function body. URL_len.H must come first since the others size their
// arrays from it.
namespace
{
    using Foam::scalar;
    #include "URL/URL_len.H"
    #include "URL/r_list.H"
    #include "URL/rplus_list.H"
    #include "URL/U.H"
    #include "URL/L.H"
    #include "URL/R.H"
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::syntheticTurbulenceInletFvPatchVectorField::
syntheticTurbulenceInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    Utau_(dict.lookupOrDefault<scalar>("Utau", 0.0)),
    UtauIsInDict_(dict.found("Utau")),
    geometryMode_(dict.lookup<word>("geometryMode")),
    wallLongOrShort_(dict.lookupOrDefault<word>("wallLongOrShort","long")),
    arbitraryModeRadius_(dict.lookupOrDefault<scalar>("arbitraryModeRadius",-1)),
    arbitraryModeCenter_(dict.lookupOrDefault<vector>("arbitraryModeCenter",vector(0,0,0))),
    arbitraryModeCenterIsInDict_(dict.found("arbitraryModeCenter")),
    U_(size(), Zero),
    massFlowMode_(dict.lookupOrDefault<bool>("massFlowMode",false)),
    inFlow_(
        Function1<scalar>::New
        (
            massFlowMode_ ? "massFlow" : "Umean",
            this->db().time().userUnits(),
            iF.dimensions(),
            dict
        )
    ),
    Umean_(0),
    Uavg_(0),
    swirlNr_(dict.lookupOrDefault<scalar>("swirlNr", 0)),
    nEddy_(0),
    delta_(1),
    d_(dict.lookupOrDefault<scalar>("d", 1)),
    kappa_(0.41),
    Lref_(dict.lookupOrDefault<scalar>("Lref", 1)),
    scale_(dict.lookupOrDefault<scalar>("scale", geometryMode_ == "arbitrary" ? 0.2 : 1)),
    nCellPerEddy_(dict.lookupOrDefault<label>("nCellPerEddy", geometryMode_ == "arbitrary" ? 5 : 1)),
    notInitialized_(true),
    seed_(dict.lookupOrDefault<label>("seed", 536)),
    randGen_(seed_),
    rhoInf_(dict.lookupOrDefault<scalar>("rhoInf",1))
{
    vectorField::operator=(U_);
}

Foam::syntheticTurbulenceInletFvPatchVectorField::
syntheticTurbulenceInletFvPatchVectorField
(
    const syntheticTurbulenceInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper, false),
    Utau_(ptf.Utau_),
    UtauIsInDict_(ptf.UtauIsInDict_),
    Sr_(ptf.Sr_),
    avNu_(ptf.avNu_),
    geometryMode_(ptf.geometryMode_),
    wallLongOrShort_(ptf.wallLongOrShort_),
    arbitraryModeRadius_(ptf.arbitraryModeRadius_),
    arbitraryModeCenter_(ptf.arbitraryModeCenter_),
    arbitraryModeCenterIsInDict_(ptf.arbitraryModeCenterIsInDict_),
    U_(ptf.U_),
    length1_(ptf.length1_),
    length2_(ptf.length2_),
    centerCoord_(ptf.centerCoord_),
    patchNormal_(ptf.patchNormal_),
    inplaneVec1_(ptf.inplaneVec1_),
    inplaneVec2_(ptf.inplaneVec2_),
    massFlowMode_(ptf.massFlowMode_),
    inFlow_(ptf.inFlow_, false),
    Umean_(ptf.Umean_),
    Uavg_(ptf.Uavg_),
    swirlNr_(ptf.swirlNr_),
    eddies_(ptf.eddies_),
    nEddy_(ptf.nEddy_),
    delta_(ptf.delta_),
    d_(ptf.d_),
    kappa_(ptf.kappa_),
    Lref_(ptf.Lref_),
    scale_(ptf.scale_),
    nCellPerEddy_(ptf.nCellPerEddy_),
    minSigma_(ptf.minSigma_),
    v0_(ptf.v0_),
    maxSigmaX_(ptf.maxSigmaX_),
    notInitialized_(ptf.notInitialized_),
    seed_(ptf.seed_),
    randGen_(ptf.randGen_),
    rhoInf_(ptf.rhoInf_)
{}


Foam::syntheticTurbulenceInletFvPatchVectorField::
syntheticTurbulenceInletFvPatchVectorField
(
    const syntheticTurbulenceInletFvPatchVectorField& mwvpvf,
    const DimensionedField<vector, fvMesh>& iF
)
:
    fixedValueFvPatchVectorField(mwvpvf, iF),
    Utau_(mwvpvf.Utau_),
    UtauIsInDict_(mwvpvf.UtauIsInDict_),
    Sr_(mwvpvf.Sr_),
    avNu_(mwvpvf.avNu_),
    geometryMode_(mwvpvf.geometryMode_),
    wallLongOrShort_(mwvpvf.wallLongOrShort_),
    arbitraryModeRadius_(mwvpvf.arbitraryModeRadius_),
    arbitraryModeCenter_(mwvpvf.arbitraryModeCenter_),
    arbitraryModeCenterIsInDict_(mwvpvf.arbitraryModeCenterIsInDict_),
    U_(mwvpvf.U_),
    length1_(mwvpvf.length1_),
    length2_(mwvpvf.length2_),
    centerCoord_(mwvpvf.centerCoord_),
    patchNormal_(mwvpvf.patchNormal_),
    inplaneVec1_(mwvpvf.inplaneVec1_),
    inplaneVec2_(mwvpvf.inplaneVec2_),
    massFlowMode_(mwvpvf.massFlowMode_),
    inFlow_(mwvpvf.inFlow_, false),
    Umean_(mwvpvf.Umean_),
    Uavg_(mwvpvf.Uavg_),
    swirlNr_(mwvpvf.swirlNr_),
    eddies_(mwvpvf.eddies_),
    nEddy_(mwvpvf.nEddy_),
    delta_(mwvpvf.delta_),
    d_(mwvpvf.d_),
    kappa_(mwvpvf.kappa_),
    Lref_(mwvpvf.Lref_),
    scale_(mwvpvf.scale_),
    nCellPerEddy_(mwvpvf.nCellPerEddy_),
    minSigma_(mwvpvf.minSigma_),
    v0_(mwvpvf.v0_),
    maxSigmaX_(mwvpvf.maxSigmaX_),
    notInitialized_(mwvpvf.notInitialized_),
    seed_(mwvpvf.seed_),
    randGen_(mwvpvf.randGen_),
    rhoInf_(mwvpvf.rhoInf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::scalarField>
Foam::syntheticTurbulenceInletFvPatchVectorField::rhoOnPatch() const
{
  const surfaceScalarField& phi =
      this->internalField().mesh().lookupObject<surfaceScalarField>("phi");
  if (phi.dimensions() == dimMass/dimTime)
  {
    return tmp<scalarField>
    (
      new scalarField
      (
        patch().lookupPatchField<volScalarField, scalar>("rho")
      )
    );
  }
  return tmp<scalarField>(new scalarField(patch().size(), rhoInf_));
}

void Foam::syntheticTurbulenceInletFvPatchVectorField::updateUmean()
{
  const scalarField rho(rhoOnPatch());
  Umean_ = inFlow_->value(db().time().value())
        / (massFlowMode_ ? gSum(rho*patch().magSf()) : 1);
  /* Info << "Umean_ = " << Umean_<< endl; */
}

void Foam::syntheticTurbulenceInletFvPatchVectorField::initializeGeometry()
{
  centerCoord_ = gAverage(patch().Cf());
  /* Info << "Center coord: " << centerCoord_ << endl; */
  const pointField& faceVertices = patch().poly().localPoints();

  // Patch normal variable here points into domain
  const vectorField nf(patch().nf());
  patchNormal_ = -gSum(nf);

  // reconstructPar creates fvPatch with empty patch
  if(mag(patchNormal_) > 0.0){patchNormal_ /= mag(patchNormal_);}

  // Check that patch is planar
  const scalar error = max(magSqr(patchNormal_ + nf));
  if (error > SMALL)
  {
      WarningInFunction
          << "Patch " << patch().name() << " is not planar"
          << endl;
  }

  //determine the two in-plane vectors and the length scales for the different geometries
  inplaneVec1_ = vector(0, 1, 0) ^ patchNormal_;
  if (mag(inplaneVec1_) < 0.00001){
    inplaneVec1_ = vector(0, 0, 1) ^ patchNormal_;
  }
  inplaneVec2_ = inplaneVec1_ ^ patchNormal_;
  inplaneVec1_ /= mag(inplaneVec1_);
  inplaneVec2_ /= mag(inplaneVec2_);
  scalar distBetweenWalls(0); //distance between the walls
  if (geometryMode_ == "annulus") 
  {
    length1_ = gMin(mag(faceVertices - centerCoord_)); // inner radius
    length2_ = gMax(mag(faceVertices - centerCoord_)); // outer radius
    Info << "\nDetected an annulus with radii = " <<  length1_ << " and " << length2_ << endl << endl;
    if (length1_ < 1e-8) {
      WarningInFunction
          << "The inner radius is very small, are you sure it is not a \"geometryMode_ disk\" rather than a \"geometryMode_ annulus\"? "
          << endl;
    }
    distBetweenWalls = length2_ - length1_;

  } else if (geometryMode_ == "disk") {
    length1_ = 0.0;
    length2_ = gMax(mag(faceVertices - centerCoord_)); 
    Info << "\nDetected a disk with radius = " <<  length2_ << endl << endl;
    distBetweenWalls = 2 * length2_;

  } else if (geometryMode_ == "rectangle" || geometryMode_ == "channel") {
    scalar D = gMax(mag(faceVertices - centerCoord_));
    /* Info << "Diagonal: " << D << endl; */ 
    scalar rect2 = gMax((faceVertices - centerCoord_) & inplaneVec2_);
    scalar rect3 = gMax((faceVertices - centerCoord_) & inplaneVec1_);
    scalar rect1 = sqrt(D*D - rect2*rect2);
    scalar rect4 = sqrt(D*D - rect3*rect3);

    const scalar minCellDx(gMin(Foam::sqrt(patch().magSf())));
    vector rect5 = rect1 * inplaneVec1_ + rect2 * inplaneVec2_;
    //if the diagonal vector rect5 is outside the patch, flip it:  
    if (gMin(mag(faceVertices - centerCoord_ - rect5)) > minCellDx) {
      rect5 = -rect1 * inplaneVec1_ + rect2 * inplaneVec2_;
    }
    vector rect6 = rect3 * inplaneVec1_ + rect4 * inplaneVec2_;
    int switchInt = -1;

    //if the diagonal vectors rect 5 & 6 are pointing in the same direction:
    if (mag(rect5 ^ rect6) < minCellDx) {
      rect6 = rect3 * inplaneVec1_ + (switchInt * rect4) * inplaneVec2_;
      switchInt *= -1;
    }

    //if the diagonal vector rect6 is outside the patch, flip it:  
    if (gMin(mag(faceVertices - centerCoord_ - rect6)) > minCellDx) {
      if (switchInt == 1) {
        FatalErrorInFunction << "Failed to find inplane vectors. This should not happen." << exit(FatalError);
      }
      rect6 = rect3 * inplaneVec1_ + (switchInt * rect4) * inplaneVec2_;
    }
    
    //rebaptize:
    inplaneVec1_ = rect6 - rect5;
    inplaneVec2_ = rect5 + rect6;
    length1_ = mag(inplaneVec1_);
    length2_ = mag(inplaneVec2_);
    if(length1_ > 0.0 && length2_ > 0.0)
    {
      //order so that length2_ > length1_:
      if (length1_ > length2_) {
        scalar temp1 = length2_;
        length2_ = length1_;
        length1_ = temp1;
        vector temp2 = inplaneVec2_;
        inplaneVec2_ = inplaneVec1_;
        inplaneVec1_ = temp2;
      }
      inplaneVec1_ /= length1_;
      inplaneVec2_ /= length2_;
      Info << "\nDetected a " << geometryMode_ << " with dimensions " <<  length1_ << " and " << length2_ << endl;// << endl;
      Info << "in directions " << inplaneVec1_ << " and " <<  inplaneVec2_ << endl << endl;
    }

    distBetweenWalls = length1_;
    if (geometryMode_ == "channel" && wallLongOrShort_ == "short" ) {
      distBetweenWalls = length2_;
    }
  } else if (geometryMode_ == "arbitrary") {
    length1_ = 0.0;
    if (arbitraryModeRadius_ < 0) {
      arbitraryModeRadius_ = gMax(mag(faceVertices - centerCoord_));
    }
    length2_ = arbitraryModeRadius_; 
    distBetweenWalls = 2 * length2_;
    if (!arbitraryModeCenterIsInDict_) {
      arbitraryModeCenter_ = centerCoord_;
    }
    Info << "\nArbitrary mode with radius = " <<  length2_ <<  " and center " << arbitraryModeCenter_ << endl << endl;
  } else {
    FatalErrorInFunction << "Invalid geometryMode_, choose between: {annulus, disk, rectangle, channel, arbitrary}" << exit(FatalError);
  }

  const surfaceScalarField& phi = this->internalField().mesh().lookupObject<surfaceScalarField>("phi");
  if (phi.dimensions() == dimMass/dimTime)
  { 
    // if phi is mass flow rate
    const fvPatchField<scalar>& mu = patch().lookupPatchField<volScalarField, scalar>("mu");
    const fvPatchField<scalar>& rho = patch().lookupPatchField<volScalarField, scalar>("rho");
    avNu_ = gAverage(mu/rho);
  } else {
    // if phi is volumetric flow rate
    const fvPatchField<scalar>& nu = patch().lookupPatchField<volScalarField, scalar>("nu");
    avNu_ = gAverage(nu);
  }

  //calculate the minimum allowed friction velocity Utau_ from Reynolds number Re
  if (!UtauIsInDict_) {
    //instead of simply updateUmean(), integrate over 10 s:
    const scalarField rho(rhoOnPatch());
    Umean_ = inFlow_->integral(db().time().value(), db().time().value() + 10.0)/10.0
          / (massFlowMode_ ? gSum(rho*patch().magSf()) : 1);

    scalar Re = max(Umean_, 1e-4) * distBetweenWalls / avNu_;
    scalar sqf; //square root of the Darcy friction factor
    Info << "Re = " << Re << endl;
    if (Re < 3000) {
      sqf = 8/sqrt(Re);
    } else {
      sqf = 0.2;
      for (int i = 0; i < 10; ++i) //Newton's method with a few iterations
      {//smooth pipe expression
        sqf -= ( 1.93*log10(Re*sqf)*sqf*sqf -  sqf - 0.537*sqf*sqf ) / ( 1.93*sqf/log(10.0) + 1 );
        /* Info << "Residual = " << 1.93*log10(Re*sqf) -  1/sqf - 0.537 << endl; */
      }
    }
    Utau_ = max(Umean_, 1e-4) * sqf / sqrt(8.0);
  }

  scalar Utau_min = rplus_list[URL_len-1] * 2* avNu_ / distBetweenWalls / r_list[URL_len-1];

  //use the user defined Utau_ if it is bigger
  Utau_ = max(Utau_ , Utau_min);
  Sr_ = Utau_ / avNu_;
  Info << "Utau_ = " << Utau_ << endl;
  Info << "Scaling factor for yplus = Sr_/(Utau_min/avNu_) = " << Sr_/(Utau_min/avNu_) << endl << endl;

  //- Characteristic length scale
  delta_ = distBetweenWalls;
  if (geometryMode_ == "disk") {
    delta_ /= 2.0;
  }
}

void Foam::syntheticTurbulenceInletFvPatchVectorField::getRandomPositionAndIntensities(Foam::vector& pos, Foam::vector& intensities)
{
  //generate the random position and the intensity for one new eddy

  //random number between 0 and 1:
  std::uniform_real_distribution<> dist1(0.0, 1.0);
  //random number that is either -1 or 1:
  std::bernoulli_distribution dist2(0.5); // 50% chance for true or false  

  scalar rp1 = dist1(randGen_);
  scalar rp2 = dist1(randGen_);
  intensities[0] = dist2(randGen_) ? 1 : -1;
  intensities[1] = dist2(randGen_) ? 1 : -1;
  intensities[2] = dist2(randGen_) ? 1 : -1;

  if (geometryMode_ == "annulus") {
    rp1 *= Foam::constant::mathematical::twoPi; //theta
    // map positions between 0 and 1 to between length1_ and length2_, nonuniformly: 
    rp2 = sqrt(rp2*(length2_*length2_ - length1_*length1_) + length1_*length1_); //radius
    pos = centerCoord_ + (rp2*cos(rp1))*inplaneVec1_ + (rp2*sin(rp1))*inplaneVec2_;

  } else if (geometryMode_ == "disk") {
    rp1 *= Foam::constant::mathematical::twoPi; //theta
    // map positions between 0 and 1 to between 0 to length2_, nonuniformly: 
    rp2 = length2_*sqrt(rp2); //radius
    pos = centerCoord_ + (rp2*cos(rp1))*inplaneVec1_ + (rp2*sin(rp1))*inplaneVec2_;

  } else if (geometryMode_ == "rectangle" || geometryMode_ == "channel") {
    // map positions between 0 and 1 to between -length1_/2 to length1_/2 and -length2_/2 to length2_/2, uniformly: 
    rp1 = (rp1 - 0.5)*length1_;
    rp2 = (rp2 - 0.5)*length2_;
    pos = centerCoord_ + rp1*inplaneVec1_ + rp2*inplaneVec2_;
  } else if (geometryMode_ == "arbitrary") {
    rp1 *= Foam::constant::mathematical::twoPi; //theta
    // map positions between 0 and 1 to between 0 to length2_, nonuniformly: 
    rp2 = length2_*sqrt(rp2); //radius
    pos = arbitraryModeCenter_ + (rp2*cos(rp1))*inplaneVec1_ + (rp2*sin(rp1))*inplaneVec2_;
  }
}

void Foam::syntheticTurbulenceInletFvPatchVectorField::getRandSigmaFromPos(const vector& pos, symmTensor& R, scalar& sigmax) const
{
  //get R (Reynolds stress tensor) and Sigma (dimension of the eddy) from a specific position


  scalar r = mag(pos - centerCoord_); //position from the wall
  // u - along flow, v - away from wall, w - tangential to wall
  vector u(patchNormal_);
  vector v((centerCoord_ - pos)/r); //points from pos to centerCoord_. 

  if (geometryMode_ == "annulus") {
    if((r - length1_) < (length2_ - r)){
      r = r - length1_;
      v = -v;
    } else {
      r = length2_ - r;
    }

  } else if (geometryMode_ == "disk") {
      r = length2_ - r;

  } else if (geometryMode_ == "rectangle") {
    scalar l1 = (pos - centerCoord_) & inplaneVec1_;
    scalar l2 = (pos - centerCoord_) & inplaneVec2_;
    if ((length1_/2.0 - mag(l1)) > (length2_/2.0 - mag(l2))) {
      r = length2_/2.0 - mag(l2);
      if (l2 > 0) {
        v = -inplaneVec2_;
      } else {  
        v = inplaneVec2_;
      }
    } else {
      r = length1_/2.0 - mag(l1);
      if (l1 > 0) {
        v = -inplaneVec1_;
      } else {  
        v = inplaneVec1_;
      }
    }
  } else if (geometryMode_ == "channel") {
    if (wallLongOrShort_ == "long") {
      scalar l1 = (pos - centerCoord_) & inplaneVec1_;
      r = length1_/2.0 - mag(l1);
      if (l1 > 0) {
        v = -inplaneVec1_;
      } else {  
        v = inplaneVec1_;
      }
    } else if (wallLongOrShort_ == "short"){
      scalar l2 = (pos - centerCoord_) & inplaneVec2_;
      r = length2_/2.0 - mag(l2);
      if (l2 > 0) {
        v = -inplaneVec2_;
      } else {  
        v = inplaneVec2_;
      }
    } else {
      FatalErrorInFunction << "Invalid \"wallLongOrShort_\", choose between: {long, short}, which indicates if the walls are on the long or the short edge. The other direction is assumed to have PBC." << exit(FatalError);
    }
  }
  r *= Sr_; // y to yplus

  vector w = u^v;
  const tensor T(u,v,w); // Transformation matrix, T, for the Reynolds stress tensor

  scalar L;
  if (geometryMode_ == "arbitrary") {
    L = L_list[URL_len -1];
    R = symmTensor(Ruu_list[URL_len-1], Ruv_list[URL_len-1], Ruw_list[URL_len-1],
                                        Rvv_list[URL_len-1], Rvw_list[URL_len-1],
                                                             Rww_list[URL_len-1]);
  } else {
    const int j = std::upper_bound(rplus_list, rplus_list + URL_len, r) - rplus_list;
    if (j == 0) {
      L = L_list[0];
      R = symmTensor(Ruu_list[0], Ruv_list[0], Ruw_list[0],
                                      Rvv_list[0], Rvw_list[0],
                                                   Rww_list[0]);
    } else if (j == URL_len) {
      L = L_list[URL_len -1];
      R = symmTensor(Ruu_list[URL_len-1], Ruv_list[URL_len-1], Ruw_list[URL_len-1],
                                      Rvv_list[URL_len-1], Rvw_list[URL_len-1],
                                                   Rww_list[URL_len-1]);
    } else {
      //interpolate between the values at j and j-1
      L = -(( L_list[j-1]*(rplus_list[j]-r) + L_list[j]*(r-rplus_list[j-1]) ) / (rplus_list[j]-rplus_list[j-1]));
      R = (symmTensor(Ruu_list[j], Ruv_list[j], Ruw_list[j],
                          Rvv_list[j], Rvw_list[j],
                          Rww_list[j])*(rplus_list[j]-r) +
               symmTensor(Ruu_list[j-1], Ruv_list[j-1], Ruw_list[j-1],
                         Rvv_list[j-1], Rvw_list[j-1],
                         Rww_list[j-1])*(r-rplus_list[j-1]))/(rplus_list[j]-rplus_list[j-1]);
    }
  }

  // Transform R to local coordinate system
  const tensor Rtmp = T.T() & tensor(R) & T;
  R = symmTensor(Rtmp.xx(), Rtmp.xy(), Rtmp.xz(), Rtmp.yy(), Rtmp.yz(), Rtmp.zz());
  scalar Uscale = max(Umean_, 1e-4);
  R *= Uscale * Uscale / (Uavg_ * Uavg_) ;
  L /= -Lref_ * (Uscale / Uavg_ / avNu_);

  sigmax = min(mag(L), kappa_*delta_);
  sigmax = max(sigmax, minSigma_);
}

Foam::eddy Foam::syntheticTurbulenceInletFvPatchVectorField::oneEddy(scalar depth)
{
  vector pos;
  vector intensities;
  symmTensor R;
  scalar sigmax;
  getRandomPositionAndIntensities(pos, intensities);
  getRandSigmaFromPos(pos, R, sigmax);

  eddy e(pos,depth,sigmax,R,intensities);
  if (debug){ Pout << "Eddy pos: " << pos << endl;}
  return e;
}

void Foam::syntheticTurbulenceInletFvPatchVectorField::initializeUProfile()
{
  //initialize U and use L to compute maxSigmaX.


  const scalarField& magSf = patch().magSf(); // magnitude of face area vector
  const vectorField& patchPos = patch().Cf();

  //- Integral-length scale field
  scalarField L(size(), Zero);

  if (geometryMode_ == "arbitrary") {
    forAll(patchPos, i)
    {
      U_[i] = patchNormal_*1.0;
      L[i] = L_list[URL_len -1];
    }
  } else {
    forAll(patchPos, i)
    {
      scalar r = mag(patchPos[i] - centerCoord_); //distance from wall

      if (geometryMode_ == "annulus") {
        if((r - length1_) < (length2_ - r)){
          r = r - length1_;
        } else {
          r = length2_ - r;
        }

      } else if (geometryMode_ == "disk") {
          r = length2_ - r;

      } else if (geometryMode_ == "rectangle") {
        scalar l1 = (patchPos[i] - centerCoord_) & inplaneVec1_;
        scalar l2 = (patchPos[i] - centerCoord_) & inplaneVec2_;
        if ((length1_/2.0 - mag(l1)) > (length2_/2.0 - mag(l2))) {
          r = length2_/2.0 - mag(l2);
        } else {
          r = length1_/2.0 - mag(l1);
        }

      } else if (geometryMode_ == "channel") {
        if (wallLongOrShort_ == "long") {
          scalar l1 = (patchPos[i] - centerCoord_) & inplaneVec1_;
          r = length1_/2.0 - mag(l1);
        } else if (wallLongOrShort_ == "short"){
          scalar l2 = (patchPos[i] - centerCoord_) & inplaneVec2_;
          r = length2_/2.0 - mag(l2);
        } else {
          FatalErrorInFunction << "Invalid \"wallLongOrShort_\", choose between: {long, short}, which indicates if the walls are on the long or the short edge. The other direction is assumed to have PBC." << exit(FatalError);
        }
      }

      r *= Sr_; // y to yplus

      const int j = std::upper_bound(rplus_list, rplus_list + URL_len, r) - rplus_list;
      if (j == 0) {
        U_[i] = patchNormal_*0.0;
        L[i] = L_list[0];
      } else if (j >= URL_len) {
        U_[i] = patchNormal_*U_list[URL_len -1];
        L[i] = L_list[URL_len -1];
      } else {
        //interpolate between the values at j and j-1
        U_[i] = patchNormal_*(( U_list[j-1]*(rplus_list[j]-r) + U_list[j]*(r-rplus_list[j-1]) ) / (rplus_list[j]-rplus_list[j-1]));
        L[i] = -(( L_list[j-1]*(rplus_list[j]-r) + L_list[j]*(r-rplus_list[j-1]) ) / (rplus_list[j]-rplus_list[j-1]));
      }

    }
  }
  // Normalize U
  Uavg_ = gSum(U_ & -patch().Sf())/ gSum(magSf);
  U_ *= 1.0/Uavg_; //scaled to correct mean later in updateCoeff
  updateUmean();
  L /= -Lref_ * (max(Umean_, 1e-4) / Uavg_ / avNu_);

  //add velocity in the angular direction
  if (geometryMode_ == "annulus" || geometryMode_ == "disk") {
    scalar Rm = (length1_+length2_)/2; //mean radius

    const scalarField rho(rhoOnPatch());
    const scalar G1 = gSum(rho * magSqr(U_) * magSf * magSqr(patchPos - centerCoord_));
    const scalar G2 = gSum(rho * magSqr(U_) * magSf);

    scalar swirlFrac = swirlNr_ * Rm * G2 / G1;

    forAll(patchPos, i)
    {
      scalar rFromC = mag(patchPos[i] - centerCoord_); //distance from center

      vector wSwirl =  patchNormal_ ^ (patchPos[i] - centerCoord_);
      wSwirl /= mag(wSwirl);
      // Utheta = swirlFrac *radius* axial-speed
      U_[i] += mag(U_[i]) * swirlFrac * rFromC * wSwirl;
    }

  }

  
  //for testing:
  /* Pout << "Testing stuff; remember to uncomment Upatch afterwards!!"<<endl; */
  /* U_ = L * vector(1.0,0.0,0.0); //kommenterte ut upatch nede også */
  /*
  symmTensor R;
  scalar sigmahei;
  forAll(patchPos, i)
  {
    getRandSigmaFromPos(patchPos[i], R, sigmahei);
    U_[i] =  vector(R.xx(),R.xy(),R.xz()); //kommenterte ut upatch nede også
    //U_[i] =  vector(R.yy(),R.yz(),R.zz()); //kommenterte ut upatch nede også
  }
  */

  // From DFSEM:
  // (PCF:Eq. 14)
  const scalarField cellDx(max(Foam::sqrt(magSf), 2/patch().deltaCoeffs())); //longest dimension of each cell
  minSigma_ = gMax(nCellPerEddy_*cellDx); //make sure all eddies_ are longer than determined by minimum nr cells per eddy

  scalarField sigmax_(size(), Zero);

  // Inialise eddy box extents
  forAll(*this, faceI)
  {
      scalar& s = sigmax_[faceI];

      // Average length scale (SST:Eq. 24)
      // Personal communication regarding (PCR:Eq. 14)
      //  - the min operator in Eq. 14 is a typo, and should be a max operator
      s = min(mag(L[faceI]), kappa_*delta_);
      s = max(s, minSigma_); 
  }

  // Maximum extent across all processors
  maxSigmaX_ = gMax(sigmax_);

  // Eddy box volume
  if (geometryMode_ == "arbitrary") {
    //to make the density on the turbulent part of the patch the same:
    v0_ = 2*Foam::constant::mathematical::pi*arbitraryModeRadius_*arbitraryModeRadius_*maxSigmaX_;
  } else {
    v0_ = 2*gSum(magSf)*maxSigmaX_;
  }

  if (debug)
  {
      Info<< "Patch: " << patch().poly().name() << " eddy box:" << nl
          << "    volume    : " << v0_ << nl
          << "    maxSigmaX : " << maxSigmaX_ << nl
          << endl;
  }
}

void Foam::syntheticTurbulenceInletFvPatchVectorField::initializeEddies()
{
    /* #include "randomNumbers/initDepth.H" */

    //random number between 0 and 1:
    std::uniform_real_distribution<> dist3(-maxSigmaX_, maxSigmaX_);

    DynamicList<eddy> eddiesLocal(size());

    // Initialise eddy properties
    scalar sumVolEddy = 0;

    scalar depthVal;
    while (sumVolEddy/v0_ < d_)
    {
        depthVal = dist3(randGen_);

        eddy e = oneEddy(depthVal);
        eddiesLocal.append(e);
        sumVolEddy += e.volume();
    }
    eddies_.transfer(eddiesLocal);

    nEddy_ = eddies_.size();
    Info << "Initialized patch " << patch().name() << " with " << nEddy_<< " eddies. "<<endl<<endl;
}

void Foam::syntheticTurbulenceInletFvPatchVectorField::convectEddies()
{
    const scalar deltaT = db().time().deltaTValue();

    updateUmean();
    forAll(eddies_, eddyI)
    {
        eddy& e = eddies_[eddyI];
        e.move(deltaT*Umean_);

        const scalar position0 = e.x();

        // Check to see if eddy has exited downstream box plane
        if (position0 > maxSigmaX_)
        {
            // DFSEM moves the eddy by only maxSigmaX_, not 2*maxSigmaX_. This means that although the original eddy distribution is in (-maxSigmaX_, maxSigmaX_), it will evolve into (0, maxSigmaX_).
            e = oneEddy(position0 - 2*floor(position0/maxSigmaX_)*maxSigmaX_);
        }
    }

}

Foam::vector Foam::syntheticTurbulenceInletFvPatchVectorField::uPrimeEddy(const vector& point) const
{
    vector uPrime(Zero);

    forAll(eddies_, k)
    {
        const eddy& e = eddies_[k];
        uPrime += e.uPrime(point, patchNormal_);

        // Adding ghost eddies_ in order to have correct R near walls
        if (geometryMode_ == "channel" || geometryMode_ == "rectangle") {

            if (mag(point + length2_*inplaneVec2_ - e.position0()) < maxSigmaX_) {
              uPrime += e.uPrime(point + length2_*inplaneVec2_, patchNormal_);
            } else if (mag(point - length2_*inplaneVec2_ - e.position0()) < maxSigmaX_) {
              uPrime += e.uPrime(point - length2_*inplaneVec2_, patchNormal_);
            } 
            
            if (mag(point + length1_*inplaneVec1_ - e.position0()) < maxSigmaX_) {
              uPrime += e.uPrime(point + length1_*inplaneVec1_, patchNormal_);
            } else if (mag(point - length1_*inplaneVec1_ - e.position0()) < maxSigmaX_) {
              uPrime += e.uPrime(point - length1_*inplaneVec1_, patchNormal_);
            } 

        }
        else if (geometryMode_ == "disk") {
            const scalar magPC = mag(point - centerCoord_);
            if (magPC > SMALL) {
                vector uv((point - centerCoord_)/magPC); // unit vector pointing from center to point
                if (mag(point - 2*length2_*uv - e.position0()) < maxSigmaX_) {
                  uPrime += e.uPrime(point - 2*length2_*uv, patchNormal_);
                }
            }
        }
        else if (geometryMode_ == "annulus") {
            const scalar magPC = mag(point - centerCoord_);
            if (magPC > SMALL) {
                vector uv((point - centerCoord_)/magPC); // unit vector pointing from center to point
                if (mag(point - 2*length1_*uv - e.position0()) < maxSigmaX_) {
                  uPrime += e.uPrime(point - 2*length1_*uv, patchNormal_);
                } else if (mag(point - 2*length2_*uv - e.position0()) < maxSigmaX_) {
                  uPrime += e.uPrime(point - 2*length2_*uv, patchNormal_);
                }
            }
        }

    }
    return uPrime;
}

void Foam::syntheticTurbulenceInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    if (notInitialized_) {
      initializeGeometry();
      initializeUProfile();
      initializeEddies();
      notInitialized_ = false;
    }

    convectEddies();

    // scale with mean, then add eddies_, then rescale with mean again:
    
    updateUmean();
    vectorField Upatch = U_*Umean_;

    const scalar c2 = scale_*Foam::sqrt((315.0/16.0)*v0_/scalar(nEddy_));

    const pointField& Cf = patch().Cf();

    forAll(Upatch, faceI)
    {
        Upatch[faceI] += c2*uPrimeEddy(Cf[faceI]);
    }
    // Re-scale to ensure correct flow rate
    scalar fCorr = 0.0;
    if (massFlowMode_) {
      //fCorr is inverse mass flow rate here
      const scalarField rho(rhoOnPatch());
      fCorr = 1.0/gSum(rho*(Upatch & -patch().Sf())); //total mass flow rate
    } else {
      //fCorr is area average U normal to the patch
      fCorr = 
        gSum(patch().magSf())
       /gSum(Upatch & -patch().Sf()); 
    }

    Upatch *= inFlow_->value(db().time().value())*fCorr; //* if massFlowMode_: mass flow rate / sum mass flow over patch, if not then Umean_/ U_average

    vectorField::operator=(Upatch);

    fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::syntheticTurbulenceInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry
    (
        os,
        this->db().time().userUnits(),
        this->internalField().dimensions(),
        inFlow_()
    );
    writeEntry(os, "massFlowMode", massFlowMode_);
    writeEntry(os, "Utau", Utau_);
    writeEntry(os, "geometryMode", geometryMode_);
    if (geometryMode_ == "channel") {
      writeEntry(os, "wallLongOrShort", wallLongOrShort_);
    }
    if (geometryMode_ == "arbitrary") {
      writeEntry(os, "arbitraryModeRadius", arbitraryModeRadius_);
      if (arbitraryModeCenterIsInDict_) {
        writeEntry(os, "arbitraryModeCenter", arbitraryModeCenter_);
      }
    }
    writeEntry(os, "swirlNr", swirlNr_);
    writeEntry(os, "d", d_);
    writeEntry(os, "Lref", Lref_);
    writeEntry(os, "scale", scale_);
    writeEntry(os, "nCellPerEddy", nCellPerEddy_);
    writeEntry(os, "rhoInf", rhoInf_);
    writeEntry(os, "seed", seed_);
}

 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        syntheticTurbulenceInletFvPatchVectorField
    );
}
 
// ************************************************************************* //
