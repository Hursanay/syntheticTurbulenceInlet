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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::syntheticTurbulenceInletFvPatchVectorField::
syntheticTurbulenceInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    Utau(dict.lookupOrDefault<scalar>("Utau", 0.0)),
    UtauIsInDict(dict.found("Utau")),
    geometryMode(dict.lookup<word>("geometryMode")),
    wallLongOrShort(dict.lookupOrDefault<word>("wallLongOrShort","long")),
    arbitraryModeRadius(dict.lookupOrDefault<scalar>("arbitraryModeRadius",-1)),
    arbitraryModeCenter(dict.lookupOrDefault<vector>("arbitraryModeCenter",vector(0,0,0))),
    arbitraryModeCenterIsInDict(dict.found("arbitraryModeCenter")),
    U_(size(), Zero),
    massFlowMode(dict.lookupOrDefault<bool>("massFlowMode",false)),
    inFlow(
        Function1<scalar>::New
        (
            massFlowMode ? "massFlow" : "Umean",
            this->db().time().userUnits(),
            iF.dimensions(),
            dict
        )
    ),
    Umean(0),
    Uavg(0),
    swirlNr(dict.lookupOrDefault<scalar>("swirlNr", 0)),
    nEddy_(0),
    delta_(1),
    d_(dict.lookupOrDefault<scalar>("d", 1)),
    kappa_(0.41),
    Lref(dict.lookupOrDefault<scalar>("Lref", 1)),
    scale_(dict.lookupOrDefault<scalar>("scale", 1)),
    nCellPerEddy_(dict.lookupOrDefault<label>("nCellPerEddy", 1)),
    notInitialized_(true),
    seed_(dict.lookupOrDefault<label>("seed", 536)),
    randGen(seed_),
    rhoInf(dict.lookupOrDefault<scalar>("rhoInf",1))
{
    vectorField::operator=(U_);
}

Foam::syntheticTurbulenceInletFvPatchVectorField::
syntheticTurbulenceInletFvPatchVectorField
(
    const syntheticTurbulenceInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper, false),
    Utau(ptf.Utau),
    Sr(ptf.Sr),
    avNu(ptf.avNu),
    geometryMode(ptf.geometryMode),
    wallLongOrShort(ptf.wallLongOrShort),
    arbitraryModeRadius(ptf.arbitraryModeRadius),
    arbitraryModeCenter(ptf.arbitraryModeCenter),
    arbitraryModeCenterIsInDict(ptf.arbitraryModeCenterIsInDict),
    U_(ptf.U_),
    length1(ptf.length1),
    length2(ptf.length2),
    centerCoord(ptf.centerCoord),
    massFlowMode(ptf.massFlowMode),
    inFlow(ptf.inFlow, false),
    Umean(ptf.Umean),
    Uavg(ptf.Uavg),
    swirlNr(ptf.swirlNr),
    nEddy_(ptf.nEddy_),
    delta_(ptf.delta_),
    d_(ptf.d_),
    kappa_(ptf.kappa_),
    Lref(ptf.Lref),
    scale_(ptf.scale_),
    nCellPerEddy_(ptf.nCellPerEddy_),
    minSigma_(ptf.minSigma_),
    v0_(ptf.v0_),
    maxSigmaX_(ptf.maxSigmaX_),
    notInitialized_(ptf.notInitialized_),
    seed_(ptf.seed_),
    randGen(ptf.randGen),
    rhoInf(ptf.rhoInf)
{}


Foam::syntheticTurbulenceInletFvPatchVectorField::
syntheticTurbulenceInletFvPatchVectorField
(
    const syntheticTurbulenceInletFvPatchVectorField& mwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(mwvpvf, iF),
    Utau(mwvpvf.Utau),
    Sr(mwvpvf.Sr),
    avNu(mwvpvf.avNu),
    geometryMode(mwvpvf.geometryMode),
    wallLongOrShort(mwvpvf.wallLongOrShort),
    arbitraryModeRadius(mwvpvf.arbitraryModeRadius),
    arbitraryModeCenter(mwvpvf.arbitraryModeCenter),
    arbitraryModeCenterIsInDict(mwvpvf.arbitraryModeCenterIsInDict),
    U_(mwvpvf.U_),
    length1(mwvpvf.length1),
    length2(mwvpvf.length2),
    centerCoord(mwvpvf.centerCoord),
    massFlowMode(mwvpvf.massFlowMode),
    inFlow(mwvpvf.inFlow, false),
    Umean(mwvpvf.Umean),
    Uavg(mwvpvf.Uavg),
    swirlNr(mwvpvf.swirlNr),
    nEddy_(mwvpvf.nEddy_),
    delta_(mwvpvf.delta_),
    d_(mwvpvf.d_),
    kappa_(mwvpvf.kappa_),
    Lref(mwvpvf.Lref),
    scale_(mwvpvf.scale_),
    nCellPerEddy_(mwvpvf.nCellPerEddy_),
    minSigma_(mwvpvf.minSigma_),
    v0_(mwvpvf.v0_),
    maxSigmaX_(mwvpvf.maxSigmaX_),
    notInitialized_(mwvpvf.notInitialized_),
    seed_(mwvpvf.seed_),
    randGen(mwvpvf.randGen),
    rhoInf(mwvpvf.rhoInf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::syntheticTurbulenceInletFvPatchVectorField::updateUmean()
{
  const surfaceScalarField& phi = this->internalField().mesh().lookupObject<surfaceScalarField>("phi");
  if (phi.dimensions() == dimMass/dimTime)
  { 
    // if phi is mass flow rate
    const fvPatchField<scalar>& rho = patch().lookupPatchField<volScalarField, scalar>("rho");
    Umean = inFlow->value(db().time().value())/(massFlowMode ? gSum(rho*patch().magSf()) : 1);
  } else {
    // if phi is volumetric flow rate
    Umean = inFlow->value(db().time().value())/(massFlowMode ? gSum(rhoInf*patch().magSf()) : 1);
  }
  /* Info << "Umean = " << Umean<< endl; */
}

void Foam::syntheticTurbulenceInletFvPatchVectorField::initializeGeometry()
{
  centerCoord = gAverage(patch().Cf());
  /* Info << "Center coord: " << centerCoord << endl; */
  const pointField& faceVertices = patch().patch().localPoints();

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
  if (geometryMode == "annulus") 
  {
    length1 = gMin(mag(faceVertices - centerCoord)); // inner radius
    length2 = gMax(mag(faceVertices - centerCoord)); // outer radius
    Info << "\nDetected an annulus with radii = " <<  length1 << " and " << length2 << endl << endl;
    if (length1 < 1e-8) {
      WarningInFunction
          << "The inner radius is very smaller, are you sure it is not a \"geometryMode disk\" but an \"geometryMode annulus\"? " 
          << endl;
    }
    distBetweenWalls = length2 - length1;

  } else if (geometryMode == "disk") {
    length1 = 0.0;
    length2 = gMax(mag(faceVertices - centerCoord)); 
    Info << "\nDetected a disk with radius = " <<  length2 << endl << endl;
    distBetweenWalls = 2 * length2;
    Info << "distBetweenWalls = " <<distBetweenWalls<<endl;

  } else if (geometryMode == "rectangle" || geometryMode == "channel") {
    scalar D = gMax(mag(faceVertices - centerCoord));
    /* Info << "Diagonal: " << D << endl; */ 
    scalar rect2 = gMax((faceVertices - centerCoord) & inplaneVec2_);
    scalar rect3 = gMax((faceVertices - centerCoord) & inplaneVec1_);
    scalar rect1 = sqrt(D*D - rect2*rect2);
    scalar rect4 = sqrt(D*D - rect3*rect3);

    const scalar minCellDx(gMin(Foam::sqrt(patch().magSf())));
    vector rect5 = rect1 * inplaneVec1_ + rect2 * inplaneVec2_;
    //if the diagonal vector rect5 is outside the patch, flip it:  
    if (gMin(mag(faceVertices - centerCoord - rect5)) > minCellDx) {
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
    if (gMin(mag(faceVertices - centerCoord - rect6)) > minCellDx) {
      if (switchInt == 1) {
        FatalErrorInFunction << "Failed to find inplane vectors. This should not happen." << exit(FatalError);
      }
      rect6 = rect3 * inplaneVec1_ + (switchInt * rect4) * inplaneVec2_;
    }
    
    //rebaptize:
    inplaneVec1_ = rect6 - rect5;
    inplaneVec2_ = rect5 + rect6;
    length1 = mag(inplaneVec1_);
    length2 = mag(inplaneVec2_);
    if(length1 > 0.0 && length2 > 0.0)
    {
      //order so that length2 > length1:
      if (length1 > length2) {
        scalar temp1 = length2;
        length2 = length1;
        length1 = temp1;
        vector temp2 = inplaneVec2_;
        inplaneVec2_ = inplaneVec1_;
        inplaneVec1_ = temp2;
      }
      inplaneVec1_ /= length1;
      inplaneVec2_ /= length2;
      Info << "\nDetected a " << geometryMode << " with dimensions " <<  length1 << " and " << length2 << endl;// << endl;
      Info << "in directions " << inplaneVec1_ << " and " <<  inplaneVec2_ << endl << endl;
    }

    distBetweenWalls = length1;
    if (geometryMode == "channel" && wallLongOrShort == "short" ) {
      distBetweenWalls = length2;
    }
  } else if (geometryMode == "arbitrary") {
    length1 = 0.0;
    if (arbitraryModeRadius < 0) {
      arbitraryModeRadius = gMax(mag(faceVertices - centerCoord));
    }
    length2 = arbitraryModeRadius; 
    distBetweenWalls = 2 * length2;
    if (!arbitraryModeCenterIsInDict) {
      arbitraryModeCenter = centerCoord;
    }
    Info << "\nArbitrary mode with radius = " <<  length2 <<  " and center " << arbitraryModeCenter << endl << endl;
  } else {
    FatalErrorInFunction << "Invalid geometryMode, choose between: {annulus, disk, rectangle, channel, arbitrary}" << exit(FatalError);
  }

  const surfaceScalarField& phi = this->internalField().mesh().lookupObject<surfaceScalarField>("phi");
  if (phi.dimensions() == dimMass/dimTime)
  { 
    // if phi is mass flow rate
    const fvPatchField<scalar>& mu = patch().lookupPatchField<volScalarField, scalar>("mu");
    const fvPatchField<scalar>& rho = patch().lookupPatchField<volScalarField, scalar>("rho");
    avNu = gAverage(mu/rho);
  } else {
    // if phi is volumetric flow rate
    const fvPatchField<scalar>& nu = patch().lookupPatchField<volScalarField, scalar>("nu");
    avNu = gAverage(nu);
  }

  //calculate the minimum allowed friction velocity Utau from Reynolds number Re
  if (!UtauIsInDict) {
    //instead of simply updateUmean(), integrate over 10 s:
    const surfaceScalarField& phi = this->internalField().mesh().lookupObject<surfaceScalarField>("phi");
    if (phi.dimensions() == dimMass/dimTime)
    { 
      // if phi is mass flow rate
      const fvPatchField<scalar>& rho = patch().lookupPatchField<volScalarField, scalar>("rho");
      Umean = inFlow->integral(db().time().value(), db().time().value() + 10.0)/10.0/(massFlowMode ? gSum(rho*patch().magSf()) : 1);
    } else {
      // if phi is volumetric flow rate
      Umean = inFlow->integral(db().time().value(), db().time().value() + 10.0)/10.0/(massFlowMode ? gSum(rhoInf*patch().magSf()) : 1);
    }
    /* Info << "Umean for Re = " << Umean<< endl; */

    scalar Re = max(Umean, 1e-4) * distBetweenWalls / avNu;
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
    Utau = max(Umean, 1e-4) * sqf / sqrt(8.0);
  }

  #include "URL/URL_len.H"
  #include "URL/r_list.H"
  #include "URL/rplus_list.H"
  scalar Utau_min = rplus_list[URL_len-1] * 2* avNu / distBetweenWalls / r_list[URL_len-1];

  //use the user defined Utau if it is bigger
  Utau = max(Utau , Utau_min);
  Sr = Utau / avNu;
  Info << "Utau = " << Utau << endl;
  Info << "Scaling factor for yplus = Sr/(Utau_min/avNu) = " << Sr/(Utau_min/avNu) << endl << endl;

  //- Characteristic length scale
  delta_ = distBetweenWalls;
  if (geometryMode == "disk") {
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

  scalar rp1 = dist1(randGen);
  scalar rp2 = dist1(randGen);
  intensities[0] = dist2(randGen) ? 1 : -1;
  intensities[1] = dist2(randGen) ? 1 : -1;
  intensities[2] = dist2(randGen) ? 1 : -1;

  if (geometryMode == "annulus") {
    rp1 *= 6.28318530718; //theta
    // map positions between 0 and 1 to between length1 and length2, nonuniformly: 
    rp2 = sqrt(rp2*(length2*length2 - length1*length1) + length1*length1); //radius
    pos = centerCoord + (rp2*cos(rp1))*inplaneVec1_ + (rp2*sin(rp1))*inplaneVec2_;

  } else if (geometryMode == "disk") {
    rp1 *= 6.28318530718; //theta
    // map positions between 0 and 1 to between 0 to length2, nonuniformly: 
    rp2 = length2*sqrt(rp2); //radius
    pos = centerCoord + (rp2*cos(rp1))*inplaneVec1_ + (rp2*sin(rp1))*inplaneVec2_;

  } else if (geometryMode == "rectangle" || geometryMode == "channel") {
    // map positions between 0 and 1 to between -length1/2 to length1/2 and -length2/2 to length2/2, uniformly: 
    rp1 = (rp1 - 0.5)*length1;
    rp2 = (rp2 - 0.5)*length2;
    pos = centerCoord + rp1*inplaneVec1_ + rp2*inplaneVec2_;
  } else if (geometryMode == "arbitrary") {
    rp1 *= 6.28318530718; //theta
    // map positions between 0 and 1 to between 0 to length2, nonuniformly: 
    rp2 = length2*sqrt(rp2); //radius
    pos = arbitraryModeCenter + (rp2*cos(rp1))*inplaneVec1_ + (rp2*sin(rp1))*inplaneVec2_;
  }
}

void Foam::syntheticTurbulenceInletFvPatchVectorField::getRandSigmaFromPos(vector pos, symmTensor& R, scalar& sigmax)
{
  //get R (Reynolds stress tensor) and Sigma (dimension of the eddy) from a specific position

  #include "URL/URL_len.H"
  #include "URL/rplus_list.H"
  #include "URL/R.H" //- Reynolds stress tensor field
  #include "URL/L.H" //- Integral-length scale field

  scalar r = mag(pos - centerCoord); //position from the wall
  // u - along flow, v - away from wall, w - tangential to wall
  vector u(patchNormal_);
  vector v((centerCoord - pos)/r); //points from pos to centerCoord. 

  if (geometryMode == "annulus") {
    if((r - length1) < (length2 - r)){
      r = r - length1;
      v = -v;
    } else {
      r = length2 - r;
    }

  } else if (geometryMode == "disk") {
      r = length2 - r;

  } else if (geometryMode == "rectangle") {
    scalar l1 = (pos - centerCoord) & inplaneVec1_;
    scalar l2 = (pos - centerCoord) & inplaneVec2_;
    if ((length1/2.0 - mag(l1)) > (length2/2.0 - mag(l2))) {
      r = length2/2.0 - mag(l2);
      if (l2 > 0) {
        v = -inplaneVec2_;
      } else {  
        v = inplaneVec2_;
      }
    } else {
      r = length1/2.0 - mag(l1);
      if (l1 > 0) {
        v = -inplaneVec1_;
      } else {  
        v = inplaneVec1_;
      }
    }
  } else if (geometryMode == "channel") {
    if (wallLongOrShort == "long") {
      scalar l1 = (pos - centerCoord) & inplaneVec1_;
      r = length1/2.0 - mag(l1);
      if (l1 > 0) {
        v = -inplaneVec1_;
      } else {  
        v = inplaneVec1_;
      }
    } else if (wallLongOrShort == "short"){
      scalar l2 = (pos - centerCoord) & inplaneVec2_;
      r = length2/2.0 - mag(l2);
      if (l2 > 0) {
        v = -inplaneVec2_;
      } else {  
        v = inplaneVec2_;
      }
    } else {
      FatalErrorInFunction << "Invalid \"wallLongOrShort\", choose between: {long, short}, which indicates if the walls are on the long or the short edge. The other direction is assumed to have PBC." << exit(FatalError);
    }
  }
  r *= Sr; // y to yplus

  vector w = u^v;
  const tensor T(u,v,w); // Transformation matrix, T, for the Reynolds stress tensor

  scalar L;
  if (geometryMode == "arbitrary") {
    L = L_list[URL_len -1];
    R = symmTensor(Ruu_list[URL_len-1], Ruv_list[URL_len-1], Ruw_list[URL_len-1],
                                        Rvv_list[URL_len-1], Rvw_list[URL_len-1],
                                                             Rww_list[URL_len-1]);
  } else {
    int j = 0;
    while (j < URL_len && r > rplus_list[j]){
      j++;
    }
    if (j == 0) {
      L = L_list[0];
      R = symmTensor(Ruu_list[0], Ruv_list[0], Ruw_list[0],
                                      Rvv_list[0], Rvw_list[0],
                                                   Rww_list[0]);
    } else if (r > rplus_list[URL_len-1]) {
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
  updateUmean();
  scalar Uscale = max(Umean, 1e-4);
  R *= Uscale * Uscale / (Uavg * Uavg) ;
  L /= -Lref * (Uscale / Uavg / avNu);

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
  
  #include "URL/URL_len.H" 
  #include "URL/rplus_list.H" 
  #include "URL/U.H"
  #include "URL/L.H" 


  const scalarField& magSf = patch().magSf(); // magnitude of face area vector
  const vectorField& patchPos = patch().Cf();

  //- Integral-length scale field
  scalarField L(size(), Zero);

  if (geometryMode == "arbitrary") {
    forAll(patchPos, i)
    {
      U_[i] = patchNormal_*1.0;
      L[i] = L_list[URL_len -1];
    }
  } else {
    forAll(patchPos, i)
    {
      scalar r = mag(patchPos[i] - centerCoord); //distance from wall

      if (geometryMode == "annulus") {
        if((r - length1) < (length2 - r)){
          r = r - length1;
        } else {
          r = length2 - r;
        }

      } else if (geometryMode == "disk") {
          r = length2 - r;

      } else if (geometryMode == "rectangle") {
        scalar l1 = (patchPos[i] - centerCoord) & inplaneVec1_;
        scalar l2 = (patchPos[i] - centerCoord) & inplaneVec2_;
        if ((length1/2.0 - mag(l1)) > (length2/2.0 - mag(l2))) {
          r = length2/2.0 - mag(l2);
        } else {
          r = length1/2.0 - mag(l1);
        }

      } else if (geometryMode == "channel") {
        if (wallLongOrShort == "long") {
          scalar l1 = (patchPos[i] - centerCoord) & inplaneVec1_;
          r = length1/2.0 - mag(l1);
        } else if (wallLongOrShort == "short"){
          scalar l2 = (patchPos[i] - centerCoord) & inplaneVec2_;
          r = length2/2.0 - mag(l2);
        } else {
          FatalErrorInFunction << "Invalid \"wallLongOrShort\", choose between: {long, short}, which indicates if the walls are on the long or the short edge. The other direction is assumed to have PBC." << exit(FatalError);
        }
      }

      r *= Sr; // y to yplus

      int j = 0;
      while (j < URL_len && r > rplus_list[j])
      {
        j++;
      }
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
  Uavg = gSum(U_ & -patch().Sf())/ gSum(magSf);
  U_ *= 1.0/Uavg; //scaled to correct mean later in updateCoeff
  updateUmean();
  L /= -Lref * (max(Umean, 1e-4) / Uavg / avNu);

  //add velocity in the angular direction
  if (geometryMode == "annulus" || geometryMode == "disk") {
    scalar Rm = (length1+length2)/2; //mean radius

    const surfaceScalarField& phi = this->internalField().mesh().lookupObject<surfaceScalarField>("phi");
    scalar G1 = 0.0;
    scalar G2 = 0.0;
    if (phi.dimensions() == dimMass/dimTime)
    { 
      // if phi is mass flow rate
      const fvPatchField<scalar>& rho = patch().lookupPatchField<volScalarField, scalar>("rho");
      G1 = gSum(rho * magSqr(U_) * magSf * magSqr(patchPos - centerCoord));
      G2 = gSum(rho * magSqr(U_) * magSf);
    } else {
      // if phi is volumetric flow rate
      G1 = gSum(rhoInf * magSqr(U_) * magSf * magSqr(patchPos - centerCoord));
      G2 = gSum(rhoInf * magSqr(U_) * magSf);
    }

    scalar swirlFrac = swirlNr * Rm * G2 / G1; 

    forAll(patchPos, i)
    {
      scalar rFromC = mag(patchPos[i] - centerCoord); //distance from center

      vector wSwirl =  patchNormal_ ^ (patchPos[i] - centerCoord);
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
  minSigma_ = gMax(nCellPerEddy_*cellDx); //make sure all eddies are longer than determined by minimum nr cells per eddy

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
  if (geometryMode == "arbitrary") {
    //to make the density on the turbulent part of the patch the same:
    v0_ = 2*3.14*arbitraryModeRadius*arbitraryModeRadius*maxSigmaX_;
  } else {
    v0_ = 2*gSum(magSf)*maxSigmaX_;
  }

  if (debug)
  {
      Info<< "Patch: " << patch().patch().name() << " eddy box:" << nl
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

    DynamicList<eddy> eddies_(size());

    // Initialise eddy properties
    scalar sumVolEddy = 0;

    scalar depthVal;
    while (sumVolEddy/v0_ < d_)
    {
        depthVal = dist3(randGen);
        
        eddy e = oneEddy(depthVal);
        eddies_.append(e);
        sumVolEddy += e.volume();
    }
    eddies.transfer(eddies_);

    nEddy_ = eddies.size();
    Info << "Initialized patch " << patch().name() << " with " << nEddy_<< " eddies. "<<endl<<endl;
}

void Foam::syntheticTurbulenceInletFvPatchVectorField::convectEddies()
{
    const scalar deltaT = db().time().deltaTValue();

    updateUmean();
    forAll(eddies, eddyI)
    {
        eddy& e = eddies[eddyI];
        e.move(deltaT*Umean);

        const scalar position0 = e.x();

        // Check to see if eddy has exited downstream box plane
        if (position0 > maxSigmaX_)
        {
            // DFSEM moves the eddy by only maxSigmaX_, not 2*maxSigmaX_. This means that although the original eddy distribution is in (-maxSigmaX_, maxSigmaX_), it will evolve into (0, maxSigmaX_).
            e = oneEddy(position0 - 2*floor(position0/maxSigmaX_)*maxSigmaX_);
        }
    }

}

Foam::vector Foam::syntheticTurbulenceInletFvPatchVectorField::uPrimeEddy(vector point)
{
    vector uPrime(Zero);

    forAll(eddies, k)
    {
        const eddy& e = eddies[k];
        uPrime += e.uPrime(point, patchNormal_);

        // Adding ghost eddies in order to have correct R near walls
        if (geometryMode == "channel" or geometryMode == "rectangle") {

            if (mag(point + length2*inplaneVec2_ - e.position0()) < maxSigmaX_) {
              uPrime += e.uPrime(point + length2*inplaneVec2_, patchNormal_);
            } else if (mag(point - length2*inplaneVec2_ - e.position0()) < maxSigmaX_) {
              uPrime += e.uPrime(point - length2*inplaneVec2_, patchNormal_);
            } 
            
            if (mag(point + length1*inplaneVec1_ - e.position0()) < maxSigmaX_) {
              uPrime += e.uPrime(point + length1*inplaneVec1_, patchNormal_);
            } else if (mag(point - length1*inplaneVec1_ - e.position0()) < maxSigmaX_) {
              uPrime += e.uPrime(point - length1*inplaneVec1_, patchNormal_);
            } 

        }
        else if (geometryMode == "disk") {
            vector uv((point-centerCoord)/mag(point-centerCoord)); // unit vector pointing from center to point
            if (mag(point - 2*length2*uv - e.position0()) < maxSigmaX_) {
              uPrime += e.uPrime(point - 2*length2*uv, patchNormal_);
            }
        }
        else if (geometryMode == "annulus") {
            vector uv((point-centerCoord)/mag(point-centerCoord)); // unit vector pointing from center to point
            if (mag(point - 2*length1*uv - e.position0()) < maxSigmaX_) {
              uPrime += e.uPrime(point - 2*length1*uv, patchNormal_);
            } else if (mag(point - 2*length2*uv - e.position0()) < maxSigmaX_) {
              uPrime += e.uPrime(point - 2*length2*uv, patchNormal_);
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

    // scale with mean, then add eddies, then rescale with mean again:
    
    updateUmean();
    vectorField Upatch = U_*Umean;

    const scalar c2 = scale_*Foam::sqrt((315.0/16.0)*v0_/scalar(nEddy_));

    const pointField& Cf = patch().Cf();

    forAll(Upatch, faceI)
    {
        Upatch[faceI] += c2*uPrimeEddy(Cf[faceI]);
    }
    // Re-scale to ensure correct flow rate
    scalar fCorr = 0.0;
    if (massFlowMode) {
      //fCorr is inverse ass flow rate here
      const surfaceScalarField& phi = this->internalField().mesh().lookupObject<surfaceScalarField>("phi");
      if (phi.dimensions() == dimMass/dimTime)
      { 
        // if phi is mass flow rate
        const fvPatchField<scalar>& rho = patch().lookupPatchField<volScalarField, scalar>("rho");
        fCorr = 1.0/gSum(rho*(Upatch & -patch().Sf())); //total mass flow rate
      } else {
        // if phi is volumetric flow rate
        fCorr = 1.0/gSum(rhoInf*(Upatch & -patch().Sf())); //total mass flow rate
      }
    } else {
      //fCorr is area average U normal to the patch
      fCorr = 
        gSum(patch().magSf())
       /gSum(Upatch & -patch().Sf()); 
    }

    Upatch *= inFlow->value(db().time().value())*fCorr; //* if massFlowMode: mass flow rate / sum mass flow over patch, if not then Umean/ U_average

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
        inFlow()
    );
    writeEntry(os, "massFlowMode", massFlowMode);
    writeEntry(os, "Utau", Utau);
    writeEntry(os, "geometryMode", geometryMode);
    if (geometryMode == "channel") {
      writeEntry(os, "wallLongOrShort", wallLongOrShort);
    }
    if (geometryMode == "arbitrary") {
      writeEntry(os, "arbitraryModeRadius", arbitraryModeRadius);
      if (arbitraryModeCenterIsInDict) {
        writeEntry(os, "arbitraryModeCenter", arbitraryModeCenter);
      }
    }
    writeEntry(os, "swirlNr", swirlNr);
    writeEntry(os, "delta", delta_);
    writeEntry(os, "d", d_);
    writeEntry(os, "kappa", kappa_);
    writeEntry(os, "Lref", Lref);
    writeEntry(os, "scale", scale_);
    writeEntry(os, "nCellPerEddy", nCellPerEddy_);
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
