Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
 IOobject
 (
  "p_rgh",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );

Info<< "Reading field U\n" << endl;
volVectorField U
(
 IOobject
 (
  "U",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );
  




#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());


const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();



IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Reading latent heat of fusion
dimensionedScalar L(transportProperties.lookup("L"));

// Reading melting temperature range
dimensionedScalar Ts(transportProperties.lookup("Ts"));
dimensionedScalar Tl(transportProperties.lookup("Tl"));

dimensionedScalar Tmelt
(
    "Tmelt",
    (Tl+Ts)/scalar(2.0)
);

// Reading large D'arcy-type source term constant DCl
dimensionedScalar DCl(transportProperties.lookup("DCl"));
dimensionedScalar DCs(transportProperties.lookup("DCs"));

dimensionedScalar rhos(transportProperties.lookup("rhos"));

dimensionedScalar cps(transportProperties.lookup("cps"));
dimensionedScalar cpl(transportProperties.lookup("cpl"));
dimensionedScalar cpg(transportProperties.lookup("cpg"));

dimensionedScalar lambdas(transportProperties.lookup("lambdas"));
dimensionedScalar lambdal(transportProperties.lookup("lambdal"));
dimensionedScalar lambdag(transportProperties.lookup("lambdag"));



// Need to store rho for ddt(rho, U)
volScalarField rho
(
 IOobject
 (
  "rho",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT
  ),
 alpha1*rho1 + alpha2*rho2,
 alpha1.boundaryField().types()
 );
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
 IOobject
 (
  "rhoPhi",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 fvc::interpolate(rho)*phi
 );


// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
 incompressible::turbulenceModel::New(U, phi, mixture)
 );


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
 IOobject
 (
  "p",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 p_rgh + rho*gh
 );

volScalarField T
(
 IOobject
 (
  "T",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
  mesh
);




volScalarField alphal
(
   IOobject
 (
  "alphal",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 max(min((T-Ts)/(Tl-Ts),scalar(1)),scalar(0))
);
alphal.write();
volScalarField alphas ("alphas", scalar(1)-alphal);
alphas.write();

volScalarField Tcorr("Tcorr", (Tl-Ts)*alphal + Ts);
//volScalarField alpha1fl("alpha1fl",alpha1*fl);


volScalarField DC
   (
       IOobject
       (
           "DC",
           runTime.timeName(),
           mesh
       ),
       DCl*pow(alpha2,2)/(pow(alpha1,3)+DCs)
   );

   volScalarField lambda
   (
     IOobject
       (
           "lambda",
           runTime.timeName(),
           mesh
       ),
       alpha1*lambdal + alpha2*lambdas
   );

/*   dimensionedScalar DT
(
transportProperties.lookup("DT")
);
 */
   volScalarField cp
   (
       IOobject
       (
           "cp",
           runTime.timeName(),
           mesh
       ),
       alpha1*cpl + alpha2*cps
   );
   cp.oldTime();

   surfaceScalarField phiCp
   (
       IOobject
       (
           "phi*cp",
           runTime.timeName(),
           mesh
       ),
       phi*fvc::interpolate(cp)
   );


   //Calculating enthalpy and corrected temperature
   Info<< "Calculating enthalpy and corrected temperature\n" << endl;

   volScalarField h
   (
       IOobject
       (
           "h",
           runTime.timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
       ),
       cp*(T-Tmelt) + alpha1*alphal*L
   );
   h.write();


   //Creating fields for gravity and pressure
   Info<< "Calculating field g.h\n" << endl;

 //  volScalarField gh("gh", g & mesh.C());
  // surfaceScalarField ghf("ghf", g & mesh.Cf());


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
 p,
 p_rgh,
 pimple.dict(),
 pRefCell,
 pRefValue
 );

if (p_rgh.needReference())
 {
     p += dimensionedScalar
         (
          "p",
          p.dimensions(),
          pRefValue - getRefCellValue(p, pRefCell)
          );
     p_rgh = p - rho*gh;
 }

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// MULES flux from previous time-step
surfaceScalarField alphaPhi
(
 IOobject
 (
  "alphaPhi",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
  ),
 phi*fvc::interpolate(alpha1)
 );

// MULES Correction
tmp<surfaceScalarField> talphaPhiCorr0;
