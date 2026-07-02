#include "createRDeltaT.H"

//----------------- ADD IO dict -----------------
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
//-------------------------------------------------

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


//----------------------------------------
// Concetration field
//----------------------------------------
Info<< "Reading field c\n" << endl;
volScalarField c
(
    IOobject
    (
        "c",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Creating field cif\n" << endl;
volScalarField cif
(
    IOobject
    (
        "cif",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE    // written explicitly by rv_cEqn/pn_cEqn only
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

//----------------------------------------
// Temperature field
//----------------------------------------
Info<< "Reading field T\n" << endl;
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
//--------------------------------------------

#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi, T, c);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();
const dimensionedScalar& D1 = mixture.D1();
const dimensionedScalar& D2 = mixture.D2();
const dimensionedScalar& cp1 = mixture.cp1();
const dimensionedScalar& cp2 = mixture.cp2();
const dimensionedScalar& kappa1 = mixture.kappa1();
const dimensionedScalar& kappa2 = mixture.kappa2();


//-----------------------------------------
Info << "Reading D, cSat from transportProperties\n" << endl;
dimensionedScalar cSat0("cSat0", dimensionSet(0, 0, 0, 0, 0, 0, 0), readScalar(transportProperties.lookup("cSat0")));

// ---- cSat as volScalarField (updated each step by updateCSat_T.H) ----------
Info<< "Creating volScalarField cSat\n" << endl;
volScalarField cSat
(
    IOobject
    (
        "cSat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    cSat0          // uniform initialisation from transportProperties
);
 
//-------------------------------------------
// Surface tension
//   sigma0   : reference value read from transportProperties (same key "sigma"
//              that immiscibleIncompressibleTwoPhaseMixture already uses)
//   sigma    : volScalarField updated in cEqn.H via Marangoni law
//              sigma = sigma0 + dsdT*(T - Tr) + dsdc*(c_eff - cr)
//-------------------------------------------
Info<< "Reading Marangoni / surface-tension parameters\n" << endl;
 
dimensionedScalar sigma0
(
    "sigma0", dimensionSet(1, 0, -2, 0, 0, 0, 0),
    readScalar(transportProperties.lookup("sigma0"))
);
 
// dσ/dT  [N/(m K)]  — typically negative for aqueous solutions
dimensionedScalar dsdT
(
    "dsdT", dimensionSet(1, 0, -2, -1, 0, 0, 0),
    readScalar(transportProperties.lookup("dsdT"))
);
 
// dσ/dc  [N/m per unit mass fraction]
dimensionedScalar dsdc
(
    "dsdc", dimensionSet(1, 0, -2, 0, 0, 0, 0),
    readScalar(transportProperties.lookup("dsdc"))
);
 
// Reference temperature and concentration for Marangoni / buoyancy
dimensionedScalar Tr
(
    "Tr", dimTemperature,
    readScalar(transportProperties.lookup("Tr"))
);
 
dimensionedScalar cr
(
    "cr", dimless,
    readScalar(transportProperties.lookup("cr"))
);

// sigma as volScalarField — initially uniform at sigma0, updated in cEqn.H
Info<< "Creating volScalarField sigma\n" << endl;
volScalarField sigma
(
    IOobject
    (
        "sigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    sigma0
);

//-------------------------------------------
// Boussinesq buoyancy parameters
//   betaT1/2 : thermal expansion coefficients  [1/K]
//   betaC1/2 : solutal expansion coefficients  [1/(unit c)]
//   rho1/2Ref: reference densities for Boussinesq  [kg/m³]
//-------------------------------------------
Info<< "Reading buoyancy parameters\n" << endl;
 
dimensionedScalar rho1Ref
(
    "rho1Ref", dimDensity,
    readScalar(transportProperties.lookup("rho1Ref"))
);
dimensionedScalar rho2Ref
(
    "rho2Ref", dimDensity,
    readScalar(transportProperties.lookup("rho2Ref"))
);
dimensionedScalar betaT1
(
    "betaT1", dimensionSet(0, 0, 0, -1, 0, 0, 0),
    readScalar(transportProperties.lookup("betaT1"))
);
dimensionedScalar betaT2
(
    "betaT2", dimensionSet(0, 0, 0, -1, 0, 0, 0),
    readScalar(transportProperties.lookup("betaT2"))
);
dimensionedScalar betaC1
(
    "betaC1", dimless,
    readScalar(transportProperties.lookup("betaC1"))
);
dimensionedScalar betaC2
(
    "betaC2", dimless,
    readScalar(transportProperties.lookup("betaC2"))
);
//-------------------------------------------
// Enthalpy of absorption [J/kg]
//-------------------------------------------
dimensionedScalar Ha
(
    "Ha", dimensionSet(0, 2, -2, 0, 0, 0, 0),
    readScalar(transportProperties.lookup("Ha"))
);
//------------------------------------------
// Number of inner energy loops (read from PIMPLE dict or default to 1)
label nEnergyLoops
(
    pimple.dict().lookupOrDefault<label>("nEnergyLoops", 1)
);
 
//-------------------------------------------
// Interface normal fields (persistent — updated in UEqn.H each step)
//   gradAlpha : volume-of-fluid gradient [1/m]
//   nHat      : unit interface normal (gas→liquid)
//   deltaN    : small stabilisation for nHat denominator [1/m]
//-------------------------------------------
Info<< "Creating interface-normal fields\n" << endl;
 
const dimensionedScalar deltaN
(
    "deltaN",
    dimless/dimLength,
    scalar(1e-8) / std::pow(average(mesh.V()).value(), 1.0/3.0)
);
 
volVectorField gradAlpha
(
    IOobject
    (
        "gradAlpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad(alpha1)
);
 
volVectorField nHat
(
    IOobject
    (
        "nHat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    gradAlpha / (mag(gradAlpha) + deltaN)
);
//-------------------------------------------
// Diffusivity face field
//-------------------------------------------
surfaceScalarField Df
(
    IOobject
    (
        "Df",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(alpha1) * D1 // + (scalar(1.0) - alphaf1) *D2
);

volScalarField m_abs
(
    IOobject
    (
        "m_abs",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0) 
);
//-------------------------------------------

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
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

//------------------- add thermal diffusivity kappa for moniter Fourier number--
volScalarField kappa
(
    IOobject
    (
	"kappa",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT
    ),
    fvc::average(mixture.kappaf())

);
//------------------------------------------------------------------------------



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

//////////////modification: cal rhocp rhocpphi ///////////////////////////////////
volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1*cp1 + alpha2*rho2*cp2,
    alpha1.boundaryField().types()
);
rhoCp.oldTime();        //what does it do? to calc first strp of ddt(rhoCp)?

// heat flux — initial construction only (alphaPhi10 not yet declared).
// Updated CONSISTENTLY every time step in cEqn.H / cl_cEqn.H using:
//   rhoCpPhi = alphaPhi10*(rho1*cp1-rho2*cp2) + phi*(rho2*cp2)
// to guarantee discrete conservation ddt(rhoCp)+div(rhoCpPhi)≈0.
surfaceScalarField rhoCpPhi
(
    IOobject
    (
        "rhoCpPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rhoCp)*phi
);
//////////////end/////////////////////////////////////////////////////////////////
typedef incompressibleInterPhaseTransportModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    > transportModelType;

autoPtr<transportModelType> turbulence
(
    new transportModelType
    (
        rho, U, phi, rhoPhi, mixture
    )
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

#include "createMRF.H"
#include "createFvOptions.H"



///////////// absorption model selector ////////////////////
// Only two production models are kept in interAbsFoam:
//   blocked       = face-based blocked-laplacian + implicit source (cEqn.H)
//   liquidCentred = cell-based Castro-style source + invadec       (cl_cEqn.H)
//
// bk_cEqn.H / rv_cEqn.H / pn_cEqn.H were exploratory variants with
// worse performance and are not wired into this solver.
word absorptionModel
(
    transportProperties.lookupOrDefault<word>
    (
        "absorptionModel",
        "blocked"
    )
);

Info<< "Absorption model: " << absorptionModel << nl
    << "  blocked       = face-based blocked-laplacian   (cEqn.H)"    << nl
    << "  liquidCentred = cell-based Castro src + invadec (cl_cEqn.H)" << nl
    << endl;

if
(
    absorptionModel != "blocked"
 && absorptionModel != "liquidCentred"
)
{
    FatalErrorInFunction
        << "Unknown absorptionModel: " << absorptionModel
        << "\nValid options: blocked  liquidCentred"
        << exit(FatalError);
}

///////////// initialise Mabs accumulator & output file ////////////////////
scalar Mabs_acc = 0.0;
label  absWriteInterval = 10;

fileName absDir(runTime.path()/"postProcessing"/"absed");
mkDir(absDir);

autoPtr<OFstream> absFilePtr
(
    new OFstream(absDir/"absorptionAccum.dat")
);

// Header — matches the 3 columns actually written by cEqn.H / cl_cEqn.H:
//   runTime.value(), Mabs_acc, Mabs_acc/dz
(*absFilePtr) << "# Time[s]  Mabs_acc[kg]  Mabs_acc_dz[kg/m]" << nl;
absFilePtr->flush();
