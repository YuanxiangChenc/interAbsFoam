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
const dimensionedScalar& Pr1 = mixture.Pr1();
const dimensionedScalar& Pr2 = mixture.Pr2();

// Thermal conductivity kappa is NOT stored directly in the mixture
// class (there is no mixture.kappa1() accessor). Derive it from
// the definition of the Prandtl number:
//
//   Pr = mu * cp / kappa   →   kappa = mu * cp / Pr = rho * nu * cp / Pr
//
// nuModel->nu() returns a volScalarField (spatially varying because
// viscosity models like Arrhenius depend on T). For Newtonian phases
// nu is uniform, so we take the domain average to reduce it to a
// scalar for TEqn.H's face-field construction below.
//
// If a spatially varying kappa is ever needed, build kappa1_field
// and kappa2_field as volScalarFields directly using nuModel->nu().
Info<< "Deriving kappa1, kappa2 from rho * nu * cp / Pr\n" << endl;

const dimensionedScalar nu1_avg
(
    "nu1_avg",
    dimensionSet(0, 2, -1, 0, 0, 0, 0),
    gAverage(mixture.nuModel1().nu()())
);
const dimensionedScalar nu2_avg
(
    "nu2_avg",
    dimensionSet(0, 2, -1, 0, 0, 0, 0),
    gAverage(mixture.nuModel2().nu()())
);

const dimensionedScalar kappa1
(
    "kappa1",
    dimensionSet(1, 1, -3, -1, 0, 0, 0),
    rho1.value() * nu1_avg.value() * cp1.value() / Pr1.value()
);
const dimensionedScalar kappa2
(
    "kappa2",
    dimensionSet(1, 1, -3, -1, 0, 0, 0),
    rho2.value() * nu2_avg.value() * cp2.value() / Pr2.value()
);

Info<< "  kappa1 = " << kappa1.value() << " W/(m K)" << nl
    << "  kappa2 = " << kappa2.value() << " W/(m K)" << endl;


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

// ============================================================
// createFields.H additions
//
// Add these blocks to createFields.H so that the interface
// geometry fields updated by updateInterfaceGeometry.H are
// persistent (they live between timesteps and inner iterations).
//
// Place these near the existing gradAlpha / nHat declarations
// (which stay as-is — those are still cell-centred and used
// by UEqn for Marangoni).
// ============================================================

const scalar gasAlphaThresh
(
    pimple.dict().lookupOrDefault<scalar>("gasAlphaThresh", 0.5)
);



// ------------------------------------------------------------
// Sfhat — face-normal unit vector.  MESH-DEPENDENT ONLY, not
// alpha-dependent.  Build once and never update.
// ------------------------------------------------------------
surfaceVectorField Sfhat
(
    IOobject
    (
        "Sfhat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh.Sf() / mesh.magSf()
);


// ------------------------------------------------------------
// gradAlphaf — face-native interface gradient
// nHatf     — face-native unit interface normal
// nDotSf    — nHatf . Sfhat  [-]
// Updated by updateInterfaceGeometry.H once per timestep.
// ------------------------------------------------------------
surfaceVectorField gradAlphaf
(
    IOobject
    (
        "gradAlphaf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(gradAlpha)         // gradAlpha already exists above
);

surfaceVectorField nHatf
(
    IOobject
    (
        "nHatf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    gradAlphaf
  / ( mag(gradAlphaf)
    + dimensionedScalar("small", dimless/dimLength, 1e-8)
    )
);

surfaceScalarField nDotSf
(
    IOobject
    (
        "nDotSf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    nHatf & Sfhat
);


// ------------------------------------------------------------
// isInterfaceFace — 1 where owner-neighbour alphas straddle 0.5
// D_masked        — alpha*D1 with interface faces zeroed
// Updated by updateInterfaceGeometry.H once per timestep.
// ------------------------------------------------------------
surfaceScalarField isInterfaceFace
(
    IOobject
    (
        "isInterfaceFace",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

surfaceScalarField D_masked
(
    IOobject
    (
        "D_masked",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(alpha1) * D1        // initial value
);


// ------------------------------------------------------------
// Q_abs — absorption heat source [W/m3]
// Computed each PIMPLE iteration by computeQabs.H (before TEqn),
// independently of m_abs (which is only produced later in cEqn.H).
// Persistent so TEqn.H and diagnostics can read it.
// ------------------------------------------------------------
volScalarField Q_abs
(
    IOobject
    (
        "Q_abs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
);


//------------------------------------------------------------
// Persistent Marangoni and buoyancy force fields.
//
// These are declared here (not in UEqn.H) for two reasons:
//   1. They must exist in timestep 0/ so paraFoam picks them up
//      as available volume fields.  Locally-constructed IOobjects
//      in UEqn.H are written to disk, but paraFoam builds its
//      field list from timestep 0 and skips fields that aren't
//      there.
//   2. UEqn.H is called every PIMPLE inner iteration — declaring
//      them there means constructor/destructor pair per iteration,
//      which is wasteful and can leave stale references.
//
// UEqn.H assigns to these fields; runTime.write() handles the
// AUTO_WRITE output automatically.
// ------------------------------------------------------------
volVectorField FMarangoni
(
    IOobject
    (
        "FMarangoni", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimForce/dimVolume, vector::zero)
);
 
volVectorField Fsvt
(
    IOobject
    (
        "Fsvt", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimForce/dimVolume, vector::zero)
);
 
volVectorField Fsvc
(
    IOobject
    (
        "Fsvc", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimForce/dimVolume, vector::zero)
);
 
volVectorField Fbuoy
(
    IOobject
    (
        "Fbuoy", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimForce/dimVolume, vector::zero)
);

// ------------------------------------------------------------
// m_abs tracking scalars for inner-loop residual.
//
// mAbs_T:  the m_abs value that TEqn used for its Q_abs source.
//          Sampled after computeMabs #1 (cSat = f(T_prev)),
//          BEFORE the T solve.
//
// mAbs_c:  the m_abs value that cEqn will use for its source.
//          Sampled after computeMabs #2 (cSat = f(T_new)),
//          BEFORE the c solve.
//
// Residual = |mAbs_c - mAbs_T| / max(mAbs_c, eps)
//
// Physical meaning: when the T-solve barely changes cSat, the
// m_abs computed with T_prev matches the m_abs computed with
// T_new, and the T-c coupling has converged.  
// ------------------------------------------------------------
scalar mAbs_T             = 0.0;
scalar mAbs_c             = 0.0;
scalar mAbs_residual_iter = 1.0;    // start at 1 = "not yet converged"
 

#include "createMRF.H"
#include "createFvOptions.H"

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
