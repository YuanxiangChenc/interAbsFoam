// ============================================================
// cEqn.H  —  LiBr solution concentration equation
//
// Physics (LiBr absorption):
//   c     = LiBr mass fraction in solution            [-]
//   cSat  = equilibrium LiBr mass fraction (volScalarField, T-dependent)
//   c > cSat : concentrated solution absorbs water
//   driving force = (c_liq - cSat) > 0
//
// cSat is a volScalarField updated by updateCSat_T.H before this file.
// sigma is updated at the END of this file using:
//   sigma = sigma0 + dsdT*(T - Tr) + dsdc*(c_eff - cr)
// where c_eff blends c → cSat near the interface to avoid spurious
// Marangoni spikes from diffusion-blocked cells.
// ============================================================

alpha1.max(0.0);
alpha1.min(1.0);

// Face-interpolated cSat for gradient calculations
surfaceScalarField cSatf(fvc::interpolate(cSat));

const labelUList& own = mesh.owner();
const labelUList& nei = mesh.neighbour();

// ------------------------------------------------------------
// 1. Interface face mask
// ------------------------------------------------------------
surfaceScalarField isInterfaceFace
(
    IOobject
    (
        "isInterfaceFace",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

forAll(isInterfaceFace, facei)
{
    scalar alphaO = alpha1[own[facei]];
    scalar alphaN = alpha1[nei[facei]];
    if (   (alphaO > 0.5 && alphaN < 0.5)
        || (alphaO < 0.5 && alphaN > 0.5) )
    {
        isInterfaceFace[facei] = 1.0;
    }
}

// ------------------------------------------------------------
// 2. Masked diffusivity: alpha*D1 in bulk, zero at interface faces
// ------------------------------------------------------------
surfaceScalarField D_masked
(
    IOobject
    (
        "D_masked",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (scalar(1.0) - isInterfaceFace) * fvc::interpolate(alpha1) * D1
);

// ------------------------------------------------------------
// 3. Interface normal at faces
// ------------------------------------------------------------
surfaceVectorField gradAlphaf
(
    IOobject
    (
        "gradAlphaf",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(fvc::grad(alpha1))
);

surfaceVectorField nHatf
(
    IOobject
    (
        "nHatf",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    gradAlphaf
  / (mag(gradAlphaf) + dimensionedScalar("small", dimless/dimLength, 1e-8))
);

surfaceVectorField Sfhat
(
    IOobject
    (
        "Sfhat",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh.Sf() / mesh.magSf()
);

surfaceScalarField nDotSf
(
    IOobject
    (
        "nDotSf",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    nHatf & Sfhat
);

// ------------------------------------------------------------
// 4. One-sided dc/dn at interface faces  [1/m]
//    cSat is now a volScalarField — use per-cell values
// ------------------------------------------------------------
surfaceScalarField dc_f
(
    IOobject
    (
        "dc_f",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless/dimLength, 0.0)
);

forAll(dc_f, facei)
{
    if (isInterfaceFace[facei] > 0.5)
    {
        scalar alphaO = alpha1[own[facei]];
        scalar alphaN = alpha1[nei[facei]];
        scalar delta  = 1.0 / mesh.deltaCoeffs()[facei];   // [m]

        // Use local cSat at each cell
        scalar cSatO = cSat[own[facei]];
        scalar cSatN = cSat[nei[facei]];

        scalar cO = alphaO*c[own[facei]] + (scalar(1)-alphaO)*cSatO;
        scalar cN = alphaN*c[nei[facei]] + (scalar(1)-alphaN)*cSatN;

        if (alphaO >= alphaN)
            dc_f[facei] = max((cO - cN)/delta, scalar(0));
        else
            dc_f[facei] = max((cN - cO)/delta, scalar(0));
    }
}

// ------------------------------------------------------------
// 5. Implicit absorption coefficient  [1/s]
// ------------------------------------------------------------
volScalarField coeffImpl
(
    IOobject
    (
        "coeffImpl",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless/dimTime, 0.0)
);

coeffImpl.primitiveFieldRef() =
    fvc::surfaceSum
    (
        isInterfaceFace
      * D1
      * mesh.deltaCoeffs()
      * mag(nDotSf)
      * mesh.magSf()
    )().primitiveField()
  / mesh.V().field();

coeffImpl.correctBoundaryConditions();

// Restrict to liquid cells only
coeffImpl *=
    pos(alpha1 - dimensionedScalar("half", dimless, 0.5));

// ------------------------------------------------------------
// 6. Delta diagnostics
// ------------------------------------------------------------
scalar deltaMin =  GREAT;
scalar deltaMax = -GREAT;

forAll(isInterfaceFace, facei)
{
    if (isInterfaceFace[facei] > 0.5)
    {
        scalar delta = 1.0 / mesh.deltaCoeffs()[facei];
        deltaMin = min(deltaMin, delta);
        deltaMax = max(deltaMax, delta);
    }
}

reduce(deltaMin, minOp<scalar>());
reduce(deltaMax, maxOp<scalar>());

Info<< "delta at iface [m]    = " << deltaMin
    << " / " << deltaMax << endl;

// ------------------------------------------------------------
// 7. Save c before solve (matrix-consistent m_abs)
// ------------------------------------------------------------
volScalarField c_pre
(
    IOobject
    (
        "c_pre",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    c
);

// ------------------------------------------------------------
// 8. Concentration equation
//    cSat enters as a volScalarField source term
// ------------------------------------------------------------
fvScalarMatrix cEqnM
(
    fvm::ddt(c)
  + fvm::div(phi, c)
  - fvm::laplacian(D_masked, c)
  + fvm::Sp(coeffImpl, c)
 ==
    coeffImpl * cSat     // volScalarField * volScalarField → OK
);

fvOptions.constrain(cEqnM);
cEqnM.relax();
cEqnM.solve();
fvOptions.correct(c);

// Force c = cSat in gas cells (use local cell value)
forAll(c, celli)
{
    if (alpha1[celli] < 0.01)
    {
        c[celli] = cSat[celli];
    }
}
c.correctBoundaryConditions();

c.max(0.0);
c.min(1.0);

// ------------------------------------------------------------
// 9. m_abs  [kg/m3/s]  — face flux, liquid cell only
//    Use pre-solve c and local cSat per cell
// ------------------------------------------------------------
volScalarField m_abs_field
(
    IOobject
    (
        "m_abs_field",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
);

scalarField& mabsRef = m_abs_field.primitiveFieldRef();

forAll(dc_f, facei)
{
    if (isInterfaceFace[facei] > 0.5)
    {
        scalar alphaO = alpha1[own[facei]];
        scalar alphaN = alpha1[nei[facei]];
        scalar delta  = 1.0 / mesh.deltaCoeffs()[facei];

        // Local cSat per cell
        scalar cSatO = cSat[own[facei]];
        scalar cSatN = cSat[nei[facei]];

        // Recompute dc_f with pre-solve c and local cSat
        scalar cO_pre = alphaO*c_pre[own[facei]]
                      + (scalar(1)-alphaO)*cSatO;
        scalar cN_pre = alphaN*c_pre[nei[facei]]
                      + (scalar(1)-alphaN)*cSatN;

        scalar dcPre = (alphaO >= alphaN)
                     ? max((cO_pre - cN_pre)/delta, scalar(0))
                     : max((cN_pre - cO_pre)/delta, scalar(0));

        // Use local cSat of the liquid cell as denominator
        label liqCell = (alphaO >= alphaN) ? own[facei] : nei[facei];
        scalar cSatLiq = max(cSat[liqCell], scalar(1e-6)); // guard /0

        scalar jFace = (D1.value() / cSatLiq)
                     * dcPre
                     * mag(nDotSf[facei])
                     * mesh.magSf()[facei];

        mabsRef[liqCell] +=
            rho1.value() * jFace / mesh.V()[liqCell];
    }
}

m_abs_field.correctBoundaryConditions();
m_abs = m_abs_field;

// ------------------------------------------------------------
// 10. Write m_abs at write times
// ------------------------------------------------------------
if (runTime.writeTime())
{
    m_abs.write();
}

// ------------------------------------------------------------
// 11. Update sigma volScalarField using Marangoni law
//
//     c_eff blends c toward cSat near the interface to prevent
//     spurious Marangoni spikes in diffusion-blocked cells.
//     c_eff = alpha1*c + (1-alpha1)*cSat
//           → equals c in bulk liquid, cSat in gas/near interface
//
//     sigma = sigma0 + dsdT*(T - Tr) + dsdc*(c_eff - cr)
// ------------------------------------------------------------
{
    volScalarField c_eff
    (
        IOobject
        (
            "c_eff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1*c + (scalar(1) - alpha1)*cSat
    );

    sigma = sigma0 + dsdT*(T - Tr) + dsdc*(c_eff - cr);

    // Bound: sigma must stay positive
    sigma.max(dimensionedScalar("sigMin", sigma.dimensions(), 1e-4));
}

// ------------------------------------------------------------
// 12. Update rhoCp and rhoCpPhi after alpha / rho update
// ------------------------------------------------------------
rhoCp = alpha1*rho1*cp1 + alpha2*rho2*cp2;
rhoCpPhi = fvc::interpolate(rhoCp)*phi;

// ------------------------------------------------------------
// 13. Diagnostics
// ------------------------------------------------------------
surfaceScalarField magNDotSf
(
    IOobject
    (
        "magNDotSf",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mag(nDotSf)
);

const scalar domainVol = gSum(mesh.V().field());
const scalar lx        = 0.076;
const scalar ly        = 0.020;
const scalar dz        = domainVol / (lx * ly);

const scalar dt    = runTime.deltaTValue();
const scalar GmAbs = fvc::domainIntegrate(m_abs).value();
Mabs_acc += GmAbs * dt;

Info<< "min/max(c)            = " << min(c).value()
    << " / " << max(c).value() << endl;
Info<< "min/max(cSat)         = " << gMin(cSat) << " / " << gMax(cSat) << endl;
Info<< "min/max(c-cSat)       = " << min(c - cSat).value()
    << " / " << max(c - cSat).value() << endl;
Info<< "min/max(coeffImpl)    = " << min(coeffImpl).value()
    << " / " << max(coeffImpl).value() << endl;
Info<< "min/max(sigma)        = " << gMin(sigma) << " / " << gMax(sigma) << endl;
Info<< "delta at iface [m]    = " << deltaMin
    << " / " << deltaMax << endl;
Info<< "int(isInterfaceFace)  = "
    << gSum(isInterfaceFace.primitiveField()) << " faces" << endl;
Info<< "int|nDotSf|*iface     = "
    << gSum
       (
           isInterfaceFace.primitiveField()
         * magNDotSf.primitiveField()
         * mesh.magSf().primitiveField()
       ) << " [m2]" << endl;
Info<< "domain volume         = " << domainVol << " [m3]" << endl;
Info<< "dz (mesh depth)       = " << dz << " [m]" << endl;
Info<< "m_abs      [kg/s]     = " << GmAbs << endl;
Info<< "m_abs/dz   [kg/m/s]   = " << GmAbs / dz << endl;
Info<< "Mabs_acc   [kg]       = " << Mabs_acc << endl;
Info<< "Mabs/dz    [kg/m]     = " << Mabs_acc / dz << endl;

if ((runTime.timeIndex() % absWriteInterval) == 0)
{
    (*absFilePtr) << runTime.value() << " "
                  << Mabs_acc << " "
                  << Mabs_acc / dz
                  << nl;
    absFilePtr->flush();
}
