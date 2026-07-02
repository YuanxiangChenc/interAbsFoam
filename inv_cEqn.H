// ============================================================
// iv_cEqn.H  —  cEqn.H + Castro invadec pre-solve replacement
//
// Based exactly on cEqn.H with two changes only:
//
//   ADDED:   invadec2 + invadec pre-solve blocks
//            threshold matched to cEqn.H: alpha < 0.5 for gas
//            (not 1e-6 — that was too tight and missed interface
//             gas cells which are the ones that matter for dc_f)
//
//   REMOVED: post-solve hard reset
//            invadec2 pre-solve replaces it with correct timing
//
// Root cause of m_abs collapse in previous version:
//   invadec2 threshold was alpha < 1e-6 → only truly pure gas
//   cells were reset. Interface gas cells (0 < alpha < 0.5)
//   kept stale post-solve c values → dc_f = (c_liq - c_gas_stale)
//   instead of (c_liq - cSat) → wrong gradient → wrong m_abs.
//
//   Fix: reset all cells with alpha < 0.5 (the gas side of the
//   interface) to cSat if they border liquid, 0 if isolated.
//   This matches cEqn.H's alpha < 0.01 reset in spirit but
//   uses the correct Castro two-case logic and correct timing.
// ============================================================

alpha1.max(0.0);
alpha1.min(1.0);

const labelUList& own = mesh.owner();
const labelUList& nei = mesh.neighbour();

// ============================================================
// PRE-SOLVE A: invadec2
//
// Reset ALL gas-side cells (alpha < 0.5) before dc_f is built.
// This includes interface gas cells (0 < alpha < 0.5) which are
// the ones that feed into dc_f and dominate m_abs.
//
// Without this, post-solve c in gas-side interface cells is
// whatever the matrix produced — not cSat — and dc_f is wrong.
//
// Castro threshold: FF <= MINF where MINF = 0.
// In VOF with MULES, alpha is not binary — interface cells have
// 0 < alpha < 1. We use alpha < 0.5 to catch all gas-side cells
// including partial interface cells.
//
//   alpha < 0.5 AND borders liquid  → c = cSat
//   alpha < 0.5 AND no liquid nbr   → c = 0
// ============================================================
{
    scalarField hasLiquidNbr(mesh.nCells(), 0.0);
    forAll(own, facei)
    {
        if (alpha1[own[facei]] > 0.5) hasLiquidNbr[nei[facei]] = 1.0;
        if (alpha1[nei[facei]] > 0.5) hasLiquidNbr[own[facei]] = 1.0;
    }
    forAll(c, celli)
    {
        if (alpha1[celli] < 0.5)
        {
            c[celli] = (hasLiquidNbr[celli] > 0.5)
                     ? cSat[celli]    // gas-side interface cell → cSat
                     : 0.0;           // isolated gas cell       → 0
        }
    }
    c.correctBoundaryConditions();
}

// ============================================================
// PRE-SOLVE B: invadec
//
// Newly-flooded cells (alpha now > 0.5, was <= 0.5 last step).
// Set c = VOF-flux-weighted average of liquid neighbours.
// ============================================================
{
    const scalarField& alphaOld = alpha1.oldTime().primitiveField();
    forAll(c, celli)
    {
        if (alpha1[celli] > 0.5 && alphaOld[celli] <= 0.5)
        {
            scalar wSum = 0.0;
            scalar cSum = 0.0;
            const labelList& cellFaces = mesh.cells()[celli];
            forAll(cellFaces, fi)
            {
                label facei = cellFaces[fi];
                if (facei >= mesh.nInternalFaces()) continue;
                label o   = own[facei];
                label n   = nei[facei];
                label nbr = (o == celli) ? n : o;
                if (alpha1[nbr] <= 0.5) continue;
                scalar influx = (o == celli)
                              ? max(-alphaPhi10[facei], 0.0)
                              : max( alphaPhi10[facei], 0.0);
                if (influx > SMALL)
                {
                    wSum += influx;
                    cSum += influx * c[nbr];
                }
            }
            if (wSum > SMALL)
                c[celli] = cSum / wSum;
        }
    }
    c.correctBoundaryConditions();
}

// ============================================================
// From here: identical to cEqn.H
// ============================================================

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
// 2. Blocked diffusivity
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
    nHatf & (mesh.Sf() / mesh.magSf())
);

// ------------------------------------------------------------
// 4. One-sided dc/dn at interface faces [1/m]
//
//    With invadec2 above: gas-side cells have c = cSat, so
//      cN = alphaN*cSat + (1-alphaN)*cSat = cSat  (gas side)
//      cO = alphaO*c_liq + (1-alphaO)*cSat        (liquid side)
//    dc_f = (cO - cSat)/delta = alphaO*(c_liq-cSat)/delta
//    Clean driving force, no stale gas-side contamination.
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
        scalar delta  = 1.0 / mesh.deltaCoeffs()[facei];

        scalar cO = alphaO*c[own[facei]]
                  + (scalar(1)-alphaO)*cSat[own[facei]];
        scalar cN = alphaN*c[nei[facei]]
                  + (scalar(1)-alphaN)*cSat[nei[facei]];

        if (alphaO >= alphaN)
            dc_f[facei] = max((cO - cN)/delta, scalar(0));
        else
            dc_f[facei] = max((cN - cO)/delta, scalar(0));
    }
}

// ------------------------------------------------------------
// 5. Implicit coefficient [1/s]
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
coeffImpl *= pos(alpha1 - dimensionedScalar("half", dimless, 0.5));

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
Info<< "delta at iface [m]    = " << deltaMin << " / " << deltaMax << endl;

// ------------------------------------------------------------
// 7. Save c before solve
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
// ------------------------------------------------------------
fvScalarMatrix cEqnM
(
    fvm::ddt(c)
  + fvm::div(phi, c)
  - fvm::laplacian(D_masked, c)
  + fvm::Sp(coeffImpl, c)
 ==
    coeffImpl * cSat
);

fvOptions.constrain(cEqnM);
cEqnM.relax();
cEqnM.solve();
fvOptions.correct(c);

// NO post-solve reset — invadec2 handles gas cells next time step
c.max(0.0);
c.min(1.0);

// ------------------------------------------------------------
// 9. m_abs [kg/m3/s] — face flux, liquid cell only, uses c_pre
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

        scalar cSatO = cSat[own[facei]];
        scalar cSatN = cSat[nei[facei]];

        scalar cO_pre = alphaO*c_pre[own[facei]]
                      + (scalar(1)-alphaO)*cSatO;
        scalar cN_pre = alphaN*c_pre[nei[facei]]
                      + (scalar(1)-alphaN)*cSatN;

        scalar dcPre = (alphaO >= alphaN)
                     ? max((cO_pre - cN_pre)/delta, scalar(0))
                     : max((cN_pre - cO_pre)/delta, scalar(0));

        label liqCell = (alphaO >= alphaN) ? own[facei] : nei[facei];
        scalar cSatLiq = max(cSat[liqCell], scalar(1e-6));

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
// 10. Update cif, sigma, rhoCp
// ------------------------------------------------------------
cif = alpha1*c + (scalar(1) - alpha1)*cSat;

sigma = sigma0 + dsdT*(T - Tr) + dsdc*(cif - cr);
sigma.max(dimensionedScalar("sigMin", sigma.dimensions(), 1e-4));

rhoCp = alpha1*rho1*cp1 + alpha2*rho2*cp2;
rhoCpPhi = fvc::interpolate(rhoCp)*phi;

// ------------------------------------------------------------
// 11. Write and diagnostics
// ------------------------------------------------------------
if (runTime.writeTime())
{
    m_abs.write();
    cif.write();
}

surfaceScalarField magNDotSf
(
    IOobject("magNDotSf", mesh.time().timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
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
                  << Mabs_acc        << " "
                  << Mabs_acc / dz
                  << nl;
    absFilePtr->flush();
}
