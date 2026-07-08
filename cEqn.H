// ============================================================
// cEqn.H  —  face-based ("blocked") LiBr concentration equation
// Castro-aligned form; runs INSIDE the PIMPLE inner loop.
//
// Interface geometry is pre-computed once per timestep by
// updateInterfaceGeometry.H before the inner loop starts.
// This file consumes those persistent fields:
//   isInterfaceFace, D_masked, nDotSf, nHatf
// and does not recompute them.
//
// cSat is updated by the second updateCSat_T.H call immediately
// before this file (uses T_new, matching Castro absTCv2).
// ============================================================

alpha1.max(0.0);
alpha1.min(1.0);

const labelUList& own = mesh.owner();
const labelUList& nei = mesh.neighbour();

// ------------------------------------------------------------
// 1. One-sided dc/dn at interface faces  [1/m]
//
//    Castro absTCv2 (pure form, no cif blend):
//      alphaO > alphaN  → owner liquid, neighbour gas
//        dc = (c[own] - cSat[nei]) / delta
//      alphaO < alphaN  → neighbour liquid, owner gas
//        dc = (c[nei] - cSat[own]) / delta
//
//    Clamped >= 0 to reject spurious desorption from numerical
//    noise in near-equilibrium cells.
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

        scalar dc;
        if (alphaO > alphaN)
        {
            dc = (c[own[facei]] - cSat[nei[facei]]) / delta;
        }
        else
        {
            dc = (c[nei[facei]] - cSat[own[facei]]) / delta;
        }
        dc_f[facei] = max(dc, scalar(0));
    }
}

// ------------------------------------------------------------
// 2. Implicit diagonal and per-cell RHS field
//
//    Castro absTCv2 discrete source per liquid cell:
//      S[i] = -D/delta * (c[i] - cSat_gas) * |nDotSf| * Sf / V
//
//    Split into implicit LHS + RHS:
//      + fvm::Sp(coeffImpl, c)                  (diagonal)
//      == coeffImpl_B                           (RHS)
//    where
//      coeffImpl   = D/delta * |nDotSf| * Sf / V           [1/s]
//      coeffImpl_B = D/delta * cSat_gas * |nDotSf| * Sf / V [1/s]
//
//    Net source in liquid cell:
//      coeffImpl_B - coeffImpl*c
//    = D/delta * (cSat_gas - c) * |nDotSf| * Sf / V
//    < 0 when c > cSat_gas  → absorption (drives c down).
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
coeffImpl *=
    pos(alpha1 - dimensionedScalar("half", dimless, 0.5));

// RHS field: face loop accumulating cSat_gas into the liquid cell
volScalarField coeffImpl_B
(
    IOobject
    (
        "coeffImpl_B",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless/dimTime, 0.0)
);

{
    scalarField& bRef = coeffImpl_B.primitiveFieldRef();

    forAll(isInterfaceFace, facei)
    {
        if (isInterfaceFace[facei] > 0.5)
        {
            scalar faceCoeff =
                D1.value()
              * mesh.deltaCoeffs()[facei]
              * mag(nDotSf[facei])
              * mesh.magSf()[facei];

            if (alpha1[own[facei]] > alpha1[nei[facei]])
            {
                // owner = liquid, neighbour = gas
                bRef[own[facei]] +=
                    faceCoeff * cSat[nei[facei]] / mesh.V()[own[facei]];
            }
            else
            {
                // neighbour = liquid, owner = gas
                bRef[nei[facei]] +=
                    faceCoeff * cSat[own[facei]] / mesh.V()[nei[facei]];
            }
        }
    }
}

coeffImpl_B.correctBoundaryConditions();
coeffImpl_B *=
    pos(alpha1 - dimensionedScalar("half", dimless, 0.5));

// ------------------------------------------------------------
// 3. Save c before solve — matrix-consistent m_abs
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
// 4. Concentration equation
//
//    ddt(c) + div(u,c) - laplacian(D_masked,c)
//    + fvm::Sp(coeffImpl, c) == coeffImpl_B
//
//    Interface diffusion is REPLACED by the coeffImpl / coeffImpl_B
//    absorption pair — D_masked is zero on interface faces.
// ------------------------------------------------------------
fvScalarMatrix cEqnM
(
    fvm::ddt(c)
  + fvm::div(phi, c)
  - fvm::laplacian(D_masked, c)
  + fvm::Sp(coeffImpl, c)
 ==
    coeffImpl_B
);

fvOptions.constrain(cEqnM);
cEqnM.relax();
cEqnM.solve();
fvOptions.correct(c);

c.max(0.0);
c.min(1.0);

// No post-solve gas reset — invadec2 at the top of the next
// inner iteration handles it, using the fresh cSat.

// ------------------------------------------------------------
// 5. m_abs  [kg/m3/s]  — pure gradient, c_pre-consistent
//
//    m_abs[liq] = rho1/cSat_liq * D * dcPre * |nDotSf| * Sf / V
//
//    Castro absTCv2: mx = D*(c-cSat)/delta*N;
//    Castro m_abs = mx*rho/cSat (mass balance identity).
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

{
    scalarField& mabsRef = m_abs_field.primitiveFieldRef();

    forAll(isInterfaceFace, facei)
    {
        if (isInterfaceFace[facei] > 0.5)
        {
            scalar alphaO = alpha1[own[facei]];
            scalar alphaN = alpha1[nei[facei]];
            scalar delta  = 1.0 / mesh.deltaCoeffs()[facei];

            label  liqCell;
            scalar cSatGas;
            if (alphaO > alphaN)
            {
                liqCell = own[facei];
                cSatGas = cSat[nei[facei]];
            }
            else
            {
                liqCell = nei[facei];
                cSatGas = cSat[own[facei]];
            }

            scalar dcPre = max(
                (c_pre[liqCell] - cSatGas) / delta, scalar(0));
            scalar cSatLiq = max(cSat[liqCell], scalar(1e-6));

            scalar jFace = (D1.value() / cSatLiq)
                         * dcPre
                         * mag(nDotSf[facei])
                         * mesh.magSf()[facei];

            mabsRef[liqCell] +=
                rho1.value() * jFace / mesh.V()[liqCell];
        }
    }
}

m_abs_field.correctBoundaryConditions();
m_abs = m_abs_field;

if (runTime.writeTime())
{
    m_abs.write();
}

// ------------------------------------------------------------
// 6. Diagnostics
//
//    cif and sigma are updated in UEqn.H (Marangoni only).
//    rhoCp is updated in alphaEqnSubCycle.H after alpha advance.
// ------------------------------------------------------------
const scalar domainVol = gSum(mesh.V().field());
const scalar lx        = 0.076;
const scalar ly        = 0.020;
const scalar dz        = domainVol / (lx * ly);

const scalar dt    = runTime.deltaTValue();
const scalar GmAbs = fvc::domainIntegrate(m_abs).value();
Mabs_acc += GmAbs * dt;

Info<< "min/max(c)            = " << min(c).value()
    << " / " << max(c).value() << endl;
Info<< "min/max(cSat)         = " << gMin(cSat)
    << " / " << gMax(cSat) << endl;
Info<< "min/max(c-cSat)       = " << min(c - cSat).value()
    << " / " << max(c - cSat).value() << endl;
Info<< "min/max(coeffImpl)    = " << min(coeffImpl).value()
    << " / " << max(coeffImpl).value() << endl;
Info<< "min/max(coeffImpl_B)  = " << min(coeffImpl_B).value()
    << " / " << max(coeffImpl_B).value() << endl;
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
