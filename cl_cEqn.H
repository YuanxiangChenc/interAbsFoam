// ============================================================
// cl_cEqn.H
//
// Liquid-concentrated cell-centre source with invadec pre-solve.
// Follows ca_cEqn.H exactly. Only two changes from ca_cEqn.H:
//   1. cSat: dimensionedScalar -> volScalarField
//            cSat.value() -> cSat[celli] / cSat (field algebra)
//   2. Output: write format/columns matched to interAbsFoam's
//      absFilePtr header (Time, Mabs_acc, Mabs_acc/dz)
//
// Source : 2 * D1 * (grad(cif) . grad(alpha)) * pos(alpha-0.5)
//          Only liquid cells (alpha > 0.5) receive the source.
//          Factor 2 rationale: with invadec2 in place, gas cells
//          have c = cSat so grad(cif)_gas ≈ 0 — the interface
//          signal sits almost entirely on the liquid side already.
//          Factor 2 is therefore a mild overcount here, unlike
//          ca_cEqn.H where gas cells had stale c and both sides
//          carried real signal. Watch the m_abs magnitude vs:1
//
//          ca_cEqn.H to confirm: if cl gives ~2x ca, remove the 2.
//
// m_abs  : source * rho1 / cSat  [kg/m3/s]
//          Consistent: source [1/s] * rho1 [kg/m3] / cSat [-]
//          = absorbed water mass per unit volume per unit time.
//
// Diffusion: full laplacian(alpha*D1, c), no blocking.
// ============================================================

alpha1.max(0.0);
alpha1.min(1.0);

const labelUList& own = mesh.owner();
const labelUList& nei = mesh.neighbour();

// ============================================================
// PRE-SOLVE A: invadec2
//
// Must run BEFORE cif is built so gas-side c is correct.
//   gas cell bordering liquid → c = cSat
//   isolated gas cell         → c = 0
// ============================================================
{
    scalarField hasLiquidNbr(mesh.nCells(), 0.0);
    forAll(own, facei)
    {
        if (alpha1[own[facei]] > 1e-3) hasLiquidNbr[nei[facei]] = 1.0;
        if (alpha1[nei[facei]] > 1e-3) hasLiquidNbr[own[facei]] = 1.0;
    }
    forAll(c, celli)
    {
        if (alpha1[celli] < 1e-3)
        {
            c[celli] = (hasLiquidNbr[celli] > 0.5)
                     ? cSat[celli]   // CHANGED: cSat.value() -> cSat[celli]
                     : 0.0;
        }
    }
    c.correctBoundaryConditions();
}

// ============================================================
// PRE-SOLVE B: invadec
//
// Newly-flooded cells (alpha now > 0.5, was <= 0.5 last step):
// c = VOF-flux-weighted average of liquid neighbours (inflow only).
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
// ------------------------------------------------------------
// 1. Diffusivity: alpha*D1, no blocking
// ------------------------------------------------------------
Df = fvc::interpolate(alpha1) * D1;

// ------------------------------------------------------------
// 2. cif — built AFTER invadec2 so gas side = cSat
//
//    cif = alpha*c + (1-alpha)*cSat
//    liquid bulk (alpha=1):   cif = c
//    gas bordering (alpha~0): cif = alpha*cSat + (1-alpha)*cSat = cSat
//    isolated gas  (alpha=0): cif = 0*(0) + 1*cSat = cSat
//    Wait — isolated gas has c=0 from invadec2:
//      cif = 0*0 + 1*cSat = cSat  (still cSat, correct)
//    Interface cell (alpha=0.7): cif = 0.7*c + 0.3*cSat
// ------------------------------------------------------------
cif = alpha1 * c + (scalar(1.0) - alpha1) * cSat;

// ------------------------------------------------------------
// 3. Gradients at cell centres
// ------------------------------------------------------------
volVectorField gradAlpha
(
    IOobject("gradAlpha_cl", mesh.time().timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    fvc::grad(alpha1)
);

volVectorField gradCif
(
    IOobject("gradCif_cl", mesh.time().timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    fvc::grad(cif)
);

// ------------------------------------------------------------
// 4. Liquid-side absorption source [1/s]
//
//    absSource = 2 * D1 * (grad(cif) . grad(alpha)) * pos(alpha-0.5)
//
//    pos(alpha - 0.5): = 1 if alpha > 0.5, = 0 otherwise
//    Restricts source to liquid cells only.
//
//    Factor 2: compensates for restricting to one side.
//    With invadec2, gas cell grad(cif) ≈ 0, so liquid cell
//    already carries ~full signal. Factor 2 may overcount.
//    → watch m_abs vs ca_cEqn.H:
//      if cl ≈ 2*ca: remove factor 2
//      if cl ≈ ca:   factor 2 is appropriate (split was ~50/50)
//
//    Clamped >= 0: absorption only.
// ------------------------------------------------------------
volScalarField liquidMask
(
    IOobject("liquidMask_cl", mesh.time().timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    pos(alpha1 - dimensionedScalar("half", dimless, 0.5))
);

volScalarField absSource
(
    IOobject("absSource_cl", mesh.time().timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    max
    (
        scalar(2.0) * D1 * (gradCif & gradAlpha) * liquidMask,
        dimensionedScalar("zero", dimless/dimTime, 0.0)
    )
);

// ------------------------------------------------------------
// 5. m_abs [kg/m3/s]
//
//    m_abs = absSource * rho1 / cSat
//
//    Dimensional check:
//      absSource [1/s] * rho1 [kg/m3] / cSat [-] = [kg/m3/s] ✓
//
//    Physical meaning:
//      absSource drives c (LiBr mass fraction) down at rate [1/s].
//      Water absorbed per unit volume = rate of c decrease
//      * rho1 / cSat, where 1/cSat converts from LiBr fraction
//      change to water mass absorbed.
//      (When c drops by cSat, one unit of water mass per rho1
//       has been absorbed — consistent with Castro's prefactor.)
//
//    CHANGED: cSat is now volScalarField, so (rho1/cSat) is
//    plain field algebra — no .value() needed, no forAll loop
//    needed. Division is per-cell automatically.
// ------------------------------------------------------------
m_abs = (rho1 / cSat) * absSource;

// ------------------------------------------------------------
// 6. Concentration equation
//
//    ddt(c) + div(u,c) - laplacian(alpha*D1, c) = -absSource
//
//    absSource only non-zero in liquid cells (liquidMask).
//    Gas cells evolve freely via convection after solve.
//    invadec2 resets them before next time step.
// ------------------------------------------------------------
fvScalarMatrix cEqnM
(
    fvm::ddt(c)
  + fvm::div(phi, c)
  - fvm::laplacian(Df, c)
 ==
  - absSource
);

fvOptions.constrain(cEqnM);
cEqnM.relax();
cEqnM.solve();
fvOptions.correct(c);

c.max(0.0);
c.min(1.0);

// ------------------------------------------------------------
// 7. Diagnostics
// ------------------------------------------------------------
const scalar domainVol = gSum(mesh.V().field());
const scalar lx        = 0.076;
const scalar ly        = 0.020;
const scalar dz        = domainVol / (lx * ly);

const scalar dt    = runTime.deltaTValue();
const scalar GmAbs = fvc::domainIntegrate(m_abs).value();
Mabs_acc += GmAbs * dt;

// interface area proxy: integral of |grad(alpha)| * liquidMask
const scalar intNorm = fvc::domainIntegrate
(
    mag(gradAlpha) * liquidMask
).value();

Info<< "min/max(c)          = " << min(c).value()
    << " / " << max(c).value() << endl;
Info<< "min/max(cif)        = " << min(cif).value()
    << " / " << max(cif).value() << endl;
Info<< "min/max(absSource)  = " << min(absSource).value()
    << " / " << max(absSource).value() << endl;
Info<< "int|N|*liquid [m2]  = " << intNorm << endl;
Info<< "m_abs [kg/s]        = " << GmAbs << endl;
Info<< "m_abs/dz [kg/m/s]   = " << GmAbs/dz << endl;
Info<< "Mabs_acc [kg]       = " << Mabs_acc << endl;
Info<< "Mabs/dz  [kg/m]     = " << Mabs_acc/dz << endl;

if (runTime.writeTime())
{
    m_abs.write();
    cif.write();
}

// CHANGED: output format matched to interAbsFoam's absFilePtr —
// same 3-column write as cEqn.H (Time, Mabs_acc, Mabs_acc/dz)
if ((runTime.timeIndex() % absWriteInterval) == 0)
{
    (*absFilePtr) << runTime.value() << " "
                  << Mabs_acc        << " "
                  << Mabs_acc/dz
                  << nl;
    absFilePtr->flush();
}
