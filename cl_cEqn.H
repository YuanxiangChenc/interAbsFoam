// ============================================================
// cl_cEqn.H
//
// Liquid-concentrated cell-centre source with invadec pre-solve.
// Ported into interAbsFoam: cSat is now a volScalarField
// (updated by updateCSat_T.H), so all cSat.value() references
// become per-cell cSat[celli] lookups.
//
// Source : 2 * D1 * (grad(cif) . grad(alpha)) * pos(alpha-0.5)
//          Only liquid cells (alpha > 0.5) receive the source.
//          Factor 2 rationale: with invadec2 in place, gas cells
//          have c = cSat so grad(cif)_gas ~ 0 -- the interface
//          signal sits almost entirely on the liquid side already.
//          Factor 2 is therefore a mild overcount here, unlike
//          ca_cEqn.H where gas cells had stale c and both sides
//          carried real signal. Compare m_abs vs ca_cEqn.H to
//          confirm: if cl gives ~2x ca, remove the 2.
//
// m_abs  : source * rho1 / cSat  [kg/m3/s]
//          Consistent: source [1/s] * rho1 [kg/m3] / cSat [-]
//          = absorbed water mass per unit volume per unit time.
//
// Diffusion: full laplacian(alpha*D1, c), no blocking.
//
// invadecThreshold (CHANGED 1e-6 -> 1e-3):
//   Original cl_cEqn.H used a hasLiquidNeighbour threshold of
//   alpha > 1e-6 to decide if a gas cell borders liquid. This was
//   later found (in iv_cEqn.H debugging) to be too tight: with
//   alpha < 1e-6 reserved for "pure gas", many genuine interface
//   gas cells (0.01 < alpha < 0.5) fell through the gap and were
//   NOT reset to cSat before the solve, contaminating grad(cif)
//   and collapsing m_abs. Widening to 1e-3 catches more of the
//   true interface band while still excluding only the deepest
//   bulk-gas cells from the liquid-neighbour test.
// ============================================================

alpha1.max(0.0);
alpha1.min(1.0);

const labelUList& own = mesh.owner();
const labelUList& nei = mesh.neighbour();

// Threshold for "this cell counts as liquid" when building the
// hasLiquidNeighbour mask in invadec2. Raised from 1e-6 to 1e-3
// per request -- catches more genuine interface-band gas cells.
const scalar invadecThreshold = 1e-3;

// ============================================================
// PRE-SOLVE A: invadec2
//
// Must run BEFORE cif is built so gas-side c is correct.
//   gas cell bordering liquid (alpha of neighbour > threshold)
//     -> c = cSat (local cell value, volScalarField)
//   isolated gas cell (no liquid neighbour within threshold)
//     -> c = 0
// ============================================================
{
    scalarField hasLiquidNbr(mesh.nCells(), 0.0);
    forAll(own, facei)
    {
        if (alpha1[own[facei]] > invadecThreshold)
        {
            hasLiquidNbr[nei[facei]] = 1.0;
        }
        if (alpha1[nei[facei]] > invadecThreshold)
        {
            hasLiquidNbr[own[facei]] = 1.0;
        }
    }

    forAll(c, celli)
    {
        if (alpha1[celli] < invadecThreshold)
        {
            c[celli] = (hasLiquidNbr[celli] > 0.5)
                     ? cSat[celli]   // interface gas cell -> local cSat
                     : 0.0;          // isolated gas cell  -> 0
        }
    }
    c.correctBoundaryConditions();
}

// ============================================================
// PRE-SOLVE B: invadec
//
// Newly flooded cells (alpha now > 0.5, was <= 0.5 last step)
// get c from flux-weighted average of liquid neighbours.
// Uses alphaPhi10 as the VOF flux (volume of liquid per face).
// Only inflow faces contribute (liquid arriving into cell).
// Castro: cant[i,j] = sum(|flux_in| * c_neighbour) / sum(|flux_in|)
// ============================================================
{
    const scalarField& alphaOld = alpha1.oldTime().primitiveField();

    forAll(c, celli)
    {
        // newly flooded: was gas, now liquid
        if (alpha1[celli] > 0.5 && alphaOld[celli] <= 0.5)
        {
            scalar wSum = 0.0;
            scalar cSum = 0.0;

            const labelList& faces = mesh.cells()[celli];
            forAll(faces, fi)
            {
                label facei = faces[fi];
                if (facei >= mesh.nInternalFaces()) continue; // skip boundary

                label o = own[facei];
                label n = nei[facei];
                label nbr = (o == celli) ? n : o;

                // only liquid neighbours contribute
                if (alpha1[nbr] <= 0.5) continue;

                // VOF flux: alphaPhi10 is volume of liquid per face per dt
                // positive = owner->neighbour direction
                scalar flux = alphaPhi10[facei];

                // inflow into celli:
                //   celli is owner     -> inflow when flux < 0
                //   celli is neighbour -> inflow when flux > 0
                scalar influx = (o == celli) ? max(-flux, 0.0)
                                              : max( flux, 0.0);

                if (influx > 0.0)
                {
                    wSum += influx;
                    cSum += influx * c[nbr];
                }
            }

            if (wSum > SMALL)
            {
                c[celli] = cSum / wSum;
            }
            // else: no inflow from liquid neighbours, keep current c
        }
    }
    c.correctBoundaryConditions();
}

// ------------------------------------------------------------
// 1. Interface gradient magnitude
//    gradAlpha_mag = |∇α| — used for diagnostics and integral
//    (gradAlphaVec built in step 2 for the dot product)
// ------------------------------------------------------------
volScalarField gradAlpha_mag(mag(fvc::grad(alpha1)));

// ------------------------------------------------------------
// 2. Blended cell-centre concentration (cif) and Castro-style
//    cell-centre absorption source
//
//    cif = alpha1*c + (1-alpha1)*cSat
//
//    source = 2 * D1 * (grad(cif) . grad(alpha)) * pos(alpha1-0.5)
//
//    liquidGate = pos(alpha1-0.5): BINARY gate, 1 if alpha1>0.5
//    else 0. This is NOT the alpha1 continuous weighting that was
//    correctly removed earlier (that one attenuated the source at
//    the interface cell itself, alpha~0.5 -> factor 0.5, which was
//    wrong). The binary gate does a different job:
//
//    |grad(alpha)| is built from a CENTRED difference, so it is
//    SYMMETRIC across the interface — it picks up the same jump
//    magnitude on the liquid-side cell AND the gas-side cell that
//    straddle the interface. Without a liquid/gas selector, the *2
//    factor deposits the full doubled source on BOTH sides of the
//    interface band, double-counting the same physical jump.
//
//    liquidGate keeps the full magnitude on the liquid side (no
//    attenuation there) and zeroes only the mirror copy that
//    |grad(alpha)| would otherwise deposit in the gas-side cell.
//    This matches the binary liquidGate used in UEqn.H for
//    FMarangoni/Fsvt/Fsvc — both equations now use the same
//    liquid-side localisation convention.
//
//    Why Castro's FF doesn't transfer directly: Castro's grid
//    refines the interface-normal calculation onto a 2n x 2m
//    sub-grid (calcN), and FF in his formula plays a role specific
//    to that staggered/refined bookkeeping. In OpenFOAM's collocated
//    cell-centred grid, fvc::grad(alpha1) already IS the localised
//    interface signal — the binary liquidGate is what restricts
//    deposition to the liquid side, not a continuous alpha1 weight.
//
//    The factor 2 is retained as a constant (not 2*alpha1) — it
//    reflects Castro's two-sided finite-difference stencil
//    (E+W, N+S contributions) and, combined with liquidGate,
//    deposits the full interfacial source once, on the liquid side.
//
//    NOTE: absSource CAN be negative if grad(cif).grad(alpha) < 0,
//    which happens when c < cSat (desorption) or during VOF alpha
//    oscillations near the interface. The source is clamped to
//    zero (absorption-only) before being used for m_abs and the
//    c-equation RHS.
// ------------------------------------------------------------
cif = alpha1*c + (scalar(1) - alpha1)*cSat;

volVectorField gradCif(fvc::grad(cif));
volVectorField gradAlphaVec(fvc::grad(alpha1));

// Binary liquid-side gate: 1 where alpha1>0.5, 0 where alpha1<=0.5
// Same gate used in UEqn.H for FMarangoni/Fsvt/Fsvc
volScalarField liquidGate
(
    IOobject("liquidGate", mesh.time().timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    pos(alpha1 - dimensionedScalar("half", dimless, 0.5))
);

// Raw source — can be negative (desorption / alpha oscillation)
volScalarField absSourceRaw
(
    IOobject("absSourceRaw", mesh.time().timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    scalar(2) * D1 * (gradCif & gradAlphaVec) * liquidGate
);

// Clamped source for m_abs and c-equation RHS (absorption only)
volScalarField absSource
(
    IOobject
    (
        "absSource",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    max
    (
        absSourceRaw,
        dimensionedScalar("zero", absSourceRaw.dimensions(), 0.0)
    )
);

// ------------------------------------------------------------
// 3. Diffusivity field (no blocking -- full mixture diffusion)
// ------------------------------------------------------------
Df = fvc::interpolate(alpha1) * D1;

// ------------------------------------------------------------
// 4. Concentration equation
//
//    ddt(c) + div(phi,c) - laplacian(Df,c) = -absSource
//
//    Gas cells evolve freely via convection after solve.
//    invadec2 resets them again before next time step's solve.
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
// 5. m_abs  [kg/m3/s]
//
//    Castro absTCv: mabs = rho1/cp * D*(grad(cif).N)*FF*dV
//    Isothermal limit: cp -> cSat
//    mabs = (rho1/cSat) * absSource
//
//    absSource is already clamped >= 0 so m_abs >= 0 always.
//    (Negative dot product = desorption / alpha oscillation —
//    not a physical absorption event, excluded from accumulation.)
// ------------------------------------------------------------
{
    scalarField& mabsRef = m_abs.primitiveFieldRef();
    const scalarField& srcRef = absSource.primitiveField();  // clamped
    const scalarField& cSatRef = cSat.primitiveField();

    forAll(mabsRef, celli)
    {
        scalar cSatLocal = max(cSatRef[celli], scalar(1e-6));
        mabsRef[celli] = srcRef[celli] * rho1.value() / cSatLocal;
    }
    m_abs.correctBoundaryConditions();
}

// ------------------------------------------------------------
// 6. Update sigma using cif  (mirrors cEqn.H step 11)
//
//    cif was updated at step 2 above:
//      cif = alpha1*c + (1-alpha1)*cSat
//    sigma = sigma0 + dsdT*(T - Tr) + dsdc*(cif - cr)
//
//    UNIT NOTE (transportProperties):
//      dsdc [N/m per unit mass fraction 0-1]
//           Castro -6.4e-4 is per wt% -> enter -6.4e-2
//      cr   must equal c0 = 0.5522 (NOT cSat0)
// ------------------------------------------------------------
sigma = sigma0 + dsdT*(T - Tr) + dsdc*(cif - cr);
sigma.max(dimensionedScalar("sigMin", sigma.dimensions(), 1e-4));

Info<< "sigma range:  " << gMin(sigma.primitiveField())
    << " – " << gMax(sigma.primitiveField())
    << "  (sigma0 = " << sigma0.value() << " N/m)" << endl;

// ------------------------------------------------------------
// 7. Update rhoCp and rhoCpPhi
// ------------------------------------------------------------
rhoCp    = alpha1*rho1*cp1 + alpha2*rho2*cp2;
rhoCpPhi = fvc::interpolate(rhoCp)*phi;

// ------------------------------------------------------------
// 8. Diagnostics
// ------------------------------------------------------------
const scalar domainVol = gSum(mesh.V().field());
const scalar lx        = 0.076;
const scalar ly        = 0.020;
const scalar dz        = domainVol / (lx * ly);

const scalar dt    = runTime.deltaTValue();
const scalar GmAbs = fvc::domainIntegrate(m_abs).value();
Mabs_acc += GmAbs * dt;

// interface area proxy: integral of |grad(alpha)|
const scalar intNorm = fvc::domainIntegrate(gradAlpha_mag).value();

Info<< "min/max(c)          = " << min(c).value()
    << " / " << max(c).value() << endl;
Info<< "min/max(cSat)       = " << gMin(cSat) << " / " << gMax(cSat) << endl;
Info<< "min/max(cif)        = " << min(cif).value()
    << " / " << max(cif).value() << endl;
Info<< "min/max(absSource)  = " << min(absSource).value()
    << " / " << max(absSource).value() << endl;
Info<< "invadecThreshold    = " << invadecThreshold << endl;
Info<< "int|grad(alpha)| [m2] = " << intNorm << endl;
Info<< "m_abs [kg/s]        = " << GmAbs << endl;
Info<< "m_abs/dz [kg/m/s]   = " << GmAbs/dz << endl;
Info<< "Mabs_acc [kg]       = " << Mabs_acc << endl;
Info<< "Mabs/dz  [kg/m]     = " << Mabs_acc/dz << endl;

if (runTime.writeTime())
{
    m_abs.write();
    cif.write();
    absSource.write();
}

if ((runTime.timeIndex() % absWriteInterval) == 0)
{
    (*absFilePtr) << runTime.value() << " "
                  << Mabs_acc        << " "
                  << Mabs_acc/dz
                  << nl;
    absFilePtr->flush();
}
