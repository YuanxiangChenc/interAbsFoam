// TEqn.H
// ============================================================
// Temperature equation for the two-phase LiBr/H₂O system.
//
//   ρcₚ * dT/dt + ρcₚ * ∇·(U T) - ∇·(κ ∇T) = Q_abs
//
// FIXES applied vs original:
//
//  1. Checkerboard in gas cells:
//     Root cause: rhoCp2 = rho2*cp2 ≈ 32 J/(m³ K) is ~110,000×
//     smaller than rhoCp1 ≈ 3.55e6.  The off-diagonal from
//     fvm::div(rhoCpPhi, T) at interface faces carries
//     interp(rhoCp) ≈ rhoCp1/2 — these coefficients are
//     ~8000× larger than the gas cell diagonal (rhoCp2*V/dt).
//     The matrix is violently ill-conditioned in gas rows,
//     causing checkerboard EVEN WITH Q_abs = 0.
//
//     Fix A — implicit gas anchor:
//       Add fvm::Sp(gasAnchorCoeff, T) == gasAnchorCoeff*T_gas_ref
//       in gas cells only.  gasAnchorCoeff uses rho1*cp1 (liquid
//       thermal inertia scale), NOT rho2*cp2 — because the
//       off-diagonal from rhoCpPhi at interface faces carries
//       rhoCp1-scale coefficients.  An anchor based on rhoCp2
//       is ~100,000× too weak and leaves the gas rows ill-
//       conditioned (off-diag/diag ≈ 8000), causing checkerboard
//       even with Q_abs = 0.  Using rhoCp1 gives off-diag/diag
//       ≈ 0.15 — well-conditioned.
//
//     Fix B — hard reset after solve:
//       For cells with alpha < alphaGasCutoff, overwrite T with
//       T_gas_ref (a reference gas temperature read from dict or
//       defaulting to Tr).  Mirrors the c = cSat reset in cEqn.H.
//       T_gas_ref should be the saturation temperature of the
//       vapour at system pressure, or simply the initial gas T.
//
//  2. kappaf masking:
//     mixture.kappaf() interpolates kappa to faces using alpha.
//     Near-interface faces between a liquid and a gas cell get
//     kappaf ≈ 0.5*kappa1 + 0.5*kappa2 ≈ kappa1/2 — acceptable.
//     But in pure-gas regions kappaf ≈ kappa2 ≈ 0, and combined
//     with near-zero rhoCp2 the laplacian adds nothing to the
//     matrix diagonal.  Use the mixture kappaf as-is (it already
//     weights by alpha); the anchor term is the real fix.
// ============================================================

T.correctBoundaryConditions();

// Face-interpolated thermal conductivity — HARMONIC MEAN, not arithmetic.
//
// At an interface face heat conducts through both phases in series:
//   R_total = R_liquid + R_gas
// The effective conductivity is the harmonic mean:
//   kappaf = 2*k1*k2 / (k1+k2)  (weighted by liquid fraction on each side)
//
// For our case: k1=0.44, k2=0.018 → harmonic mean = 0.035 W/(m K)
// vs arithmetic mean = 0.228 W/(m K) (6.5× larger — physically wrong).
//
// Harmonic mean:
//   (a) Reflects the correct physics — the low-k gas side is the
//       bottleneck for cross-interface conduction (Patankar §4).
//   (b) Prevents heat from over-conducting into gas cells.
//   (c) Slightly improves matrix conditioning (smaller laplacian
//       off-diagonal at interface faces).
//
// Implementation: use OpenFOAM's harmonic scheme via manual construction,
// since mixture.kappaf() default is arithmetic (linear interp of kappa).
surfaceScalarField alpha1f(fvc::interpolate(alpha1));
surfaceScalarField kappaf
(
    "kappaf",
    scalar(1.0) /
    (
        alpha1f/kappa1 + (scalar(1.0) - alpha1f)/kappa2
    )
);

// Read T bounds and gas reference temperature from PIMPLE dict
scalar TMin(pimple.dict().lookupOrDefault<scalar>("TMin", 250.0));
scalar TMax(pimple.dict().lookupOrDefault<scalar>("TMax", 400.0));
// T_gas_ref: temperature to anchor gas cells to (K).
// Set to initial gas temperature or vapour saturation temperature.
scalar T_gas_ref(pimple.dict().lookupOrDefault<scalar>("T_gas_ref", Tr.value()));
scalar alphaGasCutoff(pimple.dict().lookupOrDefault<scalar>("alphaGasCutoff", 0.01));

// Absorption heat source [W/m³]:  Ha [J/kg] * m_abs [kg/(m³ s)]
volScalarField Q_abs
(
    IOobject
    (
        "Q_abs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE     // write for diagnosis
    ),
    Ha * m_abs
);

// ------------------------------------------------------------
// Gas-cell implicit anchor coefficient  [J/(m³ K s)]
//
//   BINARY gate: nonzero only in PURE GAS cells (alpha < alphaGasCutoff).
//   Previous version used continuous alpha2 weighting, which was a bug:
//   an interface cell with alpha=0.7 has alpha2=0.3, giving anchor
//   coeff = 0.3 * rhoCp1/dt ≈ 1e10. That's enormous — it fights the
//   Q_abs heat source and pulls interface cell T back to Tr, so the
//   liquid interface never actually heats up. That's why T rises in
//   gas (from residual leaks / anchor overshoot) but NOT in the
//   liquid interface cells where the absorption physically happens.
//
//   Binary gate: pos(alphaGasCutoff - alpha1)
//     = 1 if alpha1 < alphaGasCutoff (pure gas — anchor active)
//     = 0 if alpha1 > alphaGasCutoff (interface + liquid — no anchor)
//   → Q_abs freely deposits heat in interface cells
//   → pure gas cells still get strong anchor to prevent checkerboard
// ------------------------------------------------------------
const dimensionedScalar rhoCp1_dim("rhoCp1", dimEnergy/dimVolume/dimTemperature,
    rho1.value()*cp1.value());
const scalar dt_inv = 1.0 / runTime.deltaTValue();
const dimensionedScalar alphaGasCutoff_dim("alphaGasCutoff", dimless, alphaGasCutoff);

volScalarField gasAnchorCoeff
(
    IOobject("gasAnchorCoeff", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    pos(alphaGasCutoff_dim - alpha1)
      * rhoCp1_dim
      * dimensionedScalar("dtinv", dimless/dimTime, dt_inv)
);

const dimensionedScalar T_gas_ref_dim("T_gas_ref", dimTemperature, T_gas_ref);

for (int TEqnCount = 0; TEqnCount < nEnergyLoops; TEqnCount++)
{
    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp, T)
      + fvm::div(rhoCpPhi, T)
      - fvm::laplacian(kappaf, T)
      + fvm::Sp(gasAnchorCoeff, T)        // implicit gas anchor (lhs)
     ==
        Q_abs
      + gasAnchorCoeff * T_gas_ref_dim    // matching explicit rhs
    );

    TEqn.relax();
    TEqn.solve();

    // Hard clamp: prevents runaway in poorly initialised cells
    forAll(T, cellI)
    {
        T[cellI] = max(min(T[cellI], TMax), TMin);
    }

    const scalar TMax_field = gMax(T.primitiveField());
    const scalar TMin_field = gMin(T.primitiveField());
    Info<< "Energy loop " << TEqnCount
        << "  T range: (" << TMin_field << ", " << TMax_field << ")" << endl;
}

// ------------------------------------------------------------
// Hard reset of T in pure-gas cells after solve
//
// Gas cells have no physical T evolution (no LiBr, no absorption
// heat source). Resetting them to T_gas_ref:
//   a) eliminates any residual checkerboard
//   b) provides a consistent initial condition for the next step
//   c) keeps the boundary condition at the gas/interface face
//      physically meaningful for the kappaf computation
//
// Only applied below alphaGasCutoff (default 0.01) — not at the
// interface itself where mixed cells carry real temperature info.
// ------------------------------------------------------------
forAll(T, cellI)
{
    if (alpha1[cellI] < alphaGasCutoff)
    {
        T[cellI] = T_gas_ref;
    }
}

fvOptions.correct(T);
T.correctBoundaryConditions();

if (runTime.writeTime())
{
    Q_abs.write();
}

Info<< "min/max(T) = " << gMin(T.primitiveField())
    << " / " << gMax(T.primitiveField()) << " [K]" << endl;
