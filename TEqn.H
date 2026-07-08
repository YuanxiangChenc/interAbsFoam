// ============================================================
// TEqn.H  —  Castro-aligned temperature equation
//
//   rho*cp * dT/dt + rho*cp * div(U T) - div(kappa grad T) = Q_abs
//
// Q_abs is now computed BEFORE TEqn assembly by computeQabs.H,
// using the current c and cSat (fresh from updateCSat_T call 1,
// T_prev — matches Castro absTCv timing).  Q_abs is a persistent
// volScalarField (declared in createFields.H) so both TEqn.H and
// the diagnostics can read it.
//
// This decouples Q_abs from m_abs entirely:
//   - Q_abs  is the heat source seen by TEqn      (built from c, cSat_prev)
//   - m_abs  is the water mass source in cEqn.H  (built from c_pre, cSat_new)
// They are physically consistent (Q_abs = Ha * m_abs at convergence),
// but computed independently at their respective solve points,
// exactly as Castro does with absTCv and absTCv2.
// ============================================================

T.correctBoundaryConditions();

// ------------------------------------------------------------
// Face-interpolated thermal conductivity — HARMONIC MEAN.
// At an interface face heat conducts through both phases in series:
//   R_total = R_liquid + R_gas → effective k = 2*k1*k2/(k1+k2)
// For our case k1=0.44, k2=0.018 → harmonic ≈ 0.035 W/(m K)
// vs arithmetic 0.228 — the low-k gas is the bottleneck.
// ------------------------------------------------------------
surfaceScalarField alpha1f(fvc::interpolate(alpha1));
surfaceScalarField kappaf
(
    "kappaf",
    scalar(1.0) /
    (
        alpha1f/kappa1 + (scalar(1.0) - alpha1f)/kappa2
    )
);

// ------------------------------------------------------------
// T bounds (read from PIMPLE dict)
// ------------------------------------------------------------
const scalar TMin(pimple.dict().lookupOrDefault<scalar>("TMin", 250.0));
const scalar TMax(pimple.dict().lookupOrDefault<scalar>("TMax", 400.0));

// ------------------------------------------------------------
// Temperature equation
//
//   ddt(rhoCp T) + div(rhoCpPhi, T) - laplacian(kappaf, T) == Q_abs
//
// Note: no gas-cell anchor, no post-solve T reset.  If the gas
// cells show checkerboard, revisit those separately — first
// verify Q_abs magnitudes are sensible from computeQabs.H output.
// ------------------------------------------------------------
fvScalarMatrix TEqn
(
    fvm::ddt(rhoCp, T)
  + fvm::div(rhoCpPhi, T)
  - fvm::laplacian(kappaf, T)
 ==
    Q_abs                              // built by computeQabs.H
);

TEqn.relax();
TEqn.solve();

// Hard clamp against runaway from bad initial conditions
forAll(T, cellI)
{
    T[cellI] = max(min(T[cellI], TMax), TMin);
}

fvOptions.correct(T);
T.correctBoundaryConditions();

if (runTime.writeTime())
{
    Q_abs.write();
}

Info<< "min/max(T) = " << gMin(T.primitiveField())
    << " / " << gMax(T.primitiveField()) << " [K]" << endl;
