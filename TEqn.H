// ============================================================
// TEqn.H  —  Castro-aligned temperature equation
//
//   rho*cp * dT/dt + rho*cp * div(U T) - div(kappa grad T) = Q_abs
//
// Q_abs is built by computeQabs.H BEFORE this file, using c and
// cSat(T_prev) — matches Castro absTCv timing.
// ============================================================

T.correctBoundaryConditions();

// ------------------------------------------------------------
// Harmonic face conductivity (matches Castro mediarm)
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
// T bounds
// ------------------------------------------------------------
const scalar TMin(pimple.dict().lookupOrDefault<scalar>("TMin", 250.0));
const scalar TMax(pimple.dict().lookupOrDefault<scalar>("TMax", 400.0));

// ------------------------------------------------------------
// Temperature equation
// ------------------------------------------------------------
fvScalarMatrix TEqn
(
    fvm::ddt(rhoCp, T)
  + fvm::div(rhoCpPhi, T)
  - fvm::laplacian(kappaf, T)
 ==
    Q_abs
);

TEqn.relax();
TEqn.solve();

forAll(T, cellI)
{
    T[cellI] = max(min(T[cellI], TMax), TMin);
}

fvOptions.correct(T);
T.correctBoundaryConditions();

// ------------------------------------------------------------
// Diagnostic: m_abs_from_Q  [kg/s]
//
//   m_abs = Q_abs / Ha    (by definition Q_abs = Ha * m_abs)
//
// This is the mass rate implied by the heat source that was
// just used in TEqn.  It should match m_abs_from_c (computed
// in cEqn.H) at PIMPLE convergence — small mismatch in early
// iterations is normal (TEqn used cSat(T_prev), cEqn will use
// cSat(T_new)).  Consistency check: if these two diverge by
// more than ~10% at PIMPLE convergence, look for a bug.
// ------------------------------------------------------------
const scalar mAbsFromQ =
    fvc::domainIntegrate(Q_abs).value() / Ha.value();

Info<< "T  min/max      [K]     = "
    << gMin(T.primitiveField())
    << " / " << gMax(T.primitiveField()) << endl;
Info<< "Q_abs min/max   [W/m3]  = "
    << gMin(Q_abs.primitiveField())
    << " / " << gMax(Q_abs.primitiveField()) << endl;
Info<< "m_abs_from_Q    [kg/s]  = " << mAbsFromQ << endl;

// ------------------------------------------------------------
// Interface temperature diagnostic — physical validation.
//
// Absorption releases heat at the interface (Q_abs > 0 there),
// so T at the interface should rise above the initial value.
// If T_interface stays flat while Q_abs is nonzero, either the
// heat source is not being applied to interface cells, or the
// diffusion is spreading the heat away too fast (bulk T rises
// but interface doesn't stand out).
//
// Face T = 0.5 * (T[owner] + T[neighbour]) — since T is
// continuous across the interface in this formulation, this is
// well-defined and equals the interface T to leading order.
// ------------------------------------------------------------
{
    scalar Tif_min =  GREAT;
    scalar Tif_max = -GREAT;
    scalar Tif_sum = 0.0;
    label  Tif_cnt = 0;

    const labelUList& own_t = mesh.owner();
    const labelUList& nei_t = mesh.neighbour();

    forAll(isInterfaceFace, facei)
    {
        if (isInterfaceFace[facei] > 0.5)
        {
            const scalar Tf = 0.5*(T[own_t[facei]] + T[nei_t[facei]]);
            Tif_min = min(Tif_min, Tf);
            Tif_max = max(Tif_max, Tf);
            Tif_sum += Tf;
            ++Tif_cnt;
        }
    }

    reduce(Tif_min, minOp<scalar>());
    reduce(Tif_max, maxOp<scalar>());
    reduce(Tif_sum, sumOp<scalar>());
    reduce(Tif_cnt, sumOp<label>());

    const scalar Tif_mean =
        (Tif_cnt > 0) ? (Tif_sum / Tif_cnt) : 0.0;

    Info<< "T_interface min/mean/max [K] = "
        << Tif_min << " / " << Tif_mean << " / "
        << Tif_max
        << "   (n = " << Tif_cnt << " faces)" << endl;
}
