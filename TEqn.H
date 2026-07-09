// ============================================================
// TEqn.H  —  Castro-aligned temperature equation
//
//   rho*cp * dT/dt + rho*cp * div(U T) - div(kappa grad T) = Q_abs
//
// Q_abs = Ha * m_abs, where m_abs was built by computeMabs.H
// BEFORE this file, using cSat(T_prev) — matches Castro absTCv.
//
// Q_abs >= 0 always by construction (m_abs is bounded >= 0).
// ============================================================

T.correctBoundaryConditions();

// ------------------------------------------------------------
// Q_abs from m_abs — one-line consistency identity
// ------------------------------------------------------------
Q_abs = Ha * m_abs;

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
// Diagnostics
// ------------------------------------------------------------
const scalar mAbsFromQ_TEqn =
    fvc::domainIntegrate(Q_abs).value() / Ha.value();

Info<< "T  min/max      [K]     = "
    << gMin(T.primitiveField()) << " / "
    << gMax(T.primitiveField()) << endl;
Info<< "Q_abs min/max   [W/m3]  = "
    << gMin(Q_abs.primitiveField()) << " / "
    << gMax(Q_abs.primitiveField()) << endl;
Info<< "m_abs (from Q)  [kg/s]  = " << mAbsFromQ_TEqn << endl;

// Interface temperature diagnostic
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
        << Tif_min << " / " << Tif_mean << " / " << Tif_max
        << "   (n = " << Tif_cnt << " faces)" << endl;
}
