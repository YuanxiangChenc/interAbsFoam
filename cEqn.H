// ============================================================
// cEqn.H  —  Castro-aligned LiBr concentration equation
//
// Uses m_abs (already computed by computeMabs.H just before
// this file, with cSat(T_new)) as the explicit source:
//
//   dc/dt + div(u,c) - laplacian(D,c) = -m_abs * cSat / rho1
//
// Derivation: absorbing water at rate m_abs [kg/(m3 s)] into
// solution of density rho1 dilutes the LiBr fraction at rate
//   dc/dt|_source = -m_abs * cSat / rho1
// (Castro absTCv2 identity, since m_abs = rho1*mx/cSat).
//
// Source <= 0 by construction (m_abs >= 0 bounded in
// computeMabs.H).  Consistent with the Q_abs used in TEqn.H,
// modulo cSat(T_prev) vs cSat(T_new) — the difference between
// these two is the inner-loop residual computed in
// interAbsFoam.C after this file returns.
// ============================================================

alpha1.max(0.0);
alpha1.min(1.0);

// ------------------------------------------------------------
// Explicit source
// ------------------------------------------------------------
volScalarField cSource
(
    "cSource",
    -m_abs * cSat / rho1                   // [1/s], <= 0
);

// ------------------------------------------------------------
// c equation
// ------------------------------------------------------------
fvScalarMatrix cEqnM
(
    fvm::ddt(c)
  + fvm::div(phi, c)
  - fvm::laplacian(D_masked, c)
 ==
    cSource
);

fvOptions.constrain(cEqnM);
cEqnM.relax();
cEqnM.solve();
fvOptions.correct(c);

c.max(0.0);
c.min(1.0);

// ------------------------------------------------------------
// Post-solve gas-cell reset — Castro-faithful.
// ------------------------------------------------------------
#include "invadec2.H"

// ------------------------------------------------------------
// Accumulated mass (uses the mAbs_c value the source was built
// from — that is the rate this iteration actually applied).
// ------------------------------------------------------------
const scalar dt = runTime.deltaTValue();
Mabs_acc += mAbs_c * dt;

// ------------------------------------------------------------
// Split c min/max diagnostic (catches overshoots)
// ------------------------------------------------------------
scalar cLiqMin =  GREAT;
scalar cLiqMax = -GREAT;
scalar cGasMin =  GREAT;
scalar cGasMax = -GREAT;
label  nLiq = 0;
label  nGas = 0;

const scalar diagThresh
(
    pimple.dict().lookupOrDefault<scalar>("liqAlphaThresh", 0.5)
);

forAll(c, celli)
{
    if (alpha1[celli] >= diagThresh)
    {
        cLiqMin = min(cLiqMin, c[celli]);
        cLiqMax = max(cLiqMax, c[celli]);
        ++nLiq;
    }
    else
    {
        cGasMin = min(cGasMin, c[celli]);
        cGasMax = max(cGasMax, c[celli]);
        ++nGas;
    }
}
reduce(cLiqMin, minOp<scalar>());
reduce(cLiqMax, maxOp<scalar>());
reduce(cGasMin, minOp<scalar>());
reduce(cGasMax, maxOp<scalar>());
reduce(nLiq,    sumOp<label>());
reduce(nGas,    sumOp<label>());

Info<< "c liquid min/max [-]    = "
    << cLiqMin << " / " << cLiqMax
    << "   (n = " << nLiq << ")" << endl;
Info<< "c gas    min/max [-]    = "
    << cGasMin << " / " << cGasMax
    << "   (n = " << nGas << ")" << endl;
Info<< "Mabs_acc         [kg]   = " << Mabs_acc << endl;

if (runTime.writeTime())
{
    m_abs.write();
    Q_abs.write();
    cSource.write();
}

if ((runTime.timeIndex() % absWriteInterval) == 0)
{
    (*absFilePtr) << runTime.value() << " "
                  << Mabs_acc         << " "
                  << mAbs_T           << " "
                  << mAbs_c           << " "
                  << mAbs_residual_iter
                  << nl;
    absFilePtr->flush();
}
