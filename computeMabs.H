// ============================================================
// computeMabs.H
//
// Compute m_abs [kg/m3/s] directly from the face-gradient of c
// at every interface face.  This is the SINGLE source of truth
// for the absorption rate — used by:
//   - TEqn.H  as Q_abs = Ha * m_abs
//   - cEqn.H  as cSource = -m_abs * cSat / rho1
//   - diagnostics for consistency checks
//
// Physics (Castro absTCv/absTCv2 identity):
//   m_abs[liq] = rho1/cSat_liq * D * dc * |nDotSf| * Sf / V
//   where dc = (c_liq - cSat_gas) / delta
//
// BOUNDED: dc = max(dc, 0) → m_abs >= 0 everywhere.
//   Absorption-only model.  When c drops below cSat_gas
//   (which shouldn't happen at steady state, but can in
//   transients), we freeze at zero rather than desorb.
//   This is physically defensible for an absorber operating
//   above the saturation limit, and prevents oscillation
//   between absorption and desorption within a PIMPLE iter.
//
// Called MULTIPLE times per PIMPLE iteration:
//   1. after updateCSat_T call 1 + invadec2  → m_abs_1 (drives Q_abs in TEqn)
//   2. after updateCSat_T call 2 + invadec2  → m_abs_2 (drives cEqn source)
//   3. after cEqn solve                       → m_abs_3 (post-solve check)
//
// Residual for inner-loop convergence:
//   |integrate(m_abs_3) - integrate(m_abs_2)| / integrate(m_abs_2)
// ============================================================

// Reset field to zero, then accumulate over interface faces
m_abs = dimensionedScalar("zero", m_abs.dimensions(), 0.0);

{
    const labelUList& own_m = mesh.owner();
    const labelUList& nei_m = mesh.neighbour();

    scalarField& mRef = m_abs.primitiveFieldRef();
    const scalarField& c_       = c.primitiveField();
    const scalarField& cSat_    = cSat.primitiveField();
    const scalarField& alpha_   = alpha1.primitiveField();
    const scalarField& V_       = mesh.V().field();
    const scalarField& deltaCf_ = mesh.deltaCoeffs().primitiveField();
    const scalarField& magSf_   = mesh.magSf().primitiveField();

    forAll(isInterfaceFace, facei)
    {
        if (isInterfaceFace[facei] > 0.5)
        {
            const label o = own_m[facei];
            const label n = nei_m[facei];
            const scalar delta = 1.0 / deltaCf_[facei];

            label  liqCell;
            scalar cSatGas;
            if (alpha_[o] > alpha_[n])
            {
                liqCell = o;
                cSatGas = cSat_[n];
            }
            else
            {
                liqCell = n;
                cSatGas = cSat_[o];
            }

            // BOUNDED gradient: dc >= 0 → absorption only, no desorption
            const scalar dc =
                max((c_[liqCell] - cSatGas) / delta, scalar(0));

            const scalar cSatLiq =
                max(cSat_[liqCell], scalar(1e-6));

            // Face contribution to m_abs [kg/(m3 s)]
            //   = rho1/cSat_liq * D * dc * |nDotSf| * Sf  ... divided by V
            const scalar mFace =
                rho1.value() / cSatLiq
              * D1.value() * dc
              * mag(nDotSf[facei])
              * magSf_[facei];

            mRef[liqCell] += mFace / V_[liqCell];
        }
    }
}

m_abs.correctBoundaryConditions();
