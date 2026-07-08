// ============================================================
// computeQabs.H
//
// Compute the absorption heat source Q_abs directly from the
// face-gradient of concentration — exactly Castro absTCv.
//
// Runs BEFORE TEqn assembly.  Consumes the CURRENT c field and
// the CURRENT cSat field.  Does NOT depend on m_abs (which is
// only produced later, in cEqn.H).  This decouples TEqn's Q_abs
// from cEqn's m_abs, mirroring Castro's two independent face
// loops (absTCv sees T_prev, absTCv2 sees T_new).
//
// Castro absTCv (lines 583-596):
//   Qx = habs * rho1/cSat * D * (dce*NXe + dcw*NXw)/2 * 2*FF
//   Qy = habs * rho1/cSat * D * (dcn*NYn + dcs*NYs)/2 * 2*FF
//   Q  = (Qx + Qy)                                    [W/m3]
//
// Discrete OpenFOAM analog per liquid cell:
//   Q_abs[i] = Ha * rho1 * D / cSat[i]
//            * SUM over interface faces touching cell i of
//              (c_liq - cSat_gas) / delta * |nDotSf| * Sf / V
//
// Units:
//   Ha            [J/kg]
//   rho1          [kg/m3]
//   D             [m2/s]
//   dc            [1/m]
//   |nDotSf|*Sf/V [1/m]     → face-summed [1/m2]
//   1/cSat        [-]
//   → Q_abs       [J/(m3 s)] = [W/m3]     ✓
//
// The FF factor (Castro's 2*FF) is implicit via nDotSf: for a
// cell deep in liquid, nDotSf = 0 on all faces → Q_abs = 0
// automatically. No explicit FF multiplication needed.
// ============================================================

// Reset Q_abs to zero, then accumulate over interface faces
Q_abs = dimensionedScalar("zero", Q_abs.dimensions(), 0.0);

{
    const labelUList& ownQ = mesh.owner();
    const labelUList& neiQ = mesh.neighbour();

    scalarField& qRef = Q_abs.primitiveFieldRef();
    const scalarField& c_    = c.primitiveField();
    const scalarField& cSat_ = cSat.primitiveField();
    const scalarField& V_    = mesh.V().field();
    const scalarField& deltaCoeffs_ = mesh.deltaCoeffs().primitiveField();
    const scalarField& magSf_       = mesh.magSf().primitiveField();

    forAll(isInterfaceFace, facei)
    {
        if (isInterfaceFace[facei] > 0.5)
        {
            const label o = ownQ[facei];
            const label n = neiQ[facei];
            const scalar delta = 1.0 / deltaCoeffs_[facei];

            label  liqCell;
            scalar cSatGas;
            if (alpha1[o] > alpha1[n])
            {
                liqCell = o;
                cSatGas = cSat_[n];
            }
            else
            {
                liqCell = n;
                cSatGas = cSat_[o];
            }

            const scalar dc = max(
                (c_[liqCell] - cSatGas) / delta, scalar(0));
            const scalar cSatLiq = max(cSat_[liqCell], scalar(1e-6));

            // Absorption heat rate for this face into the liquid cell
            //   = Ha * rho1/cSatLiq * D * dc * |nDotSf| * Sf / V
            const scalar qFace = Ha.value() * rho1.value() / cSatLiq
                               * D1.value() * dc
                               * mag(nDotSf[facei])
                               * magSf_[facei];

            qRef[liqCell] += qFace / V_[liqCell];
        }
    }
}

Q_abs.correctBoundaryConditions();

Info<< "Q_abs (interface heat) min/max [W/m3] = "
    << gMin(Q_abs) << " / " << gMax(Q_abs) << endl;
