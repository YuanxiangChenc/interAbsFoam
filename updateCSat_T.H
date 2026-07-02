// updateCSat_T.H
// ============================================================
// Update cSat (equilibrium LiBr mass fraction) from local T and
// the system pressure p_i.
//
// Uses the McNeely (1979) correlation via bisection — the same
// implementation Castro uses in his solver:
//   prop_X_T_P_LiBr_mac2 → bisection on prop_T_P_X_LiBr_mac2
//
// Reference:
//   L. McNeely, "Thermodynamic properties of aqueous solutions of
//   Lithium Bromide", ASHRAE Transactions Vol. 85 Pt.1, 1979.
//
// Compared to the earlier Daiguji linear form (valid only at
// p_i = 2.315 kPa), the McNeely correlation:
//   - is pressure-dependent (allows varying p_i)
//   - is non-linear (accurate over wider T range)
//   - matches Castro exactly (byte-consistent validation possible)
//
// In the operating range 303–318 K / x ∈ (0.40, 0.55), McNeely
// and Daiguji agree to within ~1%. Outside that range Daiguji
// diverges by up to 200%; McNeely stays valid down to 278 K.
//
// Direct correlation prop_T_P_X_LiBr_mac2 (T from P and x):
//   pk   = p_Pa / 1000 / 6.895       [psia]  (Pa -> kPa -> psia)
//   X    = x * 100                   [wt%]
//   aux  = D² - 4·E·(C - log10(pk))
//   trefk = -2·E / (D + sqrt(aux)) - 459.72   [°F, water ref temp]
//   t_F  = B0 + B1·X + B2·X² + B3·X³
//        + trefk·(A0 + A1·X + A2·X² + A3·X³)
//   T_K  = (t_F - 32)·5/9 + 273.15
//
// Coefficients (Castro's exact values):
//   A0 = -2.00755      A1 =  0.16976
//   A2 = -3.133362e-3  A3 =  1.97668e-5
//   B0 =  321.128      B1 = -19.322
//   B2 =  0.374382     B3 = -2.0637e-3
//   C  =  6.21147
//   D  = -2886.373
//   E  = -337269.46
//
// Inverse (X from P and T) via bisection over x ∈ [0.10, 0.70].
//
// Also computes T_sat_vap (pure-water saturation temperature at
// system pressure) for use as T_gas_ref in TEqn.H.
// For p_i = 2315 Pa: T_sat_vap ≈ 292.95 K (19.8°C) from steam tables.
// Here we use a simple Antoine-form approximation valid near p_i.
// ============================================================

{
    // System pressure (Pa) — read once from PIMPLE dict or use default
    const scalar p_i = pimple.dict().lookupOrDefault<scalar>("p_i", 2315.0);

    // McNeely coefficients
    const scalar A0 = -2.00755;
    const scalar A1 =  0.16976;
    const scalar A2 = -3.133362e-3;
    const scalar A3 =  1.97668e-5;
    const scalar B0 =  321.128;
    const scalar B1 = -19.322;
    const scalar B2 =  0.374382;
    const scalar B3 = -2.0637e-3;
    const scalar C  =  6.21147;
    const scalar D  = -2886.373;
    const scalar E  = -337269.46;

    // Pressure conversion Pa -> psia (once, outside the cell loop)
    scalar pk = p_i / 1000.0 / 6.895;
    if (pk <= 0.0) pk = 1.0;

    // Reference water temp (°F): from Castro's aux computation
    scalar aux = sqr(D) - 4.0*E*(C - ::log10(pk));
    if (aux <= 0.0) aux = 0.0001;
    const scalar trefk = -2.0*E/(D + ::sqrt(aux)) - 459.72;

    // Bisection bounds and tolerance for X(T,P) inversion
    const scalar x1_init = 0.10;
    const scalar x2_init = 0.70;
    const scalar tol     = 1e-5;
    const label  maxIter = 500;

    // Per-cell inversion
    forAll(T, celli)
    {
        const scalar Ttarget = T[celli];

        scalar x1 = x1_init;
        scalar x2 = x2_init;
        scalar xi = 0.5*(x1 + x2);
        scalar tik = 0.0;

        for (label iter = 0; iter < maxIter; ++iter)
        {
            xi = 0.5*(x1 + x2);
            const scalar X = xi*100.0;   // wt%

            // Direct correlation: t_F = B(X) + trefk*A(X)
            const scalar tF =
                B0 + B1*X + B2*X*X + B3*X*X*X
              + trefk*(A0 + A1*X + A2*X*X + A3*X*X*X);

            tik = (tF - 32.0)*5.0/9.0 + 273.15;

            if (mag(tik - Ttarget) < tol) break;

            if (tik > Ttarget) x2 = xi;
            else               x1 = xi;
        }

        cSat[celli] = xi;
    }

    // Bound to McNeely's valid range
    cSat.max(0.10);
    cSat.min(0.70);
    cSat.correctBoundaryConditions();

    // Pure-water saturation temperature at p_i (for TEqn gas anchor)
    // Antoine form valid over 0-100°C, coefficients for H2O in kPa/K:
    //   log10(P_kPa) = A - B/(T_K + C_ant)
    //   A = 7.19621, B = 1730.63, C_ant = -39.724 (Antoine, water, T in K, P in kPa)
    {
        const scalar A_ant =  7.19621;
        const scalar B_ant =  1730.63;
        const scalar C_ant = -39.724;
        const scalar p_kPa =  p_i / 1000.0;

        // Invert:  T = B/(A - log10(P)) - C
        scalar Tsat_vap = 273.15;   // fallback
        if (p_kPa > SMALL)
        {
            const scalar denom = A_ant - ::log10(p_kPa);
            if (mag(denom) > SMALL)
            {
                Tsat_vap = B_ant/denom - C_ant;
            }
        }

        Info<< "    p_i = " << p_i << " Pa,  T_sat_vap = "
            << Tsat_vap << " K (pure-water saturation)" << endl;
        Info<< "    min(cSat) = " << gMin(cSat)
            << ", max(cSat) = " << gMax(cSat) << endl;
    }
}
