// ============================================================
// Momentum equation with variable-σ Marangoni and buoyancy.
//
// Tangential Marangoni follows Castro tension.c surftension():
//   F_st = [(I-n̂n̂)·∇σ] * |∇α| * 2*α₁
// Normal capillary from mixture.sigmaK() scaled by sigma/sigma0.
// Buoyancy: Boussinesq with betaT and betaC per phase.
// ============================================================

MRF.correctBoundaryVelocity(U);

// cif and sigma updated here (formerly in cEqn.H).
// cif uses freshest cSat (call 1 of updateCSat_T, T_prev).
cif   = alpha1*c + (scalar(1) - alpha1)*cSat;
cif.max(0.0);
cif.min(1.0);            // guard against any c overshoo
sigma = sigma0 + dsdT*(T - Tr) + dsdc*(cif - cr);
sigma.max(dimensionedScalar("sigMin", sigma.dimensions(), 1e-4));

// ------------------------------------------------------------
// 2. Normal surface-tension / capillary force  [Pa/m on face]
//
//    sigmaK_var = sigma(T,cif) * K  scaled from mixture.sigmaK()
//    Fsvn enters fvc::reconstruct() consistent with interFoam.
// ------------------------------------------------------------
volScalarField sigmaK_var
(
    IOobject("sigmaK_var", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mixture.sigmaK() * (sigma / sigma0)
);

surfaceScalarField Fsvn
(
    IOobject("Fsvn", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    fvc::interpolate(sigmaK_var) * fvc::snGrad(alpha1)
);

// ------------------------------------------------------------
// 3. Tangential Marangoni force  [N/m³]
//
//    F_st = [(I-n̂n̂)·∇σ] * |∇α| * 2 * pos(α₁-0.5)
//
//    |∇α| : CSF delta function [1/m] — interface area per unit
//           volume. SYMMETRIC across the interface: a centred
//           difference picks up the alpha jump equally on the
//           liquid-side cell AND the gas-side cell straddling it.
//
//    pos(α₁-0.5) : binary liquid-side selector (1 if α₁>0.5, else 0).
//           This is NOT the same mistake as the removed alpha1
//           weighting. alpha1 is a CONTINUOUS attenuator that
//           weakens the force at the interface cell itself
//           (alpha~0.5 → factor 0.5). pos(alpha1-0.5) is a BINARY
//           gate that does not touch the magnitude on the liquid
//           side — it only zeroes the mirror contribution that
//           |∇α| deposits on the gas-side cell of the same
//           interface pair. Without this gate, *2 double-counts:
//           both the liquid cell and its gas neighbour receive the
//           full doubled force from the same interface jump.
//
//    factor 2: constant, from Castro's two-sided finite-difference
//    stencil (E+W, N+S averaging) — combined with the liquid-side
//    gate, deposits the full interfacial force once, on the liquid
//    side only. Matches the same pattern used in cl_cEqn.H for
//    absSource: 2*D1*(gradCif & gradAlpha) * pos(alpha1-0.5).
//
//    sigma is already built from cif = alpha1*c + (1-alpha1)*cSat
//    in cEqn.H / cl_cEqn.H. In gas cells cif = cSat and T is anchored
//    to Tr, so sigma naturally collapses to sigma0 there — no extra
//    sigma_liq masking is needed. Use grad(sigma) directly.
// ------------------------------------------------------------
volVectorField gradSigma(fvc::grad(sigma));

// Tangential projector: (I - n̂n̂)·∇σ  [N/m²]
volVectorField gradSigma_tang
(
    IOobject("gradSigma_tang", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    gradSigma - nHat*(nHat & gradSigma)
);

// Binary liquid-side gate: 1 where alpha1>0.5, 0 where alpha1<=0.5
// Same gate used in cl_cEqn.H for absSource — keeps both equations
// consistent in how they localise to the liquid side of the interface.
volScalarField liquidGate
(
    IOobject("liquidGate", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    pos(alpha1 - dimensionedScalar("half", dimless, 0.5))
);

// F_Marangoni = [(I-n̂n̂)·∇σ] * |∇α| * 2 * liquidGate   [N/m³] 
// modifcation : no liquida gate needed here( aling to castro code) 
volVectorField FMarangoni
(
    IOobject("FMarangoni", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    gradSigma_tang * mag(gradAlpha) * scalar(2) * liquidGate
);

// Diagnostic decomposition into thermal (Fsvt) and solutal (Fsvc)
// Same localisation: |∇α| * 2 * liquidGate
// T_liq = alpha1*T + (1-alpha1)*Tr  — one-sided T gradient
// cif already updated in cEqn.H/cl_cEqn.H before UEqn runs
volScalarField T_liq
(
    IOobject("T_liq", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    alpha1*T + (scalar(1) - alpha1)*Tr
);

volVectorField gradT(fvc::grad(T));

// note grad(cif ) move the calcuetaion of c.blending here and since c.bledning ins not longer needed in cEqn.H 
volVectorField gradCif_v(fvc::grad(cif));

// remove liquid gate and replace with FF to alignt with the Castrso code 
volVectorField Fsvt
(
    IOobject("Fsvt", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    dsdT * (gradT - nHat*(nHat & gradT))
         * mag(gradAlpha) * scalar(2) * liquidGate
);

volVectorField Fsvc
(
    IOobject("Fsvc", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    dsdc * (gradCif_v - nHat*(nHat & gradCif_v))
         * mag(gradAlpha) * scalar(2) * liquidGate
);

// ------------------------------------------------------------
// 4. Buoyancy (Boussinesq)  [N/m³]
//
//    Fbuoy = -g * [α₁ρ₁ʳᵉᶠ(βT₁ΔT + βC₁Δc)
//                +α₂ρ₂ʳᵉᶠ(βT₂ΔT + βC₂Δc)]
//    ΔT = T - Tr,  Δc = c - cr
//    betaC1 = -11.09 per unit mass fraction (×100 from Castro wt%)
// ------------------------------------------------------------
volVectorField Fbuoy
(
    IOobject("Fbuoy", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    -g *
    (
        alpha1 * rho1Ref * (betaT1*(T - Tr) + betaC1*(c - cr))
      + alpha2 * rho2Ref * (betaT2*(T - Tr) + betaC2*(c - cr))
    )
);
 
// ------------------------------------------------------------
// 5. Momentum equation
// ------------------------------------------------------------
fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevRhoReff(rho, U)
 ==
    fvOptions(rho, U)
);

UEqn.relax();
fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
                Fsvn                      // normal capillary
              - ghf*fvc::snGrad(rho)      // gravity head
              - fvc::snGrad(p_rgh)        // modified pressure
            ) * mesh.magSf()
        )
      + Fbuoy        // Boussinesq buoyancy
      + FMarangoni   // tangential Marangoni: Castro 2*FF*|N|*(I-nn).grad(sigma)
    );

    fvOptions.correct(U);
}

if (runTime.writeTime())
{
    FMarangoni.write();
    Fsvt.write();
    Fsvc.write();
    Fbuoy.write();
}

