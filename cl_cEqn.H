// UEqn.H
// ============================================================
// Momentum equation with variable-σ Marangoni and buoyancy.
//
// Tangential Marangoni follows Castro tension.c surftension():
//   F_st = [(I-n̂n̂)·∇σ] * |∇α| * 2*α₁
// Normal capillary from mixture.sigmaK() scaled by sigma/sigma0.
// Buoyancy: Boussinesq with betaT and betaC per phase.
// ============================================================

MRF.correctBoundaryVelocity(U);

// ------------------------------------------------------------
// 1. Update interface normal from current alpha
// ------------------------------------------------------------
gradAlpha = fvc::grad(alpha1);
nHat      = gradAlpha / (mag(gradAlpha) + deltaN);

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
//    Castro tension.c surftension() exact form:
//      Ftx = 2*(Tx²*sigt*dT/dx + Tx²*sigc*dc/dx
//             +Tx*Ty*sigt*dT/dy + Tx*Ty*sigc*dc/dy) * |N| * FF
//
//    Tensor form:  F_st = [(I-n̂n̂)·∇σ] * |∇α| * 2*α₁
//
//    |∇α| = |N|  : CSF delta function [1/m], interface area/volume
//    2*α₁ = 2*FF : Castro's only localisation — zero in gas,
//                  half-weight at interface, full in bulk liquid
//                  (but |∇α|=0 in bulk so product is zero there)
//    factor 2    : from Castro's two-sided finite-difference stencil
//                  (E+W, N+S averaging), retained to match calibration
//
//    sigma_liq = alpha1*sigma + (1-alpha1)*sigma0
//    Gas cells see sigma0 (flat) → grad(sigma_liq)=0 in pure gas.
//    One-sided gradient at interface is preserved.
// ------------------------------------------------------------
volScalarField sigma_liq
(
    IOobject("sigma_liq", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    alpha1*sigma + (scalar(1) - alpha1)*sigma0
);

volVectorField gradSigma_liq(fvc::grad(sigma_liq));

// Tangential projector: (I - n̂n̂)·∇σ  [N/m²]
volVectorField gradSigma_tang
(
    IOobject("gradSigma_tang", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    gradSigma_liq - nHat*(nHat & gradSigma_liq)
);

// F_Marangoni = [(I-n̂n̂)·∇σ] * |∇α| * 2*α₁   [N/m³]
volVectorField FMarangoni
(
    IOobject("FMarangoni", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    gradSigma_tang * mag(gradAlpha) * scalar(2) * alpha1
);

// Diagnostic decomposition into thermal (Fsvt) and solutal (Fsvc)
// Same Castro localisation: |∇α| * 2*alpha1
// T_liq = alpha1*T + (1-alpha1)*Tr  — one-sided T gradient
// cif already updated in cEqn.H/cl_cEqn.H before UEqn runs
volScalarField T_liq
(
    IOobject("T_liq", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    alpha1*T + (scalar(1) - alpha1)*Tr
);

volVectorField gradT_liq(fvc::grad(T_liq));
volVectorField gradCif_v(fvc::grad(cif));

volVectorField Fsvt
(
    IOobject("Fsvt", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    dsdT * (gradT_liq - nHat*(nHat & gradT_liq))
         * mag(gradAlpha) * scalar(2) * alpha1
);

volVectorField Fsvc
(
    IOobject("Fsvc", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    dsdc * (gradCif_v - nHat*(nHat & gradCif_v))
         * mag(gradAlpha) * scalar(2) * alpha1
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
