// UEqn.H
// ============================================================
// Momentum equation with variable-σ Marangoni and buoyancy.
//
// sigma is a volScalarField updated at the END of cEqn.H using
// c_eff (which blends c → cSat near the interface to avoid
// spurious Marangoni spikes from diffusion-blocked cells).
// Here we only update the interface normal and compute forces.
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
//    sigmaK_var = sigma(T,c_eff) * K   [N/m²]
//    mixture.sigmaK() = sigma0 * K  →  scale by sigma/sigma0
//
//    Fsvn is the face-normal flux entered into fvc::reconstruct()
//    consistent with how interFoam handles surfaceTensionForce().
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
//    F_Marangoni = (I - n̂n̂) · ∇σ  *  |∇α|  *  nearInterface
//
//    Two-layer localisation:
//      |∇α|            — VOF delta function, ~2 cells wide
//      nearInterface() — pos(α(1-α) - ε), sharper binary mask
//    Product confines the force to the 1-2 cell interface band.
//    Without nearInterface(), fvc::grad(sigma_liq) carries the
//    gradient 3-4 cells into the liquid via its stencil width,
//    and |∇α| alone (2 cells) does not cancel it completely.
//
//    sigma_liq = alpha1*sigma + (1-alpha1)*sigma0
//    Since sigma was built from cif (alpha-weighted) in cEqn.H,
//    gas cells already have sigma ~ sigma0 + dsdT*(T_gas-Tr).
//    After TEqn resets T_gas = T_gas_ref, gas cells have
//    sigma ~ sigma0 + dsdT*(T_gas_ref - Tr) = const, so
//    grad(sigma_liq) = 0 in pure gas automatically.
//    The masking sigma_liq = alpha*sigma + (1-alpha)*sigma0
//    provides additional smoothing at the transition.
// ------------------------------------------------------------
volScalarField sigma_liq
(
    IOobject("sigma_liq", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    alpha1*sigma + (scalar(1) - alpha1)*sigma0
);

volVectorField gradSigma_liq(fvc::grad(sigma_liq));

// Tangential projector: (I - n̂n̂)·∇σ
volVectorField gradSigma_tang
(
    IOobject("gradSigma_tang", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    gradSigma_liq - nHat*(nHat & gradSigma_liq)
);

// nearInterface: sharper binary mask — pos(α(1-α) - ε)
// Multiplied WITH |∇α|, not replacing it:
//   |∇α|          = CSF delta function = interface area per unit volume [1/m]
//                   this is the physical localisation — must be kept
//   nearInterface  = additional sharpening: cuts the 1-2 cell stencil
//                   tail that |∇α| leaves in adjacent liquid cells
//   product        = delta_s * H(α(1-α)-ε) — sharp AND area-weighted
volScalarField nearIface(mixture.nearInterface());

// Full tangential Marangoni force [N/m³]
// F_st = [(I - n̂n̂)·∇σ_liq] * |∇α| * nearInterface
volVectorField FMarangoni
(
    IOobject("FMarangoni", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    gradSigma_tang * mag(gradAlpha) * nearIface
);

// Diagnostic decomposition: Fsvt (thermal) and Fsvc (solutal)
// Same two-mask localisation: |∇α| * nearInterface
// Both use cif for consistency with sigma update in cEqn.H / cl_cEqn.H
volScalarField T_liq
(
    IOobject("T_liq", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    alpha1*T + (scalar(1) - alpha1)*Tr
);

volVectorField gradT_liq(fvc::grad(T_liq));
volVectorField gradCif(fvc::grad(cif));

volVectorField Fsvt
(
    IOobject("Fsvt", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    dsdT * (gradT_liq - nHat*(nHat & gradT_liq))
         * mag(gradAlpha) * nearIface
);

volVectorField Fsvc
(
    IOobject("Fsvc", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    dsdc * (gradCif - nHat*(nHat & gradCif))
         * mag(gradAlpha) * nearIface
);

// ------------------------------------------------------------
// 4. Buoyancy (Boussinesq)  [N/m³]
//
//    Fbuoy = -g * [α₁ρ₁ʳᵉᶠ(βT₁ΔT + βC₁Δc)
//                +α₂ρ₂ʳᵉᶠ(βT₂ΔT + βC₂Δc)]
//    ΔT = T - Tr,  Δc = c - cr   (use raw c for bulk, c_eff unavailable
//    as persistent field — buoyancy acts on bulk c, not interface layer)
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
                Fsvn                      // normal capillary force
              - ghf*fvc::snGrad(rho)      // gravity head
              - fvc::snGrad(p_rgh)        // modified pressure
            ) * mesh.magSf()
        )
      + Fbuoy        // Boussinesq buoyancy
      + FMarangoni   // tangential Marangoni: (I-nn).grad(sigma)*|grad(alpha)|
    );

    fvOptions.correct(U);
}

if (runTime.writeTime())
{
    FMarangoni.write();
    Fsvt.write();   // diagnostic only
    Fsvc.write();   // diagnostic only
    Fbuoy.write();
}
