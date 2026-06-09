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
// 3. Tangential Marangoni forces  [N/m³]
//
//    F_tang = (∇σ - n̂(n̂·∇σ)) |∇α|
//
//    ∇σ = dσ/dT * ∇T + dσ/dc * ∇(c_eff)
//    The tangential projector (I - n̂n̂) isolates the in-plane
//    component.  |∇α| localises the force at the interface.
//
//    We compute Fsvt and Fsvc separately for diagnostics;
//    the sum enters the momentum predictor.
// ------------------------------------------------------------
volVectorField gradT_vol(fvc::grad(T));
// Use grad(sigma) directly — sigma already encodes both T and c_eff
// Decompose: dσ/dT * ∇T  and  dσ/dc * ∇(c_eff)
// (c_eff is not stored as a persistent field; sigma is.
//  For the tangential split we recompute the two contributions.)
volVectorField gradc_eff(fvc::grad(sigma - sigma0 - dsdT*(T - Tr)) / dsdc);
// gradc_eff = ∇c_eff  recovered from sigma = sigma0 + dsdT*(T-Tr) + dsdc*(c_eff-cr)

volVectorField Fsvt
(
    IOobject("Fsvt", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    dsdT * (gradT_vol - nHat*(nHat & gradT_vol)) * mag(gradAlpha)
);

volVectorField Fsvc
(
    IOobject("Fsvc", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    dsdc * (gradc_eff - nHat*(nHat & gradc_eff)) * mag(gradAlpha)
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
                Fsvn                      // normal Marangoni / capillary
              - ghf*fvc::snGrad(rho)      // gravity head
              - fvc::snGrad(p_rgh)        // modified pressure
            ) * mesh.magSf()
        )
      + Fbuoy      // Boussinesq buoyancy
      + Fsvt       // thermal Marangoni (tangential)
      + Fsvc       // solutal Marangoni (tangential)
    );

    fvOptions.correct(U);
}

if (runTime.writeTime())
{
    Fsvt.write();
    Fsvc.write();
    Fbuoy.write();
}
