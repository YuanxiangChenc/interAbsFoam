// ============================================================
// Momentum equation with variable-σ Marangoni and buoyancy.
//
// Tangential Marangoni (Castro-aligned, single formulation):
//   F_st = [(I-n̂n̂)·∇σ] * |∇α|
//
// The |∇α| localises the force to the interface layer (natural
// CSF delta function).  No extra ×2 or liquid-side gate — the
// gradient jump in ∇α already handles interface localisation
// symmetrically across the diffuse-interface band.
//
// Diagnostic decompositions Fsvt (thermal) and Fsvc (solutal)
// use the same formulation for consistency.
//
// FMarangoni, Fsvt, Fsvc, Fbuoy are PERSISTENT fields declared
// in createFields.H — this fixes the paraFoam visualisation
// problem (fields must exist in timestep 0/ for paraFoam to
// pick them up as available volume fields).
// ============================================================

MRF.correctBoundaryVelocity(U);

// ------------------------------------------------------------
// 1. Update cif and sigma (Marangoni driver fields)
//    cif uses freshest cSat (call 1 of updateCSat_T, T_prev).
// ------------------------------------------------------------
cif = alpha1*c + (scalar(1) - alpha1)*cSat;
cif.max(0.0);
cif.min(1.0);

sigma = sigma0 + dsdT*(T - Tr) + dsdc*(cif - cr);
sigma.max(dimensionedScalar("sigMin", sigma.dimensions(), 1e-4));

// ------------------------------------------------------------
// 2. Normal capillary force  [Pa/m on face]
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
//    F_st = [(I-n̂n̂)·∇σ] * |∇α|
// ------------------------------------------------------------
volVectorField gradSigma(fvc::grad(sigma));

// Tangential projector: (I - n̂n̂)·∇σ  [N/m²]
volVectorField gradSigma_tang
(
    IOobject("gradSigma_tang", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    gradSigma - nHat*(nHat & gradSigma)
);

// FMarangoni is a persistent field in createFields.H — assign to it.
FMarangoni = gradSigma_tang * mag(gradAlpha);

// ------------------------------------------------------------
// 4. Diagnostic Marangoni decomposition (thermal vs solutal)
//    Same |∇α| localisation as FMarangoni above.  These are
//    diagnostics — they sum to ~FMarangoni if T and cif gradients
//    align, otherwise the difference is a cross-coupling term.
// ------------------------------------------------------------
volVectorField gradT   (fvc::grad(T));
volVectorField gradCif (fvc::grad(cif));

Fsvt = dsdT * (gradT   - nHat*(nHat & gradT))   * mag(gradAlpha);
Fsvc = dsdc * (gradCif - nHat*(nHat & gradCif)) * mag(gradAlpha);

// ------------------------------------------------------------
// 5. Buoyancy (Boussinesq)  [N/m³]
//    Fbuoy = -g * [α₁ρ₁ʳᵉᶠ(βT₁ΔT + βC₁Δc)
//                +α₂ρ₂ʳᵉᶠ(βT₂ΔT + βC₂Δc)]
// ------------------------------------------------------------
Fbuoy =
    -g *
    (
        alpha1 * rho1Ref * (betaT1*(T - Tr) + betaC1*(c - cr))
      + alpha2 * rho2Ref * (betaT2*(T - Tr) + betaC2*(c - cr))
    );

// ------------------------------------------------------------
// 6. Momentum equation
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
      + Fbuoy         // Boussinesq buoyancy
      + FMarangoni    // tangential Marangoni: (I-n̂n̂)·∇σ · |∇α|
    );

    fvOptions.correct(U);
}

// No explicit write() calls — FMarangoni, Fsvt, Fsvc, Fbuoy are
// AUTO_WRITE persistent fields (see createFields.H), so runTime.write()
// handles them automatically at every output timestep.
