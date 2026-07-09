// ============================================================
// updateInterfaceGeometry.H
//
// Update all interface-normal / interface-mask fields ONCE per
// timestep, right after the alpha equation and invadec.  These
// fields are then frozen for the entire inner PIMPLE loop —
// UEqn.H, cEqn.H, and TEqn.H all read them without recomputing.
//
// Fields updated (all declared as persistent in createFields.H):
//   gradAlpha        [1/m]   volVectorField    — for UEqn Marangoni
//   nHat             [-]     volVectorField    — for UEqn Marangoni
//   gradAlphaf       [1/m]   surfaceVectorField — for cEqn flux
//   nHatf            [-]     surfaceVectorField — for cEqn flux
//   nDotSf           [m2]    surfaceScalarField — for cEqn coeff
//   isInterfaceFace  [-]     surfaceScalarField — for cEqn mask
//   D_masked         [m2/s]  surfaceScalarField — for cEqn laplacian
//
// Not updated here (constant throughout run):
//   Sfhat            [-]     surfaceVectorField — mesh.Sf()/mesh.magSf()
//                                                   (built once in createFields.H)
// ============================================================

const labelUList& own_ = mesh.owner();
const labelUList& nei_ = mesh.neighbour();

// --- Volumetric interface normal (cell-centred) ---
gradAlpha = fvc::grad(alpha1);
nHat      = gradAlpha / (mag(gradAlpha) + deltaN);

// --- Face-native interface normal ---
gradAlphaf = fvc::interpolate(gradAlpha);
nHatf      = gradAlphaf
           / ( mag(gradAlphaf)
             + dimensionedScalar("small", dimless/dimLength, 1e-8)
             );

nDotSf = nHatf & Sfhat;

// --- Interface face mask (owner-neighbour straddle the threshold) ---
// Threshold read from PIMPLE dict.  Default 0.5 = classical
// "sign change of alpha - 0.5" detection.  Can be tuned to
// e.g. 0.1 for a wider interface band that captures diffuse
// interfaces without missing near-flat regions.
const scalar ifaceThresh
(
    pimple.dict().lookupOrDefault<scalar>("ifaceAlphaThresh", 0.5)
);

forAll(isInterfaceFace, facei)
{
    scalar aO = alpha1[own_[facei]];
    scalar aN = alpha1[nei_[facei]];
    isInterfaceFace[facei] =
        (   (aO > ifaceThresh && aN < ifaceThresh)
         || (aO < ifaceThresh && aN > ifaceThresh) )
        ? 1.0 : 0.0;
}
isInterfaceFace.correctBoundaryConditions();

// --- Masked diffusivity: alpha*D1 in bulk, zero on interface faces ---
D_masked =
    (scalar(1.0) - isInterfaceFace)
  * fvc::interpolate(alpha1)
  * D1;

Info<< "Interface geometry updated  |  nInterface faces = "
    << gSum(isInterfaceFace.primitiveField()) << endl;
