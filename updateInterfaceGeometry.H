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

// --- Interface face mask (owner-neighbour crosses 0.5) ---
forAll(isInterfaceFace, facei)
{
    scalar aO = alpha1[own_[facei]];
    scalar aN = alpha1[nei_[facei]];
    isInterfaceFace[facei] =
        (   (aO > 0.5 && aN < 0.5)
         || (aO < 0.5 && aN > 0.5) )
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
