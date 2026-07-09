// ============================================================
// invadec2.H  —  Castro-aligned gas-cell BC for c
//
// For every cell with alpha < gasAlphaThresh:
//   - if any neighbour is liquid (alpha >= gasAlphaThresh) →  c = cSat[cell]
//     Absorption Dirichlet BC seen by the c equation.
//   - otherwise (pure gas, no liquid nearby) →  c = 0
//     Pure gas has no LiBr.  Prevents advection noise buildup
//     since D_masked ≈ 0 gives no diffusive smoothing.
//     Castro absorcion.c lines 344-364.
//
// gasAlphaThresh is read from system/fvSolution PIMPLE dict:
//   PIMPLE
//   {
//       gasAlphaThresh 1e-2;   // default 0.5 if omitted
//   }
// A small value (1e-3 to 1e-2) matches Castro's MINF=0 more
// closely — only truly gas cells get overwritten.  0.5 is a
// broader "gas-dominated" definition.  For the stagnant pool
// case a tight threshold is safer since cells with alpha in
// (1e-3, 0.5) are physically part of the diffuse-interface
// liquid and shouldn't be reset to cSat.
//
// Called TWICE per PIMPLE iteration in interAbsFoam.C — once
// after each updateCSat_T.H.
// ============================================================
{
    const labelList& own_ic = mesh.owner();
    const labelList& nei_ic = mesh.neighbour();
    const cellList&  cells  = mesh.cells();

    scalarField& cRef = c.primitiveFieldRef();
    const scalarField& cSatRef = cSat.primitiveField();

    forAll(alpha1, celli)
    {
        if (alpha1[celli] < gasAlphaThresh)
        {
            // Gas cell — check whether any face neighbour is liquid
            bool hasLiquidNbr = false;
            const labelList& cFaces = cells[celli];

            forAll(cFaces, i)
            {
                const label facei = cFaces[i];
                if (mesh.isInternalFace(facei))
                {
                    const label nbr =
                        (own_ic[facei] == celli)
                      ?  nei_ic[facei]
                      :  own_ic[facei];

                    if (alpha1[nbr] >= gasAlphaThresh)
                    {
                        hasLiquidNbr = true;
                        break;
                    }
                }
            }

            if (hasLiquidNbr)
            {
                // Interface gas cell — cSat Dirichlet BC
                cRef[celli] = cSatRef[celli];
            }
            else
            {
                // Pure gas — no LiBr, kill advection noise
                cRef[celli] = 0.0;
            }
        }
    }
}

c.correctBoundaryConditions();
