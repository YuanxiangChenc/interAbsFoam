// ============================================================
// invadec2.H  —  Castro-aligned gas-cell BC for c
//
// For every cell with alpha < 0.5:
//   - if any neighbour is liquid (alpha >= 0.5)  →  c = cSat[cell]
//     This is the absorption Dirichlet BC seen by the c equation
//     when its diffusion stencil crosses the interface.
//   - otherwise (pure gas, no liquid nearby)     →  c = 0
//     Pure gas has no LiBr.  Not resetting these leaves them
//     susceptible to div(phi, c) advection noise, which with
//     D_masked ≈ 0 in gas has no smoothing → checkerboarding.
//     Castro absorcion.c lines 344-364 do exactly this reset.
//
// Called TWICE per PIMPLE iteration in interAbsFoam.C — once
// after each updateCSat_T.H, so the gas BC always reflects the
// freshest cSat.
// ============================================================

{
    const labelList& own_ic = mesh.owner();
    const labelList& nei_ic = mesh.neighbour();
    const cellList&  cells  = mesh.cells();

    scalarField& cRef = c.primitiveFieldRef();
    const scalarField& cSatRef = cSat.primitiveField();

    forAll(alpha1, celli)
    {
        if (alpha1[celli] < 0.01)
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

                    if (alpha1[nbr] >= 0.01)
                    {
                        hasLiquidNbr = true;
                        break;
                    }
                }
            }

            if (hasLiquidNbr)
            {
                // Interface gas cell — impose cSat Dirichlet BC
                cRef[celli] = cSatRef[celli];
            }
            else
            {
                // Pure gas — no LiBr, kill any advection noise
                cRef[celli] = 0.0;
            }
        }
    }
}

c.correctBoundaryConditions();
