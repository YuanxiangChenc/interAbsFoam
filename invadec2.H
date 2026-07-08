
alpha1.max(0.0);
alpha1.min(1.0);

const labelUList& own = mesh.owner();
const labelUList& nei = mesh.neighbour();

// ============================================================
// PRE-SOLVE A: invadec2
//
// Must run BEFORE cif is built so gas-side c is correct.
//   gas cell bordering liquid → c = cSat
//   isolated gas cell         → c = 0
// ============================================================
    scalarField hasLiquidNbr(mesh.nCells(), 0.0);
    forAll(own, facei)
    {
        if (alpha1[own[facei]] > 1e-3) hasLiquidNbr[nei[facei]] = 1.0;
        if (alpha1[nei[facei]] > 1e-3) hasLiquidNbr[own[facei]] = 1.0;
    }
    forAll(c, celli)
    {
        if (alpha1[celli] < 1e-3)
        {
            c[celli] = (hasLiquidNbr[celli] > 0.5)
                     ? cSat[celli]   // CHANGED: cSat.value() -> cSat[celli]
                     : 0.0;
        }
    }
    c.correctBoundaryConditions();

