// ============================================================
// invadec.H
// Update once every time step after alpha field solved.
// Newly-flooded cells (alpha now > 0.5, was <= 0.5 last step): 
// c = VOF-flux-weighted average of liquid neighbours (inflow only).
// ============================================================
const labelUList& own = mesh.owner();     // move to createField?
const labelUList& nei = mesh.neighbour();

const scalarField& alphaOld = alpha1.oldTime().primitiveField();
forAll(c, celli)
{
	if (alpha1[celli] > 0.5 && alphaOld[celli] <= 0.5)
	{
	    scalar wSum = 0.0;
	    scalar cSum = 0.0;
	    const labelList& cellFaces = mesh.cells()[celli];
	    forAll(cellFaces, fi)
	    {
		label facei = cellFaces[fi];
		if (facei >= mesh.nInternalFaces()) continue;
		label o   = own[facei];
		label n   = nei[facei];
		label nbr = (o == celli) ? n : o;
		if (alpha1[nbr] <= 0.5) continue;
		scalar influx = (o == celli)
			      ? max(-alphaPhi10[facei], 0.0)
			      : max( alphaPhi10[facei], 0.0);
		if (influx > SMALL)
		{
		    wSum += influx;
		    cSum += influx * c[nbr];
		}
	    }
	    if (wSum > SMALL)
		c[celli] = cSum / wSum;
	}
}
c.correctBoundaryConditions();


