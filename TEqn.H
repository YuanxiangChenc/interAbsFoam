// TEqn.H
// ============================================================
// Temperature equation for the two-phase LiBr/H₂O system.
//
//   ρcₚ * dT/dt + ρcₚ * ∇·(U T) - ∇·(κ ∇T) = Q_abs
//
// Thermal conductivity κ (W/(m K)) comes from mixture.kappaf()
// which returns the face-interpolated κ using the custom mixture
// class.  rhoCp and rhoCpPhi are declared in createFields.H and
// updated in alphaEqnSubCycle.H after every alpha solve.
//
// Q_abs = Ha * m_abs  [W/m³]  — exothermic heat of absorption,
// concentrated in liquid cells at the interface where m_abs > 0.
//
// The inner nEnergyLoops loop follows the uploaded solver pattern:
// solve → clamp → repeat.  TMin/TMax are read from PIMPLE dict.
// ============================================================

// Face-interpolated thermal conductivity from mixture model
surfaceScalarField kappaf(mixture.kappaf());

T.correctBoundaryConditions();

// Read T bounds from fvSolution/PIMPLE dict
scalar TMin(readScalar(pimple.dict().lookup("TMin")));
scalar TMax(readScalar(pimple.dict().lookup("TMax")));

// Absorption heat source [W/m³]:  Ha [J/kg] * m_abs [kg/(m³ s)]
volScalarField Q_abs
(
    IOobject
    (
        "Q_abs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Ha * m_abs
);

scalar TMax_field = 0.0;
scalar TMin_field = 0.0;

for (int TEqnCount = 0; TEqnCount < nEnergyLoops; TEqnCount++)
{
    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp, T)
      + fvm::div(rhoCpPhi, T)
      - fvm::laplacian(kappaf, T)
     ==
        Q_abs
    );

    TEqn.relax();
    TEqn.solve();

    TMax_field = gMax(T.primitiveField());
    TMin_field = gMin(T.primitiveField());

    Info<< "Energy loop " << TEqnCount
        << "  T range: (" << TMin_field << ", " << TMax_field << ")" << endl;

    // Hard clamp: prevents runaway in poorly initialised cells
    forAll(T, cellI)
    {
        T[cellI] = max(min(T[cellI], TMax), TMin);
    }
}

fvOptions.correct(T);
T.correctBoundaryConditions();

Info<< "min/max(T) = " << gMin(T.primitiveField())
    << " / " << gMax(T.primitiveField()) << " [K]" << endl;
