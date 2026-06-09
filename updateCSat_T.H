// updateCSat.H
// Correlation valid at p_i = 2.315e3 Pa:
// Mechanism of absorption enhancement by surfactant. Daiguji, H. 1997. table 1
// T [K] = 1.458*cSat_percent + 242.6
//
// cSat field in solver is mass fraction (0-1), so:
// cSat = ((T - 242.6)/1.458)/100

{
    //Info<< "Updating cSat from T using linear correlation at fixed p_i = 2.315e3 Pa" << endl;

    const dimensionedScalar Tref("Tref", dimTemperature, 242.6);
    const dimensionedScalar slope("slope", dimTemperature, 1.458);

    // First compute cSat in percent from the correlation,
    // then convert percent -> mass fraction
    cSat = (T - Tref)/slope/100.0;

    // Bound to physical range for mass fraction
    cSat.max(0.0);
    cSat.min(1.0);

    cSat.correctBoundaryConditions();

    Info<< "    min(cSat) = " << gMin(cSat)
        << ", max(cSat) = " << gMax(cSat) << endl;
}
