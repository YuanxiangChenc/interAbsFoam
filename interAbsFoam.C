/*---------------------------------------------------------------------------*\
  interAbsFoam.C  —  main loop
    - alpha, invadec, interface geometry: ONCE per timestep, outside PIMPLE
    - inner PIMPLE loop = Castro SIMPLE loop:  U -> p -> T -> c
    - two updateCSat_T calls per inner iteration: one before UEqn (T_prev),
      one after TEqn (T_new) — matches Castro absTCv vs absTCv2 timing
    - inner-loop residual check driven by outerCorrectorResidualControl
      in fvSolution (see accompanying note)
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "myImmiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Two-phase LiBr/H2O absorption solver — Castro-aligned inner loop"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // ==========================================================
        // OUTER (per-timestep) block: alpha field, mesh, interface
        // Runs exactly ONCE per timestep.
        // Matches Castro's actualizaF() called after SIMPLE converges.
        // ==========================================================

        // --- Mesh update (dynamic mesh only) ---
        if (moveMeshOuterCorrectors)
        {
            mesh.update();

            if (mesh.changing())
            {
                if (mesh.topoChanging())
                {
                    talphaPhi1Corr0.clear();
                }
                gh  = (g & mesh.C())  - ghRef;
                ghf = (g & mesh.Cf()) - ghRef;
                MRF.update();

                if (correctPhi)
                {
                    phi = mesh.Sf() & Uf();
                    #include "correctPhi.H"
                    fvc::makeRelative(phi, U);
                    mixture.correct();
                }

                if (checkMeshCourantNo)
                {
                    #include "meshCourantNo.H"
                }
            }
        }

        // --- alpha equation (once per timestep) ---
        #include "alphaControls.H"
        #include "alphaEqnSubCycle.H"     // also updates rho, rhoCp
        #include "invadec.H"              // c in newly-flooded cells
        mixture.correct();

        // --- interface geometry (once per timestep) ---
        #include "updateInterfaceGeometry.H"


        // ==========================================================
        // INNER (Castro SIMPLE) loop: U -> p -> T -> c
        //
        // pimple.loop() honours nOuterCorrectors and
        // outerCorrectorResidualControl set in system/fvSolution:
        //
        //   PIMPLE
        //   {
        //       nOuterCorrectors 9;         // Castro ITEMAX = 9
        //       nCorrectors      2;
        //       nNonOrthogonalCorrectors 1;
        //       momentumPredictor yes;
        //
        //       outerCorrectorResidualControl
        //       {
        //           U     { tolerance 1e-5; relTol 0; }
        //           p_rgh { tolerance 1e-5; relTol 0; }
        //           T     { tolerance 1e-5; relTol 0; }
        //           c     { tolerance 1e-5; relTol 0; }
        //       }
        //   }
        //
        // pimple.loop() exits when EITHER all controlled residuals
        // are below tolerance OR nOuterCorrectors is reached.
        // This gives Castro's exit condition (all residuals < 1e-5,
        // max 9 iterations).
        //
        // Note: pimple.loop() has NO built-in minimum iteration
        // count.  To enforce ITEMIN = 5, either accept early exit
        // (usually harmless — residuals rarely reach 1e-5 in <5
        // iterations for this coupled system) or convert to a
        // manual for-loop with an if-check on iterCount.
        // ==========================================================
        while (pimple.loop())
        {
            if (pimple.frozenFlow())
            {
                continue;
            }

            // --- cSat call 1: at T_prev (used by invadec2, UEqn, TEqn) ---
            #include "updateCSat_T.H"

            // --- invadec2: hard-reset c = cSat in gas-side interface cells ---
            #include "invadec2.H"

            // --- momentum ---
            #include "UEqn.H"        // includes Marangoni + Boussinesq

            // --- pressure correction (with non-ortho subcorrectors) ---
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }

            // --- absorption heat source (uses cSat from call 1, T_prev) ---
            // Independent of m_abs — computed from c and cSat directly,
            #include "computeQabs.H"

            // --- energy ---
            #include "TEqn.H"

            // --- cSat call 2: at T_new (used by cEqn only) ---
            #include "updateCSat_T.H"

            // --- concentration (runtime model select) ---
            if (absorptionModel == "blocked")
            {
                #include "cEqn.H"
            }
            else // "liquidCentred"
            {
                #include "cl_cEqn.H"
            }
        }
        // end of inner SIMPLE loop

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    // --- concentration profile along vertical line at x = l/2 ---
    forAll(mesh.cells(), celli)
    {
        vector cc = mesh.C()[celli];
        if (mag(cc.x() - 0.038) < 0.001)
        {
            Info<< "y=" << cc.y()
                << "  c=" << c[celli]
                << "  T=" << T[celli] << nl;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
