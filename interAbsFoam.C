/*---------------------------------------------------------------------------*\
  interAbsFoam.C  —  main loop
  Castro-aligned structure:
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
        // Convergence criterion:
        //   residual = |mAbs_c - mAbs_T| / mAbs_c
        // where
        //   mAbs_T = m_abs sampled BEFORE T solve  (uses cSat(T_prev))
        //   mAbs_c = m_abs sampled BEFORE c solve  (uses cSat(T_new))
        //
        // When the T solve barely changes cSat, m_abs computed with
        // T_prev matches m_abs computed with T_new → coupling has
        // converged.
        //
        // Set in system/fvSolution PIMPLE dict:
        //   nOuterCorrectors 9;             // Castro ITEMAX (default 9)
        //   nMinInnerIter    5;             // Castro ITEMIN (default 5)
        //   mAbsResidualTol  1e-3;          // convergence target
        // ==========================================================
        const label  nMinInnerIter =
            pimple.dict().lookupOrDefault<label>("nMinInnerIter", 5);
        const scalar mAbsResTol    =
            pimple.dict().lookupOrDefault<scalar>("mAbsResidualTol", 1e-3);

        label pIter = 0;

        while (pimple.loop())
        {
            ++pIter;

            if (pimple.frozenFlow())
            {
                continue;
            }

            Info<< "PIMPLE inner iteration " << pIter << endl;

            // ----- cSat call 1: T_prev -----
            #include "updateCSat_T.H"
            #include "invadec2.H"

            // ----- m_abs for TEqn source (cSat(T_prev)) -----
            #include "computeMabs.H"
            mAbs_T = fvc::domainIntegrate(m_abs).value();
            Info<< "mAbs_T (source for TEqn)  [kg/s] = " << mAbs_T << endl;

            // ----- momentum -----
            #include "UEqn.H"

            // ----- pressure correction -----
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }

            // ----- energy: uses Q_abs = Ha * m_abs (m_abs = mAbs_T-consistent) -----
            #include "TEqn.H"

            // ----- cSat call 2: T_new -----
            #include "updateCSat_T.H"
            #include "invadec2.H"

            // ----- m_abs for cEqn source (cSat(T_new)) -----
            #include "computeMabs.H"
            mAbs_c = fvc::domainIntegrate(m_abs).value();
            Info<< "mAbs_c (source for cEqn)  [kg/s] = " << mAbs_c << endl;

            // ----- inner-loop residual: how much T-solve changed m_abs -----
            mAbs_residual_iter =
                mag(mAbs_c - mAbs_T)
              / max(mag(mAbs_c), scalar(1e-20));

            // ----- concentration (uses m_abs = mAbs_c-consistent as source) -----
            
            #include "cEqn.H"
            // ----- Report + early exit on m_abs convergence -----
            Info<< "----------------------------------------------------" << nl
                << "PIMPLE iter " << pIter << " summary:" << nl
                << "  mAbs_T (T-source)   [kg/s] = " << mAbs_T << nl
                << "  mAbs_c (c-source)   [kg/s] = " << mAbs_c << nl
                << "  m_abs residual             = "
                << mAbs_residual_iter
                << "   (tol = " << mAbsResTol
                << ", min iter = " << nMinInnerIter << ")" << nl
                << "----------------------------------------------------" << endl;

            if
            (
                pIter >= nMinInnerIter
             && mAbs_residual_iter < mAbsResTol
            )
            {
                Info<< "==> Inner loop converged on m_abs after "
                    << pIter << " iterations" << nl << endl;
                break;
            }
        }

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
