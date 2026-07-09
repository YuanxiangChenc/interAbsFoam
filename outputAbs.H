// ==========================================================
// End-of-timestep post-processing.
//
// The mAbs_T and mAbs_c below are the values from the FINAL
// (converged) inner iteration of this timestep.  Mabs_acc is
// integrated using mAbs_c — the rate the c equation actually
// applied in the last iteration.
//
// absFilePtr file columns:
//   time  Mabs_acc  Mabs_acc/dz  mAbs_T  mAbs_c  residual  pIter
// ==========================================================
    const scalar dt = runTime.deltaTValue();
    Mabs_acc += mAbs_c * dt;

    const scalar domainVol = gSum(mesh.V().field());
    const scalar lx        = 0.076;
    const scalar ly        = 0.020;
    const scalar dz        = domainVol / (lx * ly);

    Info<< "==== end of timestep " << runTime.timeName()
	<< "  (" << pIter << " inner iters) ====" << nl
	<< "  mAbs_T (final)        [kg/s]  = " << mAbs_T   << nl
	<< "  mAbs_c (final)        [kg/s]  = " << mAbs_c   << nl
	<< "  m_abs residual (final)        = "
	<< mAbs_residual_iter << nl
	<< "  Mabs_acc              [kg]    = " << Mabs_acc << nl
	<< "  Mabs_acc / dz         [kg/m]  = "
	<< Mabs_acc / dz << nl
	<< endl;

    if ((runTime.timeIndex() % absWriteInterval) == 0)
    {
	(*absFilePtr)
	    << runTime.value()      << " "
	    << Mabs_acc             << " "
	    << Mabs_acc / dz        << " "
	    << mAbs_T               << " "
	    << mAbs_c               << " "
	    << mAbs_residual_iter   << " "
	    << pIter
	    << nl;
	absFilePtr->flush();
    }
