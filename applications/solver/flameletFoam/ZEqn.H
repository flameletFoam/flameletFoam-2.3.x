	fvScalarMatrix ZEqn
	(
		(
		  fvm::ddt(rho, Z)
		+ fvm::div(phi, Z)
        - fvm::laplacian(turbulence->DZEff(), Z)
        )
    );

	ZEqn.relax();
    ZEqn.solve(mesh.solver("Z"));

	turbulence->correctChi();
	turbulence->correctVarZ();

    combustion->correct();
