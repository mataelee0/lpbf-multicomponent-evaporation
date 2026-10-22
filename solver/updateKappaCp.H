{
	volScalarField kappa1
    (
        IOobject("kappa1", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        kappa 
    );

    volScalarField cp1
    (
        IOobject("cp1", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        cp
    );

	forAll(kappa, celli)
	{
		kappa1[celli] = kappaPtr->value(T[celli]) * epsilon1[celli] + kappasPtr->value(T[celli]) * (1.0 - epsilon1[celli]);
		kappa[celli] = alpha1[celli] * kappa1[celli] + (1.0 - alpha1[celli]) * kappa2l.value();

		cp1[celli] = cpPtr->value(T[celli]) * epsilon1[celli] + cpsPtr->value(T[celli]) * (1.0 - epsilon1[celli]);
		cp[celli] = alpha1[celli] * cp1[celli] + (1.0 - alpha1[celli]) * Cp2l.value();

                electrical_resistivity[celli] = resistivityPtr->value(T[celli]);

	}

    kappa.correctBoundaryConditions();
    cp.correctBoundaryConditions();
}
