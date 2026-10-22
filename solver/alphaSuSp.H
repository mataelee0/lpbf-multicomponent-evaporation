tmp<volScalarField::Internal> Sp
(
    volScalarField::Internal::New
    (
        "Sp",
        mesh,
        dimensionedScalar(dimless/dimTime, 0)
    )
);

tmp<volScalarField::Internal> Su
(
    volScalarField::Internal::New
    (
        "Su",
        mesh,
        dimensionedScalar(dimless/dimTime, 0)
    )
);


// Phase change alpha1 source
Pair<tmp<volScalarField::Internal>> phaseChangeS
(
    phaseChange.Salpha(alpha1)
);

if (phaseChangeS[0].valid())
{
    Su = phaseChangeS[0];
    Sp = phaseChangeS[1];
}

volScalarField::Internal& SuRef = Su.ref();

//const Field<scalar>& mFluxField = (total_massflux * mag(fvc::grad(alpha1)) * damper)().internalField();

//forAll(SuRef, celli)
//{
//    SuRef[celli] -= mFluxField[celli];
//}



tmp<volScalarField> divU;

if (Sp.valid())
{
    divU =
        (
            mesh.moving()
            ? fvc::div(phiCN() + mesh.phi())
            : fvc::div(phiCN())
        );
}
