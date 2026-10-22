/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  10
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "constant";
    object      physicalProperties.water;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

viscosityModel  temperatureDependentViscosity;

nu table
(
    (901.0 5.98E-07)
    (1252.0 3.30E-07)
    (1856.0 2.24E-07)
    (2790.0 1.85E-07)
);

rho 2500.0;

kappa table
(
    (925.00 87.0)
    (1100.00 93.0)
    (1300.00 100.0)
    (1500.00 107.0)
    (1700.00 114.0)
    (1900.00 121.0)
);

kappas table
(
    (300.00 140.0000)
    (500.00 140.0000)
    (700.00 140.0000)
    (806.00 140.0000)
);

cp table
(
    (901.00 1130.00)
    (2790.00 891.72)
);

cps table
(
    (298.00 953.70)
    (574.00 1073.26)
    (805.00 1196.32)
);

    
Tsolidus 805.0;
Tliquidus 901.0;
LatentHeat 332000;
beta 0.000025;

//elec_resistivity 3.50e-7;

elec_resistivity table
(
    (298.00 1.658e-07)
    (544.00 1.658e-07)
    (839.00 1.658e-07)
    (919.00 3.355e-07)
    (2790.00 5.889e-07)
);

// ************************************************************************* //
