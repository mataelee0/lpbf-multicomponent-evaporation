// ============================================================
// 0. Preprocessing
// ============================================================

// set condition
scalar alpha_min = 1e-2;

volScalarField sumAlpha
(
    IOobject
    (
        "sumAlpha",
        mesh.time().timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(dimless, 0.0)
);



// ============================================================
// 1. evaporation
// ============================================================



// ============================================================
// 2. Convection: solve alphas equation
// ============================================================
for (int i = 0; i < n_comp; i++)
{
    volScalarField& alpha = alphas[i];
    fvScalarMatrix alphasEqn
    (
        fvm::ddt(alpha)
        + fvc::div(alphaPhi1, Cs[i])
        == - massfluxs[i] * mag(gradAlpha) * damper
    );

    alphasEqn.solve();
}


// ============================================================
// 3. correction
// ============================================================
sumAlpha *= 0.0;
for (int i = 0; i < n_comp; i++)
{
    sumAlpha += alphas[i];
}


for (int i = 0; i < n_comp; i++)
{
    volScalarField& alpha = alphas[i];

    forAll(alpha, celli)
    {
        if (alpha1[celli] > alpha_min)
        {
            scalar sA = sumAlpha[celli];
            if (sA > alpha_min)
            {
                // mul
                scalar correctionFactor = alpha1[celli] / sA;
                alpha[celli] *= correctionFactor;
            }
            alpha[celli] = max(min(alpha[celli], alpha1[celli]), 0.0);
        }
        else
        {
            alpha[celli] = 0.0;
        }
    }

    alpha.correctBoundaryConditions();
}



// ============================================================
// 4. post process
// ============================================================
sumAlpha *= 0.0;
for (int i = 0; i < n_comp; i++)
{
    sumAlpha += alphas[i];
}

// 4.1 set Cs
for (int i = 0; i < n_comp; i++)
{
    forAll(Cs[i], celli)
    {
        if (sumAlpha[celli] > alpha_min)
        {
            Cs[i][celli] = alphas[i][celli] / sumAlpha[celli];
        }
    }
    Cs[i].correctBoundaryConditions();
    Cs[i] = min(max(Cs[i], 0.0), 1.0);
}


// 4.2 set mole fraction
total_mol_vol *= 0.0;

for (int i = 0; i < n_comp; i++)
{
    total_mol_vol += Cs[i] * rho_list[i] / molarmass_list[i];
}

for (int i = 0; i < n_comp; i++)
{
    scalar rho1 = rho_list[i];
    scalar molarmass1 = molarmass_list[i];

    forAll(mesh.C(), celli)
    {
        if (alpha1[celli] > alpha_min && total_mol_vol[celli] > SMALL)
        {
            Xs[i][celli] = Cs[i][celli] * rho1 / molarmass1 / total_mol_vol[celli];
        }
        else
        {
            Xs[i][celli] = xf_list[i];
        }
    }

    Xs[i].correctBoundaryConditions();
}

// 4.3 set weight percent
total_weight *= 0.0;

for (int i = 0; i < n_comp; i++)
{
    total_weight += Cs[i] * rho_list[i];
}


for (int i = 0; i < n_comp; i++)
{
    scalar rho1 = rho_list[i];
    
    forAll(mesh.C(), celli)
    {
        if (alpha1[celli] > alpha_min && total_weight[celli] > SMALL)
        {
            Ws[i][celli] = Cs[i][celli] * rho1 / total_weight[celli] * 100;
        }
        else
        {
            Ws[i][celli] = wf_list[i] * 100;
        }
    }
}       


// ============================================================
// 5. final prints
// ============================================================
forAll(mesh.C(), celli)
{
    if ((melted[celli] < epsilon1[celli]) && (T[celli] > 800))
    {
        melted[celli] = epsilon1[celli];
    }
}

volScalarField liquid_vol_fraction
(
    IOobject
    (
        "liquid_vol_fraction",
        mesh.time().timeName(),
        mesh
    ),
    melted
);

volScalarField melted_vol_fraction
(
    IOobject
    (
        "melted_vol_fraction",
        mesh.time().timeName(),
        mesh
    ),
    melted
);

liquid_vol_fraction = 0.0;
melted_vol_fraction = 0.0;

forAll(mesh.C(), celli)
{
    if ((alpha_min < alpha1[celli]) && (alpha_min < epsilon1[celli]) && (T[celli] > 800))
    {
        liquid_vol_fraction[celli] = epsilon1[celli] * alpha1[celli];
    }

    if (alpha_min < alpha1[celli] && alpha_min < melted[celli])
    {
        melted_vol_fraction[celli] = melted[celli] * alpha1[celli];
    }
}

liquid = liquid_vol_fraction;
melted_metal = melted_vol_fraction;