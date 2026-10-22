{
    // 1. Calculate thermodynamic activity coefficients

    // Find the minimum vaporization temperature
    scalar Tvaps = 10000.0;
    for (int i = 0; i < n_comp; i++)
    {
        if (Tvap_list[i] < Tvaps)
            Tvaps = Tvap_list[i];
    }
    
	for (int i = 0; i < n_comp; i++)
    {
        gammas[i] *= 0.0;
    }

    double R = constant::physicoChemical::R.value();

    forAll(mesh.C(), celli)
    {
        if (T[celli] >= Tvaps && alpha1[celli] >= 0.001)
        {       
            // Initialize variables for excess Gibbs energy calculations
            double G = 0.0;
            double DGDX[4] = { 0.0, 0.0, 0.0, 0.0 };
            double GMUR[4] = { 0.0, 0.0, 0.0, 0.0 };
            double X[4];
            for (int i = 0; i < n_comp; i++)
            {
                X[i] = Xs[i][celli];
            }

            // Calculate binary contributions to excess Gibbs energy and its derivatives
            for (int i1 = 0; i1 < n_comp - 1; i1++)
                for (int i2 = i1 + 1; i2 < n_comp; i2++)
                {
                    double XX = X[i1] - X[i2];
                    double XXX = X[i1] * X[i2];
                    double FLA = 0.0;
                    double FLB = 0.0;
                    if (i1 == 0 && i2 == 1)
                    {
                        FLA = FL1[0][i1][i2] + FL1[1][i1][i2] * T[celli] + (FL1[2][i1][i2] + FL1[3][i1][i2] * T[celli]) * XX + (FL1[4][i1][i2] + FL1[5][i1][i2] * T[celli]) * XX * XX + (-7290 + 5.5 * T[celli]) * XX * XX * XX;
                        FLB = FL1[2][i1][i2] + FL1[3][i1][i2] * T[celli] + (FL1[4][i1][i2] + FL1[5][i1][i2] * T[celli]) * 2 * XX + (-7290 + 5.5 * T[celli]) * 3 * XX * XX;
                    }
                    else if (i1 == 2 && i2 == 3)
                    {
                        FLA = FL1[0][i1][i2] + FL1[1][i1][i2] * T[celli] -95 * T[celli] * Foam::log(T[celli]) + 0.04 * T[celli] * T[celli] + (FL1[2][i1][i2] + FL1[3][i1][i2] * T[celli]) * XX + (FL1[4][i1][i2] + FL1[5][i1][i2] * T[celli]) * XX * XX;
                        FLB = FL1[2][i1][i2] + FL1[3][i1][i2] * T[celli] + (FL1[4][i1][i2] + FL1[5][i1][i2] * T[celli]) * 2 * XX;
                    }
                    else
                    {
                        FLA = FL1[0][i1][i2] + FL1[1][i1][i2] * T[celli] + (FL1[2][i1][i2] + FL1[3][i1][i2] * T[celli]) * XX + (FL1[4][i1][i2] + FL1[5][i1][i2] * T[celli]) * XX * XX;
                        FLB = FL1[2][i1][i2] + FL1[3][i1][i2] * T[celli] + (FL1[4][i1][i2] + FL1[5][i1][i2] * T[celli]) * 2 * XX;
                    }
                 

                    G += XXX * FLA;
                    DGDX[i1] += X[i2] * FLA + XXX * FLB;
                    DGDX[i2] += X[i1] * FLA - XXX * FLB;
                }

            // Calculate ternary contributions to excess Gibbs energy and its derivatives
            if (n_comp >= 3)
            {
                for (int i1 = 0; i1 < n_comp - 2; i1++)
                    for (int i2 = i1 + 1; i2 < n_comp - 1; i2++)
                        for (int i3 = i2 + 1; i3 < n_comp; i3++)
                        {
                            double X1 = X[i1];
                            double X2 = X[i2];
                            double X3 = X[i3];
                            double TLL1 = TL1[2][i1][i2][i3] + TL1[5][i1][i2][i3] * T[celli];
                            double TLL2 = TL1[3][i1][i2][i3] + TL1[6][i1][i2][i3] * T[celli];
                            double TLL3 = TL1[4][i1][i2][i3] + TL1[7][i1][i2][i3] * T[celli];
                            double XXX = X1 * X2 * X3;
                            double TTT = TLL1 + TLL2 + TLL3;
                            double TTL1 = TLL1 - TTT / 3.0;
                            double TTL2 = TLL2 - TTT / 3.0;
                            double TTL3 = TLL3 - TTT / 3.0;
                            double TLA = TL1[0][i1][i2][i3] + TL1[1][i1][i2][i3] * T[celli] + (X1 + (1 - X1 - X2 - X3) / 3.0) * TLL1 + (X2 + (1 - X1 - X2 - X3) / 3.0) * TLL2 + (X3 + (1 - X1 - X2 - X3) / 3.0) * TLL3;

                            G += XXX * TLA;
                            DGDX[i1] += X2 * X3 * TLA + XXX * TTL1;
                            DGDX[i2] += X1 * X3 * TLA + XXX * TTL2;
                            DGDX[i3] += X1 * X2 * TLA + XXX * TTL3;
                        }
            }

            // Calculate partial molar excess Gibbs energies
            for (int k = 0; k < n_comp; k++)
            {
                GMUR[k] += G;
                for (int j = 1; j < n_comp; j++)
                {
                    if (k == j)
                    {
                        GMUR[k] += (1 - X[j]) * (DGDX[j] - DGDX[0]);
                    }
                    else
                    {
                        GMUR[k] += -X[j] * (DGDX[j] - DGDX[0]);
                    }
                }
            }


            // Calculate activity coefficients
            for (int i = 0; i < n_comp; i++)
            {
                gammas[i][celli] = Foam::exp(GMUR[i] / R / T[celli]);
            }
        }
    }


    // 2. Calculate vapor pressure and Knudsen-layer parameters
    Ps_total *= 0.0;
    Alpha *= 0.0;
    Phi *= 0.0;


    for (int i = 0; i < n_comp; i++)
    {
        Pss[i] *= 0.0;
    }

    volScalarField MPs_temp
    (
        IOobject("MPs_temp", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("MPs_temp", dimensionSet(2, -1, -2, 0, -1), 0.0)
    );


    // Calculate partial vapor pressures and mean vapor molar mass
    for (int i = 0; i < n_comp; i++)
    {
        scalar Latentheatgas1 = Latentheatgas_list[i];
        scalar Tvap1 = Tvap_list[i];
        scalar molarmass1 = molarmass_list[i];
        scalar R = constant::physicoChemical::R.value();

        forAll(mesh.C(), celli)
        {
            if (T[celli] >= Tvaps && alpha1[celli] >= 0.001)
            {
                Pss[i][celli] = p0.value() * Foam::exp(Latentheatgas1 * molarmass1 * ((T[celli] - Tvap1) / (R * T[celli] * Tvap1)));
                scalar Psi = Xs[i][celli] * gammas[i][celli] * Pss[i][celli];
                Ps_total[celli] += Psi;
                MPs_temp[celli] += Psi * molarmass1;
            }            
        }
    }
    
    forAll(mesh.C(), celli)
    {
        if (Ps_total[celli] == 0.0)
            continue;
        molarmass_total[celli] = MPs_temp[celli] / Ps_total[celli];
    }

    // Calculate the Knudsen-layer Mach number
    forAll(mesh.C(), celli)
    {
        if (T[celli] >= Tvaps && alpha1[celli] >= 0.001)
        {
            // Solve the nonlinear equation using the Newton-Raphson method
            scalar Pi = constant::mathematical::pi;
            auto get_m = [&](scalar M3) {
                return Foam::pow(5.0 / 3.0 / 2.0, 0.5) * M3;
                };

            auto T3_Ts = [&](scalar m) {
                scalar term = (5.0 / 3.0 - 1.0) / (5.0 / 3.0 + 1.0) * m / 2.0;
                scalar part1 = Foam::pow(1 + Pi * Foam::pow(term, 2.0), 0.5);
                scalar part2 = -Foam::pow(Pi, 0.5) * term;
                return Foam::pow(part1 + part2, 2.0);
                };

            auto rho3_rhos = [&](scalar m) {
                scalar erfc_m = Foam::erfc(m);
                scalar exp_m2 = Foam::exp(m * m);

                scalar part1 = 1.0 / 2.0 * (2.0 * m * m + 1.0) * exp_m2 * erfc_m - m / Foam::pow(Pi, 0.5);
                scalar part2 = 1.0 - Foam::pow(Pi, 0.5) * m * exp_m2 * erfc_m;
                return Foam::pow(1.0 / T3_Ts(m), 0.5) * part1 + 1.0 / 2.0 / T3_Ts(m) * part2;
                };

            auto a3_a1 = [&](scalar m) {
                return Foam::pow(air_molarmass * T3_Ts(m) * T[celli] / molarmass_total[celli] / 300.0, 0.5);
                };

            auto f = [&](scalar M3) {
                scalar m = get_m(M3);
                scalar a31 = a3_a1(m);
                scalar term1 = (5.0 / 3.0 + 1.0) / 4.0 * M3 * a31;
                scalar part1 = 1.0 + 5.0 / 3.0 * M3 * a31 * (term1 + Foam::pow(1.0 + term1 * term1, 0.5));
                scalar part2 = 1.0 / rho3_rhos(m) / T3_Ts(m);
                scalar part3 = -Ps_total[celli] / p0.value();
                return part1 * part2 + part3;
                };

            auto f_prime = [&](scalar M3) {
                const scalar h = 1.0e-6;
                return (f(M3 + h) - f(M3 - h)) / (2.0 * h);
                };

            auto newton_raphson = [&](scalar initial_guess) {
                scalar M3 = initial_guess;
                const scalar tol = 1.0e-6;
                const int max_iter = 100;
                for (int i = 0; i < max_iter; ++i) {
                    scalar fv = f(M3);
                    scalar fp = f_prime(M3);
                    if (Foam::mag(fp) < 1e-12) {
                        throw std::runtime_error("Derivative too small");
                    }
                    scalar M3_new = M3 - fv / fp;
                    if (Foam::mag(M3_new - M3) < tol)
                        return M3_new;
                    M3 = M3_new;
                }
                Info << "Warning: Knudsen-layer Newton-Raphson solver did not converge" << endl;
                return 0.0;
                };

            scalar M3_value = newton_raphson(1.0);
            if (M3_value < 0.0)
                M3_total[celli] = 0.0;
            else if (M3_value > 1.0)
                M3_total[celli] = 1.0;
            else
                M3_total[celli] = M3_value;

            // Calculate evaporation coefficient and recoil pressure factor
            Phi[celli] = Foam::pow(2.0 * Pi * 5.0 / 3.0, 0.5) * M3_total[celli] * rho3_rhos(get_m(M3_total[celli])) * Foam::pow(T3_Ts(get_m(M3_total[celli])), 0.5);
            Alpha[celli] = 1.0 / 2.0 * (1.0 + 1.0 / 2.0 * (1.0 - Phi[celli]) * (1.0 + rho3_rhos(get_m(M3_total[celli])) * T3_Ts(get_m(M3_total[celli]))));
        }
    }
    

    // 3. Calculate recoil pressure and evaporative cooling
    pVap *= 0.0;
    Qv *= 0.0;

    for (int i = 0; i < n_comp; i++)
    {
        scalar Latentheatgas1 = Latentheatgas_list[i];
        scalar Tvap1 = Tvap_list[i];
        scalar molarmass1 = molarmass_list[i];
        scalar R = constant::physicoChemical::R.value();
        scalar Pi = constant::mathematical::pi;

        volScalarField recoil_alpha(pVap);
        recoil_alpha *= 0.0;
        volScalarField Qv_alpha(Qv);
        Qv_alpha *= 0.0;

        forAll(mesh.C(), celli)
        {
            if (T[celli] >= Tvaps && alpha1[celli] >= 0.001)
            {
                // Calculate element-wise vapor pressure and evaporative cooling
                scalar evap_coeff = Phi[celli];

                recoil_alpha[celli] += Xs[i][celli] * gammas[i][celli] * p0.value() * Foam::exp(Latentheatgas1 * molarmass1 * ((T[celli] - Tvap1) / (R * T[celli] * Tvap1)));
                Qv_alpha[celli] += Xs[i][celli] * gammas[i][celli] * evap_coeff * Latentheatgas1 * molarmass1 * p0.value() * Foam::exp(Latentheatgas1 * molarmass1 * ((T[celli] - Tvap1) / (R * T[celli] * Tvap1))) / Foam::pow(2 * Pi * molarmass1 * R * T[celli], 0.5);
            }
        }
        // Sum the contributions from all elements
        pVap += recoil_alpha;
        Qv += Qv_alpha;
    }

    // Apply the recoil pressure ratio and subtract ambient pressure
    forAll(mesh.C(), celli)
    {
        if (T[celli] >= Tvaps && alpha1[celli] >= 0.001)
        {
            scalar recoil_coeff = Alpha[celli];
            pVap[celli] = recoil_coeff * pVap[celli] - p0.value();
            if (pVap[celli] < 0.0)
            {
                pVap[celli] = 0.0;
            }
        }
    }


    // 4. Calculate multicomponent evaporation flux
    total_massflux *= 0;
    for (int i = 0; i < n_comp; i++)
    {
        scalar Latentheatgas1 = Latentheatgas_list[i];
        scalar Tvap1 = Tvap_list[i];
        scalar molarmass1 = molarmass_list[i];
        scalar rho1 = rho_list[i];
        scalar R = constant::physicoChemical::R.value();
        scalar Pi = constant::mathematical::pi;

        // Reset the evaporation flux for each element
        massfluxs[i] *= 0;

        // Calculate the volumetric evaporation flux of each element
        forAll(mesh.C(), celli)
        {
            if (T[celli] >= Tvaps && alpha1[celli] >= 0.001)
            {
                scalar evap_coeff = Phi[celli];
                massfluxs[i][celli] += evap_coeff * Xs[i][celli] * gammas[i][celli] * molarmass1 * p0.value() * Foam::exp(Latentheatgas1 * molarmass1 * ((T[celli] - Tvap1) / (R * T[celli] * Tvap1))) / Foam::pow(2 * Pi * molarmass1 * R * T[celli], 0.5) / rho1;
                total_massflux[celli] += massfluxs[i][celli];
            }
        }
    }
}