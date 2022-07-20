//
// Created by Youhi Morii on 2021/06/21.
//

#include "mfr.hpp"
#include <boost/algorithm/string.hpp>

void MicroFlowReactor::WeakFlameStFlow::readTwFile(std::string fileName) {
    std::ifstream dataFile(fileName);
    std::string line;
    while (getline(dataFile, line)) {
        std::vector<std::string> list_string;
        boost::split(list_string, line, boost::is_space(), boost::algorithm::token_compress_on);
        PositionZ.push_back(stod(list_string[0]));
        TemperatureWall.push_back(stod(list_string[1]));
    }
}

double MicroFlowReactor::WeakFlameStFlow::getTwAt(double posZ) {
    return Cantera::linearInterp(posZ, PositionZ, TemperatureWall);
}

void MicroFlowReactor::WeakFlameStFlow::evalResidual(double *x, double *rsd, int *diag, double rdt, size_t jmin,
                                                     size_t jmax) {
    //----------------------------------------------------
    // evaluate the residual equations at all required
    // grid points
    //----------------------------------------------------

    // calculation of qdotRadiation

    // The simple radiation model used was established by Y. Liu and B. Rogg [Y.
    // Liu and B. Rogg, Modelling of thermally radiating diffusion flames with
    // detailed chemistry and transport, EUROTHERM Seminars, 17:114-127, 1991].
    // This model uses the optically thin limit and the gray-gas approximation
    // to simply calculate a volume specified heat flux out of the Planck
    // absorption coefficients, the boundary emissivities and the temperature.
    // The model considers only CO2 and H2O as radiating species. Polynomial
    // lines calculate the species Planck coefficients for H2O and CO2. The data
    // for the lines is taken from the RADCAL program [Grosshandler, W. L.,
    // RADCAL: A Narrow-Band Model for Radiation Calculations in a Combustion
    // Environment, NIST technical note 1402, 1993]. The coefficients for the
    // polynomials are taken from [http://www.sandia.gov/TNF/radiation.html].

    if (m_do_radiation) {
        // variable definitions for the Planck absorption coefficient and the
        // radiation calculation:
        doublereal k_P_ref = 1.0 * OneAtm;

        // polynomial coefficients:
        const doublereal c_H2O[6] = {-0.23093, -1.12390, 9.41530, -2.99880,
                                     0.51382, -1.86840e-5};
        const doublereal c_CO2[6] = {18.741, -121.310, 273.500, -194.050,
                                     56.310, -5.8169};

        // calculation of the two boundary values
        double boundary_Rad_left = m_epsilon_left * StefanBoltz * pow(T(x, 0), 4);
        double boundary_Rad_right = m_epsilon_right * StefanBoltz * pow(T(x, m_points - 1), 4);

        // loop over all grid points
        for (size_t j = jmin; j < jmax; j++) {
            // helping variable for the calculation
            double radiative_heat_loss = 0;

            // calculation of the mean Planck absorption coefficient
            double k_P = 0;
            // absorption coefficient for H2O
            if (m_kRadiating[1] != npos) {
                double k_P_H2O = 0;
                for (size_t n = 0; n <= 5; n++) {
                    k_P_H2O += c_H2O[n] * pow(1000 / T(x, j), (double) n);
                }
                k_P_H2O /= k_P_ref;
                k_P += m_press * X(x, m_kRadiating[1], j) * k_P_H2O;
            }
            // absorption coefficient for CO2
            if (m_kRadiating[0] != npos) {
                double k_P_CO2 = 0;
                for (size_t n = 0; n <= 5; n++) {
                    k_P_CO2 += c_CO2[n] * pow(1000 / T(x, j), (double) n);
                }
                k_P_CO2 /= k_P_ref;
                k_P += m_press * X(x, m_kRadiating[0], j) * k_P_CO2;
            }

            // calculation of the radiative heat loss term
            radiative_heat_loss = 2 * k_P * (2 * StefanBoltz * pow(T(x, j), 4)
                                             - boundary_Rad_left - boundary_Rad_right);

            // set the radiative heat loss vector
            m_qdotRadiation[j] = radiative_heat_loss;
        }
    }

    for (size_t j = jmin; j <= jmax; j++) {
        //----------------------------------------------
        //         left boundary
        //----------------------------------------------

        if (j == 0) {
            // these may be modified by a boundary object

            // Continuity. This propagates information right-to-left, since
            // rho_u at point 0 is dependent on rho_u at point 1, but not on
            // mdot from the inlet.
            rsd[index(c_offset_U, 0)] =
                    -(rho_u(x, 1) - rho_u(x, 0)) / m_dz[0]
                    - (density(1) * V(x, 1) + density(0) * V(x, 0));

            // the inlet (or other) object connected to this one will modify
            // these equations by subtracting its values for V, T, and mdot. As
            // a result, these residual equations will force the solution
            // variables to the values for the boundary object
            rsd[index(c_offset_V, 0)] = V(x, 0);
            if (doEnergy(0)) {
                rsd[index(c_offset_T, 0)] = T(x, 0);
            } else {
                rsd[index(c_offset_T, 0)] = T(x, 0) - T_fixed(0);
            }
            rsd[index(c_offset_L, 0)] = -rho_u(x, 0);

            // The default boundary condition for species is zero flux. However,
            // the boundary object may modify this.
            double sum = 0.0;
            for (size_t k = 0; k < m_nsp; k++) {
                sum += Y(x, k, 0);
                rsd[index(c_offset_Y + k, 0)] =
                        -(m_flux(k, 0) + rho_u(x, 0) * Y(x, k, 0));
            }
            rsd[index(c_offset_Y + leftExcessSpecies(), 0)] = 1.0 - sum;

            // set residual of poisson's equ to zero
            rsd[index(c_offset_E, 0)] = x[index(c_offset_E, j)];
        } else if (j == m_points - 1) {
            evalRightBoundary(x, rsd, diag, rdt);
            // set residual of poisson's equ to zero
            rsd[index(c_offset_E, j)] = x[index(c_offset_E, j)];
        } else { // interior points
            evalContinuity(j, x, rsd, diag, rdt);
            // set residual of poisson's equ to zero
            rsd[index(c_offset_E, j)] = x[index(c_offset_E, j)];

            //------------------------------------------------
            //    Radial momentum equation
            //
            //    \rho dV/dt + \rho u dV/dz + \rho V^2
            //       = d(\mu dV/dz)/dz - lambda
            //-------------------------------------------------
            rsd[index(c_offset_V, j)]
                    = (shear(x, j) - lambda(x, j) - rho_u(x, j) * dVdz(x, j)
                       - m_rho[j] * V(x, j) * V(x, j)) / m_rho[j]
                      - rdt * (V(x, j) - V_prev(j));
            diag[index(c_offset_V, j)] = 1;

            //-------------------------------------------------
            //    Species equations
            //
            //   \rho dY_k/dt + \rho u dY_k/dz + dJ_k/dz
            //   = M_k\omega_k
            //-------------------------------------------------
            getWdot(x, j);
            for (size_t k = 0; k < m_nsp; k++) {
                double convec = rho_u(x, j) * dYdz(x, k, j);
                double diffus = 2.0 * (m_flux(k, j) - m_flux(k, j - 1))
                                / (z(j + 1) - z(j - 1));
                rsd[index(c_offset_Y + k, j)]
                        = (m_wt[k] * (wdot(k, j))
                           - convec - diffus) / m_rho[j]
                          - rdt * (Y(x, k, j) - Y_prev(k, j));
                diag[index(c_offset_Y + k, j)] = 1;
            }

            //-----------------------------------------------
            //    energy equation
            //
            //    \rho c_p dT/dt + \rho c_p u dT/dz
            //    = d(k dT/dz)/dz
            //      - sum_k(\omega_k h_k_ref)
            //      - sum_k(J_k c_p_k / M_k) dT/dz
            //-----------------------------------------------
            if (m_do_energy[j]) {
                setGas(x, j);

                // heat release term
                const vector_fp &h_RT = m_thermo->enthalpy_RT_ref();
                const vector_fp &cp_R = m_thermo->cp_R_ref();
                double sum = 0.0;
                double sum2 = 0.0;
                for (size_t k = 0; k < m_nsp; k++) {
                    double flxk = 0.5 * (m_flux(k, j - 1) + m_flux(k, j));
                    sum += wdot(k, j) * h_RT[k];
                    sum2 += flxk * cp_R[k] / m_wt[k];
                }
                sum *= GasConstant * T(x, j);
                double dtdzj = dTdz(x, j);
                sum2 *= GasConstant * dtdzj;

                rsd[index(c_offset_T, j)] = -m_cp[j] * rho_u(x, j) * dtdzj
                                            - divHeatFlux(x, j) - sum - sum2;
                rsd[index(c_offset_T, j)] /= (m_rho[j] * m_cp[j]);
                rsd[index(c_offset_T, j)] -= rdt * (T(x, j) - T_prev(j));
                rsd[index(c_offset_T, j)] -= (m_qdotRadiation[j] / (m_rho[j] * m_cp[j]));
                diag[index(c_offset_T, j)] = 1;
                const double locZ = z(j);
                const double gasT = T(x, j);
                const double Tw = getTwAt(locZ);
                const double mfrterm = 4.0 * NusseltNumber * (Tw - gasT) / TubeDiameter / TubeDiameter;
                rsd[index(c_offset_T, j)] += mfrterm * m_tcon[j] / (m_rho[j] * m_cp[j]);
            } else {
                // residual equations if the energy equation is disabled
                rsd[index(c_offset_T, j)] = T(x, j) - T_fixed(j);
                diag[index(c_offset_T, j)] = 0;
            }

            rsd[index(c_offset_L, j)] = lambda(x, j) - lambda(x, j - 1);
            diag[index(c_offset_L, j)] = 0;
        }
    }
}
