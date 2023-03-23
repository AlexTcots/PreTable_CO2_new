#include <iostream>
#include "species.h"
#include "FHO.h"
#include "random.h"

int main() {
    const double omega_CO2_symtric = 148000;
    const double omega_CO2_bending = 52600;
    const double omega_CO2_asymtric = 256500;
    const double Morse_CO2 = 4.3e10;
    const double m_r_CO2 = 3.655e-26;
    const double m_r_o_CO2_symmetric =  2.656780857176068e-26;
    const double m_r_o_CO2_bending = 3.625330649210353e-27;
    const double m_r_o_CO2_asymmetric = 7.250661298420705e-27;
    // 2xi = m_r/m_r_o
    const double Xi_CO2_bending = 0.5*(m_r_CO2/m_r_o_CO2_bending);
    const double Xi_CO2_symmetric = 0.5 * (m_r_CO2 / m_r_o_CO2_symmetric);
    const double Xi_CO2_asymmetric = 0.5 * (m_r_CO2 / m_r_o_CO2_asymmetric);
    const double THETA_CO2_symetric = 1918.6;
    const double THETA_CO2_bending = 959;
    const double THETA_CO2_asymtric = 3382.6;

    const double h = 6.6162e-34;
    // speed of light
    const double c = 2.99792458e8;//in meter

    double E_min = 1.5 * h * c * 1e5;
    double E_max = 1.5 * h * c * 1e8;
    int Nbins = 500;
    double seed = 1991.1112;
    int MS = 1000'000;

    const int Num_core = 16; // Use Num_core theads for OpenMP
    std::cout << "XI symmetric = " << Xi_CO2_symmetric << '\n';
    std::cout << "XI bending = " << Xi_CO2_bending << '\n';
    std::cout << "XI asymmetric = " << Xi_CO2_asymmetric << '\n';

    species CO2_symeric(THETA_CO2_symetric, m_r_CO2, Morse_CO2, Xi_CO2_symmetric, omega_CO2_symtric, 0.0);
    species CO2_bending(THETA_CO2_bending, m_r_CO2, Morse_CO2, Xi_CO2_bending, omega_CO2_bending, 0.0);
    //species CO2_asymtric(THETA_CO2_asymtric, m_r_CO2, Morse_CO2, Xi_CO2_asymmetric, omega_CO2_asymtric, 0.0);
    FHO Table_co2_m1(CO2_symeric, E_min, E_max, Nbins, seed, MS, Num_core);
    Table_co2_m1.ApplyComputing("co2_symmetric");

    FHO Table_co2_m2(CO2_bending, E_min, E_max, Nbins, seed, MS,Num_core);
    Table_co2_m2.ApplyComputing("co2_bending");

    //FHO Table_co2_m3(CO2_asymtric, E_min, E_max, Nbins, seed, MS, Num_core);
    //Table_co2_m3.ApplyComputing("co2_asymmetric");
}
