#include <iostream>
#include "species.h"
#include "FHO.h"
#include "random.h"
int main() {
    const double omega_CO2_symtric = 148000;
    const double omega_CO2_bending = 52600;
    const double omega_co2_asymtric = 256500;
    const double Morse_CO2 = 4.3e10;
    const double m_r_CO2 = 3.655e-26;
    const double m_r_o_CO2 =0.725e-26;
    const double Xi_CO2 = 2.5;
    const double THETA_CO2_symetric = 1918.6;
    const double THETA_CO2_bending = 959;
    const double THETA_CO2_asymtric = 3382.6;

    const double h = 6.6162e-34;
    // speed of light
    const double c = 2.99792458e8;//in meter

    double E_min = 1.5 * h * c * 1e5;
    double E_max = 1.5 * h * c * 1e8;
    int Nbins =500;
    double seed = 1991.1112;
    int MS = 1000'000;

  species CO2(THETA_CO2_symetric,m_r_CO2,Morse_CO2,Xi_CO2,omega_CO2_symtric,0.0);
  FHO Table_co2(CO2,E_min,E_max, Nbins,seed,MS);
  Table_co2.ApplyComputing();
}
