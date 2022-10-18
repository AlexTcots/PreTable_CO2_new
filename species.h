//
// Created by alext on 10/18/2022.
//

#ifndef PRETABLE_CO2_NEW_SPECIES_H
#define PRETABLE_CO2_NEW_SPECIES_H

#include <iostream>


class species {
private:
    const double Pi = 3.14159265359;
    const double BOLTZ = 1.38064852e-23;
    //const Planck
    const double h = 6.6162e-34;
    //h_bar;
    const double h_bar = 1.0545718e-34;
    // speed of light
    const double c = 2.99792458e8;//in meter
    const double THETA_;
    const double m_r_mole_;
    const double Morse_mole_;
    const double Xi_mole_;
    const double Omega_;
    const double Omegax_;
    const double ETA_HS = 0.25;
    const double ETA = 0.25;
public:
    species();

    species(double THETA, double m_r_mole, double Morse_mole, double Xi_mole, double Omega, double Omegax);

    species(const species &other);

    //const method;
    const double &GetPi() const;

    const double &GetBOLTZ() const;

    const double &GetH() const;

    const double &GetH_bar() const;

    const double &GetC() const;

    const double &GetTHETA() const;

    const double &GetM_R_mole() const;

    const double &GetMorse_mole() const;

    const double &GEtXi_mole() const;

    const double &GetOmega() const;

    const double &GetOmegax() const;

    const double &GetETA_HS() const;

    const double &GetETA() const;


};


#endif //PRETABLE_CO2_NEW_SPECIES_H
