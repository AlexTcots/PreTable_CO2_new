//
// Created by alext on 10/18/2022.
//

#include "species.h"

species::species() : THETA_(0), m_r_mole_(0), Morse_mole_(0), Xi_mole_(0),
                     Omega_(0), Omegax_(0) {
    std::cout << "species has been created with zero values!" << std::endl;
}

species::species(double THETA, double m_r_mole, double Morse_mole, double Xi_mole,
                 double Omega, double Omegax) :
        THETA_(THETA), m_r_mole_(m_r_mole), Morse_mole_(Morse_mole), Xi_mole_(Xi_mole),
        Omega_(Omega), Omegax_(Omegax) {
    std::cout << "species has been created ! " << std::endl;

}

species::species(const species& other):
THETA_(other.THETA_),
m_r_mole_(other.m_r_mole_),
Morse_mole_(other.Morse_mole_),
Xi_mole_(other.Xi_mole_),
Omega_(other.Omega_),
Omegax_(other.Omegax_){

}

const double &species::GetPi() const {
    return Pi;
}

const double &species::GetBOLTZ() const {
    return BOLTZ;
}

const double &species::GetH() const {
    return h;
}

const double &species::GetH_bar() const {
    return h_bar;
}

const double &species::GetC() const {
    return c;
}

const double &species::GetTHETA() const {
    return THETA_;
}

const double &species::GetM_R_mole() const {
    return m_r_mole_;
}

const double &species::GetMorse_mole() const {
    return Morse_mole_;
}

const double &species::GEtXi_mole() const {
    return Xi_mole_;
}

const double &species::GetOmega() const {
    return Omega_;
}

const double &species::GetOmegax() const {
    return Omegax_;
}

const double &species::GetETA_HS() const {
    return ETA_HS;
}

const double &species::GetETA() const {
    return ETA;
}
