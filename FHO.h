//
// Created by alext on 10/18/2022.
//

#ifndef PRETABLE_CO2_NEW_FHO_H
#define PRETABLE_CO2_NEW_FHO_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <omp.h>
#include "species.h"
#include "random.h"


class FHO {

    species gas_;
    std::fstream Ebins_output;
    std::fstream E_w_output;
    std::fstream Table_output;
    const double E_min_;
    const double E_max_;
    const int Nbins_;
    std::vector<double> Ebins_;
    std::vector<double> E_w_;
    const double rd_seed_;
    RanPark rd;
    const int MS_;
    const int Num_core_;

    //private methods;
    double Enb(double e_min, double e_max, int Nbins, int nbin);

    double Vibr_eng(int vib_lev);// Kustova's book

    double s_func(int i, int f);//Adamovich

    double Ns(int i, int f);//Adamovich

    double ave_vib_qua(double e1, double e2, double s);//Adamovich

    double theta_p(double w);


    double u_from_E(double E);


    std::pair<double, double> Epsilon_Generator();

    double gam(double theta1, double theta2, double phi1, double phi2, std::pair<double, double> e12, double y);

    double Q_func(double theta_prime, double theta1, double phi1, double theta, double w, double u, double gamma);


//modify
    double Compute_pro_if_deex(int n, int i, int f);

    double Compute_pro_if_exci(int n, int i, int f);


public:
    FHO();

    FHO(const species &gas, double E_min, double E_max, int Nbins, double seed, int MS,int Num_core);

    void ApplyComputing(const std::string& name);


};


#endif //PRETABLE_CO2_NEW_FHO_H
