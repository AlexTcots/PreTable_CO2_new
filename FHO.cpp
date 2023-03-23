//
// Created by alext on 10/18/2022.
//

#include "FHO.h"

FHO::FHO() : gas_(), E_min_(), E_max_(), Nbins_(), Ebins_(), E_w_(), rd_seed_(),
             rd(rd_seed_),MS_(),Num_core_() {
    std::cout << "FHO has been created with zero values!" << std::endl;

}

FHO::FHO(const species &gas, double E_min, double E_max, int Nbins, double seed,int MS,int Num_core) :
        gas_(gas), E_min_(E_min), E_max_(E_max), Nbins_(Nbins),
        rd_seed_(seed), rd(rd_seed_),MS_(MS),Num_core_(Num_core) {
    Ebins_.reserve(Nbins_);
    E_w_.reserve(Nbins_);
    std::cout << "FHO has been created!" << std::endl;

}

double FHO::Enb(double e_min, double e_max, int Nbins, int nbin) {
    double n = static_cast<double>(nbin);
    double Nb = static_cast<double>(Nbins - 1);
    double expent = n / Nb;

    return e_min * pow(e_max / e_min, expent);
}

double FHO::Vibr_eng(int vib_lev) {
    return gas_.GetH() * gas_.GetC() * (gas_.GetOmega() * (vib_lev + 0.5) - gas_.GetOmegax() * pow(vib_lev + 0.5, 2));
}

double FHO::s_func(int i, int f) {
    return static_cast<double>(abs(i - f));
}

double FHO::Ns(int i, int f) {
    int max = std::max(i, f);
    int min = std::min(i, f);
    double s = s_func(i, f);

    return pow(tgamma(max + 1) / tgamma(min + 1), 1. / s);
}

double FHO::ave_vib_qua(double e1, double e2, double s) {
    double deltaE_if = fabs(e1 - e2);

    return deltaE_if / (s * gas_.GetH_bar());
}

double FHO::theta_p(double w) {
    return (4 * pow(gas_.GetPi() * w, 2) * gas_.GetM_R_mole())
           / (pow(gas_.GetMorse_mole(), 2) * gas_.GetBOLTZ());
}

double FHO::u_from_E(double E) {
    return sqrt(2 * E / gas_.GetM_R_mole());
}

std::pair<double, double> FHO::Epsilon_Generator() {
    double eta = gas_.GetETA_HS();
    double xm = 1 / (2 - eta);
    double Fm = xm * pow(1 - xm, 1 - eta);
    double k1;
    double cond;
    bool flag = true;
    while (flag) {
        k1 = rd.uniform();
        cond = (k1 * pow(1 - k1, 1 - eta) / Fm);
        if (rd.uniform() < cond) {
            flag = false;
        }
    }
    double epsilon1 = k1 * rd.uniform();
    double epsilon2 = k1 - epsilon1;
    return std::make_pair(epsilon1, epsilon2);
}

double FHO::gam(double theta1, double theta2, double phi1, double phi2, std::pair<double, double> e12, double y) {
    double tmp = -0.5 * sin(2 * theta1) * cos(phi1)
                 * sqrt(e12.first) - 0.5 * sin(2 * theta2) * cos(phi2)
                                     * sqrt(e12.second) + sqrt((1 - e12.first - e12.second) * (1 - y));

    return std::max(0.0, tmp);
}

double FHO::Q_func(double theta_prime, double theta1, double phi1, double theta, double w, double u, double gamma) {
    double Q;
    Q = (theta_prime * gas_.GEtXi_mole() * pow(cos(theta1), 2) * pow(cos(phi1), 2))
        / (4 * theta * pow(sinh(gas_.GetPi() * w / (gas_.GetMorse_mole() * u * gamma)), 2));
    return Q;
}

double
FHO::Compute_pro_if_deex(int n, int i, int f) {
    double E_t_r = Ebins_.at(n);
    double p;
    double u1, y1, s1, ns1, w1, gamma1, theta_prime1, theta, q1, p1;
    p = 0;
    theta = gas_.GetTHETA();
    u1 = y1 = s1 = ns1 = w1 = gamma1 = theta_prime1 = q1 = p1 = 0;
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(Num_core_); // Use Num_core_ threads for all consecutive parallel regions
#pragma omp parallel default(none) private (u1, y1, s1, ns1, w1, gamma1, theta_prime1, q1, p1) shared(E_t_r, p, rd, MS_, i, f, theta)
#pragma omp  for reduction (+:p)
    for (int j = 0; j < MS_; ++j) {

        auto e12 = Epsilon_Generator();
        double u1 = u_from_E(E_t_r);

        double e1 = Vibr_eng(i);
        double e2 = Vibr_eng(f);
        double y1 = rd.uniform();

        double theta1, theta2, phi1, phi2;
        theta1 = rd.uniform() * 2 * gas_.GetPi();
        theta2 = rd.uniform() * 2 * gas_.GetPi();
        phi1 = rd.uniform() * 2 * gas_.GetPi();
        phi2 = rd.uniform() * 2 * gas_.GetPi();

        double s1 = s_func(i, f);

        double ns1 = Ns(i, f);

        double w1 = ave_vib_qua(e1, e2, s1);

        double gamma1 = gam(theta1, theta2, phi1, phi2, e12, y1);

        double theta_prime1 = theta_p(w1);

        double q1 = Q_func(theta_prime1, theta1, phi1, theta, w1, u1, gamma1);

        double p1 = (pow(ns1, s1) / pow(tgamma(s1 + 1), 2)) * pow(q1, s1)
                    * exp((-2 * ns1 * q1) / (s1 + 1) -
                          (pow(ns1, 2) * pow(q1, 2)) / (pow(s1 + 1, 2) * (s1 + 2)));
        p += p1;

    }
    p /= MS_;
    return p;
}

double
FHO::Compute_pro_if_exci(int n, int i, int f) {
    double E_t_r = Ebins_.at(n);
    double p;
    double u1, y1, s1, ns1, w1, gamma1, theta_prime1, theta, q1, p1, E_t_f, E_prime, g_i, g_f, em1, em2, factor;
    std::pair<double, double> new_e12;
    int count = 0;
    p = 0;
    theta = gas_.GetTHETA();
    u1 = y1 = s1 = ns1 = w1 = gamma1 = theta_prime1 = q1 = p1 = E_t_f = E_prime = g_i = g_f = em1 = em2 = factor = 0;
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(Num_core_); // Use Num_core_ threads for all consecutive parallel regions
#pragma omp parallel default(none) private (u1, y1, s1, ns1, w1, gamma1, theta_prime1, q1, p1, E_t_f, E_prime, g_i, g_f, em1, em2, factor, new_e12) shared(E_t_r, p, rd, MS_, i, f, theta, count)
#pragma omp  for reduction (+:p, count)

    for (int j = 0; j < MS_; ++j) {

        auto e12 = Epsilon_Generator();


        double e1 = Vibr_eng(i);
        double e2 = Vibr_eng(f);
        double y1 = rd.uniform();

        double theta1, theta2, phi1, phi2;
        theta1 = rd.uniform() * 2 * gas_.GetPi();
        theta2 = rd.uniform() * 2 * gas_.GetPi();
        phi1 = rd.uniform() * 2 * gas_.GetPi();
        phi2 = rd.uniform() * 2 * gas_.GetPi();

        //from deex to excitation process

        double E_t_f = E_t_r * (1 - e12.first - e12.second) + e1 - e2;
        if (E_t_f > 0) {

            ++count;


            double E_prime = E_t_f + (e12.first + e12.second) * E_t_r;
            double g_i = sqrt(2 * (1 - e12.first - e12.second) * E_t_r / gas_.GetM_R_mole());
            double g_f = sqrt(2 * E_t_f / gas_.GetM_R_mole());
            double em1 = e12.first * (E_t_r / E_prime);
            double em2 = e12.second * (E_t_r / E_prime);
            new_e12 = std::make_pair(em1, em2);
            double factor = pow(g_f / g_i, 2 - 2 * gas_.GetETA());
            double u1 = u_from_E(E_prime);
            double s1 = s_func(i, f);

            double ns1 = Ns(i, f);

            double w1 = ave_vib_qua(e1, e2, s1);

            double gamma1 = gam(theta1, theta2, phi1, phi2, new_e12, y1);

            double theta_prime1 = theta_p(w1);

            double q1 = Q_func(theta_prime1, theta1, phi1, theta, w1, u1, gamma1);

            double p1 = (pow(ns1, s1) / pow(tgamma(s1 + 1), 2)) * pow(q1, s1)
                        * exp((-2 * ns1 * q1) / (s1 + 1) -
                              (pow(ns1, 2) * pow(q1, 2)) / (pow(s1 + 1, 2) * (s1 + 2)));
            p1 *= factor;
            p += p1;
        } else {

            p += 0;
        }


    }
    if (count != 0) {
        p /= count;
        std::cout << count << ' ';
    } else {
        p = 0;
    }
    return p;
}

void FHO::ApplyComputing(const std::string& name) {
    for (int i = 0; i < Nbins_; ++i) {
        double e = Enb(E_min_, E_max_, Nbins_, i);
        double w = 0.01 * e / (1.5 * gas_.GetH() * gas_.GetC());// in 1\cm
        Ebins_.push_back(std::move(e));
        E_w_.push_back(std::move(w));

    }
    Ebins_output.open("Ebins_"+name+".txt", std::fstream::out | std::fstream::trunc);
    //outfile2<<500<<'\n';


    E_w_output.open("Ebins_omega_"+name+".txt", std::fstream::out | std::fstream::trunc);

    //std::fstream outfile4;
    //outfile4.open("table_test.txt", std::fstream::out | std::fstream::trunc);

    Table_output.open("table_dsmc_"+name+".txt", std::fstream::out | std::fstream::trunc);


    //Energy bins done;
    for (const auto &E: Ebins_) {
        Ebins_output << std::setprecision(10) << E << std::endl;

    }

    for (const auto &Ew: E_w_) {
        E_w_output << std::setprecision(10) << Ew << std::endl;
    }
    int max_level = 50;
    int count = 0;
    for (int n = 0; n < Nbins_; ++n) {

        // consider max vib level is 50; i varify (0,...50)  51 levels

        //Calculating the look-up table


        for (int i = 0; i <= max_level; ++i) {

            double psum = 0;
            for (int f = 0; f <= max_level; ++f) {
                ++count;
                int s = s_func(i, f);
                if (s <= 10 && s != 0) {
                    //chose excitation or deexcitation
                    if (i > f) { //deexcitation
                        double p = Compute_pro_if_deex(n, i, f);
                        psum += p;
                        Table_output << n << ' ' << i << ' ' << f << ' '
                                     << std::setprecision(10)
                                     << p << ' '
                                     << std::endl;
                    } else if (i < f) {// excitation
                        double e1 = Vibr_eng(i);
                        double e2 = Vibr_eng(f);
                        double e_t_r = Ebins_.at(n);
                        // Trans+rot should greater than evib1 -evib2;
                        if (e_t_r > std::fabs(e1 - e2)) {
                            double p = Compute_pro_if_exci(n, i, f);
                            psum += p;
                            Table_output << n << ' ' << i << ' ' << f << ' '
                                         << std::setprecision(10)
                                         << p << ' '
                                         << std::endl;
                        } else {
                            double p = 0;
                            Table_output << n << ' ' << i << ' ' << f << ' '
                                         << std::setprecision(10)
                                         << p << ' '
                                         << std::endl;
                        }
                    }


                } else {
                    Table_output << n << ' ' << i << ' ' << f << ' '
                                 << std::setprecision(10)
                                 << 0.0 << ' '
                                 << std::endl;
                }

                std::cout << n << ' ' << i << ' ' << f << ' ' << "complete" << ' ' << count << std::endl;
            }
            Table_output << n << ' ' << i << ' ' << max_level + 1 << ' ' << std::setprecision(10) << psum << std::endl;
            std::cout << n << ' ' << i << ' ' << max_level + 1 << ' ' << "Psum complete" << std::endl;


        }



    }
}
