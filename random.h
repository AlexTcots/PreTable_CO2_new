//
// Created by alext on 10/18/2022.
//

#ifndef PRETABLE_CO2_NEW_RANDOM_H
#define PRETABLE_CO2_NEW_RANDOM_H


class RanPark{
public:
    RanPark(int);
    RanPark(double);
    ~RanPark() {}
    void reset(double, int, int);
    double uniform();
    double gaussian();

private:
    int seed,save;
    double second;
};


#endif //PRETABLE_CO2_NEW_RANDOM_H
