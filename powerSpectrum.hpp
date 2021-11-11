#ifndef _POWER_SPECTRUM_HPP_
#define _POWER_SPECTRUM_HPP_

#include <fstream>
#include <vector>
#include <string>
#include "overDensityField.hpp"

class powerSpectrum{
    std::vector<double> P, k;
    std::vector<int> Nk;
    double k_min, k_max, delta_k;

public:
    powerSpectrum(double k_min, double k_max, double delta_k);

    void calculate(overDensityField &delta, std::vector<int> N, double norm);

    void writeFile(std::string file);
};

powerSpectrum::powerSpectrum(double k_min, double k_max, double delta_k) {
    int N = (k_max - k_min)/delta_k;
    this->k_min = k_min;
    this->k_max = k_max;
    this->delta_k = delta_k;
    this->P = std::vector<double>(N);
    this->k = std::vector<double>(N);
    this->Nk = std::vector<int>(N);
}

void powerSpectrum::calculate(overDensityField &delta, std::vector<int> N, double norm) {
    std::cout << "    Binning frequencies..." << std::endl;
    for (int i = 0; i < N[0]; ++i) {
        for (int j = 0; j < N[1]; ++j) {
            for (int k = 0; k <= N[2]/2; ++k) {
                double k_mag = delta.getFreq(i,j,k);

                if (k_mag >= this->k_min && k_mag < this->k_max) {
                    int bin = (k_mag - this->k_min)/this->delta_k;
                    std::vector<double> dk = delta.getDeltaK(i,j,k);
                    this->P[bin] += (dk[0]*dk[0] + dk[1]*dk[1])*delta.gridCorrection(i,j,k) - delta.shotNoise;
                    this->k[bin] += k_mag;
                    this->Nk[bin] += 1;
                }
            }
        }
    }

    std::cout << "    Normalizing..." << std::endl;
    for (int i = 0; i < this->P.size(); ++i) {
        if (this->Nk[i] > 0) {
            P[i] /= (Nk[i]*norm);
            k[i] /= Nk[i];
        }
    }
}

void powerSpectrum::writeFile(std::string file) {
    std::ofstream fout(file);
    fout.precision(15);
    for (int i = 0; i < this->P.size(); ++i) {
        fout << this->k[i] << " " << this->P[i] << "\n";
    }
    fout.close();
}

#endif
