#ifndef _PATCHY_READ_HPP
#define _PATCHY_READ_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "galaxy.hpp"
#include "cosmology.hpp"

std::vector<galaxy> readPatchyMock(std::string file, cosmology &cosmo, std::vector<double> &r_min,
                                    std::vector<double> &r_max) {
    std::ifstream fin(file, std::ios::in);
    double ra, dec, redshift, mstar, nbar, bias, veto, fiber_col;
    std::vector<galaxy> gals;
    while (!fin.eof()) {
        fin >> ra >> dec >> redshift >> mstar >> nbar >> bias >> veto >> fiber_col;
        if (redshift >= 0.4 && redshift < 0.7) {
            double w = (veto*fiber_col)/(1.0 + 10000.0*nbar);
            galaxy gal(ra, dec, redshift, w, nbar);
            gal.astroToCart(cosmo);
            gals.push_back(gal);

            std::vector<double> pos = gal.getCartPos();
            for (int j = 0; j < 3; ++j) {
                if (pos[j] < r_min[j]) r_min[j] = pos[j];
                if (pos[j] > r_max[j]) r_max[j] = pos[j];
            }
        }
    }
    fin.close();
    return gals;
}

std::vector<galaxy> readPatchyRandoms(std::string file, cosmology &cosmo, std::vector<double> &r_min,
                                    std::vector<double> &r_max) {
    std::ifstream fin(file, std::ios::in);
    if (!fin) {
        std::cout << "Error: Randoms file did not open.\n";
        std::cout << "Check that file " << file << " exists in the current directory.\n";
    }
    double ra, dec, redshift, mstar, nbar, bias, veto, fiber_col;
    std::vector<galaxy> rans;
    // long int count = 1;
    while (!fin.eof()) {
        fin >> ra >> dec >> redshift >> nbar >> bias >> veto >> fiber_col;
        // std::cout << count << " " << redshift << "\r";
        // count++;
        if (redshift >= 0.4 && redshift < 0.7) {
            double w = (veto*fiber_col)/(1.0 + 10000.0*nbar);
            galaxy ran(ra, dec, redshift, w, nbar);
            ran.astroToCart(cosmo);
            rans.push_back(ran);
            //std::cout << "Adding random point: " << rans.size() << "\n";

            std::vector<double> pos = ran.getCartPos();
            for (int j = 0; j < 3; ++j) {
                if (pos[j] < r_min[j]) r_min[j] = pos[j];
                if (pos[j] > r_max[j]) r_max[j] = pos[j];
            }
        }
    }
    fin.close();
    return rans;
}

#endif
