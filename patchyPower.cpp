#include <iostream>
#include <sstream>
#include <iomanip>
#include "cosmology.hpp"
#include "galaxy.hpp"
#include "densityField.hpp"
#include "overDensityField.hpp"
#include "patchyRead.hpp"
#include "powerSpectrum.hpp"
#include <omp.h>

std::string fileName(std::string base, int N, std::string ext){
    std::stringstream name;
    name << base << std::setfill('0') << std::setw(4) << N << ext;
    return name.str();
}

int main() {
    double H_0 = 100.0; // Hubble constant in km/s/Mpc
    double Omega_M = 0.3111;
    double Omega_L = 0.6889;

    std::vector<int> N = {512,1024,512};
    std::vector<double> L = {1792, 3584, 1792};
    std::vector<double> r_min(3), r_max(3);

    cosmology cosmo(H_0, Omega_M, Omega_L);

    std::string ranFile = "Patchy-Mocks-Randoms-DR12-COMPSAM_V6C_x100";
    std::string mockBase = "Patchy-Mocks-";
    std::string mockExt = ".dat";

    std::cout << "Reading randoms file..." << std::endl;
    std::vector<galaxy> rans = readPatchyRandoms(ranFile, cosmo, r_min, r_max);

    std::cout << "    Minimum position: (" << r_min[0] << ", " << r_min[1] << ", " << r_min[2] << ")\n";
    std::cout << "    Maximum position: (" << r_max[0] << ", " << r_max[1] << ", " << r_max[2] << ")\n";
    std::cout << "    Extent of box:    (" << r_max[0] - r_min[0] << ", " << r_max[1] - r_min[1] << ", "
              << r_max[2] - r_min[2] << ")\n";

    r_min[0] -= (L[0] - r_max[0] + r_min[0])/2;
    r_min[1] -= (L[1] - r_max[1] + r_min[1])/2;
    r_min[2] -= (L[2] - r_max[2] + r_min[2])/2;
    std::cout << "Minimum position: (" << r_min[0] << ", " << r_min[1] << ", " << r_min[2] << ")\n";

    std::cout << "Binning randoms..." << std::endl;
    densityField n_ran(N,L);
    n_ran.binCIC(rans, r_min);

    for (int mock = 0; mock < 2048; mock++){

    }
}
