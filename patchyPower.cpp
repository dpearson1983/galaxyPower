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

    std::string ranFile = "Patchy-Mocks-Randoms-DR12-COMPSAM_V6C_x100.dat";
    std::string mockBase = "Patchy-Mocks-DR12NGC-COMPSAM_V6C_";
    std::string mockExt = ".dat";
    std::string outBase = "Patchy_Pk_";
    std::string outExt = ".dat";

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

    std::vector<double> rMinDumb(3), rMaxDumb(3);

    std::cout << "Starting processing of mock catalogs...\n";
    for (int mock = 0; mock < 2048; mock++){
        std::string mockFile = fileName(mockBase, mock, mockExt);
        std::string outFile = fileName(outBase, mock, outExt);
        std::cout << "Processing mock: " << mockFile << std::endl;
        std::cout << "    Reading in `galaxies'...\n";
        std::vector<galaxy> gals = readPatchyMock(mockFile, cosmo, rMinDumb, rMaxDumb);

        std::cout << "    Binning galaxies...\n";
        densityField n_gals(N,L);
        n_gals.binCIC(gals, r_min);

        std::cout << "    Calculating overdensity field...\n";
        overDensityField delta(N, L);
        delta.calculate(n_gals, n_ran);

        std::cout << "    Fourier transforming overdensity field...\n";
        delta.fourierTransform();

        std::cout << "    Setting up power spectrum...\n";
        powerSpectrum Pk(0.012, 0.3, 0.008);

        std::cout << "    Calculating power spectrum...\n";
        double start = omp_get_wtime();
        Pk.calculate(delta, N, n_gals.getNbw()[2]);
        double end = omp_get_wtime();
        std::cout << "        Time to calculate: " << end - start << std::endl;

        std::cout << "    Writing power spectrum to file " << outFile << "...\n";
        Pk.writeFile(outFile);
    }
}
