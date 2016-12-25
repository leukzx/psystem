#include <iostream>
#include "psystem.h"


int Particle::pCount = 0;

int main(int argc, char **argv)
{
    Eigen::initParallel();
    PSystem psystem;
    /*psystem2.timeStep = 0.01;
    psystem2.endTime = 10;
    psystem2.eps = 0;
    psystem2.G = 1;
    psystem2.methodFunction = psystem2.methodNameToMethodPtr("RK4");
    psystem2.outFile = "out2.dat";

    Eigen::Vector3d r, v;
    double mass = 1;
    r << 0.5, 0, 0;
    v << 0, 7.071067811865475e-01, 0;
    psystem2.particles.emplace_back(mass, r, v);

    r << -0.5, 0, 0;
    v << 0, -7.071067811865475e-01, 0;
    psystem2.particles.emplace_back(mass, r, v);*/
    //psystem.readParams(argv[1]);
    try {
        psystem.readConfig(argv[1]);
    }
    catch (const libconfig::ConfigException &cfcex) {
        std::cout << "Bad settings file. Terminating."
                  << std::endl;
        return(EXIT_FAILURE);
    }


    //std::cout << psystem2.EpotPtp(psystem2.particles[0], psystem2.particles[1]);
    psystem.evolve();
    std::cout << psystem.parameters() << std::endl;
    return 0;
}
