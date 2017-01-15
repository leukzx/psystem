#include <iostream>
#include "psystem.h"


int Particle::pCount = 0;

int main(int argc, char **argv)
{

    if( argc < 3) {
      std::cerr
        << "Use: settingsFile {cpu|opencl} kernelSourceFile"
        << std::endl;
      exit(EXIT_FAILURE);
    }

    Eigen::initParallel();
    PSystem psystem;

    try {
        psystem.readConfig(argv[1]);
    }
    catch (const libconfig::ConfigException &cfcex) {
        std::cout << "Bad settings file. Terminating."
                  << std::endl;
        return(EXIT_FAILURE);
    }

    const std::string platformName(argv[2]);

    if (platformName.compare("cpu")) {
        int deviceType = platformName.compare("openclgpu")?
                    CL_DEVICE_TYPE_CPU:CL_DEVICE_TYPE_GPU;
            psystem.evolveOpenCL(deviceType, argv[3]);
    } else {
        psystem.evolve();
    }

    std::cout << psystem.parameters() << std::endl;
    return 0;
}
