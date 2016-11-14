#include <iostream>
#include "psystem.h"
int Particle::pCount = 0;
int main(int argc, char **argv)
{
    if (argc == 3) { // argc should be 2 filenames for correct execution
        std::cout << std::scientific;
        std::cout << std::setprecision(DBL_DIG);
        
        PSystem psystem, psystem1, psystem2;
        psystem1.readParticlesData(argv[1]);
        psystem2.readParticlesData(argv[2]);
        
        psystem = psystem1 - psystem2;
        

        
        if (psystem.norm() != std::numeric_limits<double>::quiet_NaN()) {
            int psdimension;
            psdimension = 2 * psystem.dimXYZ.size() * psystem.particles.size();
        
            std::cout << "Phase space dimensions number = "; 
            std::cout << psdimension << std::endl;
        
            std::cout << "Distance between two systems = "; 
            std::cout << psystem.norm() << std::endl;
        } else {
            std::cout << "Can't compute distance between systems.\n";
            std::cout << "First system's dimensions number= ";
            std::cout << 2 * psystem1.dimXYZ.size() * psystem1.particles.size() << std::endl;
            std::cout << "Second system's dimensions number = ";
            std::cout << 2 * psystem2.dimXYZ.size() * psystem2.particles.size() << std::endl;
        }
    } else {
        std::cout << "Usage: " << argv[0] << " <filename1> <filename2>\n";
        std::cout << "This program computes distance between two systems of particles in 6N-dimensional phase space\n";
        std::cout << "<filename1> is the particles data of the first system\n";
        std::cout << "<filename1> is the particles data of the second system\n";
        std::cout << "Two systems should have equal numbers of particles.\n";
    }
    return 0;
}
