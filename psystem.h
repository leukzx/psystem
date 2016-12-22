#include <iostream>
#include <cmath>
//#include <eigen3/Eigen/Dense>
#include "/usr/include/eigen3/Eigen/Dense"
#include <fstream>      // Stream class to both read and write from/to files.
#include <string>       // std::string, std::to_string
#include <iomanip>      // std::setprecision
#include <cfloat>       // DBL_DIG (for output precision)
#include <vector>
#include "intersections/intersections.h"

class boundingBox { //Axis-aligned bounding box
 public:
    Eigen::Vector3d minVertex;
    Eigen::Vector3d maxVertex;


    void setBoundingBox(std::vector<Convex>); //Set minimum bounding box
    void setBoundingBox(Eigen::Vector3d, Eigen::Vector3d); //Set bounding box

    boundingBox();
};

class Particle {
 public:
    int id;

    double m;       // Mass

    Eigen::Vector3d r;  // Radius vector
    Eigen::Vector3d v;  // Velocity

    Eigen::Vector3d a;  // Acceleration

    // Needed for tsRK4 methodFunction:
    Eigen::Vector3d r0; // Radius vector at previous time step. Also used at boundaries check.
    std::vector<Eigen::Vector3d> A = std::vector<Eigen::Vector3d>(2); //Acceleration at intermediate time steps

    static int pCount;

    Particle();// Default constructor
    Particle(double, Eigen::Vector3d, Eigen::Vector3d);
    Particle(int, double, Eigen::Vector3d, Eigen::Vector3d);
    //~Particle(); // Destructor

    Particle(const Particle &x): id(x.id), m(x.m), r(x.r), v(x.v), a(x.a), r0(x.r0), A(x.A) {} // Copy constructor
    Particle & operator=(const Particle &); // Assignment operator
    Particle operator-(const Particle &); // Minus operator

    std::string state(int); // Particle's state string with specified precision
    std::string state(); // Particle's state string, max precision
    void info(); // Full state to std::cout
};

class PSystem {
 public:
    //Eigen::Vector3d dimXYZ; // PSystem dimensions
    std::vector<Particle> particles;

    std::vector<Convex> boundaries;

    boundingBox bBox; //Axis-aligned bounding box

    double EtotInit; //Initial total energy

    //double G = 6.6740831E-20; // (kg*km/s)*(km/kg)^2 , Gravitational constant
    double G;
    double eps; // Softening term for force calculation
    double timeStep;
    double time;
    double endTime;

    //External force fields
    Eigen::Vector3d g; //Uniform field coeff F= mg
    Eigen::Vector3d centralForceCoef; // Central force coeff F = m * centralForceCoef/r^2
    Eigen::Vector3d centralForceCenter;

    typedef  void (PSystem::*methodFunctionPTR)(double);
    typedef  Eigen::Vector3d (PSystem::*forceFunctionPTR)(const Particle& , const Particle&);
    typedef  double (PSystem::*energyFunctionPTR)(const Particle& , const Particle&);

    std::string outFile;
    int writePr; // Write precision
    double writeInterval; // Write interval
    
    std::string methodName;
    std::string potentialName;
    
    
    forceFunctionPTR forceFunction; //Pointer to force function

    // Estimate deltaT
    double estimateDeltaT();

    //Return particles to system boundaries
    void checkBoundsCyclic(); // Deprecated
    void checkBoundsHard(); // Deprecated
    void checkBoundaries();
    void checkBoundingBox();
    
    //Function takes integration methodFunction name and returns pointer to function
    void (PSystem::*methodNameToMethodPtr(std::string))(double);
    
    //Parameters reading from files
    void readParams(std::string fileName);
    void readParticlesData(std::string fileName);
    void readBoundariesData(std::string fileName);

    //add particle
    void addParticle(double, Eigen::Vector3d, Eigen::Vector3d);
    void addParticle(Particle&);
    void addParticle();


    //particle initialization with random parameters
    void setRandom(Particle&);
    void setRandomAll();

    //particle2 to particle1 force
    Eigen::Vector3d ptpForce(const Particle& , const Particle&);
    Eigen::Vector3d gravityForce(const Particle& , const Particle&);
    Eigen::Vector3d ljForce(const Particle& , const Particle&);

    //External field forces
    Eigen::Vector3d extFieldUniform(Particle&);
    Eigen::Vector3d extFieldCentral(Particle&);

    // Calculates accelerations of all particles
    void calcAccel();

    //Energy of the system calculation
    
    double gravityEpot(const Particle& p1, const Particle& p2);
    double ljEpot(const Particle& p1, const Particle& p2);
    double EpotPtp(Particle&, Particle&);
    double Epot();
    double Ekin();
    double Etot();
    void EtotInitSet();

    // Integration methods
    void tsForwardEuler(double);
    void tsLeapfrog(double);
    void tsRK4(double);

    // Advance in time
    void evolve();
    std::string stateE(int); // energy state of the system string
    std::string stateP(int); // all particles state string
    std::string parameters();
    void stateToFile(int, std::string); // writing system's state to file
    void info(int); //writing system's state to std::cout



    PSystem();
    PSystem(const PSystem &); //  Assignment copy
    PSystem(PSystem &&); // Move constructor
    PSystem & operator=(const PSystem &); // Assignment copy
    PSystem & operator=(PSystem &&); // Assignment move
    // Minus operator returns empty system if number of particles
    // in both systems are not equal.
    PSystem operator-(const PSystem &);
    
    double norm(); // Norm of system's vector in phase space
};
