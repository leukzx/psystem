#include "psystem.h"


Particle::Particle()
{
    m = 1;
    r << 0, 0, 0;
    v << 0, 0, 0;
    a << 0, 0, 0;
    r0 << 0, 0, 0;
    for (auto& acc : A) { acc << 0, 0, 0;}
    id = Particle::pCount++;
}

Particle::Particle(double mass, Eigen::Vector3d rad, Eigen::Vector3d vel)
{
    m = mass;
    r = rad;
    v = vel;
    a << 0, 0, 0;
    r0 << 0, 0, 0;
    for (auto& acc : A) {acc << 0, 0, 0;}
    id = Particle::pCount++;
}

Particle::Particle(int idNum, double mass, Eigen::Vector3d rad, Eigen::Vector3d vel)
{
    id = idNum;
    m = mass;
    r = rad;
    v = vel;
    a << 0, 0, 0;
    r0 << 0, 0, 0;
    for (auto& acc : A) 
        acc << 0, 0, 0;
}

Particle& Particle::operator=(const Particle& rhp)  // Assignment operator
{
    if (this != &rhp) {
        id = rhp.id;
        m = rhp.m;
        r = rhp.r;
        v = rhp.v;
        a = rhp.a;
        r0 = rhp.r0;
        A = rhp.A;
    }
    return *this;
}

Particle Particle::operator-(const Particle &rhp)
{
    Particle p;

    p.id = (*this).id;
    p.m = (*this).m;
    p.r = (*this).r - rhp.r;
    p.v = (*this).v - rhp.v;

    return p;
}

PSystem::PSystem() // Default constructor
{
    G = 1;
    eps = 0; 
    timeStep = DBL_MAX;
    time = 0;
    endTime = 0;
    potentialName = "LJpotential";
    methodName = "RK4";
    outFile = "";
    writePr = DBL_DIG; 
    writeInterval = 0;
    EtotInit = 0;
    g << 0, 0, 0;
    centralForceCoef << 0, 0, 0;
    centralForceCenter << 0, 0, 0;
}

PSystem::PSystem(PSystem &&rhs) // Move constructor
{
    G = rhs.G;
    eps = rhs.eps;
    timeStep = rhs.timeStep;
    time = rhs.time;
    endTime = rhs.endTime;
    outFile = rhs.outFile;
    writePr = rhs.writePr;
    writeInterval = rhs.writeInterval;
    EtotInit = rhs.EtotInit;
    potentialName = rhs.potentialName;
    methodName = rhs.methodName;
    (*this).particles.swap(rhs.particles);
}

PSystem & PSystem::operator=(const PSystem &rhs) // Copy assignment
{
    G = rhs.G;
    eps = rhs.eps;
    timeStep = rhs.timeStep;
    time = rhs.time;
    endTime = rhs.endTime;
    outFile = rhs.outFile;
    writePr = rhs.writePr;
    writeInterval = rhs.writeInterval;
    EtotInit = rhs.EtotInit;
    particles = rhs.particles;
    potentialName = rhs.potentialName;
    methodName = rhs.methodName;
    return *this;
}

PSystem & PSystem::operator=(PSystem &&rhs) // Move assignment
{
    G = rhs.G;
    eps = rhs.eps;
    timeStep = rhs.timeStep;
    time = rhs.time;
    endTime = rhs.endTime;
    outFile = rhs.outFile;
    writePr = rhs.writePr;
    writeInterval = rhs.writeInterval;
    EtotInit = rhs.EtotInit;
    particles.swap(rhs.particles);
    potentialName = rhs.potentialName;
    methodName = rhs.methodName;
    return *this;
}

PSystem PSystem::operator-(const PSystem &rhs)
{
    PSystem psystem;
    if ((*this).particles.size() == rhs.particles.size()) {
        
        for (unsigned int i = 0; i < (*this).particles.size(); i++) {
            psystem.particles.emplace_back((*this).particles.at(i) - rhs.particles.at(i));
        }
        
        
        //~ for (auto& lhsParticle : (*this).particles {
            //~ for (auto& rhsParticle : rhs.particles) {
                //~ if (lhsParticle.id == rhsParticle.id) {
                    //~ psystem.particles.emplace_back(lhsParticle - rhsParticle);
                //~ }
            //~ }
        //~ }
    
    psystem.G = (*this).G;
    psystem.eps = (*this).eps;
    psystem.timeStep = (*this).timeStep;
    psystem.time = (*this).time;
    psystem.endTime = (*this).endTime;
    psystem.writePr = (*this).writePr;
    psystem.writeInterval = (*this).writeInterval;
    psystem.EtotInit = (*this).EtotInit;
    psystem.potentialName = (*this).potentialName;
    psystem.methodName = (*this).methodName;
    }

    return psystem;
}

double PSystem::norm()
{
    double norm = 0;
    
    if (particles.size() > 0) {
        for (auto & particle : particles) {
            norm += particle.r.squaredNorm() + particle.v.squaredNorm();
        }
        norm = sqrt(norm);
    } else {
        norm = std::numeric_limits<double>::quiet_NaN();
    }
    return norm;
}

std::string Particle::state(int p)
{
    std::streamsize defaultPrecision = std::cout.precision();
    int cWidth = 8; // Minimum column width is 8
    std::string delim = " "; // data delimeter

    cWidth += p;

    std::ostringstream out;

    out << std::scientific;
    out << std::setprecision(p);
    std::cout << std::right;

    out << id;
    out << delim;
    out << std::setw(cWidth) << r(0);
    out << delim;
    out << std::setw(cWidth) << r(1);
    out << delim;
    out << std::setw(cWidth) << r(2);
    out << delim;
    out << std::setw(cWidth) << v(0);
    out << delim;
    out << std::setw(cWidth) << v(1);
    out << delim;
    out << std::setw(cWidth) << v(2);
    out << delim;
    out << m;

    std::cout.unsetf(std::ios::fixed | std::ios::scientific);
    std::cout << std::setprecision(defaultPrecision);
    return out.str();

}

std::string Particle::state()
{
    return state(DBL_DIG);
}

void Particle::info()
{
    const int cWidth = 23;
    std::streamsize defaultPrecision = std::cout.precision();

    std::cout << std::scientific;
    std::cout << std::setprecision(DBL_DIG);
    std::cout << std::right;

    std::cout << "Particle's id: "<< id << std::endl;

    std::cout << "    Mass (m): ";
    std::cout << "\t";
    std::cout << std::setw(cWidth)<< m << std::endl;

    std::cout << "    Position (r): " << std::endl;
    std::cout << "\t";
    std::cout << std::setw(cWidth) << r(0);
    std::cout << " " << std::setw(cWidth) << r(1);
    std::cout << " " << std::setw(cWidth) << r(2) << std::endl;

    std::cout << "    Velocity (v): " << std::endl;
    std::cout << "\t";
    std::cout << std::setw(cWidth) << v(0);
    std::cout << " " << std::setw(cWidth) << v(1);
    std::cout << " " << std::setw(cWidth) << v(2) << std::endl;

    std::cout << "    Acceleration (a): " << std::endl;
    std::cout << "\t";
    std::cout << std::setw(cWidth) << a(0);
    std::cout << " " << std::setw(cWidth) << a(1);
    std::cout << " " << std::setw(cWidth) << a(2) << std::endl;

    std::cout.unsetf(std::ios::fixed | std::ios::scientific);
    std::cout << std::setprecision(defaultPrecision);
}

void PSystem::addParticle(double mass, Eigen::Vector3d rad, Eigen::Vector3d vel)
{
    particles.emplace_back(mass, rad, vel);
}

void PSystem::addParticle(Particle& particle) {
    particles.push_back(particle);
}

void PSystem::addParticle() {
    particles.emplace_back();
}

void PSystem::setRandom(Particle& particle)
{
    //set random coordinates
    //particle.r.setRandom();
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);
    for (int i = 0; i < particle.r.size(); i++) {
        particle.r(i)= distribution(generator);
    }
    particle.r = particle.r.array().abs();
    particle.r = bBox.minVertex
                 + particle.r.cwiseProduct(bBox.maxVertex - bBox.minVertex);

    //set random speed
    double lengthpct = 0.1; // max percent of psystem size per second to move
    //particle.v.setRandom();
    for (int i = 0; i < particle.v.size(); i++) {
        particle.v(i)= distribution(generator);
    }
    particle.v = particle.v.cwiseProduct((bBox.maxVertex - bBox.minVertex)
                                         * lengthpct);
}

void PSystem::setRandomAll()
{
    for (auto& particle : particles)
        setRandom(particle);
}

Eigen::Vector3d PSystem::ptpForce(const Particle& p1, const Particle& p2)
{
    
    Eigen::Vector3d force;
        
    forceFunctionPTR forceFunction = &PSystem::gravityForce;
    if (potentialName == "Gravitational") {forceFunction = &PSystem::gravityForce;}
    if (potentialName == "LJpotential") {forceFunction = &PSystem::ljForce;}
        
    force = (this->*forceFunction)(p1, p2);
    return force;
}

Eigen::Vector3d PSystem::gravityForce(const Particle& p1, const Particle& p2)
{
    Eigen::Vector3d force;
    Eigen::Vector3d ptpRv; //radius-vector of p2 relative to p1
    double ptpR;

    ptpRv = p2.r - p1.r;
    ptpR = ptpRv.squaredNorm() + pow(eps, 2);
    ptpR = pow(ptpR, 1.5);
    force = G * p1.m * p2.m * ptpRv / ptpR; // Gravity force
    return force;
}

Eigen::Vector3d PSystem::ljForce(const Particle& p1, const Particle& p2)
{
    Eigen::Vector3d force;
    Eigen::Vector3d r; // Radius-vector of p2 relative to p1
    double rm = 1; // Distance at which the potential reaches its minimum (F=0)
    double epsilon = 1;
    double rNorm; // Distance from p2 to p1

    r = p2.r - p1.r;
    rNorm = r.norm();

    force = 12
            * epsilon
            * (pow(rm, 6)/pow(rNorm,7) - pow(rm, 12)/pow(rNorm, 13))
            * r/rNorm;

    return force;
}

void PSystem::calcAccel()
{
    int pNum; // Number of particles
    Eigen::Matrix<Eigen::Vector3d, Eigen::Dynamic, Eigen::Dynamic> forcesM;
    //Eigen::Vector3d netForce;
    std::vector<Eigen::Vector3d> netForces;

    pNum = particles.size();
    forcesM.resize(pNum, pNum);
    #pragma omp parallel for
    for (int i = 0; i < pNum; ++i) {
        forcesM(i, i) << 0, 0, 0;
        for (int j = i + 1; j < pNum; ++j) {
            forcesM(i, j) = ptpForce(particles[i], particles[j]);
            forcesM(j, i) = - forcesM(i, j);
        }
    }
    netForces.resize(pNum, Eigen::Vector3d(0, 0, 0));
    //netForce = forcesM.rowwise().sum()(1,1); // error???
    #pragma omp parallel for
    for (int i = 0; i < pNum; ++i) {
        //netForce << 0, 0, 0;
        for (int j = 0; j < pNum; ++j){
            //netForce += forcesM(i, j);
            netForces.at(i) += forcesM(i, j);
        }
        //particles[i].a = netForce / particles[i].m;
        particles[i].a = netForces.at(i)/ particles[i].m;
               // + g
               // + centralForceCoef / ((particles[i].r - centralForceCenter).squaredNorm() + DBL_MIN);

    }
}

double PSystem::estimateDeltaT()
{
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> collisionTimeM;
    int pNum;
    double relr, relv, rela, dt, dtFrac;
    dtFrac = 0.01;
    pNum = particles.size();
    collisionTimeM.resize(pNum, pNum);
    collisionTimeM.fill(DBL_MAX);
    //#pragma omp parallel for
    for (int i = 0; i < pNum; i++) {
        //#pragma omp parallel for
        for (int j = i + 1 ; j < pNum; j++) {
            relr = (particles[i].r - particles[j].r).norm();
            relv = (particles[i].v - particles[j].v).norm();
            //rela = (ptpForce(particles[i],particles[j])/particles[i].m - ptpForce(particles[j],particles[i])/particles[j].m).norm();
            rela = (particles[i].a - particles[j].a).norm();
            //collisionTimeM(i,j) = relr / relv;
            collisionTimeM(i,j) =  std::min(relr / relv, sqrt(relr / rela));
        }
    }

    dt = dtFrac * collisionTimeM.minCoeff();

    return dt;
}

void PSystem::calcEkin()
{
    double Ek = 0;
    for (auto& particle : particles) {
        Ek += particle.m * particle.v.squaredNorm() / 2;
    }
    Ekin = Ek;
}

double PSystem::getEkin()
{
    return Ekin;
}

double PSystem::EpotPtp(Particle& p1, Particle& p2)
{
    energyFunctionPTR energyFunction = &PSystem::gravityEpot;
    if (potentialName == "Gravitational") {energyFunction = &PSystem::gravityEpot;}
    if (potentialName == "LJpotential") {energyFunction = &PSystem::ljEpot;}
    double Ep;
    Ep = (this->*energyFunction)(p1, p2);
    return Ep;
}

double PSystem::ljEpot(const Particle& p1, const Particle& p2)
{
    double Epot;
    //Eigen::Vector3d r; // Radius-vector of p2 relative to p1
    double r, rmOverRpow6, rmOverRPow12; // Distance from p2 to p1
    double rm = 1; // The distance at which the potential reaches its minimum (F=0))
    double epsilon = 1; // Depth of the potential well

    r = (p2.r - p1.r).norm();
    rmOverRpow6 = pow(rm/r, 6);
    rmOverRPow12 = pow(rmOverRpow6 , 2);

    //Epot = epsilon * (pow(rm/r, 12) - 2 * pow(rm/r, 6));
    Epot = epsilon * (rmOverRPow12 - 2 * rmOverRpow6);
    return Epot;
}

double PSystem::gravityEpot(const Particle& p1, const Particle& p2)
{
    double Epot;
    Eigen::Vector3d ptpRv;
    double ptpR;
    ptpRv = p2.r - p1.r;
    ptpR = ptpRv.squaredNorm() + pow(eps, 2);
    ptpR = pow(ptpR, 0.5);
    
    Epot = - G * p1.m * p2.m / ptpR;
    return Epot;
}

void PSystem::calcEpot()
{
    double Ep = 0;
    int pNum; //number of particles
    pNum = particles.size();
    #pragma omp parallel for
    for (int i = 0; i < pNum; ++i) {
        for (int j = i + 1; j < pNum; ++j) {
            Ep += EpotPtp(particles[i], particles[j]);
        }
    }

    Epot = Ep;
}

double PSystem::getEpot()
{
    return Epot;
}

void PSystem::calcEtot()
{
    calcEpot();
    calcEkin();
    Etot = Epot + Ekin;
}

double PSystem::getEtot()
{
    return Etot;
}

void PSystem::EtotInitSet()
{
    calcEtot();
    EtotInit = Etot;
}

void PSystem::tsForwardEuler(double dt)
{
    //calcAccel();
    for (auto& particle : particles){
        particle.r0 = particle.r; // save previous position of particle
        particle.r += particle.v * dt;
        particle.v += particle.a * dt;
    }
    calcAccel();
}

void PSystem::tsLeapfrog(double dt)
{
    for (auto& particle : particles) {
        particle.r0 = particle.r; // save previous position of particle
        particle.v += particle.a * dt / 2;
        particle.r += particle.v * dt;
    }

    calcAccel();

    for (auto& particle : particles) {
        particle.v += particle.a * dt / 2;
    }
}

void PSystem::tsRK4(double dt)
{
    //calcAccel(); // Shoud run calcAccel before first timestep!!!

    for (auto& p : particles) {
        p.A.at(0) = p.a;
        p.r0 = p.r; // saving position at the beginning of the time step
        p.r = p.r0 + p.v * dt / 2 + p.a * dt * dt / 8;
    }
    calcAccel();

    for (auto& p : particles) {
        p.A.at(1) = p.a;
        p.r = p.r0 + p.v * dt + p.a * dt * dt / 2;
    }
    calcAccel();

    for (auto& p : particles) {
        p.r = p.r0 + p.v * dt + (p.A.at(0) + 2 * p.A.at(1)) * dt * dt / 6;
        p.v += (p.A.at(0) + 4 * p.A.at(1) + p.a) * dt / 6;
    }

    calcAccel();
}


//void (PSystem::*PSystem::methodNameToMethodPtr(std::string methodName))(double)
PSystem::methodFunctionPTR PSystem::methodNameToMethodPtr(std::string mName)
{
    methodFunctionPTR mptr = &PSystem::tsRK4;
    if (mName == "ForwardEuler") {mptr = &PSystem::tsForwardEuler;}
    if (mName == "Leapfrog") {mptr = &PSystem::tsLeapfrog;}
    if (mName == "RK4") {mptr = &PSystem::tsRK4;}
    return mptr;
}

void PSystem::evolve()
{
    int nt = 0; //number of time steps
    double dt = timeStep;
    double timeStepInit = timeStep;
    int precision = 15;
    double writeTime = 0;
    methodFunctionPTR methodFunction = methodNameToMethodPtr(methodName);
    //double EPSILON;

    //EPSILON = timeStep;

    EtotInitSet();
    calcAccel(); // Initial accelerations of particals.
    
    dt = std::min({timeStepInit, estimateDeltaT(), writeInterval});
    timeStep = dt;
    stateToFile(writePr, outFile);
    info(precision);
    while ((time + dt) <= endTime) {
        (this->*methodFunction)(dt);        
        checkBoundaries();
        time += dt;
        writeTime += dt;
        nt++;
        calcEtot();
        if (std::fabs(writeInterval - writeTime) < std::numeric_limits<double>::epsilon() * dt) {
            stateToFile(writePr, outFile);
            writeTime = 0;
        }
        //dt = std::min(timeStepInit, estimateDeltaT());
        dt = std::min({timeStepInit, estimateDeltaT(), writeInterval});
        if ((writeTime + dt) > writeInterval) // Check dt for not to overshoot write time
            dt = writeInterval - writeTime;
        timeStep = dt;
        info(precision);
    }

    dt =  endTime - time; //Extra time step to endTime
    if (dt > 0) {
        timeStep = dt;
        (this->*methodFunction)(dt);
        checkBoundaries();
        time += dt;
        nt++;
        stateToFile(writePr, outFile);
        info(precision);
    }
    
    std::cout << "Number of time steps: " << nt << std::endl;
}


int PSystem::evolveOpenCL(int& deviceType, const char* kernelSourceFile)
{
    try {
        //int deviceType = CL_DEVICE_TYPE_GPU;

        std::vector<cl::Platform> platformList;

        cl::Platform platform;
        std::string platformName;

        std::vector<cl::Device> deviceList;
        cl::Device device;
        std::string deviceName;


        std::cout << "Available platforms:" << std::endl;
        cl::Platform::get(&platformList);
        for (auto& pltfm : platformList ) {
            pltfm.getInfo((cl_platform_info)CL_PLATFORM_NAME,
                             &platformName);
            std::cout << " " << platformName << std::endl;

            try {
                pltfm.getDevices(CL_DEVICE_TYPE_ALL, &deviceList);
                for(auto& deviceItem : deviceList) {
                    deviceItem.getInfo((cl_device_info)CL_DEVICE_NAME, &deviceName);
                    std::cout << "  " << deviceName << std::endl;

                }
            } catch (const cl::Error &error) {
                std::cerr << "  Error while getting device list: "
                          << error.what()
                          << '(' << error.err() << ')' << std::endl;
            }
        }

        /*

        // Select first gpu device in device list
        for (auto& pltfm : platformList){
            try {
                pltfm.getDevices(deviceType, &deviceList);
                platform = pltfm;
                device = deviceList.at(0);
                break;
            } catch (const cl::Error error) {
                continue;
            }
        }
        device.getInfo((cl_device_info)CL_DEVICE_NAME, &deviceName);
        std::cout << "Using device: " << deviceName << std::endl;
        */


        /*cl_context_properties contextProps[] =
            {CL_CONTEXT_PLATFORM,
             (cl_context_properties)(platform()), 0};
        */

        //The constructor attempts to use the first platform that
        //has a device of the specified type.
        //cl::Context context(deviceType, contextProps);
        //cl::Context context(device);
        cl::Context context(deviceType);

        std::vector<cl::Device> devices =
              context.getInfo<CL_CONTEXT_DEVICES>();

        device = devices.at(0);
        device.getInfo((cl_device_info)CL_DEVICE_NAME, &deviceName);
        std::cout << "Selected device: " << deviceName << std::endl;

        // Kernel source
        std::ifstream file(kernelSourceFile);
        std::string kernelSourceStr(std::istreambuf_iterator<char>(file),
                (std::istreambuf_iterator<char>()));
        file.close();

        // Creates vector of kernel sources. ?There are may be more then 1?
        cl::Program::Sources kernelSourceVec( 1, std::make_pair(kernelSourceStr.c_str(),
                                  kernelSourceStr.length()+1));
        // Creates an OpenCL program object for a context
        // and loads the source code specified by the text strings in each
        // element of the vector sources into the program object
        cl::Program program(context, kernelSourceVec);
        try {
            program.build(devices);
        } catch (cl::Error& error) {
            std::cerr << "Program build from kernel sources failed! "
                      << error.what()
                      << '(' << error.err() << ')' << std::endl;
            std::cerr << "Retrieving  logs ...:" << std::endl;
            for (auto& dvce : devices) {
                std::cerr << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(dvce)
                      << std::endl;
            }
            return EXIT_FAILURE;
        }



        // Prepares data to copy to device buffer
        size_t vLen = particles.size();
        size_t vSize = vLen * sizeof(cl_float4);

        cl_float4 *h_pos, *h_vel;
        h_pos = new cl_float4[vLen];
        h_vel = new cl_float4[vLen];
        cl_float *h_est;
        h_est = new cl_float[vLen * vLen];

        for (size_t i = 0; i < vLen; ++i) {
            for (unsigned int  j = 0; j < particles.at(i).r.size(); ++j) {
                 h_pos[i].s[j] = (float) particles.at(i).r(j);
                 h_vel[i].s[j] = (float) particles.at(i).v(j);
            }
            h_pos[i].s[3] = (float) particles.at(i).m;
            h_vel[i].s[3] = 0.0f;
        }

        // Array with all vertices of triangles in all boundaries

        cl_int vNum = 0; // Number of vertices in all boundaries
        for (auto& boundary : boundaries) {
            if (boundary.vertices.size() > 2) {
                // Number of triangles multiplied by 3 = number of vertices
                vNum += 3 * (boundary.vertices.size() - 3 + 1);
            }
        }
        cl_float4 *h_tri;
        h_tri = new cl_float4[vNum];

        // Fills h_tri array with vertices of all triangles of all boundaries
        size_t c = 0; // Vertices counter
        for (auto& boundary : boundaries) {
            for (unsigned int i = 1 ; i < boundary.vertices.size() - 1; i++) {
                unsigned int l[] = {0, i, i + 1};
                for (unsigned int k : l){
                    for (unsigned int j = 0;
                         j < boundary.vertices.at(k).size(); ++j) {
                        h_tri[c].s[j] = boundary.vertices.at(k)(j);
                    }
                    h_tri[c].s[3] = 0;
                    ++c;
                }
            }
        }

        cl_float dt  = (float) timeStep;

        cl_float *h_ekin, *h_epot_lj;
        h_ekin = new cl_float[vLen];
        h_epot_lj = new cl_float[vLen];

        // Creates device buffers
        cl::Buffer d_pos, d_vel, d_newPos, d_newVel, d_acc, d_est, d_tri,
                   d_ekin, d_epot_lj;
        d_pos = cl::Buffer(context, CL_MEM_READ_WRITE, vSize);
        d_vel = cl::Buffer(context, CL_MEM_READ_WRITE, vSize);
        d_newPos = cl::Buffer(context, CL_MEM_READ_WRITE, vSize);
        d_newVel = cl::Buffer(context, CL_MEM_READ_WRITE, vSize);
        d_acc = cl::Buffer(context, CL_MEM_READ_WRITE, vSize);
        d_est = cl::Buffer(context, CL_MEM_READ_WRITE,
                           vLen * vLen * sizeof(cl_float));
        d_tri = cl::Buffer(context, CL_MEM_READ_WRITE,
                           vNum * sizeof(cl_float4));
        d_ekin = cl::Buffer(context, CL_MEM_READ_WRITE, vLen * sizeof(cl_float));
        d_epot_lj = cl::Buffer(context, CL_MEM_READ_WRITE, vLen * sizeof(cl_float));

        // Creates a kernel objects
        cl::Kernel stepInTime(program, "tsForwardEulerCL");
        cl::Kernel estimateDt(program, "estimateDtCL");
        cl::Kernel calcAccCL(program, "calcAccCL");
        cl::Kernel checkBoundariesCL(program, "checkBoundariesCL");
        cl::Kernel calcEkinCL(program, "calcEkinCL");
        cl::Kernel calcEpotLJCL(program, "calcEpotLJCL");

        // Set stepInTime kernel arguments
        stepInTime.setArg(0, d_pos);
        stepInTime.setArg(1, d_vel);
        stepInTime.setArg(2, d_acc);
        stepInTime.setArg(3, dt);
        stepInTime.setArg(4, d_newPos);

        // Set estimateDt kernel arguments
        estimateDt.setArg(0, d_pos);
        estimateDt.setArg(1, d_vel);
        estimateDt.setArg(2, d_acc);
        estimateDt.setArg(3, d_est);

        // Set calcAccCL kernel arguments
        calcAccCL.setArg(0, d_pos);
        calcAccCL.setArg(1, d_acc);

        // Set checkBoundariesCL kernel arguments
        checkBoundariesCL.setArg(0, d_pos);
        checkBoundariesCL.setArg(1, d_vel);
        checkBoundariesCL.setArg(2, d_newPos);
        checkBoundariesCL.setArg(3, d_tri);
        checkBoundariesCL.setArg(4, vNum);

        // Set calcEkinCL kernel arguments
        calcEkinCL.setArg(0, d_pos);
        calcEkinCL.setArg(1, d_vel);
        calcEkinCL.setArg(2, d_ekin);

        // Set calcEpotLJCL kernel arguments
        calcEpotLJCL.setArg(0, d_pos);
        calcEpotLJCL.setArg(1, d_epot_lj);

        // Creates queue with profiling enabled
        cl::CommandQueue queue(context, device, CL_QUEUE_PROFILING_ENABLE);

        cl::Event event;

        //Enqueue a command to write to a buffer object from host memory
        queue.enqueueWriteBuffer(d_pos, CL_TRUE, 0, vSize, h_pos, NULL, &event);
        queue.enqueueWriteBuffer(d_vel, CL_TRUE, 0, vSize, h_vel, NULL, &event);
        queue.enqueueWriteBuffer(d_tri, CL_TRUE, 0,
                                 vNum * sizeof(cl_float4), h_tri, NULL, &event);

        // Initial acceleration
        queue.enqueueNDRangeKernel(calcAccCL,
                                   cl::NullRange, cl::NDRange(vLen),
                                   cl::NullRange, NULL, &event);
        event.wait();

        // Start time
        cl_ulong start =
                event.getProfilingInfo<CL_PROFILING_COMMAND_START>();

        // Initial timestep
        queue.enqueueNDRangeKernel(estimateDt,
                                   cl::NullRange, cl::NDRange(vLen, vLen),
                                   cl::NullRange, NULL, &event);
        event.wait();
        queue.enqueueReadBuffer(d_est, CL_TRUE, 0,
                                vLen * vLen * sizeof(cl_float), h_est);
        std::cout << dt << std::endl;
        dt = 0.01 * (*std::min_element(h_est, h_est + vLen * vLen));
        dt = std::min({(double) dt, timeStep});
        std::cout << dt << std::endl;

        std::cout << stateP(15) << std::endl;

        double writeTime = 0.0; // Write time counter
        EtotInitSet();
        for (size_t nt = 0; time + dt <= endTime; ++nt) {
            time += dt;
            // Move to new positions
            stepInTime.setArg(3, dt);
            queue.enqueueNDRangeKernel(stepInTime,
                                       cl::NullRange, cl::NDRange(vLen),
                                       cl::NullRange, NULL, &event);
            event.wait();

            // Check and resolve crossing of a boundaries
            queue.enqueueNDRangeKernel(checkBoundariesCL,
                                       cl::NullRange, cl::NDRange(vLen),
                                       cl::NullRange, NULL, &event);
            event.wait();

            // New positions to current positions
            queue.enqueueCopyBuffer(d_newPos, d_pos, 0, 0, vSize, NULL, &event);
            event.wait();

            // Calculate total energy
            queue.enqueueNDRangeKernel(calcEkinCL,
                                       cl::NullRange, cl::NDRange(vLen),
                                       cl::NullRange, NULL, &event);
            event.wait();
            queue.enqueueNDRangeKernel(calcEpotLJCL,
                                       cl::NullRange, cl::NDRange(vLen),
                                       cl::NullRange, NULL, &event);
            event.wait();
            queue.enqueueReadBuffer(d_ekin, CL_TRUE, 0,
                                    vLen *  sizeof(cl_float), h_ekin);
            queue.enqueueReadBuffer(d_epot_lj, CL_TRUE, 0,
                                    vLen *  sizeof(cl_float), h_epot_lj);
            Etot = 0;
            Epot = 0;
            Ekin = 0;
            for (size_t i = 0; i < vLen; ++i) {
                Epot += h_epot_lj[i];
            }
            Epot /= 2;
            for (size_t i = 0; i < vLen; ++i) {
                Ekin += h_ekin[i];
            }
            Etot = Ekin + Epot;

            // Write data to file
            writeTime += dt;
            if (std::fabs(writeInterval - writeTime)
                    < std::numeric_limits<float>::epsilon() * dt) {
                // Copy state of particles to host
                for (size_t i = 0; i < vLen; ++i) {
                    particles.at(i).r0 = particles.at(i).r;
                }
                queue.enqueueReadBuffer(d_pos, CL_TRUE, 0, vSize, h_pos);
                queue.enqueueReadBuffer(d_vel, CL_TRUE, 0, vSize, h_vel);

                for (size_t i = 0; i < vLen; ++i) {
                    for (unsigned int j = 0; j < particles.at(i).r.size();
                                                                          ++j) {
                         particles.at(i).r(j) = h_pos[i].s[j];
                         particles.at(i).v(j) = h_vel[i].s[j];
                    }
                }
                stateToFile(writePr, outFile);
                writeTime = 0;
            }
            info(15);

            // Get accelerations at new (it is now current) position
            queue.enqueueNDRangeKernel(calcAccCL,
                                       cl::NullRange, cl::NDRange(vLen),
                                       cl::NullRange, NULL, &event);
            event.wait();

            // Estimate time step for next iteration
            queue.enqueueNDRangeKernel(estimateDt,
                                       cl::NullRange, cl::NDRange(vLen, vLen),
                                       cl::NullRange, NULL, &event);
            event.wait();
            queue.enqueueReadBuffer(d_est, CL_TRUE, 0,
                                    vLen * vLen * sizeof(cl_float), h_est);
            dt = 0.01 * (*std::min_element(h_est, h_est + vLen * vLen));
            dt = std::min({(double) dt, timeStep});


            // Check dt for not to overshoot write time
            if ((writeTime + dt) > writeInterval)
                dt = writeInterval - writeTime;
            std::cout << nt << " " << dt <<" " << time <<std::endl;
        }
        //std::cout << stateP(15) << std::endl;

        //End time
        cl_ulong end=
                event.getProfilingInfo<CL_PROFILING_COMMAND_END>();
        double time = 1.e-9 * (end-start);
        std::cout << "Time for kernel to execute: " << time << std::endl;


        queue.enqueueReadBuffer(d_pos, CL_TRUE, 0, vSize, h_pos);
        queue.enqueueReadBuffer(d_vel, CL_TRUE, 0, vSize, h_vel);

        for (size_t i = 0; i < vLen; ++i) {
            for (unsigned int j = 0; j < particles.at(i).r.size(); ++j) {
                 particles.at(i).r(j) = h_pos[i].s[j];
                 particles.at(i).v(j) = h_vel[i].s[j];
            }
            //particles.at(i).m = h_pos[i].s[3];
        }
        std::cout << stateP(15) << std::endl;
        stateToFile(writePr, outFile);

    } catch (cl::Error error) {
        std::cerr << "Caught exception: " << error.what()
                  << '(' << error.err() << ')'
                  << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

/*void PSystem::checkBoundsCyclic()
{
    for (auto& particle : particles) {
        for (int i = 0; i < 3; ++i) {
            if (particle.r(i) > dimXYZ(i)) {particle.r(i) = fmod(particle.r(i), dimXYZ(i));}
            if (particle.r(i) < 0) {particle.r(i) = dimXYZ(i) + fmod(particle.r(i), dimXYZ(i));}
        }
    }
}

void PSystem::checkBoundsHard() {
    for (auto& particle : particles) {
        for (int i = 0; i < 3; ++i ){
            if (particle.r(i) > dimXYZ(i)) {
                particle.r(i) += 2*(dimXYZ(i) - particle.r(i));
                particle.v(i) = -particle.v(i);
            }
            if (particle.r(i) < 0) {
                particle.r(i) = -particle.r(i);
                particle.v(i) = -particle.v(i);
            }
        }
    }
}
*/


std::string PSystem::parameters()
{
    //int precision = 5;
    std::ostringstream out;
    //out << std::scientific;
    //out << std::setprecision(precision);

    out << "timeStep = " << timeStep << std::endl;
    out << "endTime = " << endTime << std::endl;
    out << "eps = " << eps << std::endl;
    out << "G = " << G << std::endl;
    out << "potentialName = " << potentialName << std::endl;
    out << "outFile = " << outFile << std::endl;
    out << "writePr = " << writePr << std::endl;
    out << "boundingBox = [ (" << bBox.minVertex(0)
                               << " "
                               << bBox.minVertex(1)
                               << " "
                               << bBox.minVertex(2)
                               << ") ("
                               << bBox.maxVertex(0)
                               << " "
                               << bBox.maxVertex(1)
                               << " "
                               << bBox.maxVertex(2)
                               << ") ]"
                               << std::endl;
    return out.str();
}

void PSystem::info(int p)
{
    std::streamsize defaultPrecision = std::cout.precision();

    std::cout << std::setprecision(p);
    //std::cout << std::fixed;
    std::cout << std::scientific;
    std::cout << "Time = " << time;
    std::cout << "\tTime step = " << timeStep << std::endl;
    std::cout << std::scientific;
    std::cout << stateE(p) << std::endl;

    std::cout.unsetf(std::ios::fixed | std::ios::scientific);
    std::cout << std::setprecision(defaultPrecision);
}

std::string PSystem::stateE(int p)
{
    std::string delim = " ";
    std::ostringstream out;
    out << std::scientific;
    out << std::setprecision(p);

    out << "E_kin=";
    out << getEkin();
    out << delim << "E_pot=";
    out << getEpot();
    out << delim << "E_tot=";
    out << getEtot();
    out << delim;
    if (EtotInit) {
        out << "(E_tot-E_init)/E_init=";
        out << (getEtot() - EtotInit)/EtotInit;
    } else {
        out << "E_init=0\t(E_tot-E_init)=";
        out << (getEtot() - EtotInit);
    }

    return out.str();
}

std::string PSystem::stateP(int p)
{
    std::ostringstream out;

    for (auto& particle : particles){
        out << particle.state(p);
        out << std::endl;
    }

    return out.str();
}

void PSystem::stateToFile(int p, std::string fileName)
{
    std::ofstream osFile;

    std::string cSymb = "#";

    osFile.open(fileName, std::ofstream::out | std::ofstream::app); // Output, append
    if (osFile.is_open()){
        osFile << cSymb << "Time = " << time; // Time
        osFile << "\t Time step = " << timeStep << std::endl; 
        osFile << cSymb << stateE(p) << std::endl; // Energy
        osFile << stateP(p); // All particles state
        osFile << std::endl << std::endl; // End of the data block
        osFile.close();
    } else {
        std::cout << "Error opening output file!" << std::endl;
    }

}

void PSystem::readParticlesData(std::string fileName)
{
    std::ifstream isFile;
    std::stringstream sLine;
    std::string line;
    std::string word;
    std::vector<std::string> words;
    std::string blockMark = "#Time";
    std::istream::pos_type blockPos;

    // Particle parameters to read
    //int id;
    double mass;
    Eigen::Vector3d rad, vel;

    isFile.open(fileName);

    if (isFile.is_open()){
        // Search the beginning of the last time step data block
        while (std::getline(isFile, line)) {
            sLine = std::stringstream(line);
            sLine >> word;
            if (word == blockMark) {blockPos = isFile.tellg();}
        }

        isFile.clear();
        isFile.seekg(blockPos); // Setting position to the beginning of last time block
        std::getline(isFile, line); // Skipping lines with system's energy data

        // Reading particles positions until the end of file
        while (std::getline(isFile, line)) {
            if (line == "") {break;}
            sLine = std::stringstream(line);
            sLine >> word;
            //id = std::stoi(word);

            for (int i = 0; i < 3; i++) {
                sLine >> word;
                rad(i) = std::stod(word);
            }

            for (int i = 0; i < 3; i++) {
                sLine >> word;
                vel(i) = std::stod(word);
            }

            sLine >> word;
            mass = std::stod(word);

            particles.emplace_back(mass, rad, vel);
        }

    } else {
        std::cout << "Particles data file " << fileName << " is not found." << std::endl;
    }
}

void PSystem::readBoundariesData(std::string fileName)
{
    std::ifstream isFile;
    std::stringstream sLine;
    std::string line;
    std::string word;
    std::string blockMark = "#Boundary"; // This string should preceed the sequence of boundary points
    std::istream::pos_type blockPos; // Position of boundary block
    
    
    // Boundary parameters to read
    Eigen::Vector3d vertex;

    isFile.open(fileName);

    if (isFile.is_open()){
        // Search the beginning of the first boundary data block
        while (std::getline(isFile, line)) {
            sLine = std::stringstream(line);
            sLine >> word;
            if (word == blockMark) {
                blockPos = isFile.tellg();
                (*this).boundaries.emplace_back();
                break;
            }
        }
       
        // Reading boundaries vertices until the end of file
        while (std::getline(isFile, line)) {
            if (line != "") {
                sLine = std::stringstream(line);
                sLine >> word;
                if (word == blockMark) {
                    (*this).boundaries.emplace_back();
                } else {
                    vertex(0) = std::stod(word);
                    sLine >> word; vertex(1) = std::stod(word);
                    sLine >> word; vertex(2) = std::stod(word);
                    // Adds vertex to the last element (convex) of boundaries vector.
                    (*this).boundaries.back().vertices.emplace_back(vertex);
                }
            }
        }
        bBox.setBoundingBox(boundaries);
    } else {
        std::cout << "Boundaries data file " << fileName << " is not found." << std::endl;
    }
}

/*
void PSystem::readParams(std::string fileName)
{
    //double dimX = 0, dimY = 0, dimZ = 0; // System's dimensions

    std::ifstream isFile;
    std::stringstream ss;
    std::string parametersLine, parameter;
    std::vector<std::string> parameters;

    isFile.open(fileName);

    if (isFile.is_open()) {
        while (std::getline(isFile, parametersLine)) {
            ss = std::stringstream(parametersLine);
            while (ss >> parameter) {
                parameters.push_back(parameter);
            }

            if (parameters.at(0) == "random_particles") {
                for (int i = 0; i < std::stoi(parameters.at(1)); i++){
                    addParticle();
                    setRandom(particles.back());
                }
            }
            if (parameters.at(0) == "time_step") {
                timeStep = std::stod(parameters.at(1));
            }
            if (parameters.at(0) == "end_time") {
                endTime= std::stod(parameters.at(1));
            }
            if (parameters.at(0) == "eps") {
                eps = std::stod(parameters.at(1));
            }
            if (parameters.at(0) == "G") {
                G = std::stod(parameters.at(1));
            }
            if (parameters.at(0) == "integration_method") {
                methodName = parameters.at(1);
            }
            if (parameters.at(0) == "interaction_potential") {
                potentialName = parameters.at(1);
            }
            if (parameters.at(0) == "write_to_file") {
                outFile = parameters.at(1);
            }
            if (parameters.at(0) == "particles_input_file") {
                readParticlesData(parameters.at(1)); 
            }
            if (parameters.at(0) == "write_precision") {
                writePr = std::stoi(parameters.at(1));
            }
            if (parameters.at(0) == "write_interval") {
                writeInterval = std::stod(parameters.at(1));
            }
            if (parameters.at(0) == "boundaries_input_file") {
                readBoundariesData(parameters.at(1));
            }
            parameters.erase(parameters.begin(), parameters.end());
        }
    } else {
        std::cout << "Configuration file " << fileName << " is not found." << std::endl;
    }
}
*/

int PSystem::readConfig(const char *fileName)
{
    libconfig::Config cfg;

    // Read the file. If there is an error, report it and exit.
    try {
        cfg.readFile(fileName);
    }
    catch(const libconfig::FileIOException &fioex) {
        std::cerr << "I/O error while reading file." << std::endl;
        throw;
    }
    catch(const libconfig::ParseException &pex) {
        std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
                  << " - " << pex.getError() << std::endl;
        throw;
    }
    try{
        try {
            timeStep = cfg.lookup("time_step");
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'time_step' setting in configuration file."
                      << std::endl;
            throw;
        }

        try {
            endTime = cfg.lookup("end_time");
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'end_time' setting in configuration file."
                      << std::endl;
            throw;
        }

        try {
            std::string integration_method = cfg.lookup("integration_method");
            methodName = integration_method;
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'integration_method' setting in configuration file."
                      << std::endl;
            throw;
        }

        try {
            std::string interaction_potential = cfg.lookup("interaction_potential");
            potentialName = interaction_potential;
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'interaction_potential' setting in configuration file."
                      << std::endl;
            throw;
        }

        try {
            std::string output_file = cfg.lookup("output_file");
            outFile = output_file;
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'output_file' setting in configuration file."
                      << std::endl;
            throw;
        }

        try {
            std::string particles_data_file = cfg.lookup("particles_data_file");
            readParticlesData(particles_data_file);
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'particles_data_file' setting in configuration file."
                      << std::endl;
            throw;
        }

        try {
            writePr = cfg.lookup("write_precision");
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'write_precision' setting in configuration file."
                      << std::endl;
            throw;
        }

        try {
            writeInterval =  cfg.lookup("write_interval");
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'write_interval' setting in configuration file."
                      << std::endl;
            throw;
        }

        try {
            std::string boundaries_data_file = cfg.lookup("boundaries_data_file");
            readBoundariesData(boundaries_data_file);
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'write_interval' setting in configuration file."
                      << std::endl;
            throw;
        }

        try {
            int random_particles = cfg.lookup("random_particles");
            for (int i = 0; i < random_particles; i++){
                addParticle();
                setRandom(particles.back());
            }
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            //std::cerr << "No 'write_interval' setting in configuration file." << std::endl;
        }

        //
    }
    catch (const libconfig::SettingTypeException &stex) {
        std::cout << "Something wrong in the settings file with '"
                  << stex.getPath() << "'." << std::endl;
        throw;
        //return(EXIT_FAILURE);
    }

    return(EXIT_SUCCESS);
}

void PSystem::checkBoundaries()
{
    /* There may be situations when out of boundary displacement is so large
     so mirroring position of particle relative to boundary
     places particle out of another boundary. At least two passes of checks are 
     needed. */
    bool xFlag; // Boundary crossing flag
    std::vector<Eigen::Vector3d> xPointN; // Path and boundary cross point with normal vector at this point
    Eigen::Vector3d xPoint; // Cross point
    Eigen::Vector3d n0; // Normal vector to boundary at cross point
    LineSegment path;
    Eigen::Vector3d dr; //Out of boundary displacement of particle
    
    //std::ofstream osFile;
    //osFile.open("collisions.dat", std::ofstream::out | std::ofstream::app); // Output, append


    do {
        xFlag = false;
        for (auto& particle : particles) {
            for (auto& boundary : boundaries) {
                path = LineSegment(particle.r0, particle.r);
                xPointN = intersectionPointN(path, boundary);
                xPoint = xPointN.at(0);
                n0 = xPointN.at(1);
                if (!std::isnan(xPoint(0))) {

                    /*
                    osFile << "NEW COLLISION!" << " " << time << " " << particle.id << std::endl;
                    osFile << boundary.vertices.at(0) <<  "\n "<< std::endl;
                    osFile << xPoint << "\n n0 \n " << n0 << "\n r \n " << particle.r << "\n r0 \n " << particle.r0 << "\n v\n " << particle.v << std::endl;
                    */
                    // what is r0? should it be changed for next iteration of crossing check?
                    particle.v = particle.v - 2 * particle.v.dot(n0) * n0;
                    dr = particle.r - xPoint;
                    dr = dr - 2 * dr.dot(n0) * n0;
                    particle.r = xPoint + dr;
                    //osFile << xPoint << "\n n0" << n0 << "\n r\n " << particle.r <<"\n r0\n " << particle.r0 << "\n v\n " << particle.v << std::endl;
                    xFlag = true;
                }
            }      
        }
    } while (xFlag);
    
}

boundingBox::boundingBox()
{
    minVertex << -DBL_MAX, - DBL_MAX, - DBL_MAX;
    maxVertex << DBL_MAX, DBL_MAX, DBL_MAX;
}

void boundingBox::setBoundingBox(Eigen::Vector3d min, Eigen::Vector3d max)
{
    (*this).minVertex = min;
    (*this).maxVertex = max;
}

void boundingBox::setBoundingBox(std::vector<Convex> boundaries)
{
    Eigen::Vector3d minVertex(0, 0, 0);
    Eigen::Vector3d maxVertex(0, 0, 0);

    for (auto& boundary : boundaries) {
        for (auto& vertex : boundary.vertices) {
            for (int i = 0; i < 3; ++i) {
                if (minVertex(i) > vertex(i)) minVertex(i) = vertex(i);
                if (maxVertex(i) < vertex(i)) maxVertex(i) = vertex(i);
            }
        }
    }

    setBoundingBox(minVertex, maxVertex);
}

void PSystem::checkBoundingBox() {
    for (auto& particle : particles) {
        for (int i = 0; i < 3; ++i ){
            if (particle.r(i) > (*this).bBox.maxVertex(i)) {
                particle.r(i) += 2*((*this).bBox.maxVertex(i) - particle.r(i));
                particle.v(i) = -particle.v(i);
            }
            if (particle.r(i) < (*this).bBox.minVertex(i)) {
                particle.r(i) += 2*((*this).bBox.minVertex(i) - particle.r(i));
                particle.v(i) = -particle.v(i);
            }
        }
    }
}

Eigen::Vector3d PSystem::extFieldUniform(Particle& prtcl)
{
    Eigen::Vector3d F;
    F = prtcl.m * g;
    return F;
}
