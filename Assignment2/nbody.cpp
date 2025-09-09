#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <ctime>
#include <chrono>
#include <vector>
/*
//<TODO>
Create timestamps for how long each test takes
have ~400 steps for the solar system test to verify orbits are behaving correctly
Fix python plot script so it actually outputs the correct fucking thing
README with description




*/




//define a variable type to house three coordinates
struct Vec3{
    double x;
    double y;
    double z;

    //define operators to make writing vector add/sub/mul/div easier
    Vec3 operator+(const Vec3& other) const {
        Vec3 result;
        result.x = this->x + other.x;
        result.y = this->y + other.y;
        result.z = this->z + other.z;
        return result;
    }

    Vec3 operator-(const Vec3& other) const {
        Vec3 result;
        result.x = this->x - other.x;
        result.y = this->y - other.y;
        result.z = this->z - other.z;
        return result;
    }

    Vec3 operator*(double scalar) const{
        Vec3 result;
        result.x = this->x * scalar;
        result.y = this->y * scalar;
        result.z = this->z * scalar;
        return result;
    }
    friend Vec3 operator*(double scalar, const Vec3& v) {
        return v * scalar; 
    }

    Vec3 operator/(double scalar) const {
    Vec3 result;
    result.x = this->x / scalar;
    result.y = this->y / scalar;
    result.z = this->z / scalar;
    return result;
    }
};


//Particle class definition
// with getters/setters
class Particle{
private:
    //std::string name;
    int ID;
    double mass;
    Vec3 location;
    Vec3 velocity;
    Vec3 acceleration;
    Vec3 force;

public:
    Particle(int pID, double pMass, Vec3 pLocation, Vec3 pVelocity){
        ID = pID;
        mass = pMass;
        location = pLocation;
        velocity = pVelocity;
    }
    int getPID(){
        return ID;
    }
    void setPName(int pID){
        ID = pID;
    }
    double getPMass(){
        return mass;
    }
    void setPMass(double pMass){
        mass = pMass;
    }
    Vec3 getPLocation(){
        return location;
    }
    void setPLocation(Vec3 pLocation){
        location = pLocation;
    }
    Vec3 getPVelocity(){
        return velocity;
    }
    void setPVelocity(Vec3 pVelocity){
        velocity = pVelocity;
    }
    Vec3 getPAcceleration(){
        return acceleration;
    }
    void setPAcceleration(Vec3 pAcceleration){
        acceleration = pAcceleration;
    }
    Vec3 getPForce(){
        return force;
    }
    void setPForce(Vec3 pForce){
        force = pForce;
    }
    
};
//function protos
size_t getInput(const std::string &prompt);
double initialize();
Vec3 vectorInitialize();
Vec3 getDisplacement(Vec3 a, Vec3 b);
double getMagnitude(Vec3 a);
Vec3 getDirection(Vec3 a, double magnitude);
double getForceApplied(double gravity, double mass1, double mass2, double magnitude, double softening);
Vec3 getForce(double forceApplied, Vec3 direction);
void initializeSolarSystem(std::vector<Particle> &Particles);
int main (){

    //open file to write to
    //write header
    std::ofstream outFile("simulation_output.tsv");
    if (!outFile) {
        std::cerr << "Error: could not open output file.\n";
        return 1;
    }
    //outFile << "t\tID\tmass\tx\ty\tz\tvx\tvy\tvz\tfx\tfy\tfz\n";
    //get user input
    size_t N = getInput("How many particles in the experiment?");
    size_t timeSteps = getInput("How many time steps in the experiment?");
    size_t deltaTime = getInput("Timestep size?");
    size_t M = getInput("How often should the states be reset?");

    //initialize vectors/ const variables
    std::vector<Particle> Particles;
    Particles.reserve(N);
    std::vector<int> timeStep(N);
    const double gravity = 6.67430e-11;
    const double softening = 1e9;

    //initialize Particles w/
    // ID, mass, location, and velocity
    for (size_t i = 0; i < N; ++i) {
        Particles.push_back(Particle(i, initialize(), vectorInitialize(), vectorInitialize()));
    }
    //solar system test

    /*initializeSolarSystem(Particles);
    N = Particles.size();  // now N = 10*/

    
    //main loop
    //set Force to zero at each time step
    //start timer
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t t=0; t < timeSteps; ++t){
        for (auto &p: Particles){
            p.setPForce({0, 0, 0});
        }
        
        //compute all forces required for formula
        //i.e displacement b/tw particles, 
        for (size_t i=0; i < Particles.size(); ++i){
            for (size_t j=0; j < Particles.size(); ++j){
                if (i == j) continue;

                Vec3 displacement = getDisplacement(Particles[j].getPLocation(), Particles[i].getPLocation());
                double magnitude = getMagnitude(displacement);
                Vec3 direction = getDirection(displacement, magnitude);
                double forceApplied = getForceApplied(gravity, Particles[j].getPMass(), Particles[i].getPMass(), magnitude, softening);
                Vec3 forceVector = forceApplied * direction;
                Vec3 newForce = Particles[i].getPForce();  
                newForce = newForce + forceVector;                    
                Particles[i].setPForce(newForce);            
            }
        }
        //update acceleration
         for (auto &p : Particles) {
            Vec3 acceleration = p.getPForce() / p.getPMass();
            p.setPAcceleration(acceleration);
        }

        // update velocity
        for (auto &p : Particles) {
            Vec3 newVelocity = p.getPVelocity() + p.getPAcceleration() * deltaTime;
            p.setPVelocity(newVelocity);
        }

        // update location
        for (auto &p : Particles) {
            Vec3 newLocation = p.getPLocation() + p.getPVelocity() * deltaTime;
            p.setPLocation(newLocation);
        }

        // output state every M steps
        if (t % M == 0) {
            outFile << Particles.size();
for (auto &p : Particles) {
        outFile << "\t"
            //<< t << "\t" 
            //<< p.getPID() << "\t" 
            << p.getPMass() << "\t"
            << p.getPLocation().x << "\t"
            << p.getPLocation().y << "\t"
            << p.getPLocation().z << "\t"
            << p.getPVelocity().x << "\t"
            << p.getPVelocity().y << "\t"
            << p.getPVelocity().z << "\t"
            << p.getPForce().x << "\t"
            << p.getPForce().y << "\t"
            << p.getPForce().z;
                }
                outFile << "\n";
            }
        }
    //end timer
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Simulation completed in " << elapsed.count() << " seconds.\n";
    return 0;
}

//helper function for user input
size_t getInput(const std::string &prompt) {
    size_t value;
    std::cout << prompt << " ";
    std::cin >> value;
    return value;
}

//random number gen for initializing mass
double initialize(){
    static std::random_device rd;            
    static std::mt19937 gen(rd());           
    static std::uniform_real_distribution<double> dist(1.0, 100.0);
    return dist(gen);
}

//to initialize our Vec3s
Vec3 vectorInitialize(){
    static std::random_device rd;            
    static std::mt19937 gen(rd());           
    static std::uniform_real_distribution<double> dist(1.0, 100.0);
    return{dist(gen), dist(gen), dist(gen)};
}

//displacement = distance between particles
Vec3 getDisplacement(Vec3 a, Vec3 b){
    Vec3 displacement = a - b;
    return displacement;
}

//magnitude gets our scalar (R^2 in formula)
double getMagnitude(Vec3 a){
    double magnitude = std::sqrt((a.x*a.x)+(a.y*a.y)+(a.z*a.z));
    return magnitude;
}

//gets our direction by dividing by scalar
Vec3 getDirection(Vec3 a, double magnitude){
    Vec3 direction;
    direction.x = a.x/magnitude;
    direction.y = a.y/magnitude;
    direction.z= a.z/magnitude;
    return direction;
}

//returns Force Applied which is formula 1
double getForceApplied(double gravity, double mass1, double mass2, double magnitude, double softening){
    double forceApplied = (gravity * mass1 * mass2) / ((magnitude*magnitude) + softening);
    return forceApplied;
}

//returns Force for final calculations
Vec3 getForce(double forceApplied, Vec3 direction){
    Vec3 force;
    force.x = direction.x * forceApplied;
    force.y = direction.y * forceApplied;
    force.z = direction.z * forceApplied;
    return force;
}

//solar system initialization for testing purposes
//sun -> pluto in order
void initializeSolarSystem(std::vector<Particle> &Particles) {
    Particles.clear();
    Particles.reserve(10);
    Particles.push_back(Particle(0, 1.9891e30, {0,0,0}, {0,0,0}));  
    Particles.push_back(Particle(1, 3.285e23, {5.8344e10,0,0}, {0,47870,0})); 
    Particles.push_back(Particle(2, 4.867e24, {1.07712e11,0,0}, {0,35020,0}));   
    Particles.push_back(Particle(3, 5.972e24, {1.496e11,0,0}, {0,29780,0}));  
    Particles.push_back(Particle(4, 6.39e23, {2.27392e11,0,0}, {0,24130,0}));   
    Particles.push_back(Particle(5, 1.898e27, {7.7792e11,0,0}, {0,13070,0}));   
    Particles.push_back(Particle(6, 5.683e26, {1.43317e12,0,0}, {0,9680,0}));   
    Particles.push_back(Particle(7, 8.681e25, {2.87531e12,0,0}, {0,6800,0}));  
    Particles.push_back(Particle(8, 1.024e26, {4.49548e12,0,0}, {0,5430,0}));
    Particles.push_back(Particle(9, 7.342e22, {1.49984e11,0,0}, {0,30802,0}));
}
