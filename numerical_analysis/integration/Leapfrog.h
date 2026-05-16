#ifndef LEAPFROG_H
#define LEAPFROG_H

#include <vector>
#include <string>
#include <memory>
#include "Environment.h"
#include "Rocket.h"

/**
 * @brief Represents the state of the vehicle at a specific time step.
 */
struct State {
    double t_s;                 
    double altitude_m;          
    double velocity_m_s;        
    double acceleration_m_s2;   
    double atm_pressure_pa;     
    double mass_kg;             
    double thrust_n;            
    double drag_n;              
    double dynamic_pressure_pa; 
    double mach_number;         
};

class LeapfrogIntegrator {
public:
    /**
     * @brief Constructor for the Leapfrog Integrator using Dependency Injection.
     * @param env Pointer to the planetary environment interface.
     * @param rocket Pointer to the rocket interface.
     * @param dt_s The time step in seconds.
     */
    LeapfrogIntegrator(IEnvironment* env, IRocket* rocket, double dt_s);

    std::vector<State> simulate(double t_max_s, double initial_altitude_m, double initial_velocity_m_s);
    static void plotTrajectory(const std::vector<State>& states, const std::string& filename);

private:
    IEnvironment* env;
    IRocket* rocket;
    double dt_s;
};

#endif
