#ifndef LEAPFROG_H
#define LEAPFROG_H

#include <vector>
#include <string>

struct State {
    double t;
    double x;
    double v;
    double a;
    double atm_pressure;
};

class LeapfrogIntegrator {
public:
    // Initial configuration for ISRO LVM3
    LeapfrogIntegrator(double dt, double initial_mass);
    std::vector<State> simulate(double t_max, double initial_x, double initial_v);
    
    // Saves a plot using ROOT TGraph
    static void plotTrajectory(const std::vector<State>& states, const std::string& filename);

private:
    double dt;
    double current_mass;
    
    // Ambient pressure calculation based on U.S. Standard Atmosphere (1976)
    double get_ambient_pressure(double geometric_altitude_m, double earth_radius) const;
    
    // LVM3 stage properties logic
    void get_lvm3_stage_properties(double altitude_m, double& mass_flow_rate, double& exit_velocity, double& exit_pa, double& nozzle_exit_area);
    
    // Tracks if stages have been jettisoned
    bool s200_jettisoned = false;
    bool l110_jettisoned = false;
};

#endif
