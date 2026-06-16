#ifndef DIPOLE_SIMULATION_H
#define DIPOLE_SIMULATION_H

#include <vector>
#include <string>

namespace DipoleSimulation {

    struct State {
        double t_s;
        double pos_x;
        double pos_y;
        double pos_z;
        double vel_x;
        double vel_y;
        double vel_z;
        double speed;
        double energy;
        double drift_pct;
        double dist_pos;
        double dist_neg;
    };

    /**
     * @brief Run the entire 2D/3D CIC deposition, 3D Poisson solver, Leapfrog trajectory simulation, and diagnostics.
     */
    void runSimulation();

} // namespace DipoleSimulation

#endif // DIPOLE_SIMULATION_H
