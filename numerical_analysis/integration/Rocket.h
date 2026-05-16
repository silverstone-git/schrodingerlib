#ifndef ROCKET_H
#define ROCKET_H

#include <string>

/**
 * @brief Interface for any rocket or vehicle in the simulation.
 */
class IRocket {
public:
    virtual ~IRocket() = default;
    
    virtual std::string get_name() const = 0;
    
    // Core physical properties
    virtual double get_mass_kg() const = 0;
    
    // Updates internal rocket state (e.g., fuel consumption, staging)
    virtual void update_dt(double dt_s, double altitude_m, double velocity_m_s) = 0;
    
    // Forces acting on the rocket
    virtual double get_thrust_n(double ambient_pressure_pa) const = 0;
    virtual double get_drag_force_n(double density_kg_m3, double speed_of_sound_m_s, double velocity_m_s) const = 0;
    
    // Aerodynamic helper
    virtual double get_mach(double speed_of_sound_m_s, double velocity_m_s) const = 0;
};

/**
 * @brief Concrete implementation of ISRO's LVM3 rocket.
 */
class LVM3Rocket : public IRocket {
public:
    LVM3Rocket();
    
    std::string get_name() const override { return "ISRO LVM3"; }
    double get_mass_kg() const override { return current_mass_kg; }
    void update_dt(double dt_s, double altitude_m, double velocity_m_s) override;
    double get_thrust_n(double ambient_pressure_pa) const override;
    double get_drag_force_n(double density_kg_m3, double speed_of_sound_m_s, double velocity_m_s) const override;
    double get_mach(double speed_of_sound_m_s, double velocity_m_s) const override;

private:
    double current_mass_kg;
    
    // Stage jettison flags
    bool s200_jettisoned = false;
    bool l110_jettisoned = false;
    bool plf_jettisoned = false;

    // Propellant trackers in kg
    double prop_s200_kg = 410000.0; 
    double prop_l110_kg = 116000.0;
    double prop_c25_kg  = 28000.0;
    
    // Cached engine properties from last update_dt
    double current_mass_flow_rate_kg_s = 0.0;
    double current_exit_velocity_m_s = 0.0;
    double current_exit_pa = 101325.0;
    double current_nozzle_exit_area_m2 = 0.0;
};

#endif
