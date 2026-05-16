#include "Rocket.h"
#include <cmath>
#include <algorithm>

LVM3Rocket::LVM3Rocket() : current_mass_kg(640000.0) {}

void LVM3Rocket::update_dt(double dt_s, double altitude_m, double velocity_m_s) {
    current_mass_flow_rate_kg_s = 0.0;
    current_exit_velocity_m_s = 0.0;
    current_exit_pa = 101325.0;
    current_nozzle_exit_area_m2 = 0.0;

    const double s200_dry_mass_kg = 62000.0;
    const double l110_dry_mass_kg = 15000.0;

    // Phase 1: S200 Solid Strap-ons burning
    if (prop_s200_kg > 0.0) {
        current_mass_flow_rate_kg_s = 5460.0;
        current_exit_velocity_m_s = 2690.0;
        current_exit_pa = 60000.0;
        current_nozzle_exit_area_m2 = 11.2;
        
        double burn_kg = current_mass_flow_rate_kg_s * dt_s;
        if (burn_kg > prop_s200_kg) burn_kg = prop_s200_kg;
        prop_s200_kg -= burn_kg;
    }
    // Phase 2: S200 expended, L110 core burns
    else if (prop_l110_kg > 0.0) {
        if (!s200_jettisoned) {
            current_mass_kg -= s200_dry_mass_kg;
            s200_jettisoned = true;
        }
        current_mass_flow_rate_kg_s = 571.4;
        current_exit_velocity_m_s = 2873.0;
        current_exit_pa = 52000.0;
        current_nozzle_exit_area_m2 = 1.662;
        
        double burn_kg = current_mass_flow_rate_kg_s * dt_s;
        if (burn_kg > prop_l110_kg) burn_kg = prop_l110_kg;
        prop_l110_kg -= burn_kg;
    }
    // Phase 3: L110 expended, C25 upper stage burns
    else if (prop_c25_kg > 0.0) {
        if (!l110_jettisoned) {
            current_mass_kg -= l110_dry_mass_kg;
            l110_jettisoned = true;
        }
        current_mass_flow_rate_kg_s = 44.5;
        current_exit_velocity_m_s = 4345.0;
        current_exit_pa = 5000.0;
        current_nozzle_exit_area_m2 = 2.45;
        
        double burn_kg = current_mass_flow_rate_kg_s * dt_s;
        if (burn_kg > prop_c25_kg) burn_kg = prop_c25_kg;
        prop_c25_kg -= burn_kg;
    }

    // Payload Fairing (PLF) Jettison
    const double plf_mass_kg = 1500.0;
    const double plf_jettison_alt_m = 115000.0;

    if (!plf_jettisoned && altitude_m > plf_jettison_alt_m) {
        current_mass_kg -= plf_mass_kg;
        plf_jettisoned = true;
    }

    // Apply fuel consumption to total mass and safeguard against negative mass
    current_mass_kg -= current_mass_flow_rate_kg_s * dt_s;
    if (current_mass_kg < 5000.0) current_mass_kg = 5000.0; // Payload floor
}

double LVM3Rocket::get_thrust_n(double ambient_pressure_pa) const {
    return (current_mass_flow_rate_kg_s * current_exit_velocity_m_s) + 
           (current_exit_pa - ambient_pressure_pa) * current_nozzle_exit_area_m2;
}

double LVM3Rocket::get_mach(double speed_of_sound_m_s, double velocity_m_s) const {
    if (speed_of_sound_m_s <= 0) return 0.0;
    return std::abs(velocity_m_s) / speed_of_sound_m_s;
}

double LVM3Rocket::get_drag_force_n(double density_kg_m3, double speed_of_sound_m_s, double velocity_m_s) const {
    if (density_kg_m3 <= 0.0 || speed_of_sound_m_s <= 0.0) return 0.0;

    double mach = get_mach(speed_of_sound_m_s, velocity_m_s);
    double cd = 0.25;

    if (mach < 0.8) {
        cd = 0.20 + 0.02 * (mach / 0.8);
    } else if (mach < 1.2) {
        double t = (mach - 0.8) / 0.4;
        cd = 0.22 + 0.38 * t;
    } else if (mach < 5.0) {
        cd = std::max(0.25, 0.60 * (1.2 * 1.2) / (mach * mach));
    }

    double frontal_area_m2 = (std::abs(velocity_m_s) >= 1700.0) ? 12.6 : 25.0;
    return 0.5 * density_kg_m3 * velocity_m_s * velocity_m_s * cd * frontal_area_m2;
}
