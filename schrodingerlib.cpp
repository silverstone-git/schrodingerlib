#include <TGraph.h>
#include <TCanvas.h>
#include <cmath>
#include <array>
#include <algorithm>
#include <iostream>
#include "schrodingerlib.hpp"

/**
 * @brief Computes ambient atmospheric pressure based on the U.S. Standard Atmosphere (1976).
 * @param geometric_altitude_m The altitude of the rocket/craft in meters above sea level.
 * @return The ambient pressure in Pascals (Pa). Returns 0.0 if out of atmospheric bounds (>86km).
 */
double get_ambient_pressure(double geometric_altitude_m, double EARTH_RADIUS) {
    // 1. Convert geometric altitude to geopotential altitude (accounting for gravity drop)
    double h = (EARTH_RADIUS * geometric_altitude_m) / (EARTH_RADIUS + geometric_altitude_m);

    if (h < 0.0) return 101325.0; // Cap at Sea Level pressure if underwater/below datum
    if (h > 86000.0) return 0.0;  // Beyond the 1976 model boundary (effectively vacuum)

    // 2. Define structural data layers for the U.S. Standard Atmosphere 1976
    struct AtmosphericLayer {
        double base_h;  // Geopotential base altitude (m)
        double base_P;  // Base pressure of layer (Pa)
        double base_T;  // Base temperature of layer (K)
        double lapse_L; // Temperature lapse rate (K/m)
    };

    // Pre-calculated exact historical base pressures to prevent structural drift error
    static const std::array<AtmosphericLayer, 7> layers = {{
        {0.0,     101325.0,   288.15, -0.0065}, // Troposphere
        {11000.0, 22632.1,    216.65,  0.0000}, // Tropopause / Lower Stratosphere
        {20000.0, 5474.89,    216.65,  0.0010}, // Stratosphere
        {32000.0, 868.019,    228.65,  0.0028}, // Stratosphere Upper
        {47000.0, 110.906,    270.65,  0.0000}, // Stratopause
        {51000.0, 66.9389,    270.65, -0.0028}, // Mesosphere
        {71000.0, 3.95642,    214.65, -0.0020}  // Mesosphere Upper
    }};

    // 3. Find the current operational layer using binary search (std::upper_bound)
    auto it = std::upper_bound(layers.begin(), layers.end(), h, 
        [](double val, const AtmosphericLayer& layer) { return val < layer.base_h; });
    
    // Step back to grab the active layer lower bound
    const auto& layer = *(std::prev(it));

    // 4. Physical Constants
    constexpr double g0 = 9.80665;      // Standard gravity (m/s^2)
    constexpr double M  = 0.0289644;    // Molar mass of Earth's air (kg/mol)
    constexpr double R  = 8.31432;      // Universal gas constant (J/mol*K)

    // 5. Compute Pressure based on local thermodynamic profile
    double delta_h = h - layer.base_h;

    if (std::abs(layer.lapse_L) < 1e-9) {
        // Isothermal Layer (Lapse rate is 0) -> Governed purely by the Boltzmann/Hydrostatic variant
        return layer.base_P * std::exp((-g0 * M * delta_h) / (R * layer.base_T));
    } else {
        // Lapsing Layer (Temperature changes linearly) -> Governed by Power-Law variant
        double T_local = layer.base_T + layer.lapse_L * delta_h;
        double exponent = (-g0 * M) / (R * layer.lapse_L);
        return layer.base_P * std::pow(T_local / layer.base_T, exponent);
    }
}




/**
 * @brief Updates rocket mass and returns active engine parameters based on altitude and fuel.
 */
void get_lvm3_stage_properties(double altitude_m, double dt, double& current_mass,
                               double& mass_flow_rate, double& exit_velocity,
                               double& exit_pa, double& nozzle_exit_area, bool& s200_jettisoned, bool& l110_jettisoned,
                               bool& plf_jettisoned, double& prop_s200, double& prop_l110, double& prop_c25) {

    const double s200_dry_mass_kg = 62000;
    const double l110_dry_mass_kg = 15000;

    // Phase 1: S200 Solid Strap-ons burning
    if (prop_s200 > 0.0) {
        mass_flow_rate   = 5460.0;
        exit_velocity    = 2690.0;
        exit_pa          = 60000.0;
        nozzle_exit_area = 11.2;
        double burn = mass_flow_rate * dt;
        prop_s200     = std::max(0.0, prop_s200 - burn);
        current_mass -= burn;
    }
    // Phase 2: S200 expended, L110 core burns
    else if (prop_l110 > 0.0) {
        // Jettison S200 casings once, the first time we enter this phase
        if (!s200_jettisoned) {
            current_mass   -= s200_dry_mass_kg;
            s200_jettisoned = true;
            std::cout << "s200 jettisoned, " << current_mass / 1000.0 << "tonnes, " << altitude_m / 1000.0 << "kms, " << "\n";
        }
        mass_flow_rate   = 571.4;
        exit_velocity    = 2873.0;
        exit_pa          = 52000.0;
        nozzle_exit_area = 1.662;
        double burn = mass_flow_rate * dt;
        prop_l110     = std::max(0.0, prop_l110 - burn);
        current_mass -= burn;
    }
    // Phase 3: L110 expended, C25 upper stage burns
    else if (prop_c25 > 0.0) {
        if (!l110_jettisoned) {
            current_mass   -= l110_dry_mass_kg;
            l110_jettisoned = true;
            std::cout << "l110 jettisoned, " << current_mass / 1000.0 << "tonnes, " << altitude_m / 1000.0 << "kms, " << "\n";
        }
        mass_flow_rate   = 44.5;
        exit_velocity    = 4345.0;
        exit_pa          = 5000.0;
        nozzle_exit_area = 2.45;
        double burn = mass_flow_rate * dt;
        prop_c25      = std::max(0.0, prop_c25 - burn);
        current_mass -= burn;
    }
    // Coasting: all propellant expended
    else {
        mass_flow_rate   = 0.0;
        exit_velocity    = 0.0;
        exit_pa          = 101325.0;
        nozzle_exit_area = 0.0;
    }

    // add to your state variables
    const double plf_mass_kg = 1500.0;
    const double plf_jettison_alt = 115000.0; // m

    if (!plf_jettisoned && altitude_m > plf_jettison_alt) {
      current_mass -= plf_mass_kg;
      plf_jettisoned = true;
      std::cout << "plf jettisoned, " << current_mass / 1000.0 << "tonnes, " << altitude_m / 1000.0 << "kms, " << "\n";
    }
}



/**
 * @brief Calculates aerodynamic drag force acting on the LVM3.
 *
 * Uses a piecewise U.S. Standard Atmosphere model for density, a physics-motivated
 * Cd(Mach) curve, and velocity-gated booster separation.
 *
 * @param altitude_m   Height above Earth's surface in metres.
 * @param velocity_m_s Rocket velocity (signed: positive = upward).
 * @param ambient_pressure_pa Ambient atmospheric pressure in Pascals.
 * @return Drag force magnitude in Newtons (always >= 0; caller must apply direction).
 */
double get_lvm3_drag_force(double altitude_m,
                           double velocity_m_s,
                           double ambient_pressure_pa)
{
    // ── Physical constants ────────────────────────────────────────────────────
    static constexpr double R_AIR        = 287.05;  // J/(kg·K), specific gas constant for dry air
    static constexpr double GAMMA_AIR    = 1.4;     // adiabatic index
    static constexpr double KARMAN_LINE  = 100'000.0; // m — conventional atmosphere boundary

    // ── LVM3 geometry ────────────────────────────────────────────────────────
    // Frontal area with both S200 strap-on boosters attached.
    // Core diameter ~4.0 m (area ≈ 12.6 m²) + two 3.2 m boosters ≈ 25 m² total projected.
    static constexpr double AREA_WITH_BOOSTERS    = 25.0;  // m²
    // After S200 separation only the L110 core + C25 cryogenic stage remain.
    static constexpr double AREA_WITHOUT_BOOSTERS = 12.6;  // m²
    // S200 boosters burn out and separate at ~Mach 5 (≈1700 m/s); using speed as
    // the separation trigger is more physical than a fixed altitude threshold.
    static constexpr double BOOSTER_SEP_SPEED     = 1700.0; // m/s

    // ── Early exits ──────────────────────────────────────────────────────────
    if (altitude_m > KARMAN_LINE || ambient_pressure_pa <= 0.0)
        return 0.0;

    // ── 1. Atmospheric temperature — U.S. Standard Atmosphere 1976 (simplified) ──
    //   Layer boundaries: 0–11 km (troposphere), 11–20 km (tropopause),
    //   20–32 km (lower stratosphere), 32–47 km (middle stratosphere).
    //   Above 47 km we clamp to the last lapse formula; good enough below 100 km.
    double T;
    if (altitude_m < 11'000.0) {
        T = 288.15 - 0.0065 * altitude_m;                       // troposphere
    } else if (altitude_m < 20'000.0) {
        T = 216.65;                                              // tropopause (isothermal)
    } else if (altitude_m < 32'000.0) {
        T = 216.65 + 0.001  * (altitude_m - 20'000.0);          // lower stratosphere
    } else {
        T = 228.65 + 0.0028 * (altitude_m - 32'000.0);          // middle stratosphere
    }

    // Sanity floor — prevents division by zero or negative density on bad inputs.
    if (T < 1.0) T = 1.0;

    // ── 2. Air density ρ = P / (R·T) ─────────────────────────────────────────
    double rho = ambient_pressure_pa / (R_AIR * T);

    // ── 3. Mach number ────────────────────────────────────────────────────────
    double speed_of_sound = std::sqrt(GAMMA_AIR * R_AIR * T); // always > 0 given T ≥ 1
    double mach            = std::abs(velocity_m_s) / speed_of_sound;

    // ── 4. Drag coefficient Cd(Mach) ──────────────────────────────────────────
    //
    //  Regime        | Mach range | Behaviour
    //  --------------|------------|------------------------------------------
    //  Subsonic      | < 0.8      | Nearly constant, mild compressibility rise
    //  Transonic     | 0.8 – 1.2  | Sharp rise to wave-drag peak (Max-Q region)
    //  Low supersonic| 1.2 – 5.0  | Decreasing; modelled as Cd ∝ 1/Mach² per
    //                |            | thin-body supersonic wave-drag theory
    //  Hypersonic    | ≥ 5.0      | Asymptotic floor; Newtonian flow regime
    //
    double Cd;
    if (mach < 0.8) {
        // Prandtl–Glauert compressibility correction lifts Cd slightly near 0.8.
        Cd = 0.20 + 0.02 * (mach / 0.8);                         // 0.20 → 0.22
    } else if (mach < 1.2) {
        // Transonic: linear ramp to wave-drag peak at Mach 1.2.
        double t = (mach - 0.8) / 0.4;                           // 0 → 1
        Cd = 0.22 + 0.38 * t;                                    // 0.22 → 0.60
    } else if (mach < 5.0) {
        // Supersonic decay: Cd ∝ 1/Mach² anchored to peak at Mach 1.2.
        // Gives Cd(5) ≈ 0.60 * (1.2/5.0)² ≈ 0.035 — too low, so we add a floor.
        Cd = std::max(0.25, 0.60 * (1.2 * 1.2) / (mach * mach));
    } else {
        Cd = 0.25; // Hypersonic asymptote (Newtonian / blunt-body floor)
    }

    // ── 5. Reference area — velocity-gated booster separation ────────────────
    double frontal_area = (std::abs(velocity_m_s) >= BOOSTER_SEP_SPEED)
                          ? AREA_WITHOUT_BOOSTERS
                          : AREA_WITH_BOOSTERS;

    // ── 6. Drag equation: Fd = ½ · ρ · v² · Cd · A ───────────────────────────
    double drag_force = 0.5 * rho * (velocity_m_s * velocity_m_s) * Cd * frontal_area;
    return drag_force; // magnitude; caller negates along velocity direction
}
