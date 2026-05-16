#include "Environment.h"
#include <cmath>
#include <array>
#include <algorithm>

const double G_j_m_kg2 = 6.67e-11;
const double M_earth_kg = 5.9722e24;
const double R_earth_m = 6378100.0;
const double R_AIR_j_kg_k = 287.05;
const double GAMMA_AIR = 1.4;

double EarthEnv::get_radius_m() const {
    return R_earth_m;
}

double EarthEnv::get_gravity_m_s2(double altitude_m) const {
    double r_m = R_earth_m + altitude_m;
    return (G_j_m_kg2 * M_earth_kg) / (r_m * r_m);
}

double EarthEnv::get_pressure_pa(double geometric_altitude_m) const {
    double h_geopotential_m = (R_earth_m * geometric_altitude_m) / (R_earth_m + geometric_altitude_m);

    if (h_geopotential_m < 0.0) return 101325.0; 
    if (h_geopotential_m > 86000.0) return 0.0;  

    struct AtmosphericLayer {
        double base_h_m;
        double base_P_pa;
        double base_T_k;
        double lapse_L_k_m;
    };

    static const std::array<AtmosphericLayer, 7> layers = {{
        {0.0,     101325.0,   288.15, -0.0065}, 
        {11000.0, 22632.1,    216.65,  0.0000}, 
        {20000.0, 5474.89,    216.65,  0.0010}, 
        {32000.0, 868.019,    228.65,  0.0028}, 
        {47000.0, 110.906,    270.65,  0.0000}, 
        {51000.0, 66.9389,    270.65, -0.0028}, 
        {71000.0, 3.95642,    214.65, -0.0020}  
    }};

    auto it = std::upper_bound(layers.begin(), layers.end(), h_geopotential_m, 
        [](double val, const AtmosphericLayer& layer) { return val < layer.base_h_m; });
    
    const auto& layer = *(std::prev(it));

    constexpr double g0_m_s2 = 9.80665;
    constexpr double M_kg_mol  = 0.0289644;
    constexpr double R_j_mol_k  = 8.31432;

    double delta_h_m = h_geopotential_m - layer.base_h_m;

    if (std::abs(layer.lapse_L_k_m) < 1e-9) {
        return layer.base_P_pa * std::exp((-g0_m_s2 * M_kg_mol * delta_h_m) / (R_j_mol_k * layer.base_T_k));
    } else {
        double T_local_k = layer.base_T_k + layer.lapse_L_k_m * delta_h_m;
        double exponent = (-g0_m_s2 * M_kg_mol) / (R_j_mol_k * layer.lapse_L_k_m);
        return layer.base_P_pa * std::pow(T_local_k / layer.base_T_k, exponent);
    }
}

double EarthEnv::get_temperature_k(double altitude_m) const {
    if (altitude_m < 11000.0) return 288.15 - 0.0065 * altitude_m;
    if (altitude_m < 20000.0) return 216.65;
    if (altitude_m < 32000.0) return 216.65 + 0.001 * (altitude_m - 20000.0);
    return 228.65 + 0.0028 * (altitude_m - 32000.0);
}

double EarthEnv::get_density_kg_m3(double altitude_m) const {
    if (altitude_m > 100000.0) return 0.0;
    double p = get_pressure_pa(altitude_m);
    if (p <= 0) return 0.0;
    double t = get_temperature_k(altitude_m);
    if (t < 1.0) t = 1.0;
    return p / (R_AIR_j_kg_k * t);
}

double EarthEnv::get_speed_of_sound_m_s(double altitude_m) const {
    double t = get_temperature_k(altitude_m);
    if (t < 1.0) t = 1.0;
    return std::sqrt(GAMMA_AIR * R_AIR_j_kg_k * t);
}
