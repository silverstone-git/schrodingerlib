#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

/**
 * @brief Interface for any planetary environment.
 */
class IEnvironment {
public:
    virtual ~IEnvironment() = default;
    
    // Core physical properties
    virtual double get_radius_m() const = 0;
    virtual double get_gravity_m_s2(double altitude_m) const = 0;
    
    // Atmospheric properties
    virtual double get_pressure_pa(double altitude_m) const = 0;
    virtual double get_density_kg_m3(double altitude_m) const = 0;
    virtual double get_speed_of_sound_m_s(double altitude_m) const = 0;
};

/**
 * @brief Concrete implementation of Earth using U.S. Standard Atmosphere (1976).
 */
class EarthEnv : public IEnvironment {
public:
    double get_radius_m() const override;
    double get_gravity_m_s2(double altitude_m) const override;
    double get_pressure_pa(double altitude_m) const override;
    double get_density_kg_m3(double altitude_m) const override;
    double get_speed_of_sound_m_s(double altitude_m) const override;

private:
    double get_temperature_k(double altitude_m) const;
};

#endif
