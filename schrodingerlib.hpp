#ifndef schrodingerhpp

#define schrodingerhpp


double get_ambient_pressure(double geometric_altitude_m, double EARTH_RADIUS = 6378100.0);

void get_lvm3_stage_properties(double altitude_m, double dt, double& current_mass,
                               double& mass_flow_rate, double& exit_velocity,
                               double& exit_pa, double& nozzle_exit_area, bool& s200_jettisoned, bool& l110_jettisoned,
                               bool& plf_jettisoned, double& prop_s200, double& prop_l110, double& prop_c25);

void get_lvm3_stage_properties_3D(double current_time, double dt, double& current_mass,
                                  double& mass_flow_rate, double& exit_velocity,
                                  double& exit_pa, double& nozzle_exit_area, 
                                  bool& s200_jettisoned, bool& l110_jettisoned, bool& plf_jettisoned, 
                                  double& prop_s200, double& prop_l110, double& prop_c25, bool& l110_ignited);


double get_lvm3_drag_force(double altitude_m,
                           double velocity_m_s,
                           double ambient_pressure_pa);

double get_rocket_yaw_radians(double current_time_s);

#endif // !schrodingerhpp
