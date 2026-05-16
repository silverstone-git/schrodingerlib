#ifndef HEATEQUATION_H
#define HEATEQUATION_H

#include <string>

class HeatEquation {
public:
    static double solve1D(double pos_x_m, double time_s, double thermal_diffusivity_m2_s);
    static double solve3D(double pos_x_m, double pos_y_m, double time_s, double thermal_diffusivity_m2_s);
    
    // ROOT generation
    static void generate1DGIF(double thermal_diffusivity_m2_s, int num_frames, const std::string& filename);
    static void generate3DGIF(double thermal_diffusivity_m2_s, int num_frames, const std::string& filename);
};

#endif
