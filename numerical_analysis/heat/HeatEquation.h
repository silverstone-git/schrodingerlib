#ifndef HEATEQUATION_H
#define HEATEQUATION_H

#include <string>

class HeatEquation {
public:
    static double solve1D(double x, double t, double alpha);
    static double solve3D(double x, double y, double t, double alpha);
    
    // ROOT generation
    static void generate1DGIF(double alpha, int nFrames, const std::string& filename);
    static void generate3DGIF(double alpha, int nFrames, const std::string& filename);
};

#endif
