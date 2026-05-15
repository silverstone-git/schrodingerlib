#include "Leapfrog.h"
#include <cmath>
#include <array>
#include <algorithm>
#include <TGraph.h>
#include <TCanvas.h>

LeapfrogIntegrator::LeapfrogIntegrator(double dt, double initial_mass)
    : dt(dt), current_mass(initial_mass) {}

double LeapfrogIntegrator::get_ambient_pressure(double geometric_altitude_m, double earth_radius) const {
    double h = (earth_radius * geometric_altitude_m) / (earth_radius + geometric_altitude_m);

    if (h < 0.0) return 101325.0; 
    if (h > 86000.0) return 0.0;  

    struct AtmosphericLayer {
        double base_h;
        double base_P;
        double base_T;
        double lapse_L;
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

    auto it = std::upper_bound(layers.begin(), layers.end(), h, 
        [](double val, const AtmosphericLayer& layer) { return val < layer.base_h; });
    
    const auto& layer = *(std::prev(it));

    constexpr double g0 = 9.80665;
    constexpr double M  = 0.0289644;
    constexpr double R  = 8.31432;

    double delta_h = h - layer.base_h;

    if (std::abs(layer.lapse_L) < 1e-9) {
        return layer.base_P * std::exp((-g0 * M * delta_h) / (R * layer.base_T));
    } else {
        double T_local = layer.base_T + layer.lapse_L * delta_h;
        double exponent = (-g0 * M) / (R * layer.lapse_L);
        return layer.base_P * std::pow(T_local / layer.base_T, exponent);
    }
}

void LeapfrogIntegrator::get_lvm3_stage_properties(double altitude_m, double& mass_flow_rate, double& exit_velocity, double& exit_pa, double& nozzle_exit_area) {
    if (altitude_m < 43000.0) {
        mass_flow_rate   = 5460.0;
        exit_velocity    = 2690.0;
        exit_pa          = 60000.0;
        nozzle_exit_area = 11.2;
    }
    else if (altitude_m >= 43000.0 && altitude_m < 60000.0) {
        mass_flow_rate   = 6031.4;
        exit_velocity    = 2710.0;
        exit_pa          = 52000.0;
        nozzle_exit_area = 12.862;
    }
    else if (altitude_m >= 60000.0 && altitude_m < 175000.0) {
        if (!s200_jettisoned) {
            current_mass -= 62000.0;
            s200_jettisoned = true;
        }
        mass_flow_rate   = 571.4;
        exit_velocity    = 2873.0;
        exit_pa          = 52000.0;
        nozzle_exit_area = 1.662;
    }
    else {
        if (!l110_jettisoned) {
            current_mass -= 105000.0;
            l110_jettisoned = true;
        }
        mass_flow_rate   = 44.5;
        exit_velocity    = 4345.0;
        exit_pa          = 5000.0;
        nozzle_exit_area = 2.45;
    }
}

std::vector<State> LeapfrogIntegrator::simulate(double t_max, double initial_x, double initial_v) {
    std::vector<State> results;
    double t = 0.0;
    double x = initial_x;
    double v = initial_v;
    double a = 0.0;
    
    const double G = 6.67e-11;
    const double M = 5.9722e24;
    const double R = 6378100.0;
    double g = (G * M) / (R * R);
    
    double mass_flow_rate, exit_velocity, exit_pa, nozzle_exit_area;
    get_lvm3_stage_properties(0.0, mass_flow_rate, exit_velocity, exit_pa, nozzle_exit_area);
    double atm = get_ambient_pressure(0.0, R);
    
    double thrust = (mass_flow_rate * exit_velocity) + (exit_pa - atm) * nozzle_exit_area;
    double thrust_accn = thrust / current_mass;
    
    if (thrust_accn > g) {
        a = thrust_accn - g;
    } else {
        a = 0.0;
    }
    
    results.push_back({t, x, v, a, atm});
    
    double v_half = v + a * (dt / 2.0);
    
    while (t < t_max) {
        x = x + v_half * dt;
        if (x <= R) {
            x = R;
            v_half = 0.0;
        }
        
        g = (G * M) / (x * x);
        atm = get_ambient_pressure(x - R, R);
        get_lvm3_stage_properties(x - R, mass_flow_rate, exit_velocity, exit_pa, nozzle_exit_area);
        
        current_mass = current_mass - mass_flow_rate * dt;
        
        thrust = (mass_flow_rate * exit_velocity) + (exit_pa - atm) * nozzle_exit_area;
        thrust_accn = thrust / current_mass;
        
        if (thrust_accn <= g && (x - R) <= 1e-3) {
            a = 0.0;
            v_half = 0.0;
            v = 0.0;
        } else {
            a = thrust_accn - g;
            v_half = v_half + a * dt;
            v = v_half - a * (dt / 2.0);
        }
        
        t = t + dt;
        results.push_back({t, x, v, a, atm});
    }

    return results;
}

void LeapfrogIntegrator::plotTrajectory(const std::vector<State>& states, const std::string& filename) {
    int n = states.size();
    if (n == 0) return;

    std::vector<double> x(n);
    std::vector<double> t(n);
    std::vector<double> v(n);
    std::vector<double> a(n);
    std::vector<double> atm(n);

    for (int i = 0; i < n; ++i) {
        t[i] = states[i].t;
        x[i] = states[i].x;
        v[i] = states[i].v;
        a[i] = states[i].a;
        atm[i] = states[i].atm_pressure;
    }

    TCanvas* c1 = new TCanvas("c1_leapfrog", "LVM3 Telemetry", 1200, 800);
    c1->Divide(2, 2);

    c1->cd(1);
    TGraph* gr_x = new TGraph(n, t.data(), x.data());
    gr_x->SetTitle("Altitude vs Time;Time (s);Altitude (m)");
    gr_x->Draw("ALP");

    c1->cd(2);
    TGraph* gr_v = new TGraph(n, t.data(), v.data());
    gr_v->SetTitle("Velocity vs Time;Time (s);Velocity (m/s)");
    gr_v->Draw("ALP");

    c1->cd(3);
    TGraph* gr_a = new TGraph(n, t.data(), a.data());
    gr_a->SetTitle("Acceleration vs Time;Time (s);Acceleration (m/s^2)");
    gr_a->Draw("ALP");

    c1->cd(4);
    TGraph* gr_atm = new TGraph(n, t.data(), atm.data());
    gr_atm->SetTitle("Ambient Pressure vs Time;Time (s);Pressure (Pa)");
    gr_atm->Draw("ALP");

    c1->Print(filename.c_str());
    
    delete gr_x;
    delete gr_v;
    delete gr_a;
    delete gr_atm;
    delete c1;
}
