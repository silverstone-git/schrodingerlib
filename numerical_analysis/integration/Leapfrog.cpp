#include "Leapfrog.h"
#include <cmath>
#include <algorithm>
#include <TGraph.h>
#include <TCanvas.h>
#include <TMath.h>

LeapfrogIntegrator::LeapfrogIntegrator(IEnvironment* env, IRocket* rocket, double dt_s)
    : env(env), rocket(rocket), dt_s(dt_s) {}

std::vector<State> LeapfrogIntegrator::simulate(double t_max_s, double initial_altitude_m, double initial_velocity_m_s) {
    std::vector<State> results;
    double t_s = 0.0;
    
    // Absolute position from center of planet
    double x_m = initial_altitude_m + env->get_radius_m();
    double v_m_s = initial_velocity_m_s;
    double a_m_s2 = 0.0;
    
    // Initial Environment conditions
    double current_alt_m = x_m - env->get_radius_m();
    double g_m_s2 = env->get_gravity_m_s2(current_alt_m);
    double atm_pa = env->get_pressure_pa(current_alt_m);
    double rho_kg_m3 = env->get_density_kg_m3(current_alt_m);
    double sos_m_s = env->get_speed_of_sound_m_s(current_alt_m);
    
    // Init rocket
    rocket->update_dt(0.0, current_alt_m, v_m_s);
    double thrust_n = rocket->get_thrust_n(atm_pa);
    double thrust_accn_m_s2 = thrust_n / rocket->get_mass_kg();
    
    if (thrust_accn_m_s2 > g_m_s2) {
        a_m_s2 = thrust_accn_m_s2 - g_m_s2;
    } else {
        a_m_s2 = 0.0;
    }
    
    results.push_back({t_s, current_alt_m, v_m_s, a_m_s2, atm_pa, rocket->get_mass_kg(), thrust_n, 0.0, 0.0, 0.0});
    
    double v_half_m_s = v_m_s + a_m_s2 * (dt_s / 2.0);
    
    while (t_s < t_max_s) {
        x_m += v_half_m_s * dt_s;
        if (x_m <= env->get_radius_m()) {
            x_m = env->get_radius_m();
            v_half_m_s = 0.0;
        }
        
        current_alt_m = x_m - env->get_radius_m();
        
        // Update environment
        g_m_s2 = env->get_gravity_m_s2(current_alt_m);
        atm_pa = env->get_pressure_pa(current_alt_m);
        rho_kg_m3 = env->get_density_kg_m3(current_alt_m);
        sos_m_s = env->get_speed_of_sound_m_s(current_alt_m);
        
        // Update rocket
        rocket->update_dt(dt_s, current_alt_m, v_half_m_s);
        thrust_n = rocket->get_thrust_n(atm_pa);
        thrust_accn_m_s2 = thrust_n / rocket->get_mass_kg();
        
        // Aerodynamics
        double mach = rocket->get_mach(sos_m_s, v_half_m_s);
        double drag_n = rocket->get_drag_force_n(rho_kg_m3, sos_m_s, v_half_m_s);
        double q_pa = 0.5 * rho_kg_m3 * v_half_m_s * v_half_m_s;

        // Kinematics
        if (thrust_accn_m_s2 <= g_m_s2 && current_alt_m <= 1e-3) {
            a_m_s2 = 0.0;
            v_half_m_s = 0.0;
            v_m_s = 0.0;
        } else {
            a_m_s2 = thrust_accn_m_s2 - g_m_s2;
            double drag_accn_m_s2 = drag_n / rocket->get_mass_kg();
            
            if (v_half_m_s > 0) a_m_s2 -= drag_accn_m_s2;
            else a_m_s2 += drag_accn_m_s2;

            v_half_m_s += a_m_s2 * dt_s;
            v_m_s = v_half_m_s - a_m_s2 * (dt_s / 2.0);
        }
        
        t_s += dt_s;
        results.push_back({t_s, current_alt_m, v_m_s, a_m_s2, atm_pa, rocket->get_mass_kg(), thrust_n, drag_n, q_pa, mach});
    }

    return results;
}

void LeapfrogIntegrator::plotTrajectory(const std::vector<State>& states, const std::string& filename) {
    int n = states.size();
    if (n == 0) return;

    std::vector<double> time_s(n), alt_m(n), vel_m_s(n), acc_m_s2(n), press_pa(n);
    std::vector<double> mass_kg(n), thrust_n(n), drag_n(n), q_pa(n), mach(n);

    for (int i = 0; i < n; ++i) {
        time_s[i] = states[i].t_s;
        alt_m[i]  = states[i].altitude_m;
        vel_m_s[i] = states[i].velocity_m_s;
        acc_m_s2[i] = states[i].acceleration_m_s2;
        press_pa[i] = states[i].atm_pressure_pa;
        mass_kg[i] = states[i].mass_kg;
        thrust_n[i] = states[i].thrust_n;
        drag_n[i] = states[i].drag_n;
        q_pa[i] = states[i].dynamic_pressure_pa;
        mach[i] = states[i].mach_number;
    }

    TCanvas* c1 = new TCanvas("c1_telemetry", "Vehicle Telemetry", 1600, 1200);
    c1->Divide(3, 3);

    c1->cd(1);
    TGraph* gr_alt = new TGraph(n, time_s.data(), alt_m.data());
    gr_alt->SetTitle("Altitude vs Time;Time (s);Altitude (m)");
    gr_alt->SetLineColor(kBlue);
    gr_alt->Draw("ALP");

    c1->cd(2);
    TGraph* gr_vel = new TGraph(n, time_s.data(), vel_m_s.data());
    gr_vel->SetTitle("Velocity vs Time;Time (s);Velocity (m/s)");
    gr_vel->SetLineColor(kRed);
    gr_vel->Draw("ALP");

    c1->cd(3);
    TGraph* gr_acc = new TGraph(n, time_s.data(), acc_m_s2.data());
    gr_acc->SetTitle("Net Acceleration vs Time;Time (s);Acceleration (m/s^2)");
    gr_acc->SetLineColor(kBlack);
    gr_acc->Draw("ALP");

    c1->cd(4);
    TGraph* gr_mass = new TGraph(n, time_s.data(), mass_kg.data());
    gr_mass->SetTitle("Total Vehicle Mass vs Time;Time (s);Mass (kg)");
    gr_mass->SetLineColor(kGreen+2);
    gr_mass->Draw("ALP");

    c1->cd(5);
    TGraph* gr_thrust = new TGraph(n, time_s.data(), thrust_n.data());
    gr_thrust->SetTitle("Engine Thrust vs Time;Time (s);Thrust (N)");
    gr_thrust->SetLineColor(kOrange+7);
    gr_thrust->Draw("ALP");

    c1->cd(6);
    TGraph* gr_drag = new TGraph(n, time_s.data(), drag_n.data());
    gr_drag->SetTitle("Aerodynamic Drag vs Time;Time (s);Drag Force (N)");
    gr_drag->SetLineColor(kMagenta);
    gr_drag->Draw("ALP");

    c1->cd(7);
    TGraph* gr_q = new TGraph(n, time_s.data(), q_pa.data());
    gr_q->SetTitle("Dynamic Pressure (q) vs Time;Time (s);Dynamic Pressure (Pa)");
    gr_q->SetLineColor(kTeal);
    gr_q->Draw("ALP");

    c1->cd(8);
    TGraph* gr_mach = new TGraph(n, time_s.data(), mach.data());
    gr_mach->SetTitle("Mach Number vs Time;Time (s);Mach");
    gr_mach->SetLineColor(kAzure+1);
    gr_mach->Draw("ALP");

    c1->cd(9);
    TGraph* gr_press = new TGraph(n, time_s.data(), press_pa.data());
    gr_press->SetTitle("Ambient Pressure vs Time;Time (s);Pressure (Pa)");
    gr_press->SetLineColor(kGray+2);
    gr_press->Draw("ALP");

    c1->Print(filename.c_str());
    
    delete gr_alt;
    delete gr_vel;
    delete gr_acc;
    delete gr_mass;
    delete gr_thrust;
    delete gr_drag;
    delete gr_q;
    delete gr_mach;
    delete gr_press;
    delete c1;
}
