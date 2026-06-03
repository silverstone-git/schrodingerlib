#include "TVector3.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPolyLine3D.h"
#include "TAxis.h"
#include "TView.h"
#include "TH3F.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoSphere.h"
#include "schrodingerlib.cpp"
#ifdef __CLING__
#pragma link C++ class std::vector<TVector3>+;
#endif
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

void shar_flight() {
    // --- SIMULATION CONFIGURATION ---
    const int n_ascent = 10000; 
    double dt = 0.1;
    
    const int n_orbit = 20000;  
    double dt_orbit = 1.0;

    // Environmental & Planetary Constants (WGS 84 Reference Ellipsoid Model)
    const double R_a = 6378137.0; // Semi-major axis (Equatorial)
    const double R_b = 6356752.3; // Semi-minor axis (Polar)
    double lat_rad = 13.71979552320247 * M_PI / 180.0;
    double lon_rad = 80.23038514341116 * M_PI / 180.0;
    
    const double R_pad = 6376924.0; 
    const double G = 6.6743e-11; 
    const double M = 5.9722e+24;

    // Trackers for Rocket Ascent Phase
    std::vector<TVector3> pos(n_ascent);
    std::vector<TVector3> vel(n_ascent);
    std::vector<TVector3> acc(n_ascent);
    std::vector<double> time_pad(n_ascent);
    std::vector<double> thrust_pad(n_ascent);
    std::vector<double> mass_log(n_ascent);

    // Initial Positioning Setup
    TVector3 initial_pos(R_pad * std::cos(lat_rad) * std::cos(lon_rad),
                         R_pad * std::cos(lat_rad) * std::sin(lon_rad),
                         R_pad * std::sin(lat_rad));
    pos.at(0) = initial_pos;
    vel.at(0).SetXYZ(0, 0, 0);
    time_pad.at(0) = 0.0;

    double m_rocket = 640000.0; 
    double dry_mass = 32000.0; 
    double prop_s200 = 410000.0, prop_l110 = 116000.0, prop_c25 = 28000.0;
    bool s200_jettisoned = false, l110_jettisoned = false, plf_jettisoned = false, l110_ignited = false;

    double mass_flow, exit_v, exit_p, nozzle_a, atm, thrust;
    TVector3 v_half(0,0,0);

    TVector3 local_up = pos.at(0).Unit();
    TVector3 global_z(0, 0, 1);
    TVector3 local_east = global_z.Cross(local_up).Unit();

    // --- PHASE 1: POWERED ASCENT SIMULATION ---
    for (int i = 0; i < n_ascent - 1; i++) {
        double current_t = time_pad.at(i);
        pos.at(i + 1) = pos.at(i) + v_half * dt;

        double current_alt = pos.at(i + 1).Mag() - R_pad;
        if (current_alt <= 0.0) {
            pos.at(i + 1) = local_up * R_pad; 
            v_half.SetXYZ(0,0,0);
            current_alt = 0.0;
        }
        
        atm = get_ambient_pressure(current_alt, R_pad);
        get_lvm3_stage_properties_3D(current_t, dt, m_rocket, mass_flow, exit_v, exit_p, nozzle_a,
                                     s200_jettisoned, l110_jettisoned, plf_jettisoned, 
                                     prop_s200, prop_l110, prop_c25, l110_ignited);
        
        thrust = (mass_flow * exit_v) + (exit_p - atm) * nozzle_a;
        thrust_pad.at(i + 1) = thrust;
        mass_log.at(i) = m_rocket;

        local_up = pos.at(i + 1).Unit();
        local_east = global_z.Cross(local_up).Unit(); 
        double yaw = get_rocket_yaw_radians(current_t);
        double tilt = (M_PI / 2.0) - yaw; 
        TVector3 thrust_dir = (local_up * std::cos(tilt) + local_east * std::sin(tilt)).Unit();

        TVector3 g_vec = local_up * (- (G * M) / pos.at(i + 1).Mag2());
        TVector3 thrust_acc_vec = thrust_dir * (thrust / m_rocket);
        
        if (current_alt <= 1e-3 && (thrust / m_rocket) <= g_vec.Mag()) {
            acc.at(i + 1).SetXYZ(0,0,0);
            v_half.SetXYZ(0,0,0);
            vel.at(i + 1).SetXYZ(0,0,0);
            thrust_pad.at(i + 1) = 0.0; 
        } 
        else {
            if (m_rocket > dry_mass) m_rocket -= mass_flow * dt;
            
            double drag_f = get_lvm3_drag_force(current_alt, v_half.Mag(), atm);
            TVector3 drag_acc_vec = (v_half.Mag() > 0.01) ? v_half.Unit() * (-drag_f / m_rocket) : TVector3(0,0,0);
            
            acc.at(i + 1) = thrust_acc_vec + g_vec + drag_acc_vec;
            v_half = v_half + acc.at(i + 1) * dt;
            vel.at(i + 1) = v_half; 
        }
        time_pad.at(i + 1) = current_t + dt;
    }
    mass_log.at(n_ascent - 1) = m_rocket;

    // --- PHASE 2: ORBITAL SEPARATION ---
    std::vector<TVector3> payload_pos(n_orbit);
    std::vector<TVector3> payload_vel(n_orbit);
    std::vector<TVector3> spent_stage_pos(n_orbit);
    std::vector<TVector3> spent_stage_vel(n_orbit);

    TVector3 seco_pos = pos.at(n_ascent - 1);
    TVector3 seco_vel = vel.at(n_ascent - 1);
    TVector3 separation_dir = seco_vel.Unit(); 

    payload_pos.at(0) = seco_pos;
    payload_vel.at(0) = seco_vel + (separation_dir * 1.5); 

    spent_stage_pos.at(0) = seco_pos;
    spent_stage_vel.at(0) = seco_vel - (separation_dir * 0.2); 

    for (int i = 0; i < n_orbit - 1; i++) {
        TVector3 r_p = payload_pos.at(i);
        TVector3 g_p = r_p.Unit() * (- (G * M) / r_p.Mag2());
        payload_pos.at(i + 1) = r_p + payload_vel.at(i) * dt_orbit + g_p * (0.5 * dt_orbit * dt_orbit);
        TVector3 g_p_next = payload_pos.at(i + 1).Unit() * (- (G * M) / payload_pos.at(i + 1).Mag2());
        payload_vel.at(i + 1) = payload_vel.at(i) + (g_p + g_p_next) * (0.5 * dt_orbit);

        TVector3 r_s = spent_stage_pos.at(i);
        TVector3 g_s = r_s.Unit() * (- (G * M) / r_s.Mag2());
        spent_stage_pos.at(i + 1) = r_s + spent_stage_vel.at(i) * dt_orbit + g_s * (0.5 * dt_orbit * dt_orbit);
        TVector3 g_s_next = spent_stage_pos.at(i + 1).Unit() * (- (G * M) / spent_stage_pos.at(i + 1).Mag2());
        spent_stage_vel.at(i + 1) = spent_stage_vel.at(i) + (g_s + g_s_next) * (0.5 * dt_orbit);
    }

    // =========================================================================
    // GENERATE GRAPH 1: POWERED PERFORMANCE PROFILES
    // =========================================================================
    TCanvas *c1 = new TCanvas("c1", "LVM3 Powered Ascent Phase Profiles", 1200, 900);
    c1->Divide(2, 2);

    TGraph *g_alt = new TGraph(); g_alt->SetTitle("Ellipsoidal Altitude Profile;Time (s);Altitude (m)");
    TGraph *g_acc = new TGraph(); g_acc->SetTitle("Flight Acceleration Profile;Time (s);Acc (m/s^{2})");
    TGraph *g_thr = new TGraph(); g_thr->SetTitle("Engine Dynamic Thrust Profile;Time (s);Thrust (N)");

    for (int i = 0; i < n_ascent; i++) {
        g_alt->SetPoint(i, time_pad.at(i), pos.at(i).Mag() - R_pad);
        g_acc->SetPoint(i, time_pad.at(i), acc.at(i).Mag());
        g_thr->SetPoint(i, time_pad.at(i), thrust_pad.at(i));
    }

    c1->cd(1); g_alt->Draw("AL");
    c1->cd(2); g_acc->Draw("AL");
    c1->cd(3); g_thr->Draw("AL");

    c1->cd(4);
    TPolyLine3D *path3d = new TPolyLine3D(n_ascent);
    path3d->SetLineColor(kRed + 1); path3d->SetLineWidth(3);
    for (int i = 0; i < n_ascent; i++) {
        path3d->SetPoint(i, pos.at(i).X() - initial_pos.X(), pos.at(i).Y() - initial_pos.Y(), pos.at(i).Z() - initial_pos.Z());
    }
    TH3F *frame3d = new TH3F("frame3d", "3D Trajectory Relative to Sriharikota Pad;Local X (m);Local Y (m);Local Z (m)", 
                            10, -100000, 1500000, 10, -100000, 1500000, 10, -100000, 1500000);
    frame3d->Draw(); path3d->Draw("same");
    int iret; TView *view1 = gPad->GetView(); if (view1) view1->SetView(75.0, 55.0, 110.0, iret);
    c1->Print("ascent_performance.png");

    // =========================================================================
    // GENERATE GRAPH 2: GEOCENTRIC SCALED MISSION TRACK (CALIBRATED TO KM)
    // =========================================================================
    TCanvas *c2 = new TCanvas("c2", "Geocentric Multi-Object Drift Map with Planetary Reference", 1000, 1000);
    
    // Normalized Coordinate System Frame: Bounds scaled cleanly to +/- 60,000 kilometers
    double limit3D_km = 60000.0; 
    TH3F *orbit_frame = new TH3F("orbit_frame", "Geocentric Multi-Object Drift Map (Scale: km);Inertial X (km);Inertial Y (km);Inertial Z (km)", 
                                 10, -limit3D_km, limit3D_km, 10, -limit3D_km, limit3D_km, 10, -limit3D_km, limit3D_km);
    orbit_frame->Draw();

    // 1. Initialize the ROOT Geometry Manager in kilometer units
    if (gGeoManager) { delete gGeoManager; gGeoManager = nullptr; } // Fresh reset
    gGeoManager = new TGeoManager("world", "Planetary Mission Tracking World");
    TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0, 0, 0);
    TGeoMedium *vacuum = new TGeoMedium("Vacuum", 1, matVacuum);
    TGeoVolume *top = gGeoManager->MakeBox("TOP", vacuum, 1e5, 1e5, 1e5);
    gGeoManager->SetTopVolume(top);

    // 2. Build the Earth Reference Ellipsoid in clean kilometer dimensions
    double R_a_km = R_a / 1000.0; // 6378.137 km
    double R_b_km = R_b / 1000.0; // 6356.752 km

    TGeoVolume *earth_sphere = gGeoManager->MakeSphere("EARTH_BASE", vacuum, 0.0, 1.0); 
    earth_sphere->SetLineColor(kCyan - 7); 
    earth_sphere->SetFillColor(kCyan - 9);

    TGeoScale *scale_ellipsoid = new TGeoScale("wgs84_km", R_a_km, R_a_km, R_b_km);
    top->AddNode(earth_sphere, 1, scale_ellipsoid);
    gGeoManager->CloseGeometry();

    // Render the calibrated true-scale planet
    top->Draw("same");

    // 3. Down-convert trajectory paths to kilometers for seamless visual sync
    TPolyLine3D *sat_line = new TPolyLine3D(n_orbit);
    sat_line->SetLineColor(kBlue + 1); sat_line->SetLineWidth(3);
    
    TPolyLine3D *stage_line = new TPolyLine3D(n_orbit);
    stage_line->SetLineColor(kMagenta - 3); stage_line->SetLineWidth(2); stage_line->SetLineStyle(2);

    for (int i = 0; i < n_orbit; i++) {
        sat_line->SetPoint(i, payload_pos.at(i).X() / 1000.0, payload_pos.at(i).Y() / 1000.0, payload_pos.at(i).Z() / 1000.0);
        stage_line->SetPoint(i, spent_stage_pos.at(i).X() / 1000.0, spent_stage_pos.at(i).Y() / 1000.0, spent_stage_pos.at(i).Z() / 1000.0);
    }
    sat_line->Draw("same");
    stage_line->Draw("same");

    TPolyLine3D *ascent_segment = new TPolyLine3D(n_ascent);
    ascent_segment->SetLineColor(kRed); ascent_segment->SetLineWidth(2);
    for (int i = 0; i < n_ascent; i++) {
        ascent_segment->SetPoint(i, pos.at(i).X() / 1000.0, pos.at(i).Y() / 1000.0, pos.at(i).Z() / 1000.0);
    }
    ascent_segment->Draw("same");

    TView *view2 = gPad->GetView(); 
    if (view2) {
        view2->SetView(60.0, 30.0, 125.0, iret); 
    }
    
    gPad->Modified(); gPad->Update();
    c2->Print("payload_orbital_insertion.png");

    // =========================================================================
    // AUTOMATIC LOG FILE GENERATOR SYSTEM
    // =========================================================================
    std::ofstream logfile("flight_simulation.log");
    if (!logfile.is_open()) {
        std::cerr << "!!! Error: Unable to automatically initialize local flight log file.\n";
        return;
    }

    logfile << "=================================================================================\n";
    logfile << "                 LVM3 SRIHARIKOTA MISSION FLIGHT TELEMETRY LOG                   \n";
    logfile << "=================================================================================\n\n";
    
    logfile << "[MISSION METRICS BASELINE]\n";
    logfile << "  Launch Facility : Satish Dhawan Space Centre (SDSC-SHAR)\n";
    logfile << "  Coordinates     : Lat 13.72 N, Lon 80.23 E\n";
    logfile << "  Reference Model : WGS 84 Ellipsoid Ground Plane Boundary Baseline\n\n";

    logfile << "---------------------------------------------------------------------------------\n";
    logfile << " TIME(s)   ALT(km)     VEL(m/s)    ACC(m/s2)   MASS(Tonnes)   THRUST(kN)  PHASE\n";
    logfile << "---------------------------------------------------------------------------------\n";

    for (int i = 0; i < n_ascent; i += 500) {
        logfile << std::fixed << std::setprecision(1)
                << std::setw(6)  << time_pad.at(i) << "   "
                << std::setw(8)  << (pos.at(i).Mag() - R_pad) / 1000.0 << "    "
                << std::setw(8)  << vel.at(i).Mag() << "    "
                << std::setw(8)  << acc.at(i).Mag() << "    "
                << std::setw(10) << mass_log.at(i) / 1000.0 << "    "
                << std::setw(10) << thrust_pad.at(i) / 1000.0 << "    "
                << "POWERED_ASCENT\n";
    }

    logfile << "\n[T+ 1000.0 s] >>> MISSION EVENT: DETACHMENT & DUAL ORBITAL SEPARATION SECO\n\n";

    for (int i = 0; i < n_orbit; i += 1000) {
        double t_orbit = 1000.0 + (i * dt_orbit);
        logfile << std::fixed << std::setprecision(1)
                << std::setw(6)  << t_orbit << "   "
                << std::setw(8)  << (payload_pos.at(i).Mag() - R_pad) / 1000.0 << "    "
                << std::setw(8)  << payload_vel.at(i).Mag() << "    "
                << std::setw(8)  << (- (G * M) / payload_pos.at(i).Mag2()) << "    "
                << std::setw(10) << 4.0 << "    " 
                << std::setw(10) << 0.0 << "    " 
                << "COAST_PAYLOAD\n";
                
        logfile << std::fixed << std::setprecision(1)
                << std::setw(6)  << t_orbit << "   "
                << std::setw(8)  << (spent_stage_pos.at(i).Mag() - R_pad) / 1000.0 << "    "
                << std::setw(8)  << spent_stage_vel.at(i).Mag() << "    "
                << std::setw(8)  << (- (G * M) / spent_stage_pos.at(i).Mag2()) << "    "
                << std::setw(10) << dry_mass / 1000.0 << "    "
                << std::setw(10) << 0.0 << "    "
                << "COAST_SPENT_C25\n";
        logfile << "---------------------------------------------------------------------------------\n";
    }

    logfile.close();
    std::cout << ">>> Execution Complete!\n";
    std::cout << "    - Plot 1 Exported: ascent_performance.png\n";
    std::cout << "    - Plot 2 Exported: payload_orbital_insertion.png (Calibrated Kilometer Scene)\n";
    std::cout << "    - Automated Data File Saved: flight_simulation.log\n";
}
