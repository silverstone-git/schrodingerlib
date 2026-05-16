#include "HeatEquation.h"
#include <cmath>
#include <TMath.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TSystem.h>

double HeatEquation::solve1D(double pos_x_m, double time_s, double thermal_diffusivity_m2_s) {
    if (time_s <= 0) return 0.0;
    return (1.0 / std::sqrt(4.0 * TMath::Pi() * thermal_diffusivity_m2_s * time_s)) * std::exp(-(pos_x_m * pos_x_m) / (4.0 * thermal_diffusivity_m2_s * time_s));
}

double HeatEquation::solve3D(double pos_x_m, double pos_y_m, double time_s, double thermal_diffusivity_m2_s) {
    if (time_s <= 0) return 0.0;
    return (1.0 / (4.0 * TMath::Pi() * thermal_diffusivity_m2_s * time_s)) * std::exp(-(pos_x_m * pos_x_m + pos_y_m * pos_y_m) / (4.0 * thermal_diffusivity_m2_s * time_s));
}

void HeatEquation::generate1DGIF(double thermal_diffusivity_m2_s, int num_frames, const std::string& filename) {
    gSystem->Unlink(filename.c_str()); 
    auto c1_heat = new TCanvas("c1_heat", "Heat Equation 1D", 800, 600);
    auto heat_func = new TF1("heat_func", "1.0/sqrt(4*pi*[1]*[0]) * exp(-(x*x)/(4*[1]*[0]))", -10, 10);

    heat_func->SetParameter(1, thermal_diffusivity_m2_s);
    heat_func->SetLineColor(kBlue);
    heat_func->SetMinimum(0);
    heat_func->SetMaximum(1);

    for (int i = 1; i <= num_frames; ++i) {
        double time_s = i * 0.2;
        heat_func->SetParameter(0, time_s);
        heat_func->SetTitle(Form("Time t = %.1f s;x (m);u(x,t)", time_s));
        
        heat_func->Draw();
        c1_heat->Modified();
        c1_heat->Update();

        if (i == 1) {
            c1_heat->Print(filename.c_str());
        } else {
            c1_heat->Print((filename + "+").c_str());
        }
    }
    delete heat_func;
    delete c1_heat;
}

void HeatEquation::generate3DGIF(double thermal_diffusivity_m2_s, int num_frames, const std::string& filename) {
    gSystem->Unlink(filename.c_str()); 
    auto c3_heat = new TCanvas("c3_heat", "Heat Equation 3D", 800, 600);
    auto heat3D_func = new TF2("heat3D_func", "(1.0/(4*pi*[1]*[0])) * exp(-(x*x + y*y)/(4*[1]*[0]))", -5, 5, -5, 5);

    heat3D_func->SetParameter(1, thermal_diffusivity_m2_s);
    heat3D_func->SetTitle("3D Heat Diffusion;x (m);y (m);Temperature");

    gStyle->SetPalette(kBird); 
    heat3D_func->SetMinimum(0);
    heat3D_func->SetMaximum(0.5); 

    for (int i = 1; i <= num_frames; ++i) {
        double time_s = i * 0.015;
        heat3D_func->SetParameter(0, time_s);
        
        heat3D_func->Draw("SURF1");
        
        c3_heat->Modified();
        c3_heat->Update();

        if (i == 1) {
            c3_heat->Print(filename.c_str());
        } else {
            c3_heat->Print((filename + "+").c_str());
        }
    }
    delete heat3D_func;
    delete c3_heat;
}
