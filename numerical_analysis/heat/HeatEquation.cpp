#include "HeatEquation.h"
#include <cmath>
#include <TMath.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TSystem.h>

double HeatEquation::solve1D(double x, double t, double alpha) {
    if (t <= 0) return 0.0;
    return (1.0 / std::sqrt(4.0 * TMath::Pi() * alpha * t)) * std::exp(-(x * x) / (4.0 * alpha * t));
}

double HeatEquation::solve3D(double x, double y, double t, double alpha) {
    if (t <= 0) return 0.0;
    return (1.0 / (4.0 * TMath::Pi() * alpha * t)) * std::exp(-(x * x + y * y) / (4.0 * alpha * t));
}

void HeatEquation::generate1DGIF(double alpha, int nFrames, const std::string& filename) {
    gSystem->Unlink(filename.c_str()); // Remove existing file to avoid appending indefinitely
    auto c2 = new TCanvas("c1_heat", "Heat Equation GIF", 800, 600);
    auto heat = new TF1("heat", "1.0/sqrt(4*pi*[1]*[0]) * exp(-(x*x)/(4*[1]*[0]))", -10, 10);

    heat->SetParameter(1, alpha);
    heat->SetLineColor(kBlue);
    heat->SetMinimum(0);
    heat->SetMaximum(1);

    for (int i = 1; i <= nFrames; ++i) {
        double t = i * 0.2;
        heat->SetParameter(0, t);
        heat->SetTitle(Form("Time t = %.1f;x;u(x,t)", t));
        
        heat->Draw();
        c2->Modified();
        c2->Update();

        if (i == 1) {
            c2->Print(filename.c_str());
        } else {
            c2->Print((filename + "+").c_str());
        }
    }
    delete heat;
    delete c2;
}

void HeatEquation::generate3DGIF(double alpha, int nFrames, const std::string& filename) {
    gSystem->Unlink(filename.c_str()); // Remove existing file
    auto c3 = new TCanvas("c3_heat", "3D Heat Evolution", 800, 600);
    auto heat3D = new TF2("heat3D", "(1.0/(4*pi*[1]*[0])) * exp(-(x*x + y*y)/(4*[1]*[0]))", -5, 5, -5, 5);

    heat3D->SetParameter(1, alpha);
    heat3D->SetTitle("3D Heat Diffusion;X;Y;Temperature");

    gStyle->SetPalette(kBird); 
    heat3D->SetMinimum(0);
    heat3D->SetMaximum(0.5); 

    for (int i = 1; i <= nFrames; ++i) {
        double t = i * 0.015;
        heat3D->SetParameter(0, t);
        
        heat3D->Draw("SURF1");
        
        c3->Modified();
        c3->Update();

        if (i == 1) {
            c3->Print(filename.c_str());
        } else {
            c3->Print((filename + "+").c_str());
        }
    }
    delete heat3D;
    delete c3;
}
