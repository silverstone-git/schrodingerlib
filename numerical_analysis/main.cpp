#include "polynomial/Polynomial.h"
#include "integration/Leapfrog.h"
#include "heat/HeatEquation.h"
#include "audio/AudioProcessor.h"
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>

void show_menu() {
    std::cout << "\n--- Numerical Analysis CLI ---" << std::endl;
    std::cout << "1. Polynomial Root Finding" << std::endl;
    std::cout << "2. Heat Equation (1D Solution)" << std::endl;
    std::cout << "3. Heat Equation (3D Solution)" << std::endl;
    std::cout << "4. Audio FFT Filtering" << std::endl;
    std::cout << "5. Rocket Leapfrog Simulation" << std::endl;
    std::cout << "0. Exit" << std::endl;
    std::cout << "Select an option: ";
}

void handle_polynomial() {
    std::cout << "Enter size of polynomial: ";
    int size;
    std::cin >> size;
    if (size < 1) return;

    std::cout << "Enter coefficients (a_n ... a_0): ";
    std::vector<float> coeffs(size);
    for (int i = 0; i < size; ++i) std::cin >> coeffs[i];

    Polynomial p(coeffs);
    p.print();
    std::vector<float> roots = bisection(p);
    std::cout << "Roots found: ";
    for (float r : roots) std::cout << r << " ";
    std::cout << std::endl;
}

void handle_heat1d() {
    double alpha;
    int frames;
    std::string out_file;
    std::cout << "Enter alpha, number of frames, and output file (e.g. 0.5 30 heat1d.gif): ";
    std::cin >> alpha >> frames >> out_file;
    HeatEquation::generate1DGIF(alpha, frames, out_file);
    std::cout << "1D Heat Equation GIF generated: " << out_file << std::endl;
}

void handle_heat3d() {
    double alpha;
    int frames;
    std::string out_file;
    std::cout << "Enter alpha, number of frames, and output file (e.g. 0.5 400 heat3d.gif): ";
    std::cin >> alpha >> frames >> out_file;
    HeatEquation::generate3DGIF(alpha, frames, out_file);
    std::cout << "3D Heat Equation GIF generated: " << out_file << std::endl;
}

void handle_audio() {
    std::string input, output;
    double lo, hi;
    std::cout << "Enter input wav file: ";
    std::cin >> input;
    std::cout << "Enter output wav file: ";
    std::cin >> output;
    std::cout << "Enter frequency range to REJECT (lo hi): ";
    std::cin >> lo >> hi;

    int sampleRate;
    std::vector<double> samples = AudioProcessor::readWav(input, sampleRate);
    if (samples.empty()) {
        std::cout << "Failed to read " << input << std::endl;
        return;
    }

    std::vector<double> filtered = AudioProcessor::filter(samples, sampleRate, lo, hi);
    AudioProcessor::writeWav(output, filtered, sampleRate);
    std::cout << "Filtered audio saved to " << output << std::endl;
}

void handle_leapfrog() {
    double dt, t_max;
    std::string out_file;
    std::cout << "Enter dt, t_max, and output plot file (e.g. 0.1 1000 leapfrog.pdf): ";
    std::cin >> dt >> t_max >> out_file;

    double initial_mass = 640000.0; // 640 Tonnes for LVM3
    LeapfrogIntegrator integrator(dt, initial_mass);
    const double R = 6378100.0;
    auto states = integrator.simulate(t_max, R, 0.0);

    LeapfrogIntegrator::plotTrajectory(states, out_file);
    std::cout << "Trajectory plotted and saved to " << out_file << std::endl;
}

int main() {
    int choice;
    do {
        show_menu();
        std::cin >> choice;
        switch (choice) {
            case 1: handle_polynomial(); break;
            case 2: handle_heat1d(); break;
            case 3: handle_heat3d(); break;
            case 4: handle_audio(); break;
            case 5: handle_leapfrog(); break;
            case 0: std::cout << "Exiting..." << std::endl; break;
            default: std::cout << "Invalid choice!" << std::endl;
        }
    } while (choice != 0);

    return 0;
}
