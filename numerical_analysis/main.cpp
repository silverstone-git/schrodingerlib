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
    std::cout << "5. Rocket Leapfrog Simulation (LVM3)" << std::endl;
    std::cout << "0. Exit" << std::endl;
    std::cout << "Select an option: ";
}

void handle_polynomial() {
    std::cout << "Enter size of polynomial: ";
    int poly_size;
    std::cin >> poly_size;
    if (poly_size < 1) return;

    std::cout << "Enter coefficients (a_n ... a_0): ";
    std::vector<float> coeffs_vector(poly_size);
    for (int i = 0; i < poly_size; ++i) std::cin >> coeffs_vector[i];

    Polynomial poly(coeffs_vector);
    poly.print();
    std::vector<float> found_roots = bisection(poly);
    std::cout << "Roots found: ";
    for (float root_val : found_roots) std::cout << root_val << " ";
    std::cout << std::endl;
}

void handle_heat1d() {
    double alpha_m2_s;
    int num_frames;
    std::string out_file;
    std::cout << "Enter alpha (m^2/s), number of frames, and output file (e.g. 0.5 30 heat1d.gif): ";
    std::cin >> alpha_m2_s >> num_frames >> out_file;
    HeatEquation::generate1DGIF(alpha_m2_s, num_frames, out_file);
    std::cout << "1D Heat Equation GIF generated: " << out_file << std::endl;
}

void handle_heat3d() {
    double alpha_m2_s;
    int num_frames;
    std::string out_file;
    std::cout << "Enter alpha (m^2/s), number of frames, and output file (e.g. 0.5 400 heat3d.gif): ";
    std::cin >> alpha_m2_s >> num_frames >> out_file;
    HeatEquation::generate3DGIF(alpha_m2_s, num_frames, out_file);
    std::cout << "3D Heat Equation GIF generated: " << out_file << std::endl;
}

void handle_audio() {
    std::string input_wav, output_wav;
    double lo_hz, hi_hz;
    std::cout << "Enter input wav file: ";
    std::cin >> input_wav;
    std::cout << "Enter output wav file: ";
    std::cin >> output_wav;
    std::cout << "Enter frequency range to REJECT (lo hi in Hz): ";
    std::cin >> lo_hz >> hi_hz;

    int sample_rate_hz;
    std::vector<double> pcm_samples = AudioProcessor::readWav(input_wav, sample_rate_hz);
    if (pcm_samples.empty()) {
        std::cout << "Failed to read " << input_wav << std::endl;
        return;
    }

    std::vector<double> filtered_pcm = AudioProcessor::filter(pcm_samples, sample_rate_hz, lo_hz, hi_hz);
    AudioProcessor::writeWav(output_wav, filtered_pcm, sample_rate_hz);
    std::cout << "Filtered audio saved to " << output_wav << std::endl;
}

void handle_leapfrog() {
    double dt_s, t_max_s;
    std::string out_file;
    std::cout << "Enter dt (s), total time (s), and output plot file (e.g. 0.1 200 lvm3_telemetry.pdf): ";
    std::cin >> dt_s >> t_max_s >> out_file;

    // Instantiate specific Environment and Rocket via interfaces
    std::unique_ptr<IEnvironment> earth = std::make_unique<EarthEnv>();
    std::unique_ptr<IRocket> lvm3 = std::make_unique<LVM3Rocket>();

    // Inject dependencies into the integrator
    LeapfrogIntegrator integrator(earth.get(), lvm3.get(), dt_s);
    
    // Simulate from ground level with 0 initial velocity
    auto trajectory_states = integrator.simulate(t_max_s, 0.0, 0.0);

    LeapfrogIntegrator::plotTrajectory(trajectory_states, out_file);
    std::cout << "Trajectory plotted and saved to " << out_file << std::endl;
}

int main() {
    int user_choice;
    do {
        show_menu();
        if (!(std::cin >> user_choice)) break;
        switch (user_choice) {
            case 1: handle_polynomial(); break;
            case 2: handle_heat1d(); break;
            case 3: handle_heat3d(); break;
            case 4: handle_audio(); break;
            case 5: handle_leapfrog(); break;
            case 0: std::cout << "Exiting..." << std::endl; break;
            default: std::cout << "Invalid choice!" << std::endl;
        }
    } while (user_choice != 0);

    return 0;
}
