#include "DipoleSimulation.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip>

#include <TRandom3.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TVirtualFFT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TMath.h>
#include <filesystem>

namespace DipoleSimulation {

    // Helper: Fourier wave number calculation (matches np.fft.fftfreq behavior)
    static double get_k(int i, int N) {
        if (i < (N + 1) / 2) {
            return 2.0 * TMath::Pi() * i / N;
        } else {
            return 2.0 * TMath::Pi() * (i - N) / N;
        }
    }

    // Helper: Normalize histogram as a probability density
    static void normalize_hist(TH1D* hist) {
        double integral = hist->Integral();
        if (integral > 0.0) {
            double bin_width = hist->GetBinWidth(1);
            hist->Scale(1.0 / (integral * bin_width));
        }
    }

    // Structure for Welch PSD output
    struct WelchResult {
        std::vector<double> freqs;
        std::vector<double> psd;
    };

    // Helper: Welch's method for Power Spectral Density (PSD)
    static WelchResult compute_welch(const std::vector<double>& signal, double fs, int nperseg) {
        int N = signal.size();
        int L = nperseg;
        int step = L / 2; // 50% overlap

        // Create Hanning window and compute U
        std::vector<double> w(L);
        double U = 0.0;
        for (int i = 0; i < L; ++i) {
            w[i] = 0.5 * (1.0 - std::cos(2.0 * TMath::Pi() * i / (L - 1)));
            U += w[i] * w[i];
        }

        int num_segments = 0;
        std::vector<double> avg_psd(L / 2 + 1, 0.0);

        // Create FFT object for Real-to-Complex 1D transform
        Int_t n_fft[1] = {L};
        TVirtualFFT* fft_1d = TVirtualFFT::FFT(1, n_fft, "R2C M K");
        if (!fft_1d) {
            std::cerr << "Error: TVirtualFFT R2C could not be created! (Welch FFT skipped)" << std::endl;
            return {};
        }

        std::vector<double> segment_data(L);
        std::vector<double> re_out(L);
        std::vector<double> im_out(L);

        for (int offset = 0; offset + L <= N; offset += step) {
            // 1. Extract segment and compute mean
            double sum = 0.0;
            for (int i = 0; i < L; ++i) {
                sum += signal[offset + i];
            }
            double mean = sum / L;

            // 2. Apply Window
            for (int i = 0; i < L; ++i) {
                segment_data[i] = (signal[offset + i] - mean) * w[i];
            }

            // 3. Perform FFT
            fft_1d->SetPoints(segment_data.data());
            fft_1d->Transform();
            fft_1d->GetPointsComplex(re_out.data(), im_out.data());

            // 4. Accumulate one-sided PSD
            for (int k = 0; k <= L / 2; ++k) {
                double mag_sq = re_out[k] * re_out[k] + im_out[k] * im_out[k];
                double val = mag_sq / (fs * U);
                if (k > 0 && k < L / 2) {
                    val *= 2.0; // Fold negative frequencies
                }
                avg_psd[k] += val;
            }
            num_segments++;
        }

        delete fft_1d;

        WelchResult result;
        result.freqs.resize(L / 2 + 1);
        result.psd.resize(L / 2 + 1);

        for (int k = 0; k <= L / 2; ++k) {
            result.freqs[k] = (k * fs) / L;
            result.psd[k] = avg_psd[k] / num_segments;
        }

        return result;
    }

    void runSimulation() {
        std::cout << "\n=== Starting C++/ROOT Electrostatic Dipole & CIC Simulation ===" << std::endl;
        std::filesystem::create_directories("./outputs");

        // ==========================================
        // 1. 2D Cloud-in-Cell (CIC) Charge Deposition
        // ==========================================
        std::cout << "[1/6] Running 2D CIC Charge Deposition..." << std::endl;
        constexpr int NX_2D = 50;
        constexpr int NY_2D = 50;
        constexpr double X_MIN_2D = 0.0, X_MAX_2D = 50.0;
        constexpr double Y_MIN_2D = 0.0, Y_MAX_2D = 50.0;
        double dx_2d = (X_MAX_2D - X_MIN_2D) / NX_2D;
        double dy_2d = (Y_MAX_2D - Y_MIN_2D) / NY_2D;

        TRandom3 rng(42);
        int num_particles_2d = 50000;

        TH2F* grid_2d = new TH2F("charge_grid", "2D CIC Charge Density Map;X Position;Y Position",
                                 NX_2D, X_MIN_2D, X_MAX_2D, NY_2D, Y_MIN_2D, Y_MAX_2D);

        for (int k = 0; k < num_particles_2d; ++k) {
            double r = rng.Gaus(0.0, 8.0);
            double theta = rng.Uniform(0.0, 2.0 * TMath::Pi());
            double x = 25.0 + r * std::cos(theta);
            double y = 25.0 + r * std::sin(theta);
            double q = 0.1;

            if (x < X_MIN_2D || x >= X_MAX_2D || y < Y_MIN_2D || y >= Y_MAX_2D) {
                continue;
            }

            double g_x = (x - X_MIN_2D) / dx_2d;
            double g_y = (y - Y_MIN_2D) / dy_2d;

            int i = static_cast<int>(std::floor(g_x));
            int j = static_cast<int>(std::floor(g_y));

            double frac_x = g_x - i;
            double frac_y = g_y - j;

            double w_bl = (1.0 - frac_x) * (1.0 - frac_y) * q;
            double w_br = frac_x * (1.0 - frac_y) * q;
            double w_tl = (1.0 - frac_x) * frac_y * q;
            double w_tr = frac_x * frac_y * q;

            int idx_x = std::min(i + 1, NX_2D);
            int idx_y = std::min(j + 1, NY_2D);
            int idx_x_next = std::min(idx_x + 1, NX_2D);
            int idx_y_next = std::min(idx_y + 1, NY_2D);

            grid_2d->SetBinContent(idx_x,      idx_y,      grid_2d->GetBinContent(idx_x,      idx_y)      + w_bl);
            grid_2d->SetBinContent(idx_x_next, idx_y,      grid_2d->GetBinContent(idx_x_next, idx_y)      + w_br);
            grid_2d->SetBinContent(idx_x,      idx_y_next, grid_2d->GetBinContent(idx_x,      idx_y_next) + w_tl);
            grid_2d->SetBinContent(idx_x_next, idx_y_next, grid_2d->GetBinContent(idx_x_next, idx_y_next) + w_tr);
        }

        TCanvas* c_2d = new TCanvas("c_2d", "2D CIC Charge Deposition", 800, 700);
        gStyle->SetOptStat(0);
        gStyle->SetPalette(kCool); // Cool palette
        grid_2d->SetMinimum(0.0);
        grid_2d->Draw("COLZ0");
        c_2d->Update();
        c_2d->SaveAs("outputs/cic_charge_deposition.png");
        delete c_2d;

        // ==========================================
        // 2. 3D Cloud-in-Cell (CIC) Charge Deposition
        // ==========================================
        std::cout << "[2/6] Running 3D CIC Charge Deposition..." << std::endl;
        constexpr int NX_3D = 30, NY_3D = 30, NZ_3D = 30;
        constexpr double X_MIN_3D = 0.0, X_MAX_3D = 30.0;
        constexpr double Y_MIN_3D = 0.0, Y_MAX_3D = 30.0;
        constexpr double Z_MIN_3D = 0.0, Z_MAX_3D = 30.0;
        double dx_3d = 1.0, dy_3d = 1.0, dz_3d = 1.0;

        int num_particles_3d = 100000;
        TH3F* grid_3d = new TH3F("grid_3d", "3D Gaussian Charge Blob;X;Y;Z",
                                 NX_3D, X_MIN_3D, X_MAX_3D, NY_3D, Y_MIN_3D, Y_MAX_3D, NZ_3D, Z_MIN_3D, Z_MAX_3D);

        for (int k = 0; k < num_particles_3d; ++k) {
            double x = rng.Gaus(15.0, 2.0);
            double y = rng.Gaus(15.0, 2.0);
            double z = rng.Gaus(15.0, 2.0);
            double q = 0.1;

            if (x < X_MIN_3D + 1.0 || x >= X_MAX_3D - 1.0 ||
                y < Y_MIN_3D + 1.0 || y >= Y_MAX_3D - 1.0 ||
                z < Z_MIN_3D + 1.0 || z >= Z_MAX_3D - 1.0) {
                continue;
            }

            double gx = (x - X_MIN_3D) / dx_3d;
            double gy = (y - Y_MIN_3D) / dy_3d;
            double gz = (z - Z_MIN_3D) / dz_3d;

            int i = static_cast<int>(std::floor(gx));
            int j = static_cast<int>(std::floor(gy));
            int l = static_cast<int>(std::floor(gz));

            double fx = gx - i;
            double fy = gy - j;
            double fz = gz - l;

            for (int di = 0; di <= 1; ++di) {
                for (int dj = 0; dj <= 1; ++dj) {
                    for (int dl = 0; dl <= 1; ++dl) {
                        double weight = q * (di ? fx : (1.0 - fx)) * (dj ? fy : (1.0 - fy)) * (dl ? fz : (1.0 - fz));
                        grid_3d->SetBinContent(i + di + 1, j + dj + 1, l + dl + 1,
                                               grid_3d->GetBinContent(i + di + 1, j + dj + 1, l + dl + 1) + weight);
                    }
                }
            }
        }
        std::cout << "  3D Gaussian Blob Deposition Complete." << std::endl;

        gStyle->SetCanvasPreferGL(kTRUE);
        TCanvas* c_3d = new TCanvas("c_3d", "3D Gaussian Distribution", 800, 700);
        grid_3d->SetFillColor(kAzure + 1);
        grid_3d->Draw("ISO");
        c_3d->Update();
        c_3d->SaveAs("outputs/cic_charge_3d.png");
        delete c_3d;

        // ==========================================
        // 3. 3D Spectral Poisson Solver
        // ==========================================
        std::cout << "[3/6] Solving Poisson's Equation using 3D Spectral FFT..." << std::endl;
        constexpr int NX_field = 60, NY_field = 60, NZ_field = 60;
        constexpr int N_tot = NX_field * NY_field * NZ_field;
        constexpr double epsilon_0 = 1.0;
        constexpr double sigma = 3.0;
        constexpr double A = 80.0;

        // Coordinates of dipole source centers
        constexpr double x1 = 26.0, y1 = 30.0, z1 = 30.0; // Positive source
        constexpr double x2 = 34.0, y2 = 30.0, z2 = 30.0; // Negative source

        std::vector<double> rho(N_tot, 0.0);
        for (int i = 0; i < NX_field; ++i) {
            for (int j = 0; j < NY_field; ++j) {
                for (int k = 0; k < NZ_field; ++k) {
                    double dx1 = i - x1;
                    double dy1 = j - y1;
                    double dz1 = k - z1;
                    double dx2 = i - x2;
                    double dy2 = j - y2;
                    double dz2 = k - z2;
                    double r1_sq = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
                    double r2_sq = dx2 * dx2 + dy2 * dy2 + dz2 * dz2;

                    double val = A * std::exp(-r1_sq / (2.0 * sigma * sigma)) - A * std::exp(-r2_sq / (2.0 * sigma * sigma));
                    rho[i + NX_field * (j + NY_field * k)] = val;
                }
            }
        }

        // Perform 3D forward C2C FFT on rho
        Int_t n_fft[3] = {NX_field, NY_field, NZ_field};
        TVirtualFFT* fft_forward = TVirtualFFT::FFT(3, n_fft, "C2CFORWARD");
        if (!fft_forward) {
            std::cerr << "Error: TVirtualFFT C2CFORWARD could not be initialized! (Is FFTW3 missing?)" << std::endl;
            delete grid_2d;
            delete grid_3d;
            return;
        }

        std::vector<double> re_in(N_tot, 0.0);
        std::vector<double> im_in(N_tot, 0.0);
        for (int i = 0; i < N_tot; ++i) {
            re_in[i] = rho[i];
        }

        fft_forward->SetPointsComplex(re_in.data(), im_in.data());
        fft_forward->Transform();

        std::vector<double> re_out(N_tot, 0.0);
        std::vector<double> im_out(N_tot, 0.0);
        fft_forward->GetPointsComplex(re_out.data(), im_out.data());

        // Process in Fourier Space: solve potential phi and electric fields Ex, Ey, Ez
        std::vector<double> phi_re(N_tot, 0.0);
        std::vector<double> phi_im(N_tot, 0.0);
        std::vector<double> Ex_re(N_tot, 0.0);
        std::vector<double> Ex_im(N_tot, 0.0);
        std::vector<double> Ey_re(N_tot, 0.0);
        std::vector<double> Ey_im(N_tot, 0.0);
        std::vector<double> Ez_re(N_tot, 0.0);
        std::vector<double> Ez_im(N_tot, 0.0);

        for (int i = 0; i < NX_field; ++i) {
            double kx = get_k(i, NX_field);
            for (int j = 0; j < NY_field; ++j) {
                double ky = get_k(j, NY_field);
                for (int k = 0; k < NZ_field; ++k) {
                    double kz = get_k(k, NZ_field);
                    int idx = i + NX_field * (j + NY_field * k);

                    if (i == 0 && j == 0 && k == 0) {
                        continue;
                    }

                    // Discrete Laplacian eigenvalues
                    double k_sq = 2.0 * (std::cos(kx) - 1.0) + 2.0 * (std::cos(ky) - 1.0) + 2.0 * (std::cos(kz) - 1.0);

                    // Potential: \phi = - ( \rho / \epsilon_0 ) / k_sq
                    double p_re = -re_out[idx] / k_sq;
                    double p_im = -im_out[idx] / k_sq;

                    phi_re[idx] = p_re;
                    phi_im[idx] = p_im;

                    // E-field components: E = -i * k * \phi
                    // E_re + i * E_im = -i * k * (p_re + i * p_im) = k * p_im - i * k * p_re
                    Ex_re[idx] = kx * p_im;
                    Ex_im[idx] = -kx * p_re;

                    Ey_re[idx] = ky * p_im;
                    Ey_im[idx] = -ky * p_re;

                    Ez_re[idx] = kz * p_im;
                    Ez_im[idx] = -kz * p_re;
                }
            }
        }

        // Perform 3D backward C2C FFTs
        TVirtualFFT* fft_backward = TVirtualFFT::FFT(3, n_fft, "C2CBACKWARD");
        if (!fft_backward) {
            std::cerr << "Error: TVirtualFFT C2CBACKWARD could not be initialized!" << std::endl;
            delete fft_forward;
            delete grid_2d;
            delete grid_3d;
            return;
        }

        std::vector<double> phi(N_tot, 0.0);
        std::vector<double> Ex(N_tot, 0.0);
        std::vector<double> Ey(N_tot, 0.0);
        std::vector<double> Ez(N_tot, 0.0);

        auto run_backward = [&](const std::vector<double>& re_src, const std::vector<double>& im_src, std::vector<double>& dest) {
            fft_backward->SetPointsComplex(re_src.data(), im_src.data());
            fft_backward->Transform();
            std::vector<double> temp_re(N_tot, 0.0);
            std::vector<double> temp_im(N_tot, 0.0);
            fft_backward->GetPointsComplex(temp_re.data(), temp_im.data());
            for (int idx = 0; idx < N_tot; ++idx) {
                dest[idx] = temp_re[idx] / N_tot; // Normalize FFTW3 unnormalized transform
            }
        };

        run_backward(phi_re, phi_im, phi);
        run_backward(Ex_re, Ex_im, Ex);
        run_backward(Ey_re, Ey_im, Ey);
        run_backward(Ez_re, Ez_im, Ez);

        delete fft_forward;
        delete fft_backward;
        std::cout << "  3D Poisson Solver Complete." << std::endl;

        // Trilinear Interpolation Lambda
        auto interpolate3D = [&](const std::vector<double>& field_grid, double x, double y, double z) -> double {
            x = std::clamp(x, 0.0, static_cast<double>(NX_field - 1));
            y = std::clamp(y, 0.0, static_cast<double>(NY_field - 1));
            z = std::clamp(z, 0.0, static_cast<double>(NZ_field - 1));

            int x0 = static_cast<int>(std::floor(x));
            int y0 = static_cast<int>(std::floor(y));
            int z0 = static_cast<int>(std::floor(z));

            int x1 = std::min(x0 + 1, NX_field - 1);
            int y1 = std::min(y0 + 1, NY_field - 1);
            int z1 = std::min(z0 + 1, NZ_field - 1);

            double fx = x - x0;
            double fy = y - y0;
            double fz = z - z0;

            auto get_val = [&](int i, int j, int k) {
                return field_grid[i + NX_field * (j + NY_field * k)];
            };

            double c000 = get_val(x0, y0, z0);
            double c100 = get_val(x1, y0, z0);
            double c010 = get_val(x0, y1, z0);
            double c001 = get_val(x0, y0, z1);
            double c110 = get_val(x1, y1, z0);
            double c101 = get_val(x1, y0, z1);
            double c011 = get_val(x0, y1, z1);
            double c111 = get_val(x1, y1, z1);

            double c00 = c000 * (1.0 - fx) + c100 * fx;
            double c01 = c001 * (1.0 - fx) + c101 * fx;
            double c10 = c010 * (1.0 - fx) + c110 * fx;
            double c11 = c011 * (1.0 - fx) + c111 * fx;

            double c0 = c00 * (1.0 - fy) + c10 * fy;
            double c1 = c01 * (1.0 - fy) + c11 * fy;

            return c0 * (1.0 - fz) + c1 * fz;
        };

        // ==========================================
        // 4. Leapfrog Trajectory Simulation
        // ==========================================
        std::cout << "[4/6] Running Leapfrog Integration (100,000 steps)..." << std::endl;
        constexpr double charge_q = -1.0;
        constexpr double mass = 1.0;
        constexpr double dt = 0.005;
        constexpr int steps = 100000;

        double px = 10.0, py = 30.0, pz = 30.0;   // Start away from boundaries
        double vx = 0.5, vy = -0.3, vz = 0.1;

        double ex = interpolate3D(Ex, px, py, pz);
        double ey = interpolate3D(Ey, px, py, pz);
        double ez = interpolate3D(Ez, px, py, pz);

        double ax = (charge_q / mass) * ex;
        double ay = (charge_q / mass) * ey;
        double az = (charge_q / mass) * ez;

        // Initialize velocity at half step: v_{-1/2} = v_0 - 0.5 * a_0 * dt
        double vh_x = vx - 0.5 * ax * dt;
        double vh_y = vy - 0.5 * ay * dt;
        double vh_z = vz - 0.5 * az * dt;

        std::vector<State> trajectory;
        trajectory.reserve(steps);

        double initial_energy = 0.0;

        for (int i = 0; i < steps; ++i) {
            // Reconstruct velocity at integer step: v_i = v_{i-1/2} + 0.5 * a_i * dt
            double vi_x = vh_x + 0.5 * ax * dt;
            double vi_y = vh_y + 0.5 * ay * dt;
            double vi_z = vh_z + 0.5 * az * dt;
            double speed = std::sqrt(vi_x*vi_x + vi_y*vi_y + vi_z*vi_z);

            double pot = interpolate3D(phi, px, py, pz);
            double ke = 0.5 * mass * speed * speed;
            double pe = charge_q * pot;
            double total_energy = ke + pe;

            if (i == 0) {
                initial_energy = total_energy;
            }

            double drift_pct = (initial_energy != 0.0) ? (std::abs(total_energy - initial_energy) / std::abs(initial_energy) * 100.0) : 0.0;

            double dx1 = px - x1;
            double dy1 = py - y1;
            double dz1 = pz - z1;
            double dx2 = px - x2;
            double dy2 = py - y2;
            double dz2 = pz - z2;
            double dist_pos = std::sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
            double dist_neg = std::sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);

            trajectory.push_back({
                i * dt,
                px, py, pz,
                vi_x, vi_y, vi_z,
                speed,
                total_energy,
                drift_pct,
                dist_pos,
                dist_neg
            });

            // Update velocity (kick): v_{i+1/2} = v_{i-1/2} + a_i * dt
            vh_x += ax * dt;
            vh_y += ay * dt;
            vh_z += az * dt;

            // Update position (drift): x_{i+1} = x_i + v_{i+1/2} * dt
            px += vh_x * dt;
            py += vh_y * dt;
            pz += vh_z * dt;

            // Reflective Boundaries (perfect elastic mirror collisions)
            // X-axis
            if (px < 1.0) {
                px = 2.0 - px;
                vh_x = std::abs(vh_x);
            } else if (px > NX_field - 1) {
                px = 2.0 * (NX_field - 1) - px;
                vh_x = -std::abs(vh_x);
            }
            // Y-axis
            if (py < 1.0) {
                py = 2.0 - py;
                vh_y = std::abs(vh_y);
            } else if (py > NY_field - 1) {
                py = 2.0 * (NY_field - 1) - py;
                vh_y = -std::abs(vh_y);
            }
            // Z-axis
            if (pz < 1.0) {
                pz = 2.0 - pz;
                vh_z = std::abs(vh_z);
            } else if (pz > NZ_field - 1) {
                pz = 2.0 * (NZ_field - 1) - pz;
                vh_z = -std::abs(vh_z);
            }

            // Recalculate acceleration at new position
            ex = interpolate3D(Ex, px, py, pz);
            ey = interpolate3D(Ey, px, py, pz);
            ez = interpolate3D(Ez, px, py, pz);

            ax = (charge_q / mass) * ex;
            ay = (charge_q / mass) * ey;
            az = (charge_q / mass) * ez;

            if ((i + 1) % 20000 == 0) {
                std::cout << "  Step " << std::setw(6) << (i + 1)
                          << " | Drift: " << std::setw(7) << std::fixed << std::setprecision(4) << drift_pct << "%"
                          << " | Acceleration: " << std::setprecision(3) << std::sqrt(ax*ax + ay*ay + az*az) << std::endl;
            }
        }

        std::cout << "  Initial Energy: " << initial_energy << std::endl;
        std::cout << "  Final Energy: " << trajectory.back().energy << std::endl;
        std::cout << "  Final Drift: " << trajectory.back().drift_pct << "%" << std::endl;

        // ==========================================
        // 5. Diagnostics & Signal Processing
        // ==========================================
        std::cout << "[5/6] Extracting Diagnostics (Poincaré Crossings & Welch PSD)..." << std::endl;

        // Poincaré Section crossings at Y = 30
        std::vector<double> crossings_x;
        std::vector<double> crossings_vx;
        for (size_t i = 1; i < trajectory.size(); ++i) {
            double y_prev = trajectory[i - 1].pos_y;
            double y_curr = trajectory[i].pos_y;
            if ((y_prev < y1 && y1 <= y_curr) || (y_curr < y1 && y1 <= y_prev)) {
                double frac = (y1 - y_prev) / (y_curr - y_prev + 1e-15);
                double xc = trajectory[i - 1].pos_x + frac * (trajectory[i].pos_x - trajectory[i - 1].pos_x);
                double vxc = trajectory[i - 1].vel_x + frac * (trajectory[i].vel_x - trajectory[i - 1].vel_x);
                crossings_x.push_back(xc);
                crossings_vx.push_back(vxc);
            }
        }
        std::cout << "  Poincaré Section crossings recorded: " << crossings_x.size() << std::endl;

        // Welch PSD of test particle speed
        std::vector<double> speeds(steps);
        for (int i = 0; i < steps; ++i) {
            speeds[i] = trajectory[i].speed;
        }
        WelchResult psd_result = compute_welch(speeds, 1.0 / dt, 2048);
        std::cout << "  Welch Power Spectral Density computed." << std::endl;

        // ==========================================
        // 6. Report Dashboard Generation (ROOT Canvas)
        // ==========================================
        std::cout << "[6/6] Rendering Diagnostics to dipole_diagnostics.pdf..." << std::endl;

        TCanvas* c = new TCanvas("c_dipole", "Dipole Simulation Diagnostics Dashboard", 1800, 1200);
        c->Divide(3, 3);

        // Prep data arrays for standard ROOT graphs
        std::vector<double> t_arr(steps), px_arr(steps), py_arr(steps), pz_arr(steps);
        std::vector<double> vx_arr(steps), vy_arr(steps), vz_arr(steps);
        std::vector<double> e_arr(steps), drift_arr(steps), r_pos_arr(steps), r_neg_arr(steps);
        std::vector<double> path_arr(steps);

        double path_sum = 0.0;
        for (int i = 0; i < steps; ++i) {
            t_arr[i] = trajectory[i].t_s;
            px_arr[i] = trajectory[i].pos_x;
            py_arr[i] = trajectory[i].pos_y;
            pz_arr[i] = trajectory[i].pos_z;
            vx_arr[i] = trajectory[i].vel_x;
            vy_arr[i] = trajectory[i].vel_y;
            vz_arr[i] = trajectory[i].vel_z;
            e_arr[i] = trajectory[i].energy;
            drift_arr[i] = trajectory[i].drift_pct;
            r_pos_arr[i] = trajectory[i].dist_pos;
            r_neg_arr[i] = trajectory[i].dist_neg;

            if (i > 0) {
                double ds = std::sqrt(std::pow(px_arr[i] - px_arr[i - 1], 2) +
                                      std::pow(py_arr[i] - py_arr[i - 1], 2) +
                                      std::pow(pz_arr[i] - pz_arr[i - 1], 2));
                path_sum += ds;
            }
            path_arr[i] = path_sum;
        }

        // Pad 1: 3D Trajectory
        c->cd(1);
        TGraph2D* gr_3d = new TGraph2D(steps, px_arr.data(), py_arr.data(), pz_arr.data());
        gr_3d->SetTitle("3D Trajectory Path;X;Y;Z");
        gr_3d->SetLineColor(kBlack);
        gr_3d->SetLineWidth(1);
        gr_3d->Draw("LINE");

        // Add markers for Dipole Source Charges
        TGraph2D* gr_sources = new TGraph2D(2);
        gr_sources->SetPoint(0, x1, y1, z1);
        gr_sources->SetPoint(1, x2, y2, z2);
        gr_sources->SetMarkerStyle(29); // Star symbol
        gr_sources->SetMarkerSize(2.5);
        gr_sources->SetMarkerColor(kRed + 1); // Positive is red star
        gr_sources->Draw("P SAME");

        // Pad 2: Total Energy vs Time
        c->cd(2);
        TGraph* gr_energy = new TGraph(steps, t_arr.data(), e_arr.data());
        gr_energy->SetTitle("Total Energy vs Time;Time (s);Energy");
        gr_energy->SetLineColor(kAzure + 2);
        gr_energy->Draw("AL");

        // Pad 3: Energy Drift % vs Time (Log scale)
        c->cd(3);
        TGraph* gr_drift = new TGraph(steps, t_arr.data(), drift_arr.data());
        gr_drift->SetTitle("Energy Drift Percentage;Time (s);|#DeltaE/E_{0}| (%)");
        gr_drift->SetLineColor(kOrange + 7);
        gPad->SetLogy(1);
        gr_drift->Draw("AL");

        // Pad 4: Potential Phi Slice at Z = 30 + XY Path
        c->cd(4);
        TH2F* slice_phi = new TH2F("slice_phi", "Potential Slice (Z=30) & XY Path;X;Y",
                                   NX_field, 0, NX_field, NY_field, 0, NY_field);
        for (int i = 0; i < NX_field; ++i) {
            for (int j = 0; j < NY_field; ++j) {
                int z_idx = 30;
                double val = phi[i + NX_field * (j + NY_field * z_idx)];
                slice_phi->SetBinContent(i + 1, j + 1, val);
            }
        }
        slice_phi->SetMinimum(slice_phi->GetMinimum());
        slice_phi->SetMaximum(slice_phi->GetMaximum());
        slice_phi->Draw("COLZ");

        TGraph* path_xy = new TGraph(steps, px_arr.data(), py_arr.data());
        path_xy->SetLineColor(kGreen + 2);
        path_xy->SetLineWidth(1);
        path_xy->Draw("L SAME");

        // Pad 5: Phase Portrait (X vs VX)
        c->cd(5);
        TGraph* phase_x = new TGraph(steps, px_arr.data(), vx_arr.data());
        phase_x->SetTitle("X-axis Phase Portrait;X (position);Vx (velocity)");
        phase_x->SetMarkerStyle(1);
        phase_x->SetMarkerColor(kAzure + 7);
        phase_x->Draw("AP");

        // Pad 6: Radial Distribution Histograms
        c->cd(6);
        TH1D* hist_pos = new TH1D("hist_pos", "Radial Distance Distribution;Distance (grid units);Probability", 80, 0.0, 45.0);
        TH1D* hist_neg = new TH1D("hist_neg", "Radial Distance Distribution;Distance (grid units);Probability", 80, 0.0, 45.0);
        for (int i = 0; i < steps; ++i) {
            hist_pos->Fill(r_pos_arr[i]);
            hist_neg->Fill(r_neg_arr[i]);
        }
        normalize_hist(hist_pos);
        normalize_hist(hist_neg);

        hist_pos->SetLineColor(kRed + 1);
        hist_pos->SetLineWidth(2);
        hist_pos->Draw("HIST");

        hist_neg->SetLineColor(kAzure + 2);
        hist_neg->SetLineWidth(2);
        hist_neg->Draw("HIST SAME");

        // Pad 7: Cumulative Path Length vs Step
        c->cd(7);
        TGraph* gr_path = new TGraph(steps, t_arr.data(), path_arr.data());
        gr_path->SetTitle("Cumulative Path Length vs Time;Time (s);Cumulative Distance");
        gr_path->SetLineColor(kViolet + 1);
        gr_path->Draw("AL");

        // Pad 8: Poincaré Section at Y = 30
        c->cd(8);
        if (!crossings_x.empty()) {
            TGraph* poincare = new TGraph(crossings_x.size(), crossings_x.data(), crossings_vx.data());
            poincare->SetTitle("Poincar#acute{e} Section (Y=30 crossings);X;Vx");
            poincare->SetMarkerStyle(8);
            poincare->SetMarkerSize(0.6);
            poincare->SetMarkerColor(kBlue - 4);
            poincare->Draw("AP");
        } else {
            TH2F* dummy = new TH2F("dummy", "Poincar#acute{e} Section;X;Vx", 10, 0, NX_field, 10, -2, 2);
            dummy->Draw();
        }

        // Pad 9: Welch Speed PSD
        c->cd(9);
        if (!psd_result.psd.empty()) {
            TGraph* gr_psd = new TGraph(psd_result.psd.size(), psd_result.freqs.data(), psd_result.psd.data());
            gr_psd->SetTitle("Speed Power Spectral Density (Welch);Frequency (Hz);PSD");
            gr_psd->SetLineColor(kTeal + 2);
            gPad->SetLogy(1);
            gr_psd->Draw("AL");
        }

        c->SaveAs("outputs/dipole_diagnostics.pdf");
        c->SaveAs("outputs/dipole_diagnostics.png");

        // Cleanup allocated objects
        delete grid_2d;
        delete grid_3d;
        delete gr_3d;
        delete gr_sources;
        delete gr_energy;
        delete gr_drift;
        delete slice_phi;
        delete path_xy;
        delete phase_x;
        delete hist_pos;
        delete hist_neg;
        delete gr_path;
        delete c;

        std::cout << "\n=== Simulation Finished Successfully! ===" << std::endl;
        std::cout << "Saved: outputs/cic_charge_deposition.png" << std::endl;
        std::cout << "Saved: outputs/cic_charge_3d.png" << std::endl;
        std::cout << "Saved: outputs/dipole_diagnostics.pdf" << std::endl;
        std::cout << "Saved: outputs/dipole_diagnostics.png" << std::endl;
    }

} // namespace DipoleSimulation
