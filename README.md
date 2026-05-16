# SchrödingerLib: High-Performance Numerical Analysis & Physical Simulation

**SchrödingerLib** is a comprehensive C++23 toolkit designed for high-fidelity physical simulations, signal processing, and numerical modeling. Built upon the **CERN ROOT Data Analysis Framework**, the project leverages industrial-grade libraries to deliver precise analytical results and professional telemetry visualizations.

---

## 🛠 Core Modules

### 1. Aerospace Integration (ISRO LVM3 Simulator)
A high-fidelity trajectory integrator utilizing the **Leapfrog (Midpoint) Method** for superior energy conservation over long-duration simulations.
*   **Decoupled Architecture**: Implements a **Dependency Injection** pattern, separating planetary physics (`IEnvironment`) from vehicle engineering (`IRocket`).
*   **Atmospheric Model**: Integrates the **U.S. Standard Atmosphere (1976)** to compute geopotential altitude, dynamic pressure, and ambient density across 7 layers.
*   **Dynamic Staging**: Simulates the ISRO LVM3 flight profile with sequential fuel-based staging (S200 -> L110 -> C25) and automated Payload Fairing (PLF) jettison at 115km.
*   **Aerodynamics**: Features a Mach-dependent $C_d$ drag model, accounting for the transonic wave-drag peak and dynamic reference area changes.
*   **Telemetry**: Generates a professional 3x3 telemetry grid (PDF) visualizing 9 synchronized parameters including **Max-Q (Dynamic Pressure)** and Mach transitions.

### 2. Signal Processing (FFT Toolkit)
A robust frequency-domain filtering module for PCM audio data.
*   **Algorithm**: Utilizes the **Cooley-Tukey FFT** (implemented via ROOT's `TVirtualFFT`) to perform forward and inverse transforms.
*   **Filtering**: Implements a surgical band-stop filter to eliminate user-defined frequency ranges with amplitude normalization.

### 3. Thermodynamic Analytics (Heat Diffusion)
Solves the spatial-temporal state of the heat equation using analytical Green's functions.
*   **Visualization**: Evaluates Gaussian point source diffusion in 1D and 3D space, exporting high-frame-rate animated GIFs of the thermal evolution.

### 4. Mathematical Utilities
*   **Root Finding**: Implements a robust multi-interval **Bisection Scanner** combined with **Cauchy's Bound** analysis to isolate polynomial roots.

---

## 🚀 Technical Standards

*   **Language**: C++23 (ISO Standard).
*   **Frameworks**: CERN ROOT (Primary), Qt6 (GUI).
*   **Memory Management**: RAII principles with `std::unique_ptr` for dependency injection and explicit ROOT object lifecycle management.
*   **Convention**: Strict `variable_name_unit` naming convention (e.g., `mass_kg`, `altitude_m`, `sample_rate_hz`) to ensure physical dimensional consistency.

---

## 📦 Build & Execution

### Prerequisites
*   **CERN ROOT**: Ensure `ROOTSYS` is configured in your environment.
*   **CMake**: Version 3.15 or higher.

### Compilation
```bash
cd numerical_analysis
mkdir build && cd build
cmake ..
make -j$(nproc)
```

### Running the CLI
```bash
./runme
```

---

## 🗺 Roadmap
- [ ] **Sturm's Theorem**: Exact real root counting for the Polynomial module.
- [ ] **MarsEnvironment**: Extension of `IEnvironment` for extraterrestrial simulations.
- [ ] **Monte Carlo Analysis**: Sensitivity testing for launch window variations.

*Developed for engineers and physics researchers.*
