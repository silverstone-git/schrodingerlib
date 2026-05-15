# 🚀 SchrödingerLib: Physics, Rockets, and Sound! 

Welcome to **SchrödingerLib**! If you've ever wanted to launch an Indian LVM3 rocket, visualize how heat travels through a room, or use math to "delete" annoying frequencies from a song, you're in the right place. 

Think of this as a "super-powered calculator" built with **C++** (the fast language used to make games like Fortnite) and **ROOT** (the framework used by scientists at CERN to find the Higgs Boson).

---

## 🌟 What's inside?

### 🚀 1. Rocket Science (LVM3 Simulator)
Ever wondered how a rocket actually gets to space? This isn't just a simple "up" button. Our simulator uses the **exact** specs for the ISRO LVM3 rocket.
*   **Staging**: It knows when to drop the heavy empty boosters to go faster.
*   **Atmosphere**: It calculates how the air gets thinner as you go higher using the *1976 U.S. Standard Atmosphere* model.
*   **Telemetry**: It creates cool graphs showing your altitude, speed, and acceleration!

### 🔊 2. The Audio "Eraser" (FFT Filter)
You can give this tool a `.wav` file (like a song) and tell it to delete a specific frequency range.
*   It uses something called a **Fast Fourier Transform (FFT)** to turn sound waves into a "frequency map." 
*   It's like having an EQ slider that can perfectly cut out a specific sound!

### 🌡️ 3. Heat Evolution
Watch how heat spreads out! You can generate **animated GIFs** that show a "heat pulse" smoothing out over time in 1D (a wire) or 3D (a block). It's basically math turned into a movie.

### 📐 4. Polynomial Root Finder
Give it a long equation like $x^2 - 2x - 3 = 0$, and it will find exactly where the graph hits zero. It uses the **Bisection Method**, which is like a high-speed version of the "High-Low" guessing game.

---

## 🛠️ Getting Started (For the Devs)

### 🐍 "I know Python, how do I run this C++ stuff?"
Don't worry! C++ looks a bit more complex, but we've set everything up so it's easy to build. You just need to have **ROOT** installed.

#### Step 1: Go to the folder
```bash
cd numerical_analysis
```

#### Step 2: Build it (Turn the code into an app)
```bash
mkdir -p build && cd build
cmake ..
make
```

#### Step 3: Run the Magic!
```bash
./runme
```

---

## 🎮 The "Ball on a Wedge" Game
If you want to see a physical simulation with a real window and sliders, we have a **Qt6** app in the root directory. You can change gravity, friction, and the angle of a wedge to see how a ball rolls down!

**To run the GUI app:**
```bash
cmake --preset default
cmake --build build
./build/schrodingerlib  # (Or whatever your OS calls the exe!)
```

---

## 📚 Technical Stuff (The "Nerd" Corner)
*   **Language**: C++23 (The latest and greatest).
*   **Frameworks**: [CERN ROOT](https://root.cern/) and [Qt6](https://www.qt.io/).
*   **Algorithms**: Leapfrog Integration, Cooley-Tukey FFT, Bisection Method.

*Made with ❤️ for physics lovers everywhere.*
