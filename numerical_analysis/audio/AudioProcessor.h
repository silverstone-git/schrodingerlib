#ifndef AUDIOPROCESSOR_H
#define AUDIOPROCESSOR_H

#include <vector>
#include <string>

class AudioProcessor {
public:
    static std::vector<double> readWav(const std::string& filename, int& sampleRate);
    static void writeWav(const std::string& filename, const std::vector<double>& samples, int sampleRate);
    
    // Uses ROOT's TVirtualFFT
    static std::vector<double> filter(const std::vector<double>& samples, int sampleRate, double freq_lo, double freq_hi);
};

#endif
