#ifndef AUDIOPROCESSOR_H
#define AUDIOPROCESSOR_H

#include <vector>
#include <string>

class AudioProcessor {
public:
    static std::vector<double> readWav(const std::string& filename, int& sample_rate_hz);
    static void writeWav(const std::string& filename, const std::vector<double>& pcm_samples, int sample_rate_hz);
    
    // Uses ROOT's TVirtualFFT
    static std::vector<double> filter(const std::vector<double>& pcm_samples, int sample_rate_hz, double freq_lo_hz, double freq_hi_hz);
};

#endif
