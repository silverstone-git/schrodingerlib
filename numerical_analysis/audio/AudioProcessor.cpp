#include "AudioProcessor.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <TVirtualFFT.h>

std::vector<double> AudioProcessor::readWav(const std::string& filename, int& sample_rate_hz) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) return {};

    char header[44];
    file.read(header, 44);

    sample_rate_hz = *reinterpret_cast<int*>(&header[24]);
    int data_size_bytes = *reinterpret_cast<int*>(&header[40]);
    int num_samples = data_size_bytes / sizeof(short);

    std::vector<double> pcm_samples;
    pcm_samples.reserve(num_samples);

    short raw_sample;
    while (file.read(reinterpret_cast<char*>(&raw_sample), sizeof(short))) {
        pcm_samples.push_back(static_cast<double>(raw_sample));
    }

    return pcm_samples;
}

void AudioProcessor::writeWav(const std::string& filename, const std::vector<double>& pcm_samples, int sample_rate_hz) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) return;

    int num_samples = pcm_samples.size();
    int data_size_bytes = num_samples * sizeof(short);
    int chunk_size_bytes = 36 + data_size_bytes;
    short audio_format = 1; // PCM
    short num_channels = 1; // Mono
    int byte_rate = sample_rate_hz * num_channels * sizeof(short);
    short block_align = num_channels * sizeof(short);
    short bits_per_sample = 16;

    file.write("RIFF", 4);
    file.write(reinterpret_cast<char*>(&chunk_size_bytes), 4);
    file.write("WAVE", 4);
    file.write("fmt ", 4);
    int subchunk1_size = 16;
    file.write(reinterpret_cast<char*>(&subchunk1_size), 4);
    file.write(reinterpret_cast<char*>(&audio_format), 2);
    file.write(reinterpret_cast<char*>(&num_channels), 2);
    file.write(reinterpret_cast<char*>(&sample_rate_hz), 4);
    file.write(reinterpret_cast<char*>(&byte_rate), 4);
    file.write(reinterpret_cast<char*>(&block_align), 2);
    file.write(reinterpret_cast<char*>(&bits_per_sample), 2);
    file.write("data", 4);
    file.write(reinterpret_cast<char*>(&data_size_bytes), 4);

    for (double sample_val : pcm_samples) {
        short short_val = static_cast<short>(std::clamp(sample_val, -32768.0, 32767.0));
        file.write(reinterpret_cast<char*>(&short_val), sizeof(short));
    }
}

std::vector<double> AudioProcessor::filter(const std::vector<double>& pcm_samples, int sample_rate_hz, double freq_lo_hz, double freq_hi_hz) {
    int n_samples = pcm_samples.size();
    if (n_samples <= 0) return {};

    // 1. Forward FFT using ROOT
    TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &n_samples, "R2C M K");
    fft_forward->SetPoints(pcm_samples.data());
    fft_forward->Transform();

    std::vector<double> re_arr(n_samples/2 + 1), im_arr(n_samples/2 + 1);
    for (int i = 0; i <= n_samples/2; ++i) {
        fft_forward->GetPointComplex(i, re_arr[i], im_arr[i]);
    }

    // 2. Filter
    for (int i = 0; i <= n_samples/2; ++i) {
        double freq_hz = static_cast<double>(i) * sample_rate_hz / n_samples;
        if (freq_hz >= freq_lo_hz && freq_hz <= freq_hi_hz) {
            re_arr[i] = 0.0;
            im_arr[i] = 0.0;
        }
    }

    // 3. Inverse FFT using ROOT
    TVirtualFFT *fft_backward = TVirtualFFT::FFT(1, &n_samples, "C2R M K");
    fft_backward->SetPointsComplex(re_arr.data(), im_arr.data());
    fft_backward->Transform();

    std::vector<double> filtered_pcm(n_samples);
    for (int i = 0; i < n_samples; ++i) {
        filtered_pcm[i] = fft_backward->GetPointReal(i) / n_samples; 
    }

    // Normalisation to preserve loudness
    double peak_val = 0;
    for (double v : filtered_pcm) if (std::abs(v) > peak_val) peak_val = std::abs(v);
    double orig_peak_val = 0;
    for (double v : pcm_samples) if (std::abs(v) > orig_peak_val) orig_peak_val = std::abs(v);
    double scale_factor = (peak_val > 0) ? orig_peak_val / peak_val : 1.0;
    
    for(int i=0; i<n_samples; i++){
        filtered_pcm[i] *= scale_factor;
    }

    delete fft_forward;
    delete fft_backward;

    return filtered_pcm;
}
