#include "AudioProcessor.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <TVirtualFFT.h>

std::vector<double> AudioProcessor::readWav(const std::string& filename, int& sampleRate) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) return {};

    char header[44];
    file.read(header, 44);

    sampleRate = *reinterpret_cast<int*>(&header[24]);
    int dataSize = *reinterpret_cast<int*>(&header[40]);
    int numSamples = dataSize / sizeof(short);

    std::vector<double> samples;
    samples.reserve(numSamples);

    short sample;
    while (file.read(reinterpret_cast<char*>(&sample), sizeof(short))) {
        samples.push_back(static_cast<double>(sample));
    }

    return samples;
}

void AudioProcessor::writeWav(const std::string& filename, const std::vector<double>& samples, int sampleRate) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) return;

    int numSamples = samples.size();
    int dataSize = numSamples * sizeof(short);
    int chunkSize = 36 + dataSize;
    short audioFmt = 1;
    short numChan = 1;
    int byteRate = sampleRate * numChan * sizeof(short);
    short blockAlign = numChan * sizeof(short);
    short bitsPerSample = 16;

    file.write("RIFF", 4);
    file.write(reinterpret_cast<char*>(&chunkSize), 4);
    file.write("WAVE", 4);
    file.write("fmt ", 4);
    int subchunk1Size = 16;
    file.write(reinterpret_cast<char*>(&subchunk1Size), 4);
    file.write(reinterpret_cast<char*>(&audioFmt), 2);
    file.write(reinterpret_cast<char*>(&numChan), 2);
    file.write(reinterpret_cast<char*>(&sampleRate), 4);
    file.write(reinterpret_cast<char*>(&byteRate), 4);
    file.write(reinterpret_cast<char*>(&blockAlign), 2);
    file.write(reinterpret_cast<char*>(&bitsPerSample), 2);
    file.write("data", 4);
    file.write(reinterpret_cast<char*>(&dataSize), 4);

    for (double s : samples) {
        short val = static_cast<short>(std::clamp(s, -32768.0, 32767.0));
        file.write(reinterpret_cast<char*>(&val), sizeof(short));
    }
}

std::vector<double> AudioProcessor::filter(const std::vector<double>& samples, int sampleRate, double freq_lo, double freq_hi) {
    int N = samples.size();
    if (N <= 0) return {};

    // 1. Forward FFT using ROOT
    TVirtualFFT *fft = TVirtualFFT::FFT(1, &N, "R2C M K");
    fft->SetPoints(samples.data());
    fft->Transform();

    std::vector<double> re_arr(N/2 + 1), im_arr(N/2 + 1);
    for (int i = 0; i <= N/2; ++i) {
        fft->GetPointComplex(i, re_arr[i], im_arr[i]);
    }

    // 2. Filter
    for (int i = 0; i <= N/2; ++i) {
        double freq = static_cast<double>(i) * sampleRate / N;
        if (freq >= freq_lo && freq <= freq_hi) {
            re_arr[i] = 0.0;
            im_arr[i] = 0.0;
        }
    }

    // 3. Inverse FFT using ROOT
    TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2R M K");
    fft_back->SetPointsComplex(re_arr.data(), im_arr.data());
    fft_back->Transform();

    std::vector<double> result(N);
    for (int i = 0; i < N; ++i) {
        result[i] = fft_back->GetPointReal(i) / N; // FFTW output needs /N normalisation
    }

    // Normalisation to preserve loudness
    double peak = 0;
    for (double v : result) if (std::abs(v) > peak) peak = std::abs(v);
    double origPeak = 0;
    for (double v : samples) if (std::abs(v) > origPeak) origPeak = std::abs(v);
    double scale = (peak > 0) ? origPeak / peak : 1.0;
    
    for(int i=0; i<N; i++){
        result[i] *= scale;
    }

    delete fft;
    delete fft_back;

    return result;
}
