#include "ft.hpp"

#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numbers>

int main() {
    constexpr size_t N = 8;
    fourier::ComplexVector input(N);

    // Generate 1 Hz cosine wave
    for (size_t n = 0; n < N; n++) {
        double angle = -2.0 * std::numbers::pi * 1 * n / N;
        input[n] = std::complex<double>(std::cos(angle), 0.0);
    }

    // Run DFT
    fourier::ComplexVector output = fourier::dft(input);

    std::ofstream file("data/data.csv");
    if (!file.is_open()) {
        std::cerr << "Unable to open CSV file\n";
        return 1;
    }

    file << "index,time_ampl,freq_ampl\n";

    for (size_t i = 0; i < N; i++) {
        double time_amplitude = input[i].real();
        double freq_amplitude = std::abs(output[i]);

        file << i << ","
            << time_amplitude << ","
            << freq_amplitude << "\n";
    }

    file.close();

    return 0;
}
