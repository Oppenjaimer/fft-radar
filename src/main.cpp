#include "ft.hpp"

#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numbers>
#include <random>
#include <vector>

/**
 * @brief Simulated radar target.
 */
struct Target {
    double frequency;   ///< Frequency bin, proportional to the target's radial velocity.
    double amplitude;   ///< Signal strength, reflecting target size and distance.
    double phase;       ///< Initial phase offset of the wave.
};

/**
 * @brief Program entry point.
 * @return Exit code.
 */
int main() {
    constexpr size_t N = 1024;
    fourier::ComplexVector input(N);

    // Define targets
    std::vector<Target> targets = {
        {5.0,   1.2, 0.0},
        {25.0,  0.5, std::numbers::pi / 4.0},
        {50.5,  0.8, std::numbers::pi / 2.0},
        {300.0, 1.0, 1.27}
    };

    // Initialize random noise distribution
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> noise_dist(0.0, 1.0);

    // Generate composite signal from targets
    for (size_t n = 0; n < N; n++) {
        double signal_value = 0.0;

        for (const auto& target : targets) {
            double angle = 2.0 * std::numbers::pi * target.frequency * n / N + target.phase;
            signal_value += target.amplitude * std::cos(angle);
        }

        signal_value += noise_dist(gen);
        input[n] = std::complex<double>(signal_value, 0.0);
    }

    // Run FFT
    fourier::ComplexVector processed_input = input;
    fourier::apply_window(processed_input, fourier::WindowType::Hann);
    fourier::ComplexVector output = fourier::fft(processed_input);

    // Open output CSV file
    std::ofstream file("data/data.csv");
    if (!file.is_open()) {
        std::cerr << "Unable to open CSV file\n";
        return 1;
    }

    file << "index,time_ampl,freq_ampl\n";

    // Export data to CSV
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
