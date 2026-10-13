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
        {25.0,  0.5, std::numbers::pi / 4.0},
        {50.5,  0.8, std::numbers::pi / 2.0},
        {300.0, 1.0, 1.27},
        {430.2, 1.2, 0.0},
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

    // Extract magnitudes
    std::vector<double> magnitudes(N / 2);
    for (size_t k = 0; k < N / 2; k++) {
        double magnitude = std::abs(output[k]) / static_cast<double>(N);
        if (k > 0) magnitude *= 2;

        magnitudes[k] = magnitude;
    }

    // Run CA-CFAR
    constexpr size_t num_train = 15;
    constexpr size_t num_guard = 2;
    constexpr double threshold_factor = 3.5;
    std::vector<size_t> detections = fourier::ca_cfar(magnitudes, num_train, num_guard, threshold_factor);

    // Perform clustering to group contiguous detections
    std::vector<size_t> clustered_detections;
    if (!detections.empty()) {
        size_t current_peak_bin = detections[0];

        for (size_t i = 1; i < detections.size(); i++) {
            if (detections[i] == detections[i - 1] + 1) {
                // Contiguous detections, choose one with highest energy
                if (magnitudes[detections[i]] > magnitudes[current_peak_bin])
                    current_peak_bin = detections[i];
            } else {
                // No longer traversing previous cluster
                clustered_detections.push_back(current_peak_bin);
                current_peak_bin = detections[i];
            }
        }

        clustered_detections.push_back(current_peak_bin);
    }

    // Output CA-CFAR results
    if (clustered_detections.empty()) {
        std::cout << "No targets detected\n";
    } else {
        for (size_t bin : clustered_detections) {
            std::cout << "Target found at bin " << bin << " with magnitude " << magnitudes[bin] << "\n";
        }
    }

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
