#include "ft.hpp"

#include <cmath>
#include <complex>
#include <numbers>
#include <vector>

using fourier::ComplexVector;

namespace {
    /**
     * @brief Internal FFT function which assumes input size is a power of two.
     */
    ComplexVector fft_radix2(const ComplexVector& input) {
        size_t N = input.size();
        if (N <= 1) return input;

        // Split into even and odd components
        ComplexVector even(N / 2);
        ComplexVector odd(N / 2);
        for (size_t i = 0; i < N / 2; i++) {
            even[i] = input[i * 2];
            odd[i] = input[i * 2 + 1];
        }

        // Recursive calls
        ComplexVector y_even = fft_radix2(even);
        ComplexVector y_odd = fft_radix2(odd);

        // Perform butterfly operation
        ComplexVector y(N);
        for (size_t k = 0; k < N / 2; k++) {
            double angle = -2.0 * std::numbers::pi * k / N;

            // Euler: e^{i*angle} = cos(angle) + i*sin(angle)
            std::complex<double> w(std::cos(angle), std::sin(angle));

            y[k] = y_even[k] + w * y_odd[k];
            y[k + N / 2] = y_even[k] - w * y_odd[k];
        }

        return y;
    }
}

void fourier::apply_window(ComplexVector &signal, WindowType type) {
    size_t N = signal.size();
    if (N == 0) return;

    // Rectangular window is equivalent to multiplication by 1
    if (type == WindowType::Rectangular) return;

    // Compute windowed values
    for (size_t n = 0; n < N; n++) {
        double window_value = 1.0;
        double angle = 2.0 * std::numbers::pi * n / N;

        switch (type) {
            case WindowType::Hann:
                // Hann window: 0.5*(1 - cos(2πn/N))
                window_value = 0.5 - 0.5 * std::cos(angle);
                break;

            case WindowType::Hamming:
                // Hamming window: 0.54 - 0.46*cos(2πn/N)
                window_value = 0.54 - 0.46 * std::cos(angle);
                break;

            case WindowType::Blackman:
                // Blackman window: 0.42 - 0.5*cos(2πn/N) + 0.08*cos(4πn/N)
                window_value = 0.42 - 0.5 * std::cos(angle) + 0.08 * std::cos(2 * angle);
                break;

            default:
                break;
        }

        signal[n] *= window_value;
    }
}

ComplexVector fourier::dft(const ComplexVector &input) {
    size_t N = input.size();
    if (N == 0) return {};

    // Compute output values
    ComplexVector output(N, std::complex<double>(0.0, 0.0));
    for (size_t k = 0; k < N; k++) {
        for (size_t n = 0; n < N; n++) {
            double angle = -2.0 * std::numbers::pi * k * n / N;

            // Euler: e^{i*angle} = cos(angle) + i*sin(angle)
            std::complex<double> exp_term(std::cos(angle), std::sin(angle));

            output[k] += input[n] * exp_term;
        }
    }

    return output;
}

ComplexVector fourier::fft(const ComplexVector& input) {
    size_t N = input.size();
    if (N == 0) return {};

    // Compute next power of two
    size_t power_of_two = 1;
    while (power_of_two < N) {
        power_of_two *= 2;
    }

    // Zero-pad if necessary
    ComplexVector padded_input = input;
    if (N != power_of_two)
        padded_input.resize(power_of_two, std::complex<double>(0.0, 0.0));

    return fft_radix2(padded_input);
}

std::vector<size_t> fourier::ca_cfar(const std::vector<double>& magnitudes, size_t num_train, size_t num_guard, double threshold_factor) {
    std::vector<size_t> detections;

    // Check if window is out of bounds
    size_t N = magnitudes.size();
    size_t side_window = num_train + num_guard;
    if (N <= 2 * side_window) return detections;

    // Slide CUT across spectrum
    for (size_t i = side_window; i < N - side_window; i++) {
        double noise_sum = 0.0;

        // Add up noise in training cells to the left
        for (size_t j = i - side_window; j < i - num_guard; j++) {
            noise_sum += magnitudes[j];
        }

        // Add up noise in training cells to the right
        for (size_t j = i + num_guard + 1; j <= i + side_window; j++) {
            noise_sum += magnitudes[j];
        }

        // Check if magnitude is above threshold
        double noise_avg = noise_sum / (2.0 * num_train);
        double threshold = noise_avg * threshold_factor;
        if (magnitudes[i] > threshold) detections.push_back(i);
    }

    return detections;
}
