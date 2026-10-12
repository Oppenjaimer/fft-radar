#include "ft.hpp"

#include <cmath>
#include <complex>
#include <numbers>

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
