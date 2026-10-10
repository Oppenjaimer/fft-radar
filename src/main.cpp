#include "ft.hpp"

#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <numbers>

int main() {
    constexpr size_t N = 8;
    fourier::ComplexVector input(N);

    // Generate 1 Hz cosine wave
    for (size_t n = 0; n < N; n++) {
        double angle = -2.0 * std::numbers::pi * 1 * n / N;
        input[n] = std::complex<double>(std::cos(angle), 0.0);

        std::cout << "Sample " << n << ": "
                   << std::fixed << std::setprecision(3) << input[n].real() << " + "
                   << input[n].imag() << "i\n";
    }

    // Run DFT
    fourier::ComplexVector output = fourier::dft(input);
    std::cout << "\n";

    for (size_t k = 0; k < N; ++k) {
        // Multiply by 2 because input is real
        double norm_real = output[k].real() * 2.0 / static_cast<double>(N);
        double norm_imag = output[k].imag() * 2.0 / static_cast<double>(N);
        double magnitude = std::abs(output[k]) * 2.0 / static_cast<double>(N);

        std::cout << "Bin " << k << ": |"
                  << std::setprecision(3) << norm_real << " + "
                  << std::setprecision(3) << norm_imag << "i| = "
                  << magnitude << "\n";
    }

    return 0;
}
