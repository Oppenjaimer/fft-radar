#include "ft.hpp"

#include <cmath>
#include <complex>
#include <numbers>

using fourier::ComplexVector;

ComplexVector fourier::dft(const ComplexVector &input) {
    size_t N = input.size();

    // Return zero-initialized output if input is empty to avoid division by zero
    ComplexVector output(N, std::complex<double>(0.0, 0.0));
    if (N == 0) return output;

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
