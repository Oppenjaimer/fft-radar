#pragma once

#include <complex>
#include <vector>

namespace fourier {
    using ComplexVector = std::vector<std::complex<double>>;    ///< Vector of double-precision complex numbers.

    /**
     * @brief Compute the discrete Fourier transform (DFT) of a complex signal.
     *
     * The returned frequency bins are the raw, unnormalized sums. To obtain true physical amplitudes, the magnitude
     * of each bin must be divided by N during post-processing. Standard convention omits this scaling in the forward
     * pass to preserve mathematical reversibility with the inverse DFT.
     *
     * @param input Time-domain signal array of N complex samples.
     * @return ComplexVector of N frequency bins representing the spectrum.
     */
    ComplexVector dft(const ComplexVector& input);

    /**
     * @brief Compute the fast Fourier transform (FFT) of a complex signal.
     *
     * Automatically pads the input with zeros to the next power of two if necessary.
     *
     * The returned frequency bins are the raw, unnormalized sums. To obtain true physical amplitudes, the magnitude
     * of each bin must be divided by N (padded size) during post-processing. Standard convention omits this scaling
     * in the forward pass to preserve mathematical reversibility with the inverse DFT.
     *
     * @param input Time-domain signal array of complex samples.
     * @return ComplexVector of N frequency bins representing the spectrum.
     */
    ComplexVector fft(const ComplexVector& input);
};
