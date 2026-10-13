#pragma once

#include <complex>
#include <vector>

namespace fourier {
    using ComplexVector = std::vector<std::complex<double>>;    ///< Vector of double-precision complex numbers.

    /**
     * @brief Supported windowing functions to mitigate spectral leakage.
     */
    enum class WindowType {
        Rectangular,
        Hann,
        Hamming,
        Blackman
    };

    /**
     * @brief Apply a windowing function to a time-domain signal.
     *
     * Modifies the input signal by multiplying it with the chosen window curve. This prepares the signal for FFT by
     * smoothing the boundaries, reducing spectral leakage at the cost of slight frequency resolution loss.
     *
     * @param signal Time-domain signal array to be windowed.
     * @param type Windowing function type.
     */
    void apply_window(ComplexVector& signal, WindowType type);

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
