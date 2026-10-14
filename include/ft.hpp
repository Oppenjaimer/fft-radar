#pragma once

#include <complex>
#include <vector>

namespace fourier {
    using Complex = std::complex<double>;       ///< Double-precision complex number.
    using ComplexVector = std::vector<Complex>; ///< Vector of double-precision complex numbers.

    /**
     * @brief Supported window functions to mitigate spectral leakage.
     */
    enum class WindowType {
        Rectangular,
        Hann,
        Hamming,
        Blackman
    };

    /**
     * @brief Configuration parameters for 1D CA-CFAR detection.
     */
    struct CfarConfig {
        size_t num_train = 15;          ///< Number of training cells on each side of the CUT.
        size_t num_guard = 2;           ///< Number of guard cells on each side of the CUT.
        double threshold_factor = 4.5;  ///< Multiplier applied to the estimated noise floor.
    };

    /**
     * @brief Output of the CA-CFAR detection.
     */
    struct CfarResult {
        std::vector<size_t> detected_bins;  ///< Indices of bins exceeding the dynamic threshold.
        std::vector<double> thresholds;     ///< Calculated threshold per bin (NaN where window is invalid).
    };

    /**
     * @brief Return the coherent amplitude gain of a window function.
     */
    constexpr double coherent_gain(WindowType type) {
        switch (type) {
            case WindowType::Hann:     return 0.5;
            case WindowType::Hamming:  return 0.54;
            case WindowType::Blackman: return 0.42;
            default:                   return 1.0;
        }
    }

    /**
     * @brief Apply a window function in-place to a time-domain signal.
     *
     * Modifies the input signal by multiplying it with the chosen window curve. This prepares the signal for FFT by
     * smoothing the boundaries, reducing spectral leakage at the cost of slight frequency resolution loss.
     *
     * @param signal Time-domain signal array to be windowed.
     * @param type   Window function type.
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

    /**
     * @brief Compute the normalized single-sided magnitude spectrum.
     *
     * This compensates for transform size N and window coherent gain for bins 0 to N/2-1.
     *
     * @param spectrum Frequency-domain signal array to be normalized.
     * @param type     Window function type.
     *
     * @return Vector of normalized magnitudes.
     */
    std::vector<double> normalize_magnitudes(const ComplexVector& spectrum, WindowType type);

    /**
     * @brief Apply 1D cell-averaging constant false alarm rate (CA-CFAR) detection.
     *
     * Scans a magnitude spectrum to automatically detect peaks above a dynamically calculated noise floor.
     *
     * @param magnitudes Frequency spectrum magnitudes.
     * @param config     Configuration parameters.
     *
     * @return CfarResult of the detection.
     */
    CfarResult ca_cfar(const std::vector<double>& magnitudes, const CfarConfig& config);

    /**
     * @brief Perform non-maximum suppression (clustering) on CFAR detections.
     *
     * Groups contiguous detected bins into single peaks, retaining only the index of the local maximum.
     *
     * @param detections Indices of bins initially detected by the CA-CFAR algorithm.
     * @param magnitudes Magnitude spectrum used to resolve the peak within a cluster.
     *
     * @return Vector of clustered bin indices.
     */
    std::vector<size_t> cluster_detections(const std::vector<size_t>& detections, const std::vector<double>& magnitudes);
};
