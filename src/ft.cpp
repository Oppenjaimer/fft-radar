#include "ft.hpp"

#include <bit>
#include <cmath>
#include <complex>
#include <limits>
#include <numbers>
#include <utility>
#include <vector>

using fourier::Complex;
using fourier::ComplexVector;

namespace {
    /**
     * @brief In-place iterative Cooley-Tukey radix-2 FFT.
     * @param input Time-domain signal array with power-of-two size.
     */
    void fft_radix2(ComplexVector& input) {
        size_t N = input.size();
        if (N <= 1) return;

        // Bit-reversal permutation
        for (size_t i = 1, j = 0; i < N; i++) {
            size_t bit = N >> 1;
            while (j & bit) {
                j ^= bit;
                bit >>= 1;
            }

            j ^= bit;
            if (i < j) std::swap(input[i], input[j]);
        }

        // Iterative butterfly stages
        for (size_t len = 2; len <= N; len <<= 1) {
            double angle = -2.0 * std::numbers::pi / len;
            Complex w_len(std::cos(angle), std::sin(angle));
            size_t half_len = len >> 1;

            for (size_t i = 0; i < N; i += len) {
                Complex w(1.0, 0.0);

                for (size_t k = 0; k < half_len; k++) {
                    Complex u = input[i + k];
                    Complex v = input[i + k + half_len] * w;

                    input[i + k] = u + v;
                    input[i + k + half_len] = u - v;
                    w *= w_len;
                }
            }
        }
    }
}

void fourier::apply_window(ComplexVector &signal, WindowType type) {
    size_t N = signal.size();
    if (N == 0 || type == WindowType::Rectangular) return;

    double angle_step = 2.0 * std::numbers::pi / N;

    // Compute windowed values
    for (size_t n = 0; n < N; n++) {
        double window_value = 1.0;
        double angle = angle_step * n;

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
                window_value = 0.42 - 0.5 * std::cos(angle) + 0.08 * std::cos(2.0 * angle);
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

    ComplexVector output(N, Complex(0.0, 0.0));
    double angle_step = -2.0 * std::numbers::pi / N;

    // Compute output values
    for (size_t k = 0; k < N; k++) {
        for (size_t n = 0; n < N; n++) {
            double angle = angle_step * k * n;
            Complex exp_term(std::cos(angle), std::sin(angle));

            output[k] += input[n] * exp_term;
        }
    }

    return output;
}

ComplexVector fourier::fft(const ComplexVector& input) {
    size_t N = input.size();
    if (N == 0) return {};

    // Zero-pad input
    size_t padded_size = std::bit_ceil(N);
    ComplexVector buffer(padded_size, Complex(0.0, 0.0));
    std::copy(input.begin(), input.end(), buffer.begin());

    fft_radix2(buffer);
    return buffer;
}

std::vector<double> fourier::normalize_magnitudes(const ComplexVector& spectrum, WindowType type) {
    size_t N = spectrum.size();
    if (N == 0) return {};

    // Normalization factor
    double scale = 1.0 / (N * coherent_gain(type));

    // Compute normalized magnitudes
    std::vector<double> magnitudes(N / 2);
    for (size_t k = 0; k < N / 2; k++) {
        double magnitude = std::abs(spectrum[k]) * scale;
        if (k > 0) magnitude *= 2.0; // Fold negative frequency energy into positive bins

        magnitudes[k] = magnitude;
    }

    return magnitudes;
}

fourier::CfarResult fourier::ca_cfar(const std::vector<double>& magnitudes, const CfarConfig& config) {
    size_t N = magnitudes.size();
    size_t side_window = config.num_train + config.num_guard;

    CfarResult result;
    result.thresholds.assign(N, std::numeric_limits<double>::quiet_NaN());

    // Check for invalid parameter or out-of-bounds window
    if (config.num_train == 0 || N <= 2 * side_window) return result;

    // Initialize sliding window sums for the first cell under test (CUT)
    size_t i = side_window;
    double left_sum = 0.0;
    double right_sum = 0.0;

    for (size_t j = 0; j < config.num_train; j++) {
        left_sum += magnitudes[j];
        right_sum += magnitudes[i + config.num_guard + 1 + j];
    }

    double inv_total_train = 1.0 / (2.0 * config.num_train);

    // Slide window across entire spectrum
    while (i < N - side_window) {
        double noise_avg = (left_sum + right_sum) * inv_total_train;
        double threshold = noise_avg * config.threshold_factor;
        result.thresholds[i] = threshold;

        if (magnitudes[i] > threshold) result.detected_bins.push_back(i);

        // Advance sliding window if not at the last valid index
        if (i + 1 < N - side_window) {
            left_sum += magnitudes[i - config.num_guard] - magnitudes[i - side_window];
            right_sum += magnitudes[i + side_window + 1] - magnitudes[i + config.num_guard + 1];
        }

        i++;
    }

    return result;
}

std::vector<size_t> fourier::cluster_detections(const std::vector<size_t> &detections, const std::vector<double> &magnitudes) {
    std::vector<size_t> clustered;
    if (detections.empty()) return clustered;

    size_t current_peak_bin = detections[0];

    for (size_t i = 1; i < detections.size(); i++) {
        if (detections[i] == detections[i - 1] + 1) {
            // Contiguous detections, choose one with highest energy
            if (magnitudes[detections[i]] > magnitudes[current_peak_bin])
                current_peak_bin = detections[i];
        } else {
            // No longer traversing previous cluster
            clustered.push_back(current_peak_bin);
            current_peak_bin = detections[i];
        }
    }

    clustered.push_back(current_peak_bin);
    return clustered;
}
