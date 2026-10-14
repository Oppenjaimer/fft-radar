#include "ft.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <random>
#include <ratio>
#include <utility>
#include <vector>

namespace {
    /**
     * @brief Radar simulation parameters.
     */
    struct RadarConfig {
        size_t num_samples = 8192;          ///< Number of time-domain samples.
        double sample_rate = 1e4;           ///< Sampling frequency: 10 kHz.
        double radar_freq = 24e9;           ///< Transmission frequency: 24 GHz (K-band).
        double light_speed = 299792458.0;   ///< Speed of light (m/s).
        double noise_std_dev = 1.0;         ///< Noise standard deviation.

        constexpr double hz_per_bin() const {
            return sample_rate / num_samples;
        }

        constexpr double bin_to_doppler_freq(size_t bin) const {
            return static_cast<double>(bin) * hz_per_bin();
        }

        constexpr double doppler_freq_to_velocity(double doppler_freq) const {
            double velocity = (doppler_freq * light_speed) / (2.0 * radar_freq);
            return velocity * 3.6; // Convert to km/h
        }
    };

    /**
     * @brief Physical radar target injected into the time-domain simulation.
     */
    struct SimulatedTarget {
        double frequency;   ///< Cycles per observation window.
        double amplitude;   ///< Signal strength.
        double phase;       ///< Initial phase offset.
    };

    /**
     * @brief Valid radar target extracted by the detection algorithm.
     */
    struct DetectedTarget {
        size_t bin;             ///< Frequency bin index.
        double doppler_freq;    ///< Doppler frequency shift (Hz).
        double radial_velocity; ///< Relative radial velocity (km/h).
        double magnitude;       ///< Normalized magnitude of the peak.
    };

    /**
     * @brief DFT/FFT benchmark timings.
     */
    struct Benchmark {
        double dft_time;    ///< DFT timing.
        double fft_time;    ///< FFT timing.
    };

    /**
     * @brief Synthesize time-domain radar return from multiple targets plus white Gaussian noise.
     *
     * @param config  Radar configuration parameters.
     * @param targets Targets to simulate.
     * @param seed    Random number generator seed.
     *
     * @return ComplexVector with composite radar signal.
     */
    fourier::ComplexVector generate_radar_signal(
        const RadarConfig& config,
        const std::vector<SimulatedTarget>& targets,
        unsigned int seed = std::random_device{}()
    ) {
        fourier::ComplexVector signal(config.num_samples);
        std::mt19937 gen(seed);
        std::normal_distribution<double> noise_dist(0.0, config.noise_std_dev);

        double angle_step = 2.0 * std::numbers::pi / config.num_samples;

        // Generate signal values from targets
        for (size_t n = 0; n < config.num_samples; n++) {
            double signal_value = 0.0;

            for (const auto& target : targets) {
                double angle = angle_step * target.frequency * n + target.phase;
                signal_value += target.amplitude * std::cos(angle);
            }

            signal_value += noise_dist(gen);
            signal[n] = fourier::Complex(signal_value, 0.0);
        }

        return signal;
    }

    /**
     * @brief Measure execution time of a callable.
     *
     * @tparam Func Type of the callable object.
     * @param func Callable object to be executed. Must return a non-void type.
     *
     * @return std::pair containing the return value of the callable and the elapsed execution time in ms.
     */
    template <typename Func>
    auto measure_execution_time(Func&& func) {
        auto start = std::chrono::high_resolution_clock::now();
        auto result = std::forward<Func>(func)();
        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> elapsed = end - start;
        return std::pair{ std::move(result), elapsed.count() };
    }

    /**
     * @brief Convert clustered detections into physical radar target measurements.
     *
     * @param detection_bins Indices of valid targets.
     * @param magnitudes     Normalized magnitude spectrum to extract signal strength.
     * @param config         Radar configuration parameters.
     *
     * @return Vector of resolved DetectedTarget objects.
     */
    std::vector<DetectedTarget> evaluate_targets(
        const std::vector<size_t>& detection_bins,
        const std::vector<double>& magnitudes,
        const RadarConfig& config
    ) {
        std::vector<DetectedTarget> results;
        results.reserve(detection_bins.size());

        // Resolve target for each detection bin
        for (size_t bin : detection_bins) {
            double doppler_freq = config.bin_to_doppler_freq(bin);
            results.push_back(DetectedTarget{
                .bin = bin,
                .doppler_freq = doppler_freq,
                .radial_velocity = config.doppler_freq_to_velocity(doppler_freq),
                .magnitude = magnitudes[bin]
            });
        }

        return results;
    }

    /**
     * @brief Print formatted benchmark and target evaluation.
     *
     * @param timings Benchmark timings.
     * @param targets Resolved targets.
     */
    void print_results(const Benchmark& timings, const std::vector<DetectedTarget>& targets) {
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "========= BENCHMARK =========\n"
                  << "DFT:     " << timings.dft_time << " ms\n"
                  << "FFT:     " << timings.fft_time << " ms\n"
                  << "Speedup: " << timings.dft_time / timings.fft_time << "x\n\n";

        std::cout << "========= TARGETS =========\n";
        if (targets.empty()) {
            std::cout << "No targets detected\n";
            return;
        }

        for (const auto& target : targets) {
            std::cout << "Target at bin " << target.bin << ":\n"
                      << "  Doppler frequency: " << target.doppler_freq << " Hz\n"
                      << "  Radial velocity:   " << target.radial_velocity << " km/h\n"
                      << "  Magnitude:         " << target.magnitude << "\n";
        }
    }

    /**
     * @brief Export results to a CSV file.
     *
     * @param filepath   Path to output CSV file.
     * @param raw_signal Raw time-domain signal array processed by radar.
     * @param magnitudes Normalized single-sided magnitude spectrum.
     * @param thresholds CFAR detection thresholds.
     */
    void export_results(
        const std::filesystem::path& filepath,
        const fourier::ComplexVector& raw_signal,
        const std::vector<double>& magnitudes,
        const std::vector<double>& thresholds
    ) {
        std::ofstream file(filepath);
        if (!file.is_open()) {
            std::cerr << "Unable to open CSV file\n";
            return;
        }

        file << "index,time_ampl,norm_mag,cfar_threshold\n";

        size_t N = raw_signal.size();
        size_t half_N = magnitudes.size();

        for (size_t i = 0; i < N; i++) {
            file << i << "," << raw_signal[i].real() << ",";

            if (i < half_N) {
                // First half of the spectrum
                file << magnitudes[i] << ",";
                if (std::isnan(thresholds[i])) file << "\n";
                else file << thresholds[i] << "\n";
            } else {
                // Second half of the spectrum
                file << ",\n";
            }
        }
    }
}

/**
 * @brief Program entry point.
 * @return Exit code.
 */
int main() {
    constexpr RadarConfig radar_config{};
    constexpr auto window_type = fourier::WindowType::Hann;
    constexpr fourier::CfarConfig cfar_config{};

    const std::vector<SimulatedTarget> simulated_targets = {
        {45.0,   1.2, 0.0},
        {850.3,  0.7, std::numbers::pi / 4.0},
        {1500.0, 0.25, std::numbers::pi},
        {2048.0, 1.0, 0.0},
        {2058.0, 0.9, std::numbers::pi / 3.0},
        {3800.5, 0.85, 1.27}
    };

    // Signal synthesis and windowing
    fourier::ComplexVector raw_signal = generate_radar_signal(radar_config, simulated_targets);
    fourier::ComplexVector windowed_signal = raw_signal;
    fourier::apply_window(windowed_signal, window_type);

    // Benchmark
    auto [dft_spectrum, dft_time] = measure_execution_time([&] {
        return fourier::dft(windowed_signal);
    });
    auto [fft_spectrum, fft_time] = measure_execution_time([&] {
        return fourier::fft(windowed_signal);
    });

    // Detection, CA-CFAR and clustering
    std::vector<double> magnitudes = fourier::normalize_magnitudes(fft_spectrum, window_type);
    fourier::CfarResult cfar_output = fourier::ca_cfar(magnitudes, cfar_config);
    std::vector<size_t> detection_bins = fourier::cluster_detections(cfar_output.detected_bins, magnitudes);
    std::vector<DetectedTarget> detected_targets = evaluate_targets(detection_bins, magnitudes, radar_config);

    // Display and export results
    print_results({dft_time, fft_time}, detected_targets);
    export_results("data/data.csv", raw_signal, magnitudes, cfar_output.thresholds);

    return 0;
}
