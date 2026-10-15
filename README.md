# FFT Radar Target Detection

A 1D radar signal processing and target detection engine built in C++. It synthesizes time-domain radar returns from multiple targets with white Gaussian noise, implements a custom radix-2 Cooley-Tukey fast Fourier transform (FFT), and extracts physical targets using a 1D cell-averaging constant false alarm rate (CA-CFAR) algorithm. Data visualization is handled via Python using [Pandas](https://pandas.pydata.org/) and [Matplotlib](https://matplotlib.org/).

![Screenshot](https://github.com/Oppenjaimer/fft-radar/blob/master/images/screenshot.png?raw=true)

## Features

* Implementation of both a naive **discrete Fourier transform** (DFT) and an in-place iterative **radix-2 Cooley-Tukey FFT** to convert signals into the frequency domain.
* Built-in **signal windowing** supporting rectangular, Hann, Hamming, and Blackman window functions to mitigate spectral leakage prior to Fourier processing.
* One-dimensional **CA-CFAR algorithm** to dynamically estimate the noise floor via sliding windows and calculate detection thresholds.
* **Non-maximum suppression** that groups contiguous detected bins into a single cluster, retaining only the index of the local maximum.
* **Composite synthetic radar signals generation** by injecting multiple targets (defined by frequency, amplitude, and phase) alongside white Gaussian noise.
* **Automatic conversion** of raw frequency bins into Doppler frequencies and relative radial velocities based on a 24 GHz K-band radar configuration.

## Building

The project uses a `Makefile` and requires `g++` to compile. You can automatically compile and run the project with the included shell script.

```bash
chmod +x run.sh
./run.sh
```

Alternatively, to manually generate the executable binary, run `make` in the root directory. This will compile the source files and output a standalone executable in the `bin/` directory.

Object files are placed in the `build/` directory.

## Usage

To generate the radar signal, run the DFT/FFT benchmarks, and execute the CA-CFAR target detection, run the simulation binary.

```bash
./bin/radar_sim
```

The console will output a benchmark comparing DFT and FFT execution speeds, followed by a list of resolved targets including their exact Doppler frequency, radial velocity, and normalized magnitude.

The raw signal, normalized single-sided magnitude spectrum, and CFAR thresholds are exported to `data/data.csv`. Once the CSV is generated, use the included Python script (requires `pandas` and `matplotlib`) to visualize the time domain and frequency domain graphs alongside the CFAR detections.

```bash
python scripts/plot.py
```

## License

This project is available under the MIT License.
