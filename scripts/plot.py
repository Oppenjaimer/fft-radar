import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("data/data.csv")

N = len(df)
N_half = N // 2

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# Time domain
ax1.plot(df["index"], df["time_ampl"], "tab:blue")
ax1.set_title("Time Domain")
ax1.set_xlabel("Sample ($n$)")
ax1.set_ylabel("Amplitude")
ax1.grid(True)

# Normalize amplitudes
positive_freqs = df["index"][:N_half]
normalized_mag = (df["freq_ampl"][:N_half] / N) * 2

# Undo multiplication by 2 for DC bin
if len(normalized_mag) > 0:
    normalized_mag[0] = normalized_mag[0] / 2

# Frequency domain
ax2.stem(positive_freqs, normalized_mag, basefmt=" ", markerfmt="", linefmt="tab:red")
ax2.set_title("Frequency Domain")
ax2.set_xlabel("Frequency Bin ($k$)")
ax2.set_ylabel("Normalized Magnitude")
ax2.grid(True)

plt.tight_layout()
plt.show()
