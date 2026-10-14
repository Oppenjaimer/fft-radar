import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("data/data.csv")
freq_df = df.dropna(subset=["norm_mag"]).copy()

# Identify bins where magnitude exceeds CFAR threshold
detections = freq_df[freq_df["norm_mag"] > freq_df["cfar_threshold"]]

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

# Time domain
ax1.plot(df["index"], df["time_ampl"], "tab:blue", linewidth=0.8)
ax1.set_title("Time Domain Radar Return (Signal + Noise)")
ax1.set_xlabel("Sample ($n$)")
ax1.set_ylabel("Amplitude")
ax1.grid(True)

# Frequency domain
ax2.plot(freq_df["index"], freq_df["norm_mag"], color="tab:red", label="Single-Sided Spectrum")
ax2.plot(freq_df["index"], freq_df["cfar_threshold"], color="tab:green", linestyle="--", label="CA-CFAR Threshold")
ax2.scatter(detections["index"], detections["norm_mag"], color="tab:gray", zorder=5, s=30, label="CFAR Detections")
ax2.set_title("Frequency Domain & CA-CFAR Target Detection")
ax2.set_xlabel("Frequency Bin ($k$)")
ax2.set_ylabel("Normalized Magnitude")
ax2.legend(loc="upper right")
ax2.grid(True)

plt.tight_layout()
plt.show()
