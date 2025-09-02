import matplotlib.pyplot as plt

# Load data
sizes = []
times = []
with open("results.txt") as f:
    for line in f:
        n, t = line.split()
        sizes.append(int(n))
        times.append(float(t))

# Plot
plt.figure(figsize=(8,6))
plt.plot(sizes, times, marker='o')

plt.xscale("log")  # x-axis in log scale
plt.yscale("log")  # y-axis in log scale (optional, but helps show trends)

plt.xlabel("Array size (N)")
plt.ylabel("Time (ms)")
plt.title("Merge Sort Benchmark on Centaurus")
plt.grid(True, which="both", linestyle="--", linewidth=0.5)

plt.savefig("benchmark_plot.png", dpi=300)
plt.show()
