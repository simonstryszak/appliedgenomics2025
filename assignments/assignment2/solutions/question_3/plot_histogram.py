
import matplotlib.pyplot as plt
from pathlib import Path

output_dir = Path(__file__).resolve().parent

frequencies = []
kmer_counts = []

with (output_dir / "reads.histo").open() as histogram:
    for line in histogram:
        frequency, count = map(int, line.split())
        frequencies.append(frequency)
        kmer_counts.append(count)

plt.plot(frequencies, kmer_counts)
plt.xlabel("K-mer frequency")
plt.ylabel("Number of distinct 21-mers")
plt.title("21-mer frequency spectrum")
plt.xlim(0, 200)
plt.yscale("log")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(output_dir / "reads_histogram.png", dpi=200)
plt.show()
