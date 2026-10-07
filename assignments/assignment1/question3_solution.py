from math import log
import gzip 
import matplotlib.pyplot as plt
import math 



base_counts = {
    "A": 0,
    "C": 0,
    "G": 0,
    "T": 0,
    "N": 0  
}

with gzip.open("appliedgenomics2025/assignments/assignment1/chr22.fa.gz", "rt") as fasta_file:
    for line in fasta_file:
        if line.startswith(">"):
            continue
        
        sequence_line = line.strip().upper()
        for base in sequence_line:
            if base == "A":
                base_counts["A"] += 1
            if base == "C":
                base_counts["C"] += 1
            if base == "G":
                base_counts["G"] += 1
            if base == "T":
                base_counts["T"] += 1
            else:
                base_counts["N"] += 1

print(base_counts)
print("Total Chromosome Length: ", sum(base_counts.values()))

sequence_parts = []
with gzip.open("appliedgenomics2025/assignments/assignment1/chr22.fa.gz", "rt") as fasta_file:
    for line in fasta_file:
        if line.startswith(">"):
            continue
        
        sequence_parts.append(line.strip().upper())

genome_sequence = "".join(sequence_parts)

cleaned_bases = []

for base in genome_sequence:
    if base in "AGCT":
        cleaned_bases.append(base)
    else:
        cleaned_bases.append("A")

cleaned_sequence = "".join(cleaned_bases)

##2.4: Count each kmer

k = 19
kmer_frequency = {}

for start in range(len(cleaned_sequence) - k + 1):
    kmer = cleaned_sequence[start:start+k]
    if kmer in kmer_frequency:
        kmer_frequency[kmer]+=1
    else:
        kmer_frequency[kmer] = 1

print("Distinct kmers:", len(kmer_frequency))

frequency_spectrum = {}

for count in kmer_frequency.values():
    if count in frequency_spectrum:
        frequency_spectrum[count]+=1
    else:
        frequency_spectrum[count] = 1

for frequency in range(1, 21):
    print(frequency, ":", frequency_spectrum[frequency] )

##3.3: PLot frequencies
x_values = list(frequency_spectrum.keys())
y_values = list(frequency_spectrum.values())


plt.scatter(x_values, y_values)
plt.xscale("log")
plt.yscale("log")
plt.title("19mer Distribition")
plt.xlabel("19mer Frequency")
plt.ylabel("# 19mers that occur x times")
plt.show()

##3.4 uniqueness

genome_uniqueness = 100 * frequency_spectrum.get(1, 0) / sum(kmer_frequency.values())
genome_repetitiveness = 100 - genome_uniqueness

over_1000_positions = 0

for frequency, number_of_kmers in frequency_spectrum.items():
    if frequency > 1000:
        over_1000_positions += number_of_kmers * frequency

percent_1000 = over_1000_positions / sum(kmer_frequency.values())
