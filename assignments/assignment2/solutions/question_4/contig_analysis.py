from pathlib import Path

script_dir = Path(__file__).resolve().parent
fasta_path = script_dir / "spades_assembly" / "contigs.fasta"


total_length = 0 
contigs_lengths = []

with fasta_path.open() as fasta:
    for line in fasta:
        line = line.strip()
        if line.startswith(">"):
            contigs_lengths.append(0)
            continue
        total_length += len(line)
        contigs_lengths[-1] += len(line)


contigs_lengths.sort(reverse=True)

N50_current = 0
N50_length = ""

for item in contigs_lengths:
    N50_current += item
    if N50_current >= total_length / 2:
        N50_length = str(item)
        break



print("Total contig length:", total_length)
print("Number of contigs:", len(contigs_lengths))
print("Largest contig:", contigs_lengths[0])
print("N50 Length:", N50_length)
