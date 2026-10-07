import argparse
import random

parser = argparse.ArgumentParser(
    description = "Introduce random mutations into known FASTA sequence."
)

parser.add_argument("-i", "--input", required=True)
parser.add_argument("-o", "--output", required=True)
parser.add_argument("-m", "--mutation-rate", type=float, required=True)
parser.add_argument("-s", "--seed", type=int, required=True)

args = parser.parse_args()

def read_fasta(filename):
    header = ""
    sequence = ""
    for line in filename:
        line = line.strip()
        if line.startswith(">"):
            header = line[1:]
        else:
            sequence += line.upper()
    return header, sequence

with open(args.input, "r") as chrom22_file:
    header, sequence = read_fasta(chrom22_file)

random.seed(args.seed)

number_of_mutations = round(len(sequence) * args.mutation_rate)

mutation_positions = random.sample(
    range(len(sequence)),
    number_of_mutations
)

mutation_sequence = list(sequence)
for index in mutation_positions:
    base = mutation_sequence[index]
    choices = "ACGT".replace(base, "")
    mutation_sequence[index] = random.choice(choices)
mutation_sequence = "".join(mutation_sequence)

with open(args.output, "w") as output_file:
    output_file.write(f">{header}\n")
    output_file.write(mutation_sequence + "\n")
