import argparse
import math

parser = argparse.ArgumentParser(
    "Load two DNA sequences for Jaccard Calculation"
)
parser.add_argument("-a", required=True)
parser.add_argument("-b", required=True)
parser.add_argument('-k', type=int, required=True)

args = parser.parse_args()


def clean_fasta(file):
    cleaned_file = ""
    for line in file:
        line = line.strip()
        if line.startswith(">"):
            continue
        else:
            for base in line:
                if base not in ["a", "t", "c", "g", "A", "T", "C", "G"]:
                    base = "N"
                cleaned_file += base.upper()
    return cleaned_file


def get_kmers(sequence,k):
    kmers = set()

    for index in range(len(sequence) - k + 1):
        kmer = sequence[index:index+k]
        kmers.add(kmer)

    return kmers

def compute_metrics(a, b, k):
    a_kmers = get_kmers(a, k)
    b_kmers = get_kmers(b, k)

    intersection = a_kmers & b_kmers
    union = a_kmers | b_kmers

    jaccard = len(intersection) / len(union)

    ani = ((2 * jaccard) / (1 + jaccard)) ** (1/k)
    approximate_ani = 1 + math.log((2 * jaccard)/ (1 + jaccard)) / k

    return jaccard, ani, approximate_ani

with open(args.a, "r") as input1:
    a = clean_fasta(input1)
with open(args.b, "r") as input2:
    b = clean_fasta(input2)
k = args.k

jaccard, ani, approximate_ani = compute_metrics(a, b, k)

print(jaccard, ani, approximate_ani)


    
    
