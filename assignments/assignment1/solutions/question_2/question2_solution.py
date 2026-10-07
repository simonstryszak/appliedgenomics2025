import numpy as np
import matplotlib.pyplot as plt
import math




def simulate_coverage(genome, depth, read_length, seed):
    needed_reads = genome * depth / read_length
    genome_coverage = np.zeros(genome, dtype=int)
    rng = np.random.default_rng(seed=42)
    for _ in range(int(needed_reads)):
        start_pos = rng.integers(len(genome_coverage) - read_length + 1)
        end_pos = start_pos + read_length
        genome_coverage[start_pos:end_pos] +=1

    histogram = np.bincount(genome_coverage)
    coverage_values = np.arange(len(histogram))
    observed_proportion = histogram / genome
    
    poisson_probabilities = []
    for k in coverage_values:
        probability = (
            math.exp(-depth)
            * depth ** int(k)
            / math.factorial(int(k))
        )
        poisson_probabilities.append(probability)

    normal_mean = genome_coverage.mean()
    normal_sd = genome_coverage.std()

    normal_probabilities = []
    for k in coverage_values:
        normal_probabilities.append(
            1 / (normal_sd * np.sqrt(2 * np.pi))
            * np.exp(
                -0.5
                * ((k - normal_mean) / normal_sd) ** 2
            )
        )   
    
    plt.bar(coverage_values, observed_proportion, label="Observed distribution")
    plt.plot(
        coverage_values,
        poisson_probabilities,
        marker="o",
        label="Poisson distribution"
    )
    plt.plot(
        coverage_values,
        normal_probabilities,
        marker="x",
        label="Normal approximation"
    )

    plt.xlabel("Coverage")
    plt.ylabel("Proportion of genome positions")
    plt.title(f"Coverage distribution for simulated {depth}x sequencing")
    plt.legend()    
    plt.show()

    print(genome_coverage.mean())
    print(genome_coverage.min())
    print(genome_coverage.max())
    print(genome_coverage.sum())

    zero_coverage  = observed_proportion[0]
    expected_zero_proportion = math.exp(-depth)
    expected_zero_percent = expected_zero_proportion * 100
    expected_zero_bases = expected_zero_proportion * genome

    print("percentage of genome with zero coverage: ", zero_coverage * 100)
    print("Expected percentage of genome with zero coverage:", expected_zero_percent)
    print("Number of bases with zero coverage:", int(zero_coverage * genome))
    print("Expected number of bases with zero coverage:", expected_zero_bases)

    difference = zero_coverage * genome - expected_zero_bases
    print("Observed minus expected:", difference)

    return(genome_coverage)

tenx_coverage = simulate_coverage(
    genome = 1000000,
    depth = 10,
    read_length = 100,
    seed = 42
)



