from pathlib import Path
import pandas as pd 

DATA_DIR = Path(__file__).resolve().parents[1]
chrom_files = {
    "E. coli":  "ecoli.chrom.sizes",
    "Yeast": "yeast.chrom.sizes",
    "Worm": "ce10.chrom.sizes",
    "Fruit Fly": "dm6.chrom.sizes",
    "Arabidopsis thaliana": "TAIR10.chrom.sizes",
    "Tomato": "tomato.chrom.sizes",
    "Human": "hg38.chrom.sizes",
    "Wheat": "wheat.chrom.sizes"
}   

rows = []

for species, filename in chrom_files.items():
    path = DATA_DIR / filename

    df = pd.read_csv(
        path,
        sep = "\t",
        names=["chromosomes", "length"]
    )
    
    largest_row = df.loc[df["length"].idxmax()]
    smallest_row = df.loc[df["length"].idxmin()]
    
    rows.append({
        "Species": species,
        "Genome Size": df["length"].sum(),
        "Number of Chromosomes": len(df),
        "Largest Chromosome": largest_row["chromosomes"],
        "Largest Chromosome Length": largest_row["length"],
        "Smallest Chromosome": smallest_row["chromosomes"],
        "Smallest Chromosome Length": smallest_row["length"],
        "Mean Chromosome Length": df["length"].mean()
    })

summary_df = pd.DataFrame(rows)

print(summary_df.to_string(index=False))
