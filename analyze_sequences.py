from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

# Set path to your FASTA file
fasta_file = "../data/gene.fna"

# Load sequences
records = list(SeqIO.parse(fasta_file, "fasta"))
print(f"✅ Loaded {len(records)} sequences from {fasta_file}")

# Check if sequences were found
if len(records) == 0:
    print("⚠️ No sequences found. Please check the file.")
    exit()

# Analyze sequences
results = []
for record in records:
    seq = record.seq
    gc = 100 * (seq.count("G") + seq.count("C")) / len(seq)
    results.append({
        "Sequence_ID": record.id,
        "Length": len(seq),
        "GC_Content (%)": round(gc, 2)
    })

# Save summary to CSV
df = pd.DataFrame(results)
csv_path = "../results/sequence_summary.csv"
df.to_csv(csv_path, index=False)
print(f"✅ Sequence summary saved to {csv_path}")

# Plot GC content
plt.figure(figsize=(10, 5))
plt.bar(df["Sequence_ID"], df["GC_Content (%)"])
plt.title("GC Content per Gene")
plt.xlabel("Sequence ID")
plt.ylabel("GC Content (%)")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()

plot_path = "../results/gc_content_plot.png"
plt.savefig(plot_path)
plt.show()
print(f"✅ GC content plot saved to {plot_path}")

