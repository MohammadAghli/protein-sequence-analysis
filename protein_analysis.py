#!/usr/bin/env python3
"""
Basic Protein Sequence Analysis
Dataset: Breast cancer-related proteins from UniProt
"""

import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns

# =============================================================================
# STEP 1: Load and parse the FASTA file
# =============================================================================
print("📁 Loading FASTA file...")

protein_data = []
for record in SeqIO.parse("uniprotkb_breast_cancer_AND_reviewed_tr.fasta", "fasta"):
    protein_data.append({
        'protein_id': record.id,
        'description': record.description,
        'sequence': str(record.seq),
        'sequence_length': len(record.seq)
    })

# Convert to DataFrame
df = pd.DataFrame(protein_data)
print(f"✅ Loaded {len(df)} protein sequences")

# =============================================================================
# STEP 2: Basic Data Exploration
# =============================================================================
print("\n📊 Basic statistics:")
print(f"Average sequence length: {df['sequence_length'].mean():.2f}")
print(f"Shortest sequence: {df['sequence_length'].min()}")
print(f"Longest sequence: {df['sequence_length'].max()}")

# =============================================================================
# STEP 3: Calculate amino acid composition
# =============================================================================
print("\n🧬 Calculating amino acid composition...")

# Define the 20 standard amino acids (CORRECTED VERSION)
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def calculate_aa_composition(sequence):
    """Calculate amino acid composition for a sequence"""
    total_length = len(sequence)
    composition = {}
    for aa in amino_acids:
        composition[aa] = (sequence.count(aa) / total_length) * 100  # Percentage
    return composition

# Apply to all sequences
aa_composition = df['sequence'].apply(calculate_aa_composition)
aa_df = pd.DataFrame(aa_composition.tolist())

# Add protein ID to composition dataframe
aa_df['protein_id'] = df['protein_id']

# =============================================================================
# STEP 4: Create Visualizations
# =============================================================================
print("\n📈 Creating visualizations...")

# Set up the plotting style
sns.set_style("whitegrid")
plt.figure(figsize=(15, 10))

# Plot 1: Distribution of sequence lengths
plt.subplot(2, 2, 1)
sns.histplot(df['sequence_length'], bins=20, kde=True)
plt.title('Distribution of Protein Sequence Lengths')
plt.xlabel('Sequence Length')
plt.ylabel('Frequency')

# Plot 2: Top 10 most abundant amino acids (average across all proteins)
plt.subplot(2, 2, 2)
mean_aa_composition = aa_df[list(amino_acids)].mean().sort_values(ascending=False)
sns.barplot(x=mean_aa_composition.head(10).values, y=mean_aa_composition.head(10).index)
plt.title('Top 10 Most Abundant Amino Acids')
plt.xlabel('Average Percentage (%)')

# Plot 3: Scatter plot of length vs leucine content
plt.subplot(2, 2, 3)
plt.scatter(df['sequence_length'], aa_df['L'], alpha=0.7)
plt.title('Sequence Length vs Leucine Content')
plt.xlabel('Sequence Length')
plt.ylabel('Leucine Content (%)')

# Plot 4: Boxplot of hydrophobic amino acids
plt.subplot(2, 2, 4)
hydrophobic_aas = ['A', 'V', 'L', 'I', 'P', 'F', 'M', 'W']
sns.boxplot(data=aa_df[hydrophobic_aas])
plt.title('Distribution of Hydrophobic Amino Acids')
plt.xticks(rotation=45)

plt.tight_layout()
plt.savefig('protein_analysis_results.png', dpi=300, bbox_inches='tight')
plt.show()

print("✅ Analysis complete! Check 'protein_analysis_results.png'")

# =============================================================================
# STEP 5: Save results to CSV
# =============================================================================
results_df = df[['protein_id', 'sequence_length']].copy()
for aa in amino_acids:
    results_df[f'{aa}_pct'] = aa_df[aa]

results_df.to_csv('protein_analysis_results.csv', index=False)
print("✅ Results saved to 'protein_analysis_results.csv'")
