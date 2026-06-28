import itertools
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import AlignIO, SeqIO
from Bio.SeqUtils import gc_fraction

from src.logger import logger


def pre_alignment_eda(config):
    """
    Performs the pre-alignment exploratory analysis for filtered sequences and alignment.

    Args:
        config (dict): Configuration dictionary containing the input/output file paths.

    Returns:
        dict: Summary statistics, dataframe and output file paths.

    Raises:
        KeyError: If required keys are missing in the configuration.
        FileNotFoundError: If any input file does not exist.
        ValueError: If the input files do not contain valid data.
    """
    try:
        filtered_sequences_path = config["paths"]["filtered_sequences"]
        output_dir = config["paths"]["pre_alignment_eda"]
        summary_df_path = config["paths"]["pre_alignment_summary"]
        
    except KeyError as missing_key:
        logger.error(f"Configuration error. Missing key: {missing_key}")
        raise KeyError(f"Error de configuración. Falta la clave: {missing_key}")

    if not Path(filtered_sequences_path).is_file():
        logger.error(f"Input file not found: {filtered_sequences_path}")
        raise FileNotFoundError(f"Input file not found: {filtered_sequences_path}")

    Path(output_dir).parent.mkdir(parents=True, exist_ok=True)

    records = list(SeqIO.parse(filtered_sequences_path, "fasta"))
    if not records:
        logger.error("No sequences found in filtered FASTA file")
        raise ValueError("No sequences found in filtered FASTA file")


    seq_data = []
    sequences_str = []

    for record in records:
        seq_str = str(record.seq).upper()
        sequences_str.append(seq_str)
        seq_len = len(seq_str)
        ambiguous_count = seq_str.count("N") + seq_str.count("R") + seq_str.count("Y")
        ambiguous_perc = (ambiguous_count / seq_len) * 100 if seq_len > 0 else 0
        gc_content = gc_fraction(record.seq) * 100
        seq_data.append({"ID": record.id, "Length": seq_len, "GC_Content": gc_content, "Ambiguous_Perc": ambiguous_perc})

    df_seqs = pd.DataFrame(seq_data)

    q1 = df_seqs["GC_Content"].quantile(0.25)
    q3 = df_seqs["GC_Content"].quantile(0.75)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    df_seqs["GC_Outlier"] = (df_seqs["GC_Content"] < lower_bound) | (df_seqs["GC_Content"] > upper_bound)

    pairwise_coverages = []
    pairwise_identities = []
  
    for seq1, seq2 in itertools.combinations(sequences_str, 2):
        min_len = min(len(seq1), len(seq2))
        max_len = max(len(seq1), len(seq2))
        matches= 0
        for a, b in zip(seq1, seq2):
            if a==b:
                matches+=1
        identity = (matches / min_len) * 100
        coverage = (min_len / max_len) * 100
        pairwise_identities.append(identity)
        pairwise_coverages.append(coverage)


    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()

    sns.histplot(df_seqs["Length"], bins=50, kde=True, ax=axes[0], color="skyblue")
    axes[0].set_title("Sequence Length Distribution")
    axes[0].set_xlabel("Length (bp)")
    axes[0].set_ylabel("Frequency")

    sns.histplot(df_seqs["GC_Content"], bins=50, kde=True, ax=axes[1], color="lightgreen")
    axes[1].axvline(lower_bound, color="red", linestyle="dashed", linewidth=1.5, label="Q1 - 1.5*IQR")
    axes[1].axvline(upper_bound, color="red", linestyle="dashed", linewidth=1.5, label="Q3 + 1.5*IQR")
    axes[1].set_title("GC Content per Sequence (%)")
    axes[1].set_xlabel("GC Content (%)")
    axes[1].set_ylabel("Frequency")
    axes[1].legend()

    sns.histplot(df_seqs["Ambiguous_Perc"], bins=50, kde=True, ax=axes[2], color="salmon")
    axes[2].set_title("Ambiguous Bases (N, R, Y)")
    axes[2].set_xlabel("Ambiguous Bases (%)")
    axes[2].set_ylabel("Frequency")

    sns.histplot(pairwise_identities, bins=50, kde=True, ax=axes[3], color="purple")
    axes[3].set_title("Pairwise Identity (%) - Unaligned")
    axes[3].set_xlabel("Identity (%)")
    axes[3].set_ylabel("Frequency")
    
    sns.histplot(pairwise_coverages, bins=50, kde=True, ax=axes[4], color="gold")
    axes[4].set_title("Pairwise Coverage (%) - Unaligned")
    axes[4].set_xlabel("Coverage Relative to Longest (%)")
    axes[4].set_ylabel("Frequency")
    
    axes[5].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_dir, dpi=200, bbox_inches="tight")
    plt.close(fig)

    summary_data = {
        "n_sequences": len(df_seqs),
        "mean_length": float(df_seqs["Length"].mean()),
        "median_length": float(df_seqs["Length"].median()),
        "mean_gc_content": float(df_seqs["GC_Content"].mean()),
        "median_gc_content": float(df_seqs["GC_Content"].median()),
        "max_ambiguous_perc": float(df_seqs["Ambiguous_Perc"].max()),
        "gc_outliers": int(df_seqs["GC_Outlier"].sum()),
        "mean_pairwise_identity": float(np.mean(pairwise_identities)),
        "mean_pairwise_coverage": float(np.mean(pairwise_coverages))
    }
    
    pd.DataFrame([summary_data]).to_csv(summary_df_path, index=False)

    logger.info(f"Pre-alignment EDA completed. Results saved in {output_dir} and in {summary_df_path}")

    logger.info({
        "summary": summary_data,
        "dataframe": df_seqs,
        "pairwise_identities": pairwise_identities,
        "output_files": {
            "image_resume": str(output_dir),
            "summary_csv": str(summary_df_path),}})
