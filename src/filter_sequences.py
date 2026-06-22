from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
import pandas as pd
import json
from pathlib import Path
from src.logger import logger

def has_premature_stops(seq_str):
    """
    Checks if a given nucleotide sequence has at least one valid Open Reading Frame 
    (ORF) without premature stop codons across all 6 possible reading frames.
    
    Args:
        seq_str (str): The nucleotide sequence as a string.
        
    Returns:
        bool: True if finds a valid ORF without premature stop codons, False otherwise. 
    """
    reverse_seq = str(Seq(seq_str).reverse_complement())
    for seq in [seq_str, reverse_seq]:
        for frame in [0,1,2]:
            fragment = seq[frame: len(seq) - (len(seq) - frame) % 3]
            translated_seq = str(Seq(fragment).translate())
            if translated_seq[:-1].count('*') == 0:
                return False
    return True


def filter_sequences(config):
    """
    Filters sequences from a FASTA file based on length, maximum ambiguous bases, 
    and presence of a valid ORF, streaming the valid records directly to a new file.
    
    Args:
        config (dict): A configuration dictionary containing filter parameters and 
            file paths. It must adhere to the following structure:
            {
                "Filter": {
                    "min_length": int, 
                    "max_length": int, 
                    "max_ambiguous_bases": int
                },
                "paths": {
                    "raw_sequences": str, 
                    "filtered_sequences": str
                }
            }

    Returns:
        None. The function writes the filtered FASTA sequences directly to the file 
        system path specified in the `config` dictionary.

    Raises:
        KeyError: If any required nested keys are missing from the `config` dictionary.
    """
    try:
        min_length = config["Filter"]["min_length"]
        max_length = config["Filter"]["max_length"]
        max_ambiguous_bases = config["Filter"]["max_ambiguous_bases"]
        raw_seqs = config["paths"]["raw_sequences"]
        filtered_seqs = config["paths"]["filtered_sequences"]
        
    except KeyError as missing_key:
        logger.error(f"Configuration error. Missing key: {missing_key}")
        raise KeyError(f"Error de configuración. Falta la clave: {missing_key}")
    
    Path(filtered_seqs).parent.mkdir(parents=True, exist_ok=True)

    filtered_records = []
    seen_ids= set()
    
    for record in SeqIO.parse(raw_seqs, "fasta"):
        sequence_length = len(record.seq)
        ambiguous_bases_count = record.seq.count('N')
        
        if (min_length <= sequence_length <= max_length) and (ambiguous_bases_count <= max_ambiguous_bases):
            if not has_premature_stops(str(record.seq)):
                if record.id not in seen_ids:
                    filtered_records.append(record)
                    seen_ids.add(record.id)

        
    SeqIO.write(filtered_records, filtered_seqs, "fasta")
    logger.info(f"Filtered sequences saved to {filtered_seqs}: {len(filtered_records)} records")

