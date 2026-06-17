from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
import pandas as pd
import json

def has_premature_stops(seq_str):
    reverse_seq = str(Seq(seq_str).reverse_complement())
    for seq in [seq_str, reverse_seq]:
        for frame in [0,1,2]:
            fragment = seq[frame: len(seq) - (len(seq) - frame) % 3]
            translated_seq = str(Seq(fragment).translate())
            if translated_seq[:-1].count('*') > 0:
                return False
    return True

def filter_sequences(config):
    min_length = config["Filtrar"]["min_length"]
    max_length = config["Filtrar"]["max_length"]
    max_ambiguous_bases = config["Filtrar"]["max_ambiguous_bases"]

    filtered_records = []
    seen_ids= set()
    
    for record in SeqIO.parse(config["paths"]["raw_sequences"], "fasta"):
        sequence_length = len(record.seq)
        ambiguous_bases_count = record.seq.count('N')
        
        if (min_length <= sequence_length <= max_length) and (ambiguous_bases_count <= max_ambiguous_bases):
            if not has_premature_stops(str(record.seq)):
                if record.id not in seen_ids:
                    filtered_records.append(record)
                    seen_ids.add(record.id)

        
    SeqIO.write(filtered_records, config["paths"]["filtered_sequences"], "fasta")
    print(f"Filtered sequences saved to {config['paths']['filtered_sequences']}: {len(filtered_records)} records")

