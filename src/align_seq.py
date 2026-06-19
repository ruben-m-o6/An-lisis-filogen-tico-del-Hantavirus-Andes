import shutil
import subprocess
from Bio import AlignIO
from pathlib import Path

def align_seq(config):
    
    """
    Concatenates filtered sequences with an outgroup and aligns them using MAFFT.

    Args:
        config (dict): A configuration dictionary containing file paths and 
            alignment parameters. It must adhere to the following structure:
            {
                "alignment": {
                    "IterativeRefinementMethod": str, 
                    "maxiterate": int
                },
                "paths": {
                    "filtered_sequences": str, 
                    "outgroup_sequence": str,
                    "combined_sequences": str,
                    "alignment": str
                }
            }

    Returns:
        None. The combined FASTA and the final alignment are written directly to disk.

    Raises:
        KeyError: If any required nested keys are missing from the `config` dictionary.
        FileNotFoundError: If input files are missing, or if the 'mafft' executable 
            is not installed/accessible in the system's PATH.
        RuntimeError: If the MAFFT process fails or exits with an error code.
        ValueError: If the generated alignment cannot be parsed as a valid FASTA file.
    """
    
    try: 
        filtered_path = config["paths"]["filtered_sequences"]
        outgroup_path = config["paths"]["outgroup_sequence"]
        combined_path = config["paths"]["combined_sequences"]
        aligned_path = config["paths"]["alignment"]
        IterativeRefinementMethod = config["alignment"]["IterativeRefinementMethod"]
        num_maxiterate = config["alignment"]["maxiterate"]
        
    except KeyError as missing_key:
        raise KeyError(f"Missing key in configuration: {missing_key}")
    
    for input_file in [filtered_path, outgroup_path]:
        if not Path(input_file).is_file():
            raise FileNotFoundError(f"Input file not found: {input_file}")
        
    Path(combined_path).parent.mkdir(parents=True, exist_ok=True)
    Path(aligned_path).parent.mkdir(parents=True, exist_ok=True)
    
    with open(combined_path, "w") as combined:
        for path in [filtered_path, outgroup_path]:
            with open(path, "r") as f:
                combined.write(f.read())
                combined.write("\n") 
                
    print(f"Combined sequences saved to {combined_path}")
    
    mafft_path = shutil.which("mafft")

    if mafft_path is None:
        raise RuntimeError("MAFFT not found in PATH. Install it with:\n" "  conda install -c bioconda mafft\n" "Or download from: https://mafft.cbrc.jp/alignment/software/")

    result = subprocess.run(["mafft", IterativeRefinementMethod, "--maxiterate", str(num_maxiterate), combined_path], capture_output=True, text=True)
        
    if result.stderr:
        print(result.stderr)
    
    if result.returncode != 0:
        raise RuntimeError(f"MAFFT failed:\n{result.stderr}")
    
    with open(aligned_path, "w") as aligned_file:
        aligned_file.write(result.stdout)
    
    try:
        alignment = AlignIO.read(aligned_path, "fasta")
        print(f"Alignment saved to {aligned_path}: {len(alignment)} sequences, length {alignment.get_alignment_length()}")
        
    except ValueError as parse_error:
        raise ValueError(f"MAFFT finished, but output file is not a valid FASTA file: {parse_error}")
   
    
