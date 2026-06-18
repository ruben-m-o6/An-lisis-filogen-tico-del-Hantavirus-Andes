import subprocess
from Bio import AlignIO

def align_seq(config):
    
    filtered_path = config["paths"]["filtered_sequences"]
    outgroup_path = config["paths"]["outgroup_sequence"]
    combined_path = config["paths"]["combined_sequences"]
    aligned_path = config["paths"]["alignment"]
    IterativeRefinementMethod = config["alignment"]["IterativeRefinementMethod"]
    num_maxiterate = config["alignment"]["maxiterate"]
    
    with open(combined_path, "w") as combined:
        for path in [filtered_path, outgroup_path]:
            with open(path, "r") as f:
                combined.write(f.read())
                
    print(f"Combined sequences saved to {combined_path}")
    
    result= subprocess.run(["mafft", IterativeRefinementMethod, "--maxiterate", str(num_maxiterate), combined_path], capture_output=True, text=True)
    
    if result.stderr:
        print(result.stderr)
    
    if result.returncode != 0:
        raise RuntimeError(f"MAFFT failed:\n{result.stderr}")
    
    with open(aligned_path, "w") as aligned_file:
        aligned_file.write(result.stdout)
        
    alignment = AlignIO.read(aligned_path, "fasta")
    print(f"Alignment saved to {aligned_path}: {len(alignment)} sequences, length {alignment.get_alignment_length()}")
    
