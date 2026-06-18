from pathlib import Path
import subprocess
from Bio import SeqIO
from Bio import Phylo
from matplotlib import pyplot as plt
import shutil

def build_tree(config):
    
    organism = config["Query"]["organism"]
    alignments_path = config["paths"]["alignment"]
    tree_path = config["paths"]["tree"]
    treepng_path = config["paths"]["figure_tree"]
    model = config["tree"]["model"]
    bootstrap = config["tree"]["bootstrap_num"]
    threads = config["tree"]["threads"]
    seed = config["tree"]["seed"]
    
    Path(tree_path).parent.mkdir(parents=True, exist_ok=True)
     
    outgroup_id = None
    for record in SeqIO.parse(alignments_path, "fasta"):
        if record.id.lower().startswith("outgroup_"):
            outgroup_id = record.id
            break    
        
    if outgroup_id is None:
        raise ValueError("No outgroup sequence found in alignment")
    
    iqtree_path = shutil.which("iqtree")

    if iqtree_path is None:
        raise RuntimeError("IQ-TREE not found in PATH. Install it with:\n" "  conda install -c bioconda iqtree\n" "Or download from: https://github.com/iqtree/iqtree/releases")

    result = subprocess.run([iqtree_path, "-s", alignments_path, "-o", outgroup_id, "-m", model, "-bb", str(bootstrap), "-nt", threads, "-seed", str(seed), "-pre", tree_path, "-redo"], capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"IQ-TREE failed:\n{result.stderr}")
    
    tree = Phylo.read(tree_path + ".treefile", "newick")
    Path(treepng_path).parent.mkdir(parents=True, exist_ok=True)

    tree.ladderize()

    fig, ax = plt.subplots(figsize=(16, 12)) 

    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=True)
    plt.title(f"Phylogenetic Tree of {organism} (Rooted)", fontsize=16, pad=20)
    plt.tight_layout()
    plt.savefig(treepng_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
        
    print(f"Tree figure saved to: {treepng_path}")
    