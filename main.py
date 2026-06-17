from src.seqDownload import download_sequences
from src.filter_sequences import filter_sequences
from src.align_seq import align_seq
import argparse
import json




if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Phylogenetic analysis pipeline — Andes orthohantavirus")
    parser.add_argument("--step", choices=["download", "filter", "align", "tree", "all"], required=True)
    args = parser.parse_args()
    
    with open("dataConfig.json", "r") as config_file:
        config = json.load(config_file)ali
    
    if args.step == "download":
        download_sequences(config)
    elif args.step == "filter":
        filter_sequences(config)
    elif args.step == "align":
        align_seq(config)
    elif args.step == "tree":
        