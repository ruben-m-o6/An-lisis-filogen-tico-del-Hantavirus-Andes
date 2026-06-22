from src.seqDownload import download_sequences
from src.filter_sequences import filter_sequences
from src.align_seq import align_seq
from src.tree import build_tree
from src.metadata import get_metadata
from src.logger import logger
import argparse
import json




if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Phylogenetic analysis pipeline — Andes orthohantavirus")
    parser.add_argument("--step", choices=["download", "filter", "get_metadata", "align", "tree", "all"], required=True)
    parser.add_argument("--config", type=str, default="dataConfig.json", help="Path to the configuration JSON file (default: dataConfig.json)")
    args = parser.parse_args()
    
    logger.info(f"Starting pipeline with step: {args.step} and with configuration file: {args.config}")
    
    try:
        with open(args.config, "r") as config_file:
            config = json.load(config_file)
    except FileNotFoundError:
        logger.error(f"Configuration file not found: {args.config}")
        exit(1)

    if args.step == "download":
        download_sequences(config)
    elif args.step == "filter":
        filter_sequences(config)
    elif args.step == "get_metadata":
        get_metadata(config)
    elif args.step == "align":
        align_seq(config)
    elif args.step == "tree":
        build_tree(config)
    elif args.step == "all":
        download_sequences(config)
        filter_sequences(config)
        get_metadata(config)
        align_seq(config)
        build_tree(config)
        
    logger.info("Pipeline completed successfully.")