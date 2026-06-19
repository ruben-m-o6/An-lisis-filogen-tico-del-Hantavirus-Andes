from pathlib import Path
import logging

# Create results directory if it doesn't exist
Path("results").mkdir(parents=True, exist_ok=True)

# Create project logger
logger = logging.getLogger("phylo_pipeline")
logger.setLevel(logging.INFO)

if not logger.handlers:
    
    # File Handler for logging to a file
    fh = logging.FileHandler("results/queries.log", encoding="utf-8")
    fh.setLevel(logging.INFO)
    
    # Console Handler for logging to the console
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    
    # Formatter for log messages
    fmt = logging.Formatter("%(asctime)s\t%(name)s\t%(levelname)s\t%(message)s")
    fh.setFormatter(fmt)
    ch.setFormatter(fmt)
    
    # Add handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)