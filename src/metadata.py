from Bio import SeqIO
from src.logger import logger
from Bio import Entrez
from pathlib import Path
import csv
from urllib.error import HTTPError
import time

def get_metadata(config):
    """
    Modifies the filtered sequences, and adds metadata to the sequence headers. This function 
    is designed to be used in a pipeline where sequences are filtered based on certain criteria, 
    and metadata is added to the sequence headers for downstream analysis. It also generates 
    a CSV file summarizing the extracted metadata.
    
    Args:
        config (dict): A configuration dictionary containing Entrez credentials and 
            file paths. It must adhere to the following structure:
            {
                "paths": {
                    "filtered_sequences": str, 
                    "annotated_sequences": str,
                    "csv_file": str
                },
                "Entrez": {
                    "email": str,
                    "api_key": str
                },
                "Query": {
                    "database": str
                }
            }

    Returns:
        None. The function writes the annotated FASTA sequences and the metadata CSV 
        directly to the file system paths specified in the `config` dictionary.

    Raises:
        KeyError: If any required nested keys are missing from the `config` dictionary.
        FileNotFoundError: If the input FASTA file specified in paths cannot be found.
        Exception: If there are network issues with Entrez or I/O errors during file saving.
    """
    
    logger.info("Starting metadata addition process.")
    
    try:
        filtered_sequences_path = config["paths"]["filtered_sequences"]
        annotated_sequences_path= config["paths"]["annotated_sequences"]
        csv_path= config["paths"]["csv_file_metadata"]
        Entrez_email=config["Entrez"]["email"]
        Entrez_api_key=config["Entrez"]["api_key"]
        database=config["Query"]["database"]
        
    except KeyError as missing_key  :
        logger.error(f"Missing configuration key: {missing_key}")
        raise KeyError(f"Missing configuration key: {missing_key}")

    Entrez.email=Entrez_email
    Entrez.api_key=Entrez_api_key
    
    Path(annotated_sequences_path).parent.mkdir(parents=True, exist_ok=True)
    Path(csv_path).parent.mkdir(parents=True, exist_ok=True)
    
    dataCSV = [["Accession", "Organismo", "Fecha", "País", "Hospedador", "Aislado"]]
    
    try:
        sequences = list(SeqIO.parse(filtered_sequences_path, "fasta"))
        
    except FileNotFoundError:
        logger.error(f"Filtered sequences file not found: {filtered_sequences_path}")
        raise FileNotFoundError(f"Filtered sequences file not found: {filtered_sequences_path}")
    
    for record in sequences:
        id= record.id
        for attempt in range(3):
            try:
                handle = Entrez.efetch(db=database, id= id, rettype="gb", retmode="text")
                record_handle = SeqIO.read(handle, "genbank")
                handle.close()
                break  
                
            except HTTPError as e:
                if e.code in [429, 503]:
                    sleep_time = (attempt + 1) * 3  
                    logger.warning(f"NCBI rate limit hit or server busy for {id}. Retrying in {sleep_time}s (Attempt {attempt + 1}/{3})")
                    time.sleep(sleep_time)
                else:
                    logger.error(f"HTTP Error {e.code} for sequence {id}: {e}")
                    raise e
            except Exception as e:
                logger.error(f"Unexpected network error for sequence {id}: {e}")
                raise e
        
        if record_handle is None:
            logger.error(f"Failed to fetch metadata for {id} after 3 attempts. Skipping sequence.")
            continue
        
        date= "Unknown"
        country= "Unknown"
        host= "Unknown"
        isolate= "Unknown"
        organism= "Unknown"
        
        for feature in record_handle.features:
            if feature.type == "source":
                if "collection_date" in feature.qualifiers:
                    date= feature.qualifiers["collection_date"][0].replace(" ", "_")
                if "geo_loc_name" in feature.qualifiers:
                    country= feature.qualifiers["geo_loc_name"][0].replace(" ", "_")
                if "host" in feature.qualifiers:
                    host= feature.qualifiers["host"][0].replace(" ", "_")
                if "isolate" in feature.qualifiers:
                    isolate= feature.qualifiers["isolate"][0].replace(" ", "_")
                if "organism" in feature.qualifiers:
                    organism= feature.qualifiers["organism"][0].replace(" ", "_")
                    
        newID= f"{id}|{organism}|{date}|{country}|{host}|{isolate}"
        record.id=newID
        record.name = ""
        record.description = ""
        
        dataCSV.append([id, organism, date, country, host, isolate])
        logger.info(f"Metadata for: {id} -> {newID}")
        
    try:
        SeqIO.write(sequences, annotated_sequences_path, "fasta")
        logger.info(f"FASTA file with metadata successfully saved to: {annotated_sequences_path}")
        
    except Exception as e:
        logger.error(f"Error saving the FASTA file: {e}")
        raise
        
    try:
        with open(csv_path, mode="w", newline="", encoding="utf-8") as csvFile:
            writer=csv.writer(csvFile)
            writer.writerows(dataCSV)
        logger.info(f"csv file saved successfully in: {csv_path}")
            
    except Exception as e:
        logger.error(f"Error saving the csv file: {e}")
        raise e
        
        