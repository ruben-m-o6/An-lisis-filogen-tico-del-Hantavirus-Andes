from Bio import Entrez
from src.get_outgroup import get_outgroup
import json
from urllib.error import HTTPError, URLError
from pathlib import Path

def download_sequences(config):
    """
    Downloads raw target sequences and a corresponding outgroup sequence from NCBI.

    This function authenticates with NCBI Entrez, queries the specified database for 
    the target organism and gene, and safely downloads the matching records in FASTA format. 
    It also invokes the `get_outgroup` helper function to fetch a suitable outgroup sequence, 
    modifying its FASTA header for easy identification, and saves both to disk.

    Args:
        config (dict): A configuration dictionary containing credentials, query parameters, 
            and output paths. It must adhere to the following structure:
            {
                "Entrez": {"email": str, "api_key": str},
                "Query": {"database": str, "organism": str, "gene": str, "outgroup_rank": str},
                "paths": {"raw_sequences": str, "outgroup_sequence": str}
            }

    Returns:
        None. The function writes the downloaded FASTA sequences directly to the file 
        system paths specified in the `config` dictionary.

    Raises:
        KeyError: If any required nested keys are missing from the `config` dictionary.
        ConnectionError: If there is a network failure or NCBI servers are unreachable.
        RuntimeError: If Entrez fails to parse the XML response or an internal error occurs.
    """
    try:
        email = config["Entrez"]["email"]
        api_key = config["Entrez"]["api_key"]
        database = config["Query"]["database"]
        organism = config["Query"]["organism"]
        gene = config["Query"]["gene"]
        outgroup_rank = config["Query"]["outgroup_rank"]
        raw_sequences = config["paths"]["raw_sequences"]
        outgroup_sequence = config["paths"]["outgroup_sequence"]
        
    except KeyError as missing_key:
        raise KeyError(f"Error de configuración. Falta la clave: {missing_key}")

    Entrez.email = email
    Entrez.api_key = api_key

    query_term = f'"{organism}"[Organism] AND {gene}[Title]'
    
    Path(raw_sequences).parent.mkdir(parents=True, exist_ok=True)
    Path(outgroup_sequence).parent.mkdir(parents=True, exist_ok=True)
    
    try:
        search1 = Entrez.esearch(db=database, term=query_term, retmax=0)
        total_count = int(Entrez.read(search1)['Count'])
        search1.close()


        search2 = Entrez.esearch(db=database, term=query_term, retmax=total_count)
        ids = Entrez.read(search2)['IdList']
        search2.close()

        records = Entrez.efetch(db=database, id=ids, rettype="fasta", retmode="text")

        record_count = len(ids)

        with open(raw_sequences, "w") as file:
            file.write(records.read())
        records.close()

        print(f"Sequences downloaded and saved to {raw_sequences}: {record_count}")

        outgroup_id = get_outgroup(organism, gene, email, api_key, outgroup_rank)

        if outgroup_id:
            print(f"Outgroup ID found: {outgroup_id}")
            outgroup_record = Entrez.efetch(db=database, id=outgroup_id, rettype="fasta", retmode="text")
            fasta_text = outgroup_record.read()
            modified_fasta = fasta_text.replace(">", ">outgroup_", 1)
            
            with open(outgroup_sequence, "w") as outgroup_file:
                outgroup_file.write("\n" + modified_fasta)
            outgroup_record.close()
        else:
            print("Could not find a suitable outgroup for the phylogenetic analysis.")
            
    except (HTTPError, URLError) as network_error:
        raise ConnectionError(f"Failure connecting to NCBI during download: {network_error}")
    except RuntimeError as entrez_error:
        raise RuntimeError(f"Failure processing Entrez data: {entrez_error}")
