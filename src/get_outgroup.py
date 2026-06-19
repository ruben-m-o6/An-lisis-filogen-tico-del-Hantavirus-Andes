from Bio import Entrez
from urllib.error import HTTPError, URLError
from src.logger import logger

def get_outgroup(ingroup_taxon, target_gene, email, api_key, higher_rank):
    """
    Retrieves the NCBI nucleotide ID of an outgroup sequence for phylogenetic rooting.

    This function queries the NCBI Taxonomy and Nucleotide databases to find a suitable 
    outgroup. It traverses the taxonomic lineage of the provided ingroup up to a specified 
    higher rank (e.g., family, order), and then queries for a complete coding sequence (CDS) 
    of the target gene belonging to that ancestral clade, strictly excluding the ingroup itself.

    Args:
        ingroup_taxon (str): Scientific name of the target organism.
        target_gene (str): The specific gene to search for .
        email (str): Valid email address for NCBI Entrez authentication (required by NCBI).
        api_key (str): NCBI API key to increase the request rate limit (up to 10 requests/sec).
        higher_rank (str): The taxonomic rank used to define the ancestral clade (e.g., "family", "order").

    Returns:
        str: The unique NCBI nucleotide database accession ID (or UID) of the outgroup sequence.

    Raises:
        ValueError: If the `ingroup_taxon` is not found, if the lineage lacks the 
            specified `higher_rank`, or if no valid outgroup sequences are available.
        ConnectionError: If there is a network failure or NCBI servers are unreachable 
            (captures HTTPError and URLError).
        RuntimeError: If Entrez fails to parse the XML response or encounters an 
            unexpected internal error.
    """
    Entrez.email = email
    Entrez.api_key = api_key
    
    logger.info(f"Searching for outgroup for ingroup '{ingroup_taxon}' at rank '{higher_rank}' targeting gene '{target_gene}'")

    try:
        logger.info(f"Querying NCBI Taxonomy for ingroup taxon '{ingroup_taxon}'")
        search = Entrez.esearch(db="taxonomy", term=ingroup_taxon)
        search_results = Entrez.read(search)
        if search_results["Count"] == "0":
            logger.error(f"Could not find the taxon '{ingroup_taxon}' in the NCBI taxonomy database.")
            raise ValueError(f"Could not find the taxon '{ingroup_taxon}' in the NCBI taxonomy database.")
        taxid = search_results["IdList"][0]
        search.close()
        
        records = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        taxon_info = Entrez.read(records)[0]
        records.close()
        lineage = taxon_info["LineageEx"]
        
        ancestor = None
        for t in reversed(lineage):
            if t["Rank"] == higher_rank:
                ancestor = t
                break
            
        if ancestor is None:
            logger.error(f"Could not find an ancestor with rank '{higher_rank}' for taxon '{ingroup_taxon}'.")
            raise ValueError(f"Could not find an ancestor with rank '{higher_rank}' for taxon '{ingroup_taxon}'.")
            
        query_term = f'txid{ancestor["TaxId"]}[Organism] AND {target_gene}[Title] NOT {ingroup_taxon}[Organism] AND complete cds[Title]'
        search2 = Entrez.esearch(db="nucleotide", term=query_term, retmax=1)
        ids = Entrez.read(search2)["IdList"]
        search2.close()
        
        if not ids:
            logger.error(f"Could not find an outgroup in rank '{higher_rank}'.")
            raise ValueError(f"Could not find an outgroup in rank '{higher_rank}'.")
        return ids[0]
    
    except (HTTPError, URLError) as network_error:
        logger.error(f"Connection error with NCBI: {network_error}")
        raise ConnectionError(f"Connection error with NCBI: {network_error}")
    
    except RuntimeError as entrez_error:
        logger.error(f"Error processing Entrez data: {entrez_error}")
        raise RuntimeError(f"Error processing Entrez data: {entrez_error}")
    