from Bio import Entrez
from src.get_outgroup import get_outgroup
import json

def download_sequences(config):
    
    email = config["Entrez"]["email"]
    api_key = config["Entrez"]["api_key"]
    database = config["Query"]["database"]
    organism = config["Query"]["organism"]
    gene = config["Query"]["gene"]
    outgroup_rank = config["Query"]["outgroup_rank"]

    Entrez.email = email
    Entrez.api_key = api_key

    query_term = f'"{organism}"[Organism] AND {gene}[Title]'
    

    search1 = Entrez.esearch(db=database, term=query_term, retmax=0)
    total_count = int(Entrez.read(search1)['Count'])
    search1.close()


    search2 = Entrez.esearch(db=database, term=query_term, retmax=total_count)
    ids = Entrez.read(search2)['IdList']
    search2.close()

    records = Entrez.efetch(db=database, id=ids, rettype="fasta", retmode="text")

    record_count = len(ids)

    with open(config["paths"]["raw_sequences"], "w") as file:
        file.write(records.read())
    records.close()

    print(f"Sequences downloaded and saved to {config['paths']['raw_sequences']}: {record_count}")

    outgroup_id = get_outgroup(organism, gene, email, outgroup_rank)

    if outgroup_id:
        print(f"Outgroup ID found: {outgroup_id}")
        outgroup_record = Entrez.efetch(db=database, id=outgroup_id, rettype="fasta", retmode="text")
        fasta_text = outgroup_record.read()
        modified_fasta = fasta_text.replace(">", ">outgroup_", 1)
        
        with open(config["paths"]["outgroup_sequence"], "w") as outgroup_file:
            outgroup_file.write("\n" + modified_fasta)
        outgroup_record.close()
    else:
        print("Could not find a suitable outgroup for the phylogenetic analysis.")
