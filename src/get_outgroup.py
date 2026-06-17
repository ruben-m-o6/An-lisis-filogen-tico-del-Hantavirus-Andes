from Bio import Entrez


def get_outgroup(ingroup_taxon, target_gene, email, higher_rank):
    
    Entrez.email = email
    
    search = Entrez.esearch(db="taxonomy", term=ingroup_taxon)
    taxid = Entrez.read(search)["IdList"][0]
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
        raise ValueError(f"Could not find an ancestor with rank '{higher_rank}' for taxon '{ingroup_taxon}'.")
        
    query_term = f'txid{ancestor["TaxId"]}[Organism] AND {target_gene}[Title] NOT {ingroup_taxon}[Organism] AND complete cds[Title]'
    search2 = Entrez.esearch(db="nucleotide", term=query_term, retmax=1)
    ids = Entrez.read(search2)["IdList"]
    search2.close()
    
    if not ids:
        raise ValueError(f"Could not find an outgroup in rank '{higher_rank}'.")
    return ids[0]