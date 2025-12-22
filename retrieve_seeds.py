import argparse
import sys
import requests
import time
from Bio import Entrez, SeqIO
from urllib.error import HTTPError
import io

# Set your email for NCBI Entrez
Entrez.email = "your_email@example.com"  # Replace with a real email or let user provide it

def fetch_interpro_proteins(entry_id, output_file, max_records=None):
    """
    Fetches protein sequences for a given InterPro entry ID by:
    1. Querying InterPro API for protein accessions (handling pagination).
    2. Batch fetching FASTA sequences from UniProt.
    """
    print(f"Fetching InterPro data for {entry_id}...")
    
    # Setup session with retry logic
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry
    
    session = requests.Session()
    retries = Retry(total=5, backoff_factor=1, status_forcelist=[500, 502, 503, 504])
    session.mount('https://', HTTPAdapter(max_retries=retries))
    
    # Auto-detect source database
    source_db = "interpro"
    upper_id = entry_id.upper()
    if upper_id.startswith("IPR"):
        source_db = "interpro"
    elif upper_id.startswith("PF"):
        source_db = "pfam"
    elif upper_id.startswith("CD"):
        source_db = "cdd"
    elif upper_id.startswith("TIGR"):
        source_db = "tigrfams" # or ncbifam
    elif upper_id.startswith("SM"):
        source_db = "smart"
    
    print(f"Detected database type: {source_db}")
    base_url = f"https://www.ebi.ac.uk/interpro/api/protein/uniprot/entry/{source_db}/{entry_id}/"
    headers = {"Accept": "application/json"}
    
    collected_accessions = []
    next_url = base_url + "?page_size=200" # Initial page size
    
    try:
        while next_url:
            if max_records and len(collected_accessions) >= max_records:
                break
                
            print(f"Fetching page: {next_url}...")
            # Use session instead of requests.get
            req = session.get(next_url, headers=headers, timeout=60) 
            if req.status_code != 200:
                print(f"Error fetching InterPro page: {req.status_code}")
                break
                
            data = req.json()
            results = data.get("results", [])
            
            for res in results:
                # Each result typically looks like:
                # {"metadata": {"accession": "A0A...", ...}, ...}
                # But looking at previous grep, accession is top level or in metadata
                # My grep output: "accession":"A0A009HSK7"
                # Let's check typical structure. It's usually result['metadata']['accession'] for entry endpoints
                # BUT for /protein/uniprot/entry/interpro/..., the results ARE protein objects?
                # Let's safely try to find accession
                acc = res.get("metadata", {}).get("accession")
                if not acc:
                     # Fallback if structure differs
                     acc = res.get("accession")
                
                if acc:
                    collected_accessions.append(acc)
            
            next_url = data.get("next")
            # Optional: detailed logging
            print(f"Collected {len(collected_accessions)} accessions so far...")

        # Limit if needed
        if max_records and len(collected_accessions) > max_records:
            collected_accessions = collected_accessions[:max_records]
            
        if not collected_accessions:
            print("No proteins found for this InterPro entry.")
            return

        print(f"Total proteins to retreive: {len(collected_accessions)}")
        
        # Batch fetch from UniProt
        uniprot_batch_url = "https://rest.uniprot.org/uniprotkb/accessions"
        batch_size = 100
        
        with open(output_file, "w") as out_f:
            for i in range(0, len(collected_accessions), batch_size):
                batch_ids = collected_accessions[i : i + batch_size]
                params = {
                    "accessions": ",".join(batch_ids),
                    "format": "fasta"
                }
                print(f"Downloading UniProt batch {i//batch_size + 1}...")
                up_req = session.get(uniprot_batch_url, params=params, timeout=60)
                if up_req.status_code == 200:
                   out_f.write(up_req.text)
                else:
                   print(f"Failed to fetch UniProt batch: {up_req.status_code}")
                   
        print(f"Successfully saved sequences to {output_file}")
            
    except Exception as e:
        print(f"Error during InterPro retrieval: {e}")

def fetch_ncbi_proteins(search_term, output_file, max_records=1000):
    """
    Fetches protein sequences from NCBI Protein database using a search term.
    Includes strict Python-side filtering to exclude XP_ and XM_ (predicted models).
    """
    print(f"Searching NCBI for '{search_term}' (max records: {max_records})...")
    
    try:
        # Step 1: Search
        # We still use the search filter to reduce load, but we don't trust it 100%
        full_term = f"{search_term} AND srcdb_refseq[PROP] NOT srcdb_refseq_model[PROP]"
        
        handle = Entrez.esearch(db="protein", term=full_term, usehistory="y", retmax=max_records)
        record = Entrez.read(handle)
        handle.close()
        
        count = int(record["Count"])
        if count == 0:
            print("No records found on NCBI.")
            return

        print(f"Found {count} potential records. Downloading and filtering for verified entries (NP_ only)...")
        
        webenv = record["WebEnv"]
        query_key = record["QueryKey"]
        
        # Step 2: Fetch and Filter
        batch_size = 500
        saved_count = 0
        
        with open(output_file, "w") as out_handle:
            for start in range(0, min(count, max_records), batch_size):
                end = min(count, start + batch_size)
                # print(f"Processing batch {start+1} to {end}...")
                
                try:
                    fetch_handle = Entrez.efetch(
                        db="protein",
                        retstart=start,
                        retmax=batch_size,
                        webenv=webenv,
                        query_key=query_key,
                        rettype="fasta",
                        retmode="text"
                    )
                    data = fetch_handle.read()
                    fetch_handle.close()
                    
                    # Python-side filtering using SeqIO
                    # Parse the FASTA string from this batch
                    for record in SeqIO.parse(io.StringIO(data), "fasta"):
                        # STRICT FILTER: Only allow NP_ prefix
                        # Using startswith("NP_") ensures we only get verified PROTEIN entries.
                        if not record.id.startswith("NP_"):
                            continue
                        
                        SeqIO.write(record, out_handle, "fasta")
                        saved_count += 1
                        
                except HTTPError as e:
                     print(f"Network error during batch download: {e}")
                     continue

        print(f"Successfully saved {saved_count} verified RefSeq (NP_ only) records to {output_file}")

    except Exception as e:
        print(f"Error fetching NCBI data: {e}")

def main():
    parser = argparse.ArgumentParser(description="Retrieve seed protein sequences for HMM generation.")
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--interpro", help="InterPro Entry ID (e.g., IPR000001)")
    group.add_argument("--ncbi", help="NCBI Search Term")
    
    parser.add_argument("--output", required=True, help="Output FASTA file path")
    parser.add_argument("--email", help="Email for NCBI Entrez (optional but recommended)")
    parser.add_argument("--max_records", type=int, default=100000, help="Max records to search (default: 100000)")

    args = parser.parse_args()

    if args.email:
        Entrez.email = args.email

    if args.interpro:
        fetch_interpro_proteins(args.interpro, args.output, args.max_records)
    elif args.ncbi:
        fetch_ncbi_proteins(args.ncbi, args.output, args.max_records)

if __name__ == "__main__":
    main()
