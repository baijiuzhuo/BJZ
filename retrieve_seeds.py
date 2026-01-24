#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==============================================================================
#  Gene Family Identification Pipeline - Seed Retrieval (Quad-Core Edition)
# ==============================================================================

import argparse
import sys
import io
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import time
import os
import random
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import Entrez, SeqIO

# --- Session Setup ---
def get_session():
    """Creates a requests Session with robust retry logic."""
    session = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["HEAD", "GET", "POST"]
    )
    adapter = HTTPAdapter(max_retries=retries)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session

try:
    from pipeline_utils import cluster_sequences
except ImportError:
    import random
    def cluster_sequences(recs, limit):
        if len(recs) > limit: return random.sample(recs, limit)
        return recs

# --- Classification Helper ---
def classify_sequence(rec, source):
    """
    Classifies a sequence record as 'GOLD' or 'SILVER'.
    GOLD: NCBI NP/YP, UniProt Swiss-Prot, InterPro Reviewed
    SILVER: NCBI XP, UniProt TrEMBL, InterPro Unreviewed
    """
    sid = rec.id
    
    if source == "NCBI":
        if sid.startswith("NP_") or sid.startswith("YP_"):
            return "GOLD"
        elif sid.startswith("XP_"):
            return "SILVER"
        else:
            return "SILVER" # WP, etc.
            
    elif source == "INTERPRO":
        pass
    
    return "SILVER"

def get_ncbi_many(queries, email, max_records_per_query=5000, reviewed_only=False, api_key=None):
    """
    Retrieves sequences for multiple queries.
    reviewed_only: If True, restricts to srcdb_refseq_known (NP/YP) to avoid XP (Predicted).
    api_key: NCBI API Key (optional, increases rate limit from 3/s to 10/s)
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
        print("[NCBI] Using API Key (10 requests/sec limit)")
    else:
        print("[NCBI] No API Key (3 requests/sec limit - consider getting one)")
    all_recs = []
    
    for q_idx, q in enumerate(queries):
        if not q or len(q) < 2: continue
        print(f"[NCBI] Searching for '{q}'...")
        
        # Rate limiting: Wait between queries to avoid 429
        if q_idx > 0:
            wait_time = 2 + random.uniform(0, 1)
            print(f"  (Waiting {wait_time:.1f}s to avoid rate limit...)")
            time.sleep(wait_time)
        
        # Retry logic for esearch
        search_retry = 0
        max_search_retries = 5
        ids = []
        
        while search_retry < max_search_retries:
            try:
                # Search RefSeq
                # PROP "srcdb_refseq_known" excludes predicted (XP/XM) models
                prop = "srcdb_refseq_known[PROP]" if reviewed_only else "srcdb_refseq[PROP]"
                term = f"{q} AND {prop}"
                
                handle = Entrez.esearch(db="protein", term=term, retmax=max_records_per_query)
                record = Entrez.read(handle)
                handle.close()
                
                ids = record["IdList"]
                break  # Success
                
            except Exception as e:
                err_str = str(e)
                search_retry += 1
                if "429" in err_str or "Too Many Requests" in err_str:
                    wait = 10 + search_retry * 5 + random.uniform(0, 2)
                    print(f"  [Rate Limited] Waiting {wait:.0f}s before retry {search_retry}/{max_search_retries}...")
                    time.sleep(wait)
                else:
                    print(f"  [Error] {e}")
                    time.sleep(3)
        
        count = len(ids)
        print(f"  > Found {count} entries ({'Curated' if reviewed_only else 'All RefSeq'}).")
        
        if count == 0: continue

        # Download
        batch_size = 400 # Increased slightly for speed, reliance on threads
        session = get_session()
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        
        # Create batches
        batches = [ids[i:i+batch_size] for i in range(0, count, batch_size)]
        
        def fetch_ncbi_batch(batch_ids, batch_idx):
            retry_count = 0
            max_retries = 5
            while retry_count < max_retries:
                try:
                    # Add small delay per batch to be polite
                    time.sleep(0.5 + random.uniform(0, 0.5))
                    
                    payload = {
                        "db": "protein",
                        "id": ",".join(batch_ids),
                        "rettype": "fasta",
                        "retmode": "text",
                        "email": email
                    }
                    if api_key:
                        payload["api_key"] = api_key
                    # Explicit timeout
                    resp = session.post(base_url, data=payload, timeout=60)
                    
                    if resp.status_code == 200 and resp.text.strip():
                        # Parse immediately
                        local_recs = []
                        for s in SeqIO.parse(io.StringIO(resp.text), "fasta"):
                            tag = classify_sequence(s, "NCBI")
                            # Double check filter if reviewed_only
                            if reviewed_only and tag != "GOLD":
                                continue
                            s._classification = tag
                            local_recs.append(s)
                        return local_recs
                    elif resp.status_code == 429:
                        time.sleep(10 + retry_count * 5 + random.uniform(0, 2)) # Longer jitter
                        retry_count += 1
                    else:
                        time.sleep(3)
                        retry_count += 1
                except Exception as e:
                    time.sleep(3)
                    retry_count += 1
            return []

        # Concurrent Execution - reduced workers for rate limiting
        completed = 0
        with ThreadPoolExecutor(max_workers=3) as executor: # Reduced to 3 workers
            future_to_batch = {executor.submit(fetch_ncbi_batch, b, i): i for i, b in enumerate(batches)}
            
            for future in as_completed(future_to_batch):
                recs = future.result()
                if recs:
                    all_recs.extend(recs)
                
                completed += 1
                percent = (completed / len(batches)) * 100
                print(f"    progress: {completed}/{len(batches)} batches ({percent:.1f}%)...", end="\r")
                sys.stdout.flush()
        print("")
            
    return all_recs

def get_interpro_seeds(ipr_id, max_records=5000, taxid=None, reviewed_only=False):
    """
    Retrieves seeds from InterPro.
    reviewed_only: If True, skips Unreviewed (Silver) retrieval.
    """
    print(f"\n[InterPro] Retrieving for ID: {ipr_id}")
    session = get_session()
    
    # Determined DB
    db_source = "interpro"
    if ipr_id.startswith("IPR"): db_source = "interpro"
    elif ipr_id.startswith("cd") or ipr_id.startswith("sd"): db_source = "cdd"
    elif ipr_id.startswith("PF"): db_source = "pfam"

    results = []
    
    # 1. Gold (Reviewed)
    print(f"  [Gold] Searching Swiss-Prot (Reviewed)...")
    gold_seqs = fetch_interpro_endpoint(session, db_source, ipr_id, "protein/reviewed", taxid, max_records)
    for s in gold_seqs:
        srec = SeqIO.read(io.StringIO(s), "fasta")
        srec._classification = "GOLD"
        results.append(srec)
    print(f"    > Retrieved {len(gold_seqs)} Gold sequences.")

    # 2. Silver (Unreviewed) - SKIP if reviewed_only
    if not reviewed_only:
        remaining = max_records - len(gold_seqs)
        if remaining > 0:
            print(f"  [Silver] Searching TrEMBL (Unreviewed)...")
            silver_seqs = fetch_interpro_endpoint(session, db_source, ipr_id, "protein/unreviewed", taxid, remaining)
            for s in silver_seqs:
                srec = SeqIO.read(io.StringIO(s), "fasta")
                srec._classification = "SILVER"
                results.append(srec)
            print(f"    > Retrieved {len(silver_seqs)} Silver sequences.")
        
    return results

def get_uniprot_seeds_by_name(query, max_records=5000, reviewed_only=False):
    """
    Retrieves seeds from UniProtKB by text search.
    reviewed_only: If True, skips TrEMBL (Silver) search.
    """
    print(f"\n[UniProt] Text Search for: '{query}'")
    session = get_session()
    
    results = []
    
    # query_str: "(query) AND (reviewed:true)" for Gold
    # query_str: "(query) AND (reviewed:false)" for Silver
    
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    
    # 1. Gold (Swiss-Prot)
    try:
        print(f"  [Gold] Searching Swiss-Prot for '{query}'...")
        params = {
            "query": f"{query} AND (reviewed:true)",
            "format": "fasta",
            "size": 500 # Page size
        }
        # Iterate pages?
        # UniProt pagination is via Link header or cursor.
        # Simple implementation: fetch first few pages up to limit.
        # Actually standard search API returns streamed FASTA if format=fasta.
        # But for large result sets, stream is best. 
        # We will set a reasonable size limit.
        params["size"] = min(max_records, 2000) # One page should be enough for Gold usually
        
        r = session.get(base_url, params=params, timeout=60)
        if r.status_code == 200:
            count = 0
            for rec in SeqIO.parse(io.StringIO(r.text), "fasta"):
                rec._classification = "GOLD"
                results.append(rec)
                count += 1
                if count >= max_records: break
            print(f"    > Found {count} Swiss-Prot entries.")
            
    except Exception as e:
        print(f"  [Error] Swiss-Prot search failed: {e}")
        
    # 2. Silver (TrEMBL) - SKIP if reviewed_only
    if not reviewed_only:
        try:
            remaining = max_records - len(results) # Global limit per query term?
            if remaining > 0:
                print(f"  [Silver] Searching TrEMBL for '{query}'...")
                params = {
                    "query": f"{query} AND (reviewed:false)",
                    "format": "fasta",
                    "size": min(remaining, 500) # Just fetch top 500 for variety, enable more if needed
                }
                r = session.get(base_url, params=params, timeout=60)
                if r.status_code == 200:
                    count = 0
                    for rec in SeqIO.parse(io.StringIO(r.text), "fasta"):
                        rec._classification = "SILVER"
                        results.append(rec)
                        count += 1
                    print(f"    > Found {count} TrEMBL entries.")
                    
        except Exception as e:
            print(f"  [Error] TrEMBL search failed: {e}")
        
    return results

def fetch_interpro_endpoint(session, db, acc, endpoint_type, taxid, limit):
    """Helper to fetch from specific InterPro API endpoint."""
    url = f"https://www.ebi.ac.uk/interpro/api/{endpoint_type}/entry/{db}/{acc}/"
    if taxid: url += f"taxonomy/uniprot/{taxid}/"
    
    headers = {"Accept": "application/json"}
    accessions = []
    next_url = url
    
    # 1. Collect IDs
    while next_url and len(accessions) < limit:
        try:
            resp = session.get(next_url, headers=headers, timeout=30)
            if resp.status_code != 200: break
            data = resp.json()
            for res in data.get("results", []):
                m = res.get("metadata", {})
                if "accession" in m: accessions.append(m["accession"])
            next_url = data.get("next")
        except:
            break
            
    # 2. Download Fasta (UniProt)
    if not accessions: return []
    
    # 2. Download Fasta (UniProt)
    if not accessions: return []
    
    print(f"    > Fetching {len(accessions)} sequences from UniProt (Parallel)...")
    final_fasta_strings = []
    batch_size = 100 # UniProt batch size
    batches = [accessions[i:i+batch_size] for i in range(0, len(accessions), batch_size)]
    
    def fetch_uniprot_batch(batch_ids):
        try:
            u_url = "https://rest.uniprot.org/uniprotkb/accessions"
            params = {"accessions": ",".join(batch_ids), "format": "fasta"}
            # Use separate session or just requests
            r = requests.get(u_url, params=params, timeout=45) 
            if r.status_code == 200:
                # Raw text is sufficient, we split by '>' to help caller append lists
                # But here we just return the raw string block
                return r.text
        except:
            pass
        return ""

    completed = 0
    # Boosted to 20 threads to saturate bandwidth for large datasets
    with ThreadPoolExecutor(max_workers=20) as executor:
        futures = {executor.submit(fetch_uniprot_batch, b): b for b in batches}
        
        for future in as_completed(futures):
            data = future.result()
            if data:
                parts = data.strip().split(">")
                for p in parts:
                    if not p: continue
                    final_fasta_strings.append(">" + p)
            
            completed += 1
            print(f"      Progress: {completed}/{len(batches)} batches ({(completed/len(batches)*100):.1f}%)...", end="\r")
    print("")
            
    return final_fasta_strings

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--name_full", help="Full query name (e.g. GDSL lipase)")
    parser.add_argument("--name_abbr", help="Abbreviated name (e.g. GELP)")
    parser.add_argument("--interpro", help="InterPro ID")
    parser.add_argument("--email", required=True)
    parser.add_argument("--out_prefix", required=True, help="Prefix for outputs (e.g. ./tmp/gelp)")
    parser.add_argument("--max_seeds", type=int, default=10000)
    parser.add_argument("--taxid")
    
    parser.add_argument("--reviewed_only", action="store_true", help="Only download reviewed/curated seeds")
    parser.add_argument("--api_key", help="NCBI API Key (optional, increases rate limit)")
    
    args = parser.parse_args()
    
    queries = []
    if args.name_full: queries.append(args.name_full)
    if args.name_abbr: queries.append(args.name_abbr)
    
    all_records = []
    
    # 1. NCBI (Multi-Query)
    if queries:
        # Pass reviewed_only flag deeply
        all_records.extend(get_ncbi_many(queries, args.email, args.max_seeds, args.reviewed_only, args.api_key))
        
    # 2. UniProt Text Search (Mixed)
    if queries:
        for q in queries:
            if not q or len(q) < 3: continue 
            # Pass reviewed_only flag deeply
            all_records.extend(get_uniprot_seeds_by_name(q, args.max_seeds, args.reviewed_only))
        
    # 3. InterPro (Reviewed + Unreviewed)
    if args.interpro:
        # Pass reviewed_only flag deeply
        all_records.extend(get_interpro_seeds(args.interpro, args.max_seeds, args.taxid, args.reviewed_only))
        
    # 4. Deduplicate & Separate
    gold_recs = []
    silver_recs = []
    seen = set()
    
    print("\n[Merge] Deduplicating and Stratifying...")
    for rec in all_records:
        seq_str = str(rec.seq).upper()
        if seq_str in seen: continue
        seen.add(seq_str)
        
        tag = getattr(rec, "_classification", "SILVER")
        if tag == "GOLD":
            gold_recs.append(rec)
        else:
            silver_recs.append(rec)
            
    # 4. Write Outputs
    # Gold File: Pure Gold
    gold_out = args.out_prefix + "_seeds_gold.fasta"
    SeqIO.write(gold_recs, gold_out, "fasta")
    
    # Broad File: Gold + Silver
    broad_recs = gold_recs + silver_recs
    broad_out = args.out_prefix + "_seeds_broad.fasta"
    
    # Limit Broad if huge
    if len(broad_recs) > args.max_seeds:
         broad_recs = cluster_sequences(broad_recs, args.max_seeds)
         
    SeqIO.write(broad_recs, broad_out, "fasta")
    
    print("-" * 50)
    print(f"[Summary] Quad-Core Seed Retrieval")
    print(f"  - GOLD Set (Strict): {len(gold_recs)} (Saved to {os.path.basename(gold_out)})")
    if not args.reviewed_only:
        print(f"  - BROAD Set (Mixed): {len(broad_recs)} (Saved to {os.path.basename(broad_out)})")
    else:
        print(f"  - BROAD Set (Skimmed): 0 (Skipped Silver)")
    print("-" * 50)

if __name__ == "__main__":
    main()
