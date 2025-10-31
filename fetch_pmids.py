from Bio import Entrez
import os

Entrez.email = "rafael@siensmetrica.com"

def fetch_pmids(query, start, count):
    query_with_filter = f"({query}) AND free full text[sb]"
    handle = Entrez.esearch(
        db="pubmed",
        term=query_with_filter,
        retstart=start,
        retmax=count,
        sort="relevance"
    )
    record = Entrez.read(handle)
    handle.close()
    pmids = record["IdList"]
    return pmids

if __name__ == "__main__":
    query = input("Enter topic: ")
    start_index = int(input("Enter start index (e.g., 0 for beginning): "))
    max_pmids = int(input("Enter number of PMIDs to fetch: "))
    file_path = "pmids.txt"

    print(f"Fetching PMIDs {start_index} to {start_index + max_pmids} for '{query}'...")

    pmids = fetch_pmids(query, start=start_index, count=max_pmids)

    with open(file_path, "w") as f:
        f.write(",".join(pmids))

    print(f"✅ Fetched {len(pmids)} PMIDs.")
    print(f"Saved to {file_path}")
