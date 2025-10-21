import requests
import pandas as pd
import time
from Bio import Entrez

Entrez.email = "rafael@siensmetrica.com"
API_BASE = "http://3.81.212.180:8000"

def fetch_pmids(query, max_pmids):
    handle = Entrez.esearch(
        db="pubmed",
        term=f"{query} AND free full text[sb]",
        retmax=max_pmids,
        sort="relevance"
    )
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_score(url):
    try:
        r = requests.get(url)
    except Exception as e:
        print(f"Connection error to {url}: {e}")
        return "timeout"

    if r.status_code == 200:
        print(f"{url} successfully connected")
        data = r.json()
        if "task" in data:
            task_id = data["task"]
            task_url = f"{API_BASE}/task/{task_id}/wait"
            print(f"waiting on task {task_url}")
            
            try:
                r = requests.get(task_url, timeout=1800)  #30 minutes, this stops control flow
            except requests.Timeout:
                print(f"Request to {url} timed out after 30 minutes.")
                requests.get(cancel_url)
                print("task cancelled")
                return "timeout"

            if r.status_code == 200:
                print(f"{task_url} successfully connected")
                task_data = r.json()
                return task_data.get("result")
            else:
                print(f"task failed to connect (status {r.status_code})")
                cancel_url = f"{API_BASE}/task/{task_id}/cancel"
                requests.get(cancel_url)
                print("task cancelled")
                return "timeout"

        return data.get("result", data)

    else:
        print(f"failed to connect to {url} (status {r.status_code})")
        return "timeout"

def get_article_info(pmid):
    url = f"{API_BASE}/article/{pmid}"
    data = fetch_score(url)
    if data == "timeout":
        return {"title": "timeout", "journal_issn": None}
    return data or {}

def get_iso_abv(issn):
    if not issn:
        return None
    url = f"{API_BASE}/journal/{issn}"
    data = fetch_score(url)
    if data == "timeout":
        return None
    return (data or {}).get("iso_abv")

def get_author_score(pmid):
    url = f"{API_BASE}/article/{pmid}/author_scores"
    data = fetch_score(url)
    if data == "timeout" or data == None:
        return "timeout"
    return (data or {}).get("collective_all_score")

def get_journal_score(iso_abv):
    if not iso_abv:
        return "timeout"
    url = f"{API_BASE}/journal/{iso_abv}/impact"
    data = fetch_score(url)
    if data == "timeout" or data == None:
        return "timeout"
    raw_result = (data or {}).get("score")
    if raw_result is None:
        return "timeout"
    return raw_result * 3.25

def get_media_score(pmid):
    url = f"{API_BASE}/article/{pmid}/media"
    data = fetch_score(url)
    if data == "timeout" or data == None:
        return "timeout"
    return data

def get_authentic_score(pmid): #task will not connect and return timeout when not cached
    url = f"{API_BASE}/article/{pmid}/copyleaks"
    data = fetch_score(url)
    if data == "timeout" or data == None:
        return "timeout"

    result = (data or {}).get("results", data or {})
    if isinstance(result, list) and len(result) > 1:
        return result[1].get("probability")
    elif isinstance(result, list) and len(result) > 0:
        return result[0].get("probability")
    return "timeout"

def get_methods_score(pmid): #task will not connect and return timeout when not cached
    url = f"{API_BASE}/article/{pmid}/analysis"
    data = fetch_score(url)
    if data == "timeout":
        return "timeout"

    rigor_score = data.get("rigor_score") or 0
    empiricity_score = data.get("empiricity_score") or 0
    novelty_score = data.get("novelty_score") or 0
    return (rigor_score + empiricity_score + novelty_score) / 3

def get_overall_score(pmid): #when not cached, data is {'elapsed_ms': 31807} despite succesful wait
    url = f"{API_BASE}/article/{pmid}/final"
    data = fetch_score(url)
    if data == "timeout" or data == None:
        return "timeout"
    if isinstance(data, dict):
        return data.get("score")
    return data

def get_license(pmid):
    url = f"{API_BASE}/dataset/access?pmids={pmid}"
    data = fetch_score(url)
    if data == "timeout":
        return "NA"
    if data.get("hits") and len(data["hits"]) > 0:
        return data["hits"][0].get("License", "No license info")
    return "NA"

def process_pmid(pmid):
    print(f"\n{'='*50}\nFetching data for PMID {pmid}...\n")


    license = get_license(pmid)
    '''
    if license != "CC BY":
        print(f"PMID {pmid} appears to be paywalled or unavailable. Skipping.")
        return {
            "PMID": pmid,
            "Title": "paywall",
            "Author": "paywall",
            "Journal": "paywall",
            "Media": "paywall",
            "Authentic": "paywall",
            "Methods": "paywall",
            "Overall": "paywall"
        }
    '''

    article_info = get_article_info(pmid)
    title = article_info.get("title", "N/A")
    issn = article_info.get("journal_issn")
    iso_abv = get_iso_abv(issn)

    author = get_author_score(pmid)
    journal = get_journal_score(iso_abv)
    media = get_media_score(pmid)

    authentic = get_authentic_score(pmid)
    methods = get_methods_score(pmid)
    overall = get_overall_score(pmid)

    #call again because it only works when in cache
    authentic = get_authentic_score(pmid)
    methods = get_methods_score(pmid)
    overall = get_overall_score(pmid)
    print(f"license: {license}")

    return {
        "PMID": pmid,
        "Title": title,
        "Author": author,
        "Journal": journal,
        "Media": media,
        "Authentic": authentic,
        "Methods": methods,
        "Overall": overall
    }

if __name__ == "__main__":
    print("Select input mode:")
    print("1. Paste PMIDs manually")
    print("2. Search PMIDs by topic")
    mode = input("Enter 1 or 2: ").strip()

    if mode == "1":
        pasted = input("Paste PMIDs separated by commas: ").strip()
        PMIDs = [pmid.strip() for pmid in pasted.split(",") if pmid.strip()]
        print(f"Loaded {len(PMIDs)} PMIDs from input.")

        if not PMIDs:
            print("No valid PMIDs provided. Exiting.")
            exit()

    elif mode == "2":
        query = input("Enter topic: ").strip()
        max_pmids = input("Enter number of PMIDs to fetch: ").strip()
        pmids = fetch_pmids(query, max_pmids)
        print(f"Fetched {len(pmids)} PMIDs for query: {query}")

        if not pmids:
            print("Couldn't find any PMIDs for this topic.")
            exit()

        PMIDs = [pmid.strip() for pmid in pmids if pmid.strip()]

    else:
        print("Invalid option. Exiting.")
        exit()

    start_time = time.time()
    all_results = []
    total_pmids = len(PMIDs)

    for i, pmid in enumerate(PMIDs, start=1):
        print(f"\n=== [{i}/{total_pmids}] Processing PMID {pmid} ===")
        result = process_pmid(pmid)
        all_results.append(result)

    end_time = time.time()
    total_time = end_time - start_time
    hours, remainder = divmod(total_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    print(f"Time elapsed: {int(hours)}h {int(minutes)}m {int(seconds)}s")

    df = pd.DataFrame(all_results)

    # Safely convert 'Overall' to numeric; non-numeric values become NaN
    df["Overall"] = pd.to_numeric(df["Overall"], errors="coerce")

    # Sort descending by Overall score, NaN (timeouts/paywalls) go to the bottom
    df = df.sort_values(by="Overall", ascending=False, na_position="last")

    df.to_csv("article_scores.csv", index=False, encoding="utf-8")
    print("\nall results saved to 'article_scores.csv'")
    if mode == 2:
        print(f"topic was: {query}")