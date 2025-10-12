import requests
import pandas as pd

API_BASE = "http://3.81.212.180:8000"

def fetch_score(url):
    r = requests.get(url)
    if r.status_code == 200:
        print(f"{url} successfully connected")
        data = r.json()
        if "task" in data:
            task_id = data["task"]
            task_url = f"{API_BASE}/task/{task_id}/wait"
            print(f"waiting on {task_url}")
            r = requests.get(task_url)
            if r.status_code == 200:
                print(f"{task_url} successfully connected")
                task_data = r.json()
                return task_data.get("result")
        return data.get("result")
    else:
        print(f"failed to connect to {url} (status {r.status_code})")
        return {}

def get_article_info(pmid):
    url = f"{API_BASE}/article/{pmid}"
    return fetch_score(url) or {}

def get_iso_abv(issn):
    if not issn:
        return None
    url = f"{API_BASE}/journal/{issn}"
    result = fetch_score(url) or {}
    return result.get("iso_abv")

def get_author_score(pmid):
    url = f"{API_BASE}/article/{pmid}/author_scores"
    result = fetch_score(url) or {}
    return result.get("collective_score")

def get_journal_score(iso_abv):
    if not iso_abv:
        return None
    url = f"{API_BASE}/journal/{iso_abv}/impact"
    result = fetch_score(url) or {}
    return result.get("score")

def get_media_score(pmid):
    url = f"{API_BASE}/article/{pmid}/media"
    result = fetch_score(url) or {}
    return result.get("score") if isinstance(result, dict) else result

def get_authentic_score(pmid):
    url = f"{API_BASE}/article/{pmid}/copyleaks"
    result = fetch_score(url)
    if isinstance(result, list) and len(result) > 1:
        return result[1].get("probability")
    elif isinstance(result, list) and len(result) > 0:
        return result[0].get("probability")
    return None

def get_methods_score(pmid):
    url = f"{API_BASE}/article/{pmid}/analysis"
    result = fetch_score(url) or {}
    return result.get("rigor_score")

def get_overall_score(pmid):
    url = f"{API_BASE}/article/{pmid}/final"
    result = fetch_score(url) or {}
    return result.get("score") if isinstance(result, dict) else result

def process_pmid(pmid):
    print(f"\n{'='*50}\nFetching data for PMID {pmid}...\n")

    article_info = get_article_info(pmid)
    title = article_info.get("title", "N/A")
    issn = article_info.get("journal_issn")
    iso_abv = get_iso_abv(issn)

    author = get_author_score(pmid) or "timeout"
    journal = get_journal_score(iso_abv) or "timeout"
    media = get_media_score(pmid) or "timeout"
    authentic = get_authentic_score(pmid) or "timeout"
    methods = get_methods_score(pmid) or "timeout"
    overall = get_overall_score(pmid) or "timeout"

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
    pmid_input = input("Paste PMIDs separated by commas:\n> ")
    PMIDs = [pmid.strip() for pmid in pmid_input.split(",") if pmid.strip()]

    all_results = []
    for pmid in PMIDs:
        result = process_pmid(pmid)
        all_results.append(result)

    df = pd.DataFrame(all_results)
    df.to_csv("article_scores.csv", index=False, encoding="utf-8")

    print("\nall results saved to 'article_scores.csv'")
