import requests

API_BASE = "http://3.81.212.180:8000"
PMID = "36321004"

def fetch_score(url, key=None):
    """Fetch a score from an endpoint, using /task/{task_id}/wait for task IDs."""
    try:
        r = requests.get(url)
        if r.status_code == 200:
            print(f"{url} successfully connected")
            data = r.json()
            #print("this raw data : \n")
            #print(data)
            if not isinstance(data, dict):
                return None
            if "task" in data:
                print(f"waiting on {url} ")
                task_id = data["task"]
                task_url = f"{API_BASE}/task/{task_id}/wait"
                print(task_url)
                r = requests.get(task_url)
                if r.status_code == 200:
                    print(f"{task_url} successfully connected")
                    task_data = r.json()
                    print(task_data)
                    if isinstance(task_data, dict):
                        # Check if task/wait provides the result directly
                        if "result" in task_data and not task_data["result"].startswith("task_"):
                            result = task_data["result"]
                        elif task_data.get("status") == "completed":
                            # Re-fetch the original endpoint
                            r = requests.get(url)
                            if r.status_code == 200:
                                data = r.json()
                                result = data.get("result")
                            else:
                                return None
                        else:
                            return None  # Task not completed
                    else:
                        return None
                else:
                    return None
            else:
                result = data.get("result")
            # Process the final result
            if result is None:
                return None
            if isinstance(key, list):
                temp = data
                for k in key:
                    temp = temp.get(k, {})
                return temp if temp != {} else None
            return result
        return None
    except (requests.exceptions.RequestException, ValueError):
        return None

def get_article_info(pmid):
    url = f"{API_BASE}/article/{pmid}"
    result = fetch_score(url)
    return result if result is not None else {"title": "N/A", "journal_issn": None}

def get_iso_abv(issn):
    url = f"{API_BASE}/journal/{issn}"
    result = fetch_score(url)
    return result.get("iso_abv", "N/A") if result is not None else "timeout"

def get_author_score(pmid):
    url = f"{API_BASE}/article/{pmid}/author_scores"
    #return fetch_score(url, key=["result", "collective_score"])
    return fetch_score(url)

def get_journal_score(iso_abv):
    url = f"{API_BASE}/journal/{iso_abv}/impact"
    result = fetch_score(url)
    return result.get("score", "N/A") if result is not None else "timeout"

def get_media_score(pmid):
    url = f"{API_BASE}/article/{pmid}/media"
    return fetch_score(url)

def get_authentic_score(pmid):
    url = f"{API_BASE}/article/{pmid}/copyleaks"
    return fetch_score(url)

def get_methods_score(pmid):
    url = f"{API_BASE}/article/{pmid}/analysis"
    return fetch_score(url, key=["result", "rigor_score"])

def get_overall_score(pmid):
    url = f"{API_BASE}/article/{pmid}/final"
    return fetch_score(url)

if __name__ == "__main__":
    print(f"fetching for {PMID}...\n")

    #url = f"{API_BASE}/article/{PMID}/uncache"
    #r = requests.get(url)

    article_info = get_article_info(PMID)
    title = article_info.get("title", "N/A")
    issn = article_info.get("journal_issn")
    iso_abv = get_iso_abv(issn)

    journal = get_journal_score(iso_abv) or "timeout"
    print(f"journal: {journal}")

    methods = get_methods_score(PMID) or "timeout"
    print(f"methods: {methods}")

    '''
    print(f"Title: {title}")

    author = get_author_score(PMID) or "timeout"
    print(f"author: {author}")

    journal = get_journal_score(issn) or "timeout"
    print(f"journal: {journal}")

    media = get_media_score(PMID) or "timeout"
    print(f"media: {media}")

    #authentic = get_authentic_score(PMID) or "timeout"
    #print(f"authentic: {authentic}")

    methods = get_methods_score(PMID) or "timeout"
    print(f"methods: {methods}")

    overall = get_overall_score(PMID) or "timeout"
    print(f"overall: {overall}")

    '''
    print("\ndone")