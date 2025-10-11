import requests

API_BASE = "http://3.81.212.180:8000"
PMID = "38575319"

def fetch_score(url):
    r = requests.get(url)
    if r.status_code == 200:
        print(f"{url} successfully connected")
        data = r.json()
        if "task" in data:
            task_id = data["task"]
            task_url = f"{API_BASE}/task/{task_id}/wait"
            r = requests.get(task_url)
            print(f"waiting on {task_url}")
            if r.status_code == 200:
                print(f"{task_url} successfully connected")
                task_data = r.json()
                return task_data.get("result")
        return data.get("result")

def get_article_info(pmid):
    url = f"{API_BASE}/article/{pmid}"
    return fetch_score(url)

def get_iso_abv(issn):
    url = f"{API_BASE}/journal/{issn}"
    result = fetch_score(url)
    return result.get("iso_abv")

def get_author_score(pmid):
    url = f"{API_BASE}/article/{pmid}/author_scores"
    result = fetch_score(url)
    return result.get("collective_score")

def get_journal_score(iso_abv):
    url = f"{API_BASE}/journal/{iso_abv}/impact"
    result = fetch_score(url)
    return result.get("score")

def get_media_score(pmid):
    url = f"{API_BASE}/article/{pmid}/media"
    return fetch_score(url)

def get_authentic_score(pmid):
    url = f"{API_BASE}/article/{pmid}/copyleaks"
    result = fetch_score(url).get("results")[1]
    return result.get("probability")

def get_methods_score(pmid):
    url = f"{API_BASE}/article/{pmid}/analysis"
    result = fetch_score(url)
    return result.get("rigor_score")

def get_overall_score(pmid):
    url = f"{API_BASE}/article/{pmid}/final"
    return fetch_score(url)

if __name__ == "__main__":
    print(f"fetching for {PMID}...\n")
    
    article_info = get_article_info(PMID)

    title = article_info.get("title")
    print(f"title: {title}")

    issn = article_info.get("journal_issn")
    iso_abv = get_iso_abv(issn)

    author = get_author_score(PMID) or "timeout"
    print(f"author: {author}")

    journal = get_journal_score(iso_abv) or "timeout"
    print(f"journal: {journal}")

    media = get_media_score(PMID) or "timeout"
    print(f"media: {media}")

    authentic = get_authentic_score(PMID) or "timeout"
    print(f"authentic: {authentic}")

    methods = get_methods_score(PMID) or "timeout"
    print(f"methods: {methods}")

    overall = get_overall_score(PMID) or "timeout"
    print(f"overall: {overall}")

    print("\ndone")