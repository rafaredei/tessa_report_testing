import requests
import os
import pandas as pd
import time


API_BASE = "http://3.81.212.180:8000"

def task_wait(url):
    data = requests.get(url).json()
    if "task" in data:
        task_id = requests.get(url).json()["task"]
        task_url = f"{API_BASE}/task/{task_id}/wait"
        cancel_url = f"{API_BASE}/task/{task_id}/cancel"
        print(f"waiting on task {task_url}")
        try:
            r = requests.get(task_url, timeout=600)  #10 minutes, this stops control flow
        except requests.Timeout:
            print(f"Request to {url} timed out after 10 minutes.")
            requests.get(cancel_url)
            print("task cancelled")
            return "timeout"

        if r.status_code == 200:
            print(f"{task_url} successfully connected")
            task_data = r.json()
            return task_data.get("result")
        else:
            print(f"task failed to connect (status {r.status_code})")
            requests.get(cancel_url)
            print("task cancelled")
            return "timeout"
    return data.get("result")

if __name__ == "__main__":
    file_path = "pmids.txt"
    pmids = []
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            content = f.read()
        # split on commas or newlines and clean up whitespace
        pmids = [x.strip() for x in content.replace("\n", ",").split(",") if x.strip()]
        print(f"Loaded {len(pmids)} PMIDs from {file_path}")
    else:
        print("pmids.txt not found — please create it or enter PMIDs manually.")
        exit()
    
    if not pmids:
        print("No valid PMIDs provided. Exiting.")
        exit()
    
    output_file = "article_scores.csv"
    all_results = []
    start_index = 0

    '''
    if os.path.exists(output_file):
        df_existing = pd.read_csv(output_file)
        processed_pmids = set(df_existing["PMID"].astype(str))
        all_results = df_existing.to_dict(orient="records")
        pmids = [p for p in pmids if p not in processed_pmids]
        print(f"Resuming from previous progress. {len(processed_pmids)} PMIDs already done.")
        if pmids:
            print(f"Next PMID to start from: {pmids[0]}")
        else:
            print("All PMIDs already processed. Exiting.")
            exit()
    '''

    total_pmids = len(pmids)
    print(f"\nStarting processing for {total_pmids} remaining PMIDs...\n")

    start_time = time.time()

    try:
        for i, pmid in enumerate(pmids, start=1):
            print(f"\n=== [{i}/{total_pmids}] Processing PMID {pmid} ===")
            print(f"\n{'='*50}\nFetching data for PMID {pmid}...\n")
            requests.get(f"{API_BASE}/article/{pmid}/uncache")

            article_info = requests.get(f"{API_BASE}/article/{pmid}").json().get("result")
            if article_info != None:
                title = article_info.get("title", "N/A")
                authors = article_info.get("authors", [])
                author_names = [
                    f"{a.get('forename', '').strip()} {a.get('lastname', '').strip()}".strip()
                    for a in authors if isinstance(a, dict)
                ]
                authors_str = ", ".join(author_names)
                pubmed_link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                siens_link = f"https://app.siensmetrica.com/report/{pmid}"
                author_score = requests.get(f"{API_BASE}/article/{pmid}/author_scores").json().get("result").get("collective_all_score")
                print(f"author score fetched: {author_score}")

                issn = article_info.get("journal_issn")
                iso_abv = requests.get(f"{API_BASE}/journal/{issn}").json().get("result").get("iso_abv")
                journal_score = requests.get(f"{API_BASE}/journal/{iso_abv}/impact").json().get("result").get("score") * 3.25
                print(f"journal score fetched: {journal_score}")

                #we need to add them to the cache to get the scores
                task_wait(f"{API_BASE}/article/{pmid}/media")
                task_wait(f"{API_BASE}/article/{pmid}/copyleaks")
                task_wait(f"{API_BASE}/article/{pmid}/analysis")
                task_wait(f"{API_BASE}/article/{pmid}/final")

                media_data = task_wait(f"{API_BASE}/article/{pmid}/media")
                if media_data == "timeout" or media_data is None:
                    media_score = "timeout"
                elif isinstance(media_data, dict) and "elapsed_ms" in media_data:
                    media_score = "api error"
                else:
                    media_score = media_data.get("scores", {}).get("media_score", {}).get("value") if isinstance(media_data, dict) else None
                print(f"media score fetched: {media_score}")

                authentic_data = task_wait(f"{API_BASE}/article/{pmid}/copyleaks")
                if authentic_data == "timeout" or authentic_data is None:
                    authentic_score = "timeout"
                elif isinstance(authentic_data, dict) and "elapsed_ms" in authentic_data:
                    authentic_score = "api error"
                else:
                    results = authentic_data.get("results", []) if isinstance(authentic_data, dict) else []
                    authentic_score = results[0].get("probability") if results else None
                print(f"authentic score fetched: {authentic_score}")

                methods_raw = task_wait(f"{API_BASE}/article/{pmid}/analysis")
                if methods_raw == "timeout" or methods_raw is None:
                    rigor_score = empiricity_score = novelty_score = 0
                    methods_score = "timeout"
                elif isinstance(methods_raw, dict) and "elapsed_ms" in methods_raw:
                    methods_score = "api error"
                else:
                    rigor_score = methods_raw.get("rigor_score", 0)
                    empiricity_score = methods_raw.get("empiricity_score", 0)
                    novelty_score = methods_raw.get("novelty_score", 0)
                    methods_score = (rigor_score + empiricity_score + novelty_score) / 3
                print(f"methods score fetched: {methods_score}")

                overall_data = task_wait(f"{API_BASE}/article/{pmid}/final")
                if overall_data == "timeout" or overall_data is None:
                    overall_score = "timeout"
                elif isinstance(overall_data, dict) and "elapsed_ms" in overall_data:
                    overall_score = "api error"
                else:
                    overall_score = overall_data
                print(f"overall score fetched: {overall_score}")

            else:
                print("couldn't fetch any info for this pmid")
                title = "couldn't fetch"
                authors = "couldn't fetch"
                author_names = "couldn't fetch"
                authors_str = "couldn't fetch"
                pubmed_link = "couldn't fetch"
                siens_link = "couldn't fetch"
                author_score = "couldn't fetch"
                journal_score = "couldn't fetch"
                media_score = "couldn't fetch"
                authentic_score = "couldn't fetch"
                methods_score = "couldn't fetch"
                overall_score = "couldn't fetch"

            result = {
                "PMID": pmid,
                "Title": title,
                "Authors": authors_str,
                "PubMed Link": pubmed_link,
                "Siensmetrica Report": siens_link,
                "Author Score": author_score,
                "Journal Score": journal_score,
                "Media Score": media_score,
                "Authentic Score": authentic_score,
                "Methods Score": methods_score,
                "Overall Score": overall_score
            }

            all_results.append(result)

            # --- Save progress incrementally ---
            df = pd.DataFrame(all_results)
            df["NumericOverall"] = pd.to_numeric(df["Overall Score"], errors="coerce")
            df = df.sort_values(by="NumericOverall", ascending=False, na_position="last")
            df = df.drop(columns=["NumericOverall"])
            df.to_csv(output_file, index=False, encoding="utf-8")

            print(f"✅ Progress saved after PMID {pmid}")

    except KeyboardInterrupt:
        print("\n\n🟡 Script interrupted by user. Saving current progress...")
        df = pd.DataFrame(all_results)
        df.to_csv(output_file, index=False, encoding="utf-8")
        print(f"Data saved to '{output_file}'.")
        if i < total_pmids:
            print(f"You stopped at PMID {pmid}. Next time, resume from this one or skip it.")
        exit()

    end_time = time.time()
    total_time = end_time - start_time
    hours, remainder = divmod(total_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    print(f"\n✅ All done! Time elapsed: {int(hours)}h {int(minutes)}m {int(seconds)}s")
    print(f"Results saved to '{output_file}'.")
