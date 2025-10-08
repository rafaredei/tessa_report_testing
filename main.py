from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By 
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import time
import re

#different on windows
service = Service(executable_path="/Users/rafaredei/Documents/projects/siensmetrica_report_script/chromedriver")
driver = webdriver.Chrome(service=service)


#skip to siens test
'''
driver.get("https://pubmed.ncbi.nlm.nih.gov/30782404/")
driver.execute_script("window.open('https://app.siensmetrica.com/search');")
tabs = driver.window_handles
pubmed_tab = tabs[0]
siens_tab = tabs[1]
pmid = 30782404
driver.switch_to.window(siens_tab)
# wait for the input to appear
input_element = WebDriverWait(driver, 10).until(
    EC.presence_of_element_located((By.CSS_SELECTOR, "input[data-slot='input']"))
)

input_element.send_keys(str(pmid) + Keys.ENTER) #input search

# Wait until the first row in the table is loaded
pmid_cell = WebDriverWait(driver, 10).until(
    EC.element_to_be_clickable((By.CSS_SELECTOR, "tr[data-slot='table-row'] td:first-child p"))
)

# Click it
pmid_cell.click()

input("Press Enter to quit...")
driver.quit()
'''

#open all needed tabs
driver.get("https://pubmed.ncbi.nlm.nih.gov/")
driver.execute_script("window.open('https://docs.google.com/spreadsheets/d/1l11rETdG8aHaP1KrwidW338WuI5_nJkcG93VUvXU66o/edit?gid=0#gid=0');")
driver.execute_script("window.open('https://app.siensmetrica.com/search');")

tabs = driver.window_handles
pubmed_tab = tabs[0]
sheet_tab = tabs[1]
siens_tab = tabs[2]

# ask for user input on search topic
search_topic = input("Enter the PubMed search topic: ")

# instructs the user what to do before continuing
input("Login to Google Sheets, then switch to the PubMed tab and hit Enter to continue the script...")

print("Continuing...")

#wait for search bar to load on pubmed
WebDriverWait(driver, 5).until(
    EC.presence_of_element_located((By.CLASS_NAME, "term-input.tt-input"))
)

input_element = driver.find_element(By.CLASS_NAME, "term-input.tt-input") #access search input
input_element.send_keys(search_topic + Keys.ENTER) #input search

#make sure the free full text filter is selected
checkbox = WebDriverWait(driver, 10).until(
    EC.presence_of_element_located((By.ID, "id_filter_simsearch2.ffrft"))
)
driver.execute_script("arguments[0].scrollIntoView(true);", checkbox)
if not checkbox.is_selected():
    driver.execute_script("arguments[0].click();", checkbox) #real input for checkbox is hidden, so click with JS

#get first article
first_article = driver.find_element(By.CSS_SELECTOR, "article.full-docsum[data-rel-pos='1']") 
title_link = first_article.find_element(By.CSS_SELECTOR, "a.docsum-title") #locate first article's clickable element
driver.execute_script("arguments[0].scrollIntoView(true);", title_link)  # scroll to make clickable
title_link.click()

# Wait until the URL matches a PubMed article pattern (digits between slashes)
WebDriverWait(driver, 10).until(
    lambda d: re.search(r"/\d+/?$", d.current_url)
)

#save all important info
article_url = driver.current_url

match = re.search(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)", article_url)
if match:
    pmid = match.group(1)
else:
    pmid = None

title_element = driver.find_element(By.CLASS_NAME, "heading-title")
article_title = title_element.text.strip()

print("article url: " + article_url)
print("pmid: " + pmid)
print("article title: " + article_title)

#run it through tessa
driver.switch_to.window(siens_tab)
input_element = driver.find_element(By.CSS_SELECTOR, "input[data-slot='input']") 
input_element.send_keys(pmid + Keys.ENTER) #input search

# wait for the input to appear
input_element = WebDriverWait(driver, 10).until(
    EC.presence_of_element_located((By.CSS_SELECTOR, "input[data-slot='input']"))
)

input_element.send_keys(str(pmid) + Keys.ENTER) #input search

# Wait until the first row in the table is loaded
pmid_cell = WebDriverWait(driver, 10).until(
    EC.element_to_be_clickable((By.CSS_SELECTOR, "tr[data-slot='table-row'] td:first-child p"))
)

# Click it
pmid_cell.click()

input("Press Enter to quit...")
driver.quit()



