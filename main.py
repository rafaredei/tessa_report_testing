from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By 
import time


#different on windows
service = Service(executable_path="/Users/rafaredei/Documents/projects/siensmetrica_report_script/chromedriver")

driver = webdriver.Chrome(service = service)

driver.get("https://pubmed.ncbi.nlm.nih.gov/") #go to pubmed

input_element = driver.find_element(By.CLASS_NAME, "term-input.tt-input")
input_element.send_keys("vegan longevity")

time.sleep(10)

driver.quit()