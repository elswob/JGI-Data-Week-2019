import requests
import json

## Load configuration
con_file = open("config.json")
config = json.load(con_file)
con_file.close()

def get_authors():
	#https://api.elsevier.com/content/search/author?query=affil(University%20of%20Bristol)&apiKey=bd3491583a7c3dd797815dffe47e0406
	baseurl='https://api.elsevier.com/content/search/author?query=affil'
	location='(University of Bristol)'
	url = baseurl+location+'&apiKey='+config['apikey']
	print(url)
	resp = requests.get(url).json()
	totalRes = resp['search-results']['opensearch:totalResults']
	print(totalRes)
	for r in resp['search-results']['entry']:
		rData = r

get_authors()
