from Bio import Entrez
import gzip
import csv
import time
import requests
import subprocess
import re
import os
from random import randint

#https://stackoverflow.com/questions/47559098/is-there-any-way-to-get-abstracts-for-a-given-list-of-pubmed-ids

def read_existing():
	print('Read existing downloaded pubmed data')
	pubData = []
	if os.path.exists('data/pubmed.csv'):
		with open('data/pubmed.csv', newline='') as csvfile:
			reader = csv.reader(csvfile, delimiter=',')
			#next(reader, None)
			for row in reader:
				#print(row)
				pubData.append({'pmid': row[0], 'year': row[1], 'title': row[2], 'abstract': row[3]})
	print(len(pubData))
	return pubData

def pubmed_query_to_pmids(query):
	start=time.time()
	print("\n### Getting ids for "+query+" ###")
	url="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
	params = {'db': 'pubmed', 'term': query,'retmax':'1000000000','rettype':'uilist'}
	print(url,params)
	# GET with params in URL
	r = requests.get(url, params=params)

	#create random file name
	n = 10
	ran=''.join(["%s" % randint(0, 9) for num in range(0, n)])

	ranFile = '/tmp/'+ran+'.txt'
	out = open(ranFile, 'w')
	out.write(r.text)
	out.close()
	r.status_code
	end=time.time()
	print("Time taken:",round((end-start)/60,3),"minutes")

	#count the number of pmids
	cmd = "grep -c '<Id>' "+ranFile
	pCount=0
	#print(cmd)
	#check for empty searches
	try:
		pCount = int(subprocess.check_output(cmd, shell=True))
	except:
		print("No results")

	print("Total pmids: "+str(pCount))
	maxA=1000000
	counter=0
	pmidList=[]
	if 0<pCount<maxA:
		print("\n### Parsing ids ###")
		start = time.time()
		f = open('/tmp/'+ran+'.txt', 'r')
		for line in f:
			l = re.search(r'.*?<Id>(.*?)</Id>', line)
			if l:
				pmid = l.group(1)
				pmidList.append(pmid)
	return pmidList


def get_pubmed_data_entrez(pmids):
	pubData = read_existing()

	#check if already done
	pmidsToDo = []
	for p in pmids:
		if any(d['pmid'] == p for d in pubData):
			print(p,'is done')
		else:
			pmidsToDo.append(p)

	if len(pmidsToDo)>0:
		print('Processing',pmidsToDo)
		Entrez.email = 'ben.elsworth@bristol.ac.uk'
		handle = Entrez.efetch(db="pubmed", id=','.join(map(str, pmidsToDo)),
		                       rettype="xml", retmode="text")
		records = Entrez.read(handle)

		with open('data/pubmed.csv', 'a', newline='') as csvfile:
			fieldnames = ['pmid', 'year', 'title' , 'abstract']
			writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

			for pubmed_article in records['PubmedArticle']:
				print(pubmed_article)
				if 'PMID' in pubmed_article['MedlineCitation']:
					pmid = pubmed_article['MedlineCitation']['PMID']
				else:
					continue
				if 'Article' in pubmed_article['MedlineCitation']:
					if 'Abstract' in pubmed_article['MedlineCitation']['Article']:
						abstract = pubmed_article['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
					else:
						continue
					if 'ArticleTitle' in pubmed_article['MedlineCitation']['Article']:
						title = pubmed_article['MedlineCitation']['Article']['ArticleTitle']
					else:
						continue
				else:
					continue
				if 'DateCompleted' in pubmed_article['MedlineCitation']:
					year = pubmed_article['MedlineCitation']['DateCompleted']['Year'] 
				else: 
					continue
				pubData.append({'pmid':pmid,'year':int(year),'title':title,'abstract':abstract})
				writer.writerow({'pmid': pmid, 'year': int(year), 'title': title, 'abstract': abstract})
	else:
		print('Nothing to do')
	return pubData
