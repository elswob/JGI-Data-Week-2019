from Bio import Entrez
import gzip
import csv
import time
import requests
import subprocess
import re
import os
import json
import config
from random import randint

#https://stackoverflow.com/questions/47559098/is-there-any-way-to-get-abstracts-for-a-given-list-of-pubmed-ids

def read_existing():
	print('Read existing downloaded pubmed data from',config.pubmedFile)
	pubData = []
	if os.path.exists(config.pubmedFile):
		with open(config.pubmedFile, newline='') as csvfile:
			reader = csv.reader(csvfile, delimiter='\t')
			for row in reader:
				pubData.append({'pmid': row[0], 'year': row[1], 'title': row[2], 'abstract': row[3]})
	else:
		o=open(config.pubmedFile,'w')
		o.write('pmid\tyear\ttitle\tabstract\n')
		o.close()
	print((len(pubData)-1),'publication(s) already downloaded')
	return pubData

#function to convert DOIs to PubMed IDs
def doi_to_pmid(doiList):
	baseurl='https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?tool=my_tool&email=ben.elsworth@bristol.ac.uk&ids='
	url = baseurl+",".join(doiList)+'&idType=doi&format=json'
	pmidList=set()
	try:
		resp = requests.get(url).json()
		if 'records' in resp:
			for record in resp['records']:
				if 'pmid' in record:
					pmidList.add(record['pmid'])
	except:
		print('requests error')
	return list(pmidList)

def get_pubmed_data_entrez(pmids):
	#check for existing pubmed data
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
		try:
			#query entrez
			records = Entrez.read(handle)
			#open file for output
			with open(config.pubmedFile, 'a', newline='') as csvfile:
				fieldnames = ['pmid', 'year', 'title' , 'abstract']
				writer = csv.DictWriter(csvfile, fieldnames=fieldnames,delimiter='\t')
				#loop through each article
				for pubmed_article in records['PubmedArticle']:
					if 'PMID' in pubmed_article['MedlineCitation']:
						pmid = pubmed_article['MedlineCitation']['PMID']
					else:
						print('No PMID')
						continue
					if 'Article' in pubmed_article['MedlineCitation']:
						if 'Abstract' in pubmed_article['MedlineCitation']['Article']:
							abstract = pubmed_article['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
						else:
							print('No Abstract')
							continue
						if 'ArticleTitle' in pubmed_article['MedlineCitation']['Article']:
							title = pubmed_article['MedlineCitation']['Article']['ArticleTitle']
						else:
							print('No ArticleTitle')
							continue
					else:
						print('No Article')
						continue
					if 'DateCompleted' in pubmed_article['MedlineCitation']:
						year = pubmed_article['MedlineCitation']['DateCompleted']['Year']
					else:
						#DataCompleted is missing for some
						year=0
						print('No DateCompleted')
					pubData.append({'pmid':pmid,'year':int(year),'title':title,'abstract':abstract})
					writer.writerow({'pmid': pmid, 'year': int(year), 'title': title, 'abstract': abstract})
		except:
			print('bioentrez error')

	else:
		print('Nothing to do')
	pubFilter=[]
	for p in pubData:
		if p['pmid'] in pmids:
			pubFilter.append(p)
	return pubFilter
