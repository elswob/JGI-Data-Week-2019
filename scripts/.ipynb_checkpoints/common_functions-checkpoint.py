import requests
import json
import config

from scripts.pubmed_functions import *
from scripts.nltk_functions import *

def get_ids_from_orcid_public_api(orcid):
	resp = requests.get("http://pub.orcid.org/"+orcid+"/works/",
	                    headers={'Accept':'application/orcid+json'})

	results = resp.json()
	#print(json.dumps(results,indent=4))
	#PMID_List = []
	#DOI_List = []
	pubData = []
	if 'group' in results:
		for i, result in enumerate( results['group']):
			#print(i,result)
			pubDic={}
			if 'external-ids' in result:
				for e in result['external-ids']['external-id']:
					if e['external-id-type']=='pmid':
						pmid = e['external-id-value']
						#print(pmid)
						#PMID_List.append(pmid)
						pubDic['pmid']=pmid
					elif e['external-id-type']=='doi':
						doi = e['external-id-value']
						#print(pmid)
						#DOI_List.append(doi)
						pubDic['doi']=doi
			if len(pubDic)>0:
				pubData.append(pubDic)
	else:
		print('no data found')
	return pubData

def orcid_to_pubmed(orcid_ids):
	orcidPmidData=[]
	if os.path.exists(config.orcidFile):
		with open(config.orcidFile) as f:
			for line in f:
				orcid,pmid = line.split('\t')
				orcidPmidData.append({'orcid':orcid,'pmid':pmid})
		orcidFile=open(config.orcidFile,'a')
	else:
		orcidFile=open(config.orcidFile,'w')
		orcidFile.write('orcid\tpmid\n')

	pubData=[]
	counter=0
	for o in orcid_ids:
		counter+=1
		if not any(d['orcid'] == o for d in orcidPmidData):
			print(counter,'Getting ORCID data for',o)
			orcidData=get_ids_from_orcid_public_api(o)
			pubMedIDs = set()
			doiIDs = set()
			for i in orcidData:
				if 'pmid' in i:
					pubMedIDs.add(i['pmid'])
				if 'doi' in i:
					doiIDs.add(i['doi'])
			print(len(pubMedIDs),'PMIDs')
			print(len(doiIDs),'DOIs')
			doi_pmid=doi_to_pmid(list(doiIDs))
			allPMIDs = list(set(list(pubMedIDs)+list(doi_pmid)))
			pubData=get_pubmed_data_entrez(allPMIDs)
			for p in allPMIDs:
				orcidFile.write(o+'\t'+p+'\n')
	orcidFile.close()

	return pubData

def get_ngrams(text):
        #print('\nCounting unigrams')
        unigrams=get_unigrams(text)
        #print(len(unigrams))
        #print(unigrams[:10])

        #print('\nCounting bigrams')
        bigrams = get_bigrams(text)
        #print(len(bigrams))
        #print(bigrams[:10])

        #print('\nCounting trigrams')
        trigrams = get_trigrams(text)
        #print(len(trigrams))
        #print(trigrams[:10])

        return unigrams,bigrams,trigrams

def pubmed_to_ngram():
	unigramFile=open(config.unigramFile,'w')
	unigramFile.write('pmid\tt1\tcount\n')
	bigramFile=open(config.bigramFile,'w')
	bigramFile.write('pmid\tt1\tt2\t\tcount\n')
	trigramFile=open(config.trigramFile,'w')
	trigramFile.write('pmid\tt1\tt2\tt3\tcount\n')
	with open(config.pubmedFile, newline='') as csvfile:
		reader = csv.reader(csvfile, delimiter='\t')
		next(reader, None)
		for row in reader:
			#print(row)
			pmid=row[0]
			text=row[2]+' '+row[3]
			unigrams,bigrams,trigrams=get_ngrams(text)
			for i in unigrams:
				unigramFile.write(pmid+'\t'+i['t1']+'\t'+str(i['count'])+'\n')
			for i in bigrams:
				bigramFile.write(pmid+'\t'+i['t1']+'\t'+i['t2']+'\t'+str(i['count'])+'\n')
			for i in trigrams:
				trigramFile.write(pmid+'\t'+i['t1']+'\t'+i['t2']+'\t'+i['t3']+'\t'+str(i['count'])+'\n')
			#print(unigrams[0:10])
	unigramFile.close()
	bigramFile.close()
	trigramFile.close()

#load orcid to pmid data
def load_orcid():
    print('load_orcid')
    orcidToPubmedID={}
    with open(config.orcidFile) as f:
        next(f)
        for line in f:
            orcid,pmid = line.rstrip().split('\t')
            if orcid in orcidToPubmedID:
                orcidToPubmedID[orcid].append(pmid)
            else:
                orcidToPubmedID[orcid]=[pmid]
    return orcidToPubmedID

#load the publication data
def load_pubmed():
    print('load_pubmed')
    pubmedText={}
    with open(config.pubmedFile, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        next(reader, None)
        for row in reader:
            text=row[2]+' '+row[3]
            pubmedText[row[0]]=text
    return pubmedText

#create dictionary of orcid to publication text
def orcid_to_pubmed():
    print('orcid_to_pubmed')
    orcidToPubmedID=load_orcid()
    pubmedText = load_pubmed()
    orcidToPubmed={}
    for orcid in orcidToPubmedID:
        oText=''
        for p in orcidToPubmedID[orcid]:
            if p in pubmedText:
                oText+=(pubmedText[p])
        print(len(oText))
        orcidToPubmed[orcid]=oText
    return orcidToPubmed
