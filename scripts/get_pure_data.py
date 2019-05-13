import requests
import os
import json
import gzip
import re
from time import gmtime, strftime

dataDir = './output/'
fSize=1000
#fSize=3
totalPubs = 130000
#totalPubs = 10

def get_people():
	print("Running get_people...")
	pageSize=500
	#pageSize=5
	pDic={}
	if os.path.exists(dataDir+'/pure_people.txt'):
		print('Found existing data, skipping')
	else:
		o = open(dataDir+'/pure_people.txt','w')
		for i in range(0,8):
		#for i in range(0,1):
			url = 'http://research-information.bristol.ac.uk/en/persons/search.html?filter=academic&page='+str(i)+'&pageSize='+str(pageSize)
			print(url)
			res = requests.get(url)
			uuid = re.findall('persons/(.*?)\((.*?)\).html', res.text)
			#print uuid
			for u in uuid:
				name = u[0].replace('-',' ').title()
				uuid = u[1]
				pDic[uuid]=name
		print(len(pDic))
		for p in pDic:
			o.write(p+'\t'+pDic[p]+'\n')
		o.close()

def get_info():
	print("Running get_info...")
	#check for existing data
	pureDone=set()
	if os.path.exists(dataDir+'/pure_person_to_org.txt'):
		with open (dataDir+'/pure_person_to_org.txt','r') as f:
			for line in f:
				person,org=line.split('\t')
				pureDone.add(person)

	orgDic={}
	person_to_org = {}
	orgDicOut = open(dataDir+'/pure_org_to_name.txt','a')
	person_to_org_out = open(dataDir+'/pure_person_to_org.txt','a')
	person_to_orcid = open(dataDir+'/pure_person_to_orcid.txt','a')
	counter=0
	with open(dataDir+'/pure_people.txt','r') as f:
		for line in f:
			if counter % 100 == 0:
				print(counter)
			counter+=1
			uuid,name = line.rstrip().split('\t')
			if uuid not in pureDone:
				url = 'http://research-information.bristol.ac.uk/en/persons/xxx('+uuid+').html'
				print(counter,url)
				res = requests.get(url)
				#http://research-information.bristol.ac.uk/en/organisations/school-of-social-and-community-medicine(f54add52-720b-4679-8119-24ec3aaf6f63).html
				orgs = re.findall('organisations/(.*?)\((.*?)\).html', res.text)
				orcid = re.findall('orcid.org/(.*?)".*', res.text)
				if len(orcid)>0:
					person_to_orcid.write(uuid+'\t'+orcid[0]+'\n')
				print('orcid',orcid)
				for o in orgs:
					org_name = o[0].replace('-',' ').title()
					org_uuid = o[1]
					orgDic[org_uuid]=org_name
					if uuid in person_to_org:
						person_to_org[uuid].append(org_uuid)
					else:
						person_to_org[uuid]=[org_uuid]
			else:
				print(uuid,'done')
	for i in orgDic:
		orgDicOut.write(i+'\t'+orgDic[i]+'\n')
	for i in person_to_org:
		person_to_org_out.write(i+'\t'+(',').join(person_to_org[i])+'\n')

	orgDicOut.close()
	person_to_org_out.close()
	person_to_orcid.close()

get_people()
get_info()
