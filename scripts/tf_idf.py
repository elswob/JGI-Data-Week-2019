import nltk
import string
import os
import config
import csv

from sklearn.feature_extraction.text import TfidfVectorizer
from nltk.stem.porter import PorterStemmer

#based on this
#https://www.bogotobogo.com/python/NLTK/tf_idf_with_scikit-learn_NLTK.php

def tokenize(text):
    tokens = nltk.word_tokenize(text)
    stems = []
    for item in tokens:
        stems.append(PorterStemmer().stem(item))
    return stems

def read_corpus():

	token_dict = {}

	print('Reading corpus')
	with open(config.pubmedFile, newline='') as csvfile:
		reader = csv.reader(csvfile, delimiter='\t')
		next(reader, None)
		for row in reader:
			#print(row)
			text=row[2]+' '+row[3]
			token_dict[row[0]] = text.lower().translate(str.maketrans('','',string.punctuation))

	#custom tokeniser
	#tfidf = TfidfVectorizer(tokenizer=tokenize, stop_words='english',ngram_range=(1,3))

	#sklean tokeniser
	tfidf = TfidfVectorizer(stop_words='english',ngram_range=(1,3))

	tfs = tfidf.fit_transform(token_dict.values())
	#print(tfs)
	return tfidf

def tfidf_doc(tfidf='',text=''):
	#str = 'Genomics and transcriptomics across the diversity of the Nematoda.'
	text=text.lower().translate(str.maketrans('','',string.punctuation))
	response = tfidf.transform([text])
	#print(response)

	feature_names = tfidf.get_feature_names()
	res={}
	for col in response.nonzero()[1]:
		res[feature_names[col]]=response[0, col]
	    #print(feature_names[col], ' - ', response[0, col])
	sorted_res = sorted(res.items(), key=lambda kv: kv[1], reverse=True)
	return sorted_res
