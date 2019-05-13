from __future__ import division

import re
import json
import nltk
from nltk.tokenize import word_tokenize

#from nltk.book import *
from nltk.collocations import *
#import matplotlib.pyplot as plt
#import matplotlib
#matplotlib.style.use('ggplot')
from nltk.corpus import wordnet as wordnet
from textstat.textstat import textstat
from nltk.stem.snowball import SnowballStemmer
stemmer = SnowballStemmer("english")
from nltk.stem import WordNetLemmatizer
wordnet_lemmatizer = WordNetLemmatizer()

bigram_measures = nltk.collocations.BigramAssocMeasures()
trigram_measures = nltk.collocations.TrigramAssocMeasures()

#stanford NER
#from nltk.tag.stanford import StanfordNERTagger
#jar = './data/stanford-ner-2018-10-16/stanford-ner.jar'
#model = './data/stanford-ner-2018-10-16/classifiers/english.all.3class.distsim.crf.ser.gz'
#ner_tagger = StanfordNERTagger(model, jar, encoding='utf8')

#from nltk.tag.stanford import CoreNLPPOSTagger, CoreNLPNERTagger
#from nltk.parse.corenlp import CoreNLPParser
#stpos, stner = CoreNLPPOSTagger(), CoreNLPNERTagger()

#analyse output
# awk -F '\t' '$2>50' output/nltk_counts.txt | sort -t$'\t' -k4 -nr | less

#http://www.nltk.org/book/ch01.html

def content_fraction(text):
	stopwords = nltk.corpus.stopwords.words('english')
	content = [w for w in text if w.lower() not in stopwords]
	return len(content) / len(text)

def filter_stopwords_and_length(text,wordLength):
	stopwords = nltk.corpus.stopwords.words('english')
	text_no_short_and_no_stopwords = [w for w in text if len(w)>wordLength and w.lower() not in stopwords]
	return text_no_short_and_no_stopwords

def lexical_diversity(text):
	text = word_tokenize(text.lower())
	return len(set(text)) / len(text)

def get_types(text):
	#lemmatize
	#text = wordnet_lemmatizer.lemmatize(text)
	#stem
	text = stemmer.stem(text)
	aToken = word_tokenize(text.lower())
	types = set(aToken)
	types = filter_stopwords_and_length(types,5)
	return types

def wordnet_pos_code(tag):
    '''Translation from nltk tags to Wordnet code'''
    if tag.startswith('NN'):
        return wordnet.NOUN
    elif tag.startswith('VB'):
        return wordnet.VERB
    elif tag.startswith('JJ'):
        return wordnet.ADJ
    elif tag.startswith('RB'):
        return wordnet.ADV
    else:
        return None

def get_types2(text):
	#http://stackoverflow.com/questions/36188032/how-to-get-frequency-of-words-in-text-depends-of-they-part-of-speech
	stopwords = set(nltk.corpus.stopwords.words('english'))
	tagged_words = set(nltk.word_tokenize(text))

	tagged_words = nltk.pos_tag(tagged_words)

	# Remove single-character tokens (mostly punctuation)
	tagged_words = [tagged_word for tagged_word in tagged_words if len(tagged_word[0]) > 1]

	# Remove numbers
	tagged_words = [tagged_word for tagged_word in tagged_words if not tagged_word[0].isnumeric()]

	# Remove stopwords
	tagged_words = [tagged_word for tagged_word in tagged_words if tagged_word[0] not in stopwords]

	#print tagged_words

	# Dark magic
	lemmatizer = nltk.stem.WordNetLemmatizer()
	words = []
	for tagged_word in tagged_words:
		pos = wordnet_pos_code(tagged_word[1])
		# lemmatize nouns, verbs, adjectives and adverbs
		if pos is not None:
			#words.append({'word':lemmatizer.lemmatize(tagged_word[0], pos=pos), 'pos':tagged_word[1]})
			words.append(lemmatizer.lemmatize(tagged_word[0],pos=pos))
		else:
			words.append(tagged_word[0])
	#print words
	return set(words)


#http://brandonrose.org/clustering
def tokenize_and_stem(text):
	# first tokenize by sentence, then by word to ensure that punctuation is caught as it's own token
	tokens = [word for sent in nltk.sent_tokenize(text) for word in nltk.word_tokenize(sent)]
	tokens = filter_stopwords_and_length(tokens,2)
	#content = [w for w in text if w.lower() not in stopwords]

	filtered_tokens = []
	# filter out any tokens not containing letters (e.g., numeric tokens, raw punctuation)
	for token in tokens:
		if re.search('[a-zA-Z]', token):
			filtered_tokens.append(token)
	stems = [stemmer.stem(t) for t in filtered_tokens]
	return stems

def get_unigrams(text):
	#testing stanford NER
	#words = nltk.word_tokenize(text)
	#print(ner_tagger.tag(words))

	text = text.lower()
	tokens = [word for sent in nltk.sent_tokenize(text) for word in nltk.word_tokenize(sent)]
	tokens = filter_stopwords_and_length(tokens,2)
	tokens = nltk.pos_tag(tokens)

	lemmatizer = nltk.stem.WordNetLemmatizer()
	words = {}
	for token in tokens:
		if re.search('[a-zA-Z]', token[0]):
			pos = wordnet_pos_code(token[1])
			# lemmatize nouns, verbs, adjectives and adverbs
			if pos is not None:
				term=lemmatizer.lemmatize(token[0],pos=pos)
			# if not one of these types keep as it is
			else:
				term = token[0]
			if term in words:
				words[term]+=1
			else:
				words[term]=1
	#print words
	unigramData = []
	for n in words:
		unigramData.append({'t1':n,'count':words[n]})
	sorted_unigramData = sorted(unigramData, key=lambda k: k['count'],reverse=True)
	return sorted_unigramData

def get_bigrams(text):
	aToken = word_tokenize(text.lower())
	aToken = filter_stopwords_and_length(aToken,1)

	#bigrams
	#bigram_measures = nltk.collocations.BigramAssocMeasures()
	bigrams = BigramCollocationFinder.from_words(aToken)

	#freqs
	bigrams_freqs = sorted(bigrams.ngram_fd.items(), key=lambda t: (-t[1], t[0]))
	#print(bigrams_freqs)
	bigramData = []
	for n in bigrams_freqs:
		bigramData.append({'t1':n[0][0],'t2':n[0][1],'count':n[1]})
		#bigramDic[n[0]]=n[1]
	sorted_bigramData = sorted(bigramData, key=lambda k: k['count'],reverse=True)
	return sorted_bigramData

	#all
	#all_bigrams = bigrams.nbest(bigram_measures.pmi, -1)

	#above min score
	#if len(tuple(nltk.bigrams(aToken)))>0:
	#	sorted_bigrams = sorted(bigrams.above_score(bigram_measures.raw_freq,1.0 / len(tuple(nltk.bigrams(aToken)))))
		#print sorted_bigrams
	#else:
	#	sorted_bigrams = ()

	#return all_bigrams

def get_trigrams(text):
	aToken = word_tokenize(text.lower())
	aToken = filter_stopwords_and_length(aToken,1)
	#trigrams
	#trigram_measures = nltk.collocations.TrigramAssocMeasures()
	trigrams = TrigramCollocationFinder.from_words(aToken)
	#freqs
	trigrams_freqs = sorted(trigrams.ngram_fd.items(), key=lambda t: (-t[1], t[0]))
	#print(trigrams_freqs)

	trigramData = []
	for n in trigrams_freqs:
		trigramData.append({'t1':n[0][0],'t2':n[0][1],'t3':n[0][2],'count':n[1]})
	sorted_trigramData = sorted(trigramData, key=lambda k: k['count'],reverse=True)
	return sorted_trigramData

	#all
	#all_trigrams = trigrams.nbest(trigram_measures.pmi, -1)
	#above min score
	#if len(tuple(nltk.trigrams(aToken)))>0:
	#	sorted_trigrams = sorted(trigrams.above_score(trigram_measures.raw_freq,1.0 / len(tuple(nltk.trigrams(aToken)))))
		#print sorted_trigrams
	#else:
	#	sorted_trigrams = ()
	#return all_trigrams

def run_textstat(text):
	#text = """Playing games has always been thought to be important to the development of well-balanced and creative children; however, what part, if any, they should play in the lives of adults has never been researched that deeply. I believe that playing games is every bit as important for adults as for children. Not only is taking time out to play games with our children and other adults valuable to building interpersonal relationships but is also a wonderful way to release built up tension."""

	ts_flesch_reading_ease = textstat.flesch_reading_ease(text)
	ts_smog_index = textstat.smog_index(text)
	ts_flesch_kincaid_grade = textstat.flesch_kincaid_grade(text)
	ts_coleman_liau_index = textstat.coleman_liau_index(text)
	ts_automated_readability_index = textstat.automated_readability_index(text)
	ts_dale_chall_readability_score = textstat.dale_chall_readability_score(text)
	ts_difficult_words = textstat.difficult_words(text)
	ts_linsear_write_formula = textstat.linsear_write_formula(text)
	ts_gunning_fog = textstat.gunning_fog(text)
	ts_text_standard = textstat.text_standard(text)

	return(ts_flesch_reading_ease,ts_smog_index,ts_flesch_kincaid_grade,ts_coleman_liau_index,ts_automated_readability_index,ts_dale_chall_readability_score,ts_difficult_words,ts_linsear_write_formula,ts_gunning_fog,ts_text_standard)
