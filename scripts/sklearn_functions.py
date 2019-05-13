from sklearn.feature_extraction.text import TfidfVectorizer

#https://stackoverflow.com/questions/46580932/calculate-tf-idf-using-sklearn-for-n-grams-in-python
#https://towardsdatascience.com/tf-idf-for-document-ranking-from-scratch-in-python-on-real-world-dataset-796d339a4089
#https://www.bogotobogo.com/python/NLTK/tf_idf_with_scikit-learn_NLTK.php

def tfidf_corpus():
	text = ["The quick brown fox jumped over the lazy dog.",
		"The dog.",
		"The fox"]
	# create the transform
	vectorizer = TfidfVectorizer(stop_words = 'english', ngram_range=(1,3))
	# tokenize and build vocab
	vectorizer.fit(text)
	# summarize
	print(vectorizer.vocabulary_)
	print(vectorizer.idf_)
	# encode document
	vector = vectorizer.transform([text[0]])
	# summarize encoded vector
	print(vector.shape)
	print(vector.toarray())
	#get feature names
	feature_names = vectorizer.get_feature_names()
	for col in vector.nonzero()[1]:
	    print(feature_names[col], ' - ', vector[0, col])

	#return(tfidf)
