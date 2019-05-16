# jgi-data-week-workshop

## Title
Constructing an academic knowledge graph and recommendation engine

## Logistics
10am - 1pm on 20th May

## Background
Given only freely available text, we can extract sufficient data to create knowledge graphs, representing both individual components and collectives as a whole. These graphs can be used to identify key ideas, overlapping concepts and areas of missing information. From here they can recommend specific events, identify communities, derive values for missing data and predict areas of change over time. The fundamentals of this approach, however, lie in extracting reliable and representative data for each individual component.

One specific example of this can be found in academic publications. By using unique user IDs and public data we can construct a knowledge graph for a defined set of people, e.g. a department or University.

#### What to expect

This workshop will demonstrate the techniques and methods required to do this, including the use of APIs, extracting and enriching informative text, natural language processing and constructing recommendation engines.  We will show how this kind of approach can be used to recommend collaborations, automatically identify people matching a specific piece of text and identify topic areas with high and low coverage.

## Content

#### Steps/Notebooks

1. People and publications
 - How to get data from PubMed
 - How to get data from ORCID
 - How to automate the above
2. Identity enriched terms
 - TF-IDF from scratch
 - TF-IDF using sklearn
3. Putting it all together
 - Creating a recommender
 - Matching a piece of text to people

#### Flow of notebooks

Each will create data in /output directory. Full/complete copies of each are pre-computed in /data directory.

## Setup

#### Prerequisites

Using Anaconda (recommended)

Clone tutorial repo:

```
#SSH
git clone git@github.com:elswob/JGI-Data-Week-2019.git

#HTTPS
git clone https://github.com/elswob/JGI-Data-Week-2019.git
```

Activate jupyterlab environment:

```
cd JGI-Data-Week-2019
conda env create -f environment.yml
conda activate jgi-data-week-workshop
jupyter lab
```

you will see a jupyter lab in your browser

#### Alternatively

Microsoft Azure
- https://notebooks.azure.com/ben-elsworth/projects/jgi-data-week-2019
- Requires microsoft account (University of Bristol members can use standard account)
- Can use the terminal to make changes
- `cd library`
- `git pull origin master`
- https://medium.com/@mikeclymer/integrating-azure-notebooks-jupyter-notebooks-with-github-fd847e941e4


Binder
- Public to all, not that stable though
- [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/elswob/JGI-Data-Week-2019/master)


## Questions/Issues/Suggestions

- Limiting to ORCID is not ideal, but does bring an interesting bias to potential collaborations
- How to decide cutoff for TF-IDF?
	- arbitrary number is not good, too many missing term-people relationships
- Could we just use doc2vec for comparing people and publication text, e.g. create corpus treating each person as a separate document, then find most similar document (person) for each.
- Not tested CPU/Mem requirements - might break some machines
- Major issues with using PubMed, e.g. many DOIs in an ORCID not converting, i.e. not in PubMed. Means many people are underrepresented.
- Could the text matching function be modified to match people covering all terms in the text, i.e. not matching similar people, but set of people that cover all terms.
- Tokenizing, lemmatizing, etc.

## Other info

Elsevier fingerprints white paper - https://www.elsevier.com/solutions/elsevier-fingerprint-engine/elsevier-fingerprint-engine-white-paper

TF-IDF explained - https://www.quora.com/How-does-TfidfVectorizer-work-in-laymans-terms

## Issues

Have had problems with biopython. On azure this you might be able to fix this by using the terminal
```
pip install --user biopython
```
