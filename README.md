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

## Details

#### Steps

1. Create lists of people and publications
2. Get publication data
3. Process publicaton data using NLP
4. Identify enriched concepts for each person/group
5. Create network of groups, people and concepts
6. Build recommendation engine
7. Match a piece of text to people

#### Flow of notebooks

Each will create data in /output directory. Full/complete copies of each are pre-computed in /data directory. 

#### Prerequisites

Using Anaconda (recommended)

Clone tutorial repo:

```
git clone xxx
```

Activate jupyterlab environment: 

```
cd xxx
conda env create -f environment.yml
conda activate jgi-data-week-workshop
jupyter lab
```

you will see a jupyter lab in your browser


#### Questions/Issues/Suggestions

- Limiting to ORCID is not ideal, but does bring an interesting bias to potential collaborations
- How to decide cutoff for TF-IDF?
- Could we just use doc2vec for comparing people and publication text, e.g. create corpus treating each person as a separate document, then find most similar document (person) for each.
- Not tested CPU/Mem requirements - might break some machines


#### Other info

Elsevier fingerprints white paper - https://www.elsevier.com/solutions/elsevier-fingerprint-engine/elsevier-fingerprint-engine-white-paper

TF-IDF explained - https://www.quora.com/How-does-TfidfVectorizer-work-in-laymans-terms
