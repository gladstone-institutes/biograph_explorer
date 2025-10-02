# translator_multi_gene_query


Goal: Use TCT to perform a multigene pathfinder query. Given a set of DE genes in either ensembl gene ID, ensembl transcript ID or Hugo gene IDs and a single condtition/disease. Perform multiple translator queries, collect the resulting graphs and cluster them. The goal is to visualize/distill the large number of graphs into an easily digestable summary. Ideally the clusters could be sematically labeled using a simple open source llm so that the visualization is less cluttered. The idea is to integrate the resulting graphs directly into something like a neo4j database (which could potentially be augmented with outside data like the human protein atlas) that can be queried in a simple streamlit interface to "chat" with the resulting networks (https://neo4j.com/blog/developer/rag-tutorial/) using a very small model like Qwen2.5 14B that can run on a laptop. Focus on creating a minimally viable system that conforms to the TRAPI standards.


For testing use the following gene list and COVID-19 as the disease (source: https://pmc.ncbi.nlm.nih.gov/articles/PMC11255397/):

CD6
IFITM3
IFITM2
STAT5A
KLRG1
DPP4
IL32
PIK3AP1
FYN
IL4R
