## :scientist: Studying Gene Regulatory Hierarchy in Crohn’s Disease with a Bayesian Network
> Elena Kozhevnikova  
> Institute of Molecular and Cellular Biology  
> Novosibirsk  
> elena.n.kozhevnikova@gmail.com  
> tg: e_zeste  

This repo describes my Diploma project at the Bioinformatics Institute (2024-2025)
______________________________________________________________________________________________
## Summary
Inflammatory bowel diseases (IBD) are chronic inflammatory disorders of the gastrointestinal tract that include Crohn's disease and ulcerative colitis. IBD are are chronic idiopathic disorders, known to be influenced by a variety of factors, including genetic predicposition and diet. Recent transcriptome analyses have uncovered thousands of genes associated with IBD revealing its complex genetic landscape. Several studies demonstrate strong connection between transcriptomic and metabolomic data, highlighting the interplay between gene expression and metabolism ([Massimino et al., 2021](https://www.nature.com/articles/s43588-021-00114-y)). Identifying key factors among the thousands of differentially expressed genes is essential for determining potential targets in drug development and advancing personalized healthcare. Here we applied weighted gene co-expression network analysis (WGCNA, [Langerfeld and Horvath, 2008](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559)) together with Bayesian network (BN, [Puga et al., 2015](https://www.nature.com/articles/nmeth.3550)) prediction to identify key genetic and metabolic factors involved in regulation of the inflammatory state characteristic to IBD using open data ([Braun et al., 2024](https://www.nature.com/articles/s41467-024-48106-6)). Our results identify ISX gene as the most influential factors that partially explains both an IBD marker calprotectin and the disease state. However, stool metabolites did not appear to influence the CD regulatory network significantly in this experimental desing. Interestingly, *ISX* gene is a intestine-specific homeobox transcription factor critical to control the vitamin A metabolism ([Airanthi et al., 2017](https://www.pnas.org/doi/full/10.1073/pnas.1714963114)). Moreover, genetic polymorphisms in the *ISX* gene have been associated with inflammatory bowel disease in genome-wide association studies ([Dinu et al., 2012](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0043035)).




Some helpful ideas: 
Applications of Bayesian network with WGCNA eigengenes and PCA
https://www.nature.com/articles/s41598-018-24758-5#Sec14

R package for Bayesian network learning and inference (alternative to Dagma)
https://www.bnlearn.com/
