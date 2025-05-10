# :scientist: Studying Gene Regulatory Hierarchy in Crohn’s Disease with a Bayesian Network
> Elena Kozhevnikova  
> Institute of Molecular and Cellular Biology  
> Novosibirsk  
> elena.n.kozhevnikova@gmail.com  
> tg: e_zeste  

This repo describes my Diploma project at the Bioinformatics Institute (2024-2025)  
Project supervisor: Dr. Maxim Struchalin
______________________________________________________________________________________________
## Summary
Inflammatory bowel diseases (IBD) are chronic inflammatory disorders of the gastrointestinal tract that include Crohn's disease and ulcerative colitis. IBD are are chronic idiopathic disorders, known to be influenced by a variety of factors, including genetic predicposition and diet. Recent transcriptome analyses have uncovered thousands of genes associated with IBD revealing its complex genetic landscape. Several studies demonstrate strong connection between transcriptomic and metabolomic data, highlighting the interplay between gene expression and metabolism ([Massimino et al., 2021](https://www.nature.com/articles/s43588-021-00114-y)). Identifying key factors among the thousands of differentially expressed genes is essential for determining potential targets in drug development and advancing personalized healthcare. Here we applied weighted gene co-expression network analysis (WGCNA, [Langerfeld and Horvath, 2008](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559)) together with Bayesian network (BN, [Puga et al., 2015](https://www.nature.com/articles/nmeth.3550)) prediction to identify key genetic and metabolic factors involved in regulation of the inflammatory state characteristic to IBD using open data ([Braun et al., 2024](https://www.nature.com/articles/s41467-024-48106-6)). Our results identify ISX gene as the most influential factors that partially explains both an IBD marker calprotectin and the disease state. However, stool metabolites did not appear to influence the CD regulatory network significantly in this experimental desing. Interestingly, *ISX* gene is a intestine-specific homeobox transcription factor critical to control the vitamin A metabolism ([Airanthi et al., 2017](https://www.pnas.org/doi/full/10.1073/pnas.1714963114)). Moreover, genetic polymorphisms in the *ISX* gene have been associated with inflammatory bowel disease in genome-wide association studies ([Dinu et al., 2012](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0043035)). Our result suggest that *ISX* gene should be used as one of the targets to be tested in further experimantal studies involving *in vivo* systems. 

*Keywords* Crohn's disease, Bayesian Networks, WGCNA, gene regulation, metabolism

<img src="https://raw.githubusercontent.com/Elena-Kozhevnikova/Diploma_BI_25_Bayes/main/images/BN_example.jpg" alt="Bayesian Network Example" width="1000"/>

### The experimental approach described here was inspired by [Agrahari et al., 2018](https://www.nature.com/articles/s41598-018-24758-5#Sec14). 
### We used [*bnlearn*](https://www.bnlearn.com/) R package for Bayesian network learning and inference.

### Contents:
This repo contains all *R code* used to generate data and the results as markdown files:
> [1_Data_PReparation.md](https://github.com/Elena-Kozhevnikova/Diploma_BI_25_Bayes/blob/main/code/1_Data_Preparation.md) described how individual's data were downloaded and mreged into one gene expression file.  
> [2_DGE_WGCNA.md](https://github.com/Elena-Kozhevnikova/Diploma_BI_25_Bayes/blob/main/code/2_DGE_WGCNA.md) describes differential gene expression and WGCNA analysis  
> [3_BN_construction.md](https://github.com/Elena-Kozhevnikova/Diploma_BI_25_Bayes/blob/main/code/3_BN_construction.md) shows basic network that only contains genes  
> [4_Calprotectin_Disease_Prediction.md](https://github.com/Elena-Kozhevnikova/Diploma_BI_25_Bayes/blob/main/code/4_Calprotectin_Disease_Prediction.md) demonstrates how calprotectin and disease state fit into the network and how well our BN predicts these parameters  
> [5_Metabolites.md](https://github.com/Elena-Kozhevnikova/Diploma_BI_25_Bayes/blob/main/code/5_Metabolites.md) shows the impact of fecal metaolites on the BN network  

*Data files* used in this project include:  
1.```Israel_metadata.csv``` contain Israeli patients' metadata  
2.```china_metadata.csv``` contain Chinese patient metadata  
3.```Metabolites_Israel.csv``` contail metabolomics data from Israeli patients  
4.```WGCNA_hub_genes_for_Bayesian_network_4_samples.csv``` contain the result of WGCNA on all samples whith no exclusion od outliers  
5.```datExpr.rds``` contains selected gene expression data  

*Bigger data files* fc_data.csv containing Isreali patients' gene expression dataset and modified_output.csv containing Chinese gene expression dataset are available by request. Alternatively, these data can be uploaded from its origin for [Israeli](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199906) and [Chinese](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233900) patients originally described by ([Braun et al., 2024](https://www.nature.com/articles/s41467-024-48106-6)) and processes using [1_Data_PReparation.md](https://github.com/Elena-Kozhevnikova/Diploma_BI_25_Bayes/blob/main/code/1_Data_Preparation.md).

*Results* of this work are presented in ```Results.md``` in this repo.
