## This file represents the major findings ang quality checks along the line of investigation.
Two datasets were used in this study.  
:white_check_mark: Israeli patients' data contain 42 samples with gene expression, disease status and calprotectin levels, out of which 37 samples also have metabolomic data. This data was used to build the Bayesian Network and prediction model.  
:white_check_mark: Chinese patient's data contain 40 samples with gene expression and disease status (no calprotectin level!), while metabolomic data has no attribution to the disease status. This group was used to test the model.  
:round_pushpin: Differential gene expression was performed using age, sex, and disease status.  
:round_pushpin: WGCNA was performed using age, sex, and disease status plus calprotectin and c-reactive protein levels (continious data).


### WGCNA reveald clear gene clusters within about 1100 differentially expressed genes:
For this analysis, a size of the module was set as 200.  

![image](https://github.com/user-attachments/assets/038a9178-d5e7-4f02-833a-d071802891e1)

### These clusters correlated well with Calprotectin (indicates blood in stool -> inflammation) and diagnosis:
Interestingly, green cluster correlate also with C-rective protein - a non-specific indicator of systemic inflammation. Given that green cluster contains mostly immunoglubulin genes (199 out of 200), it seems biologically relevant.  

![image](https://github.com/user-attachments/assets/726e182f-a0c6-43ba-b461-acf5afbb045e)

### Gene modules were used to select key node genes based on their impact in the module (20 genes per module).
All immunoglobulin genes were removed after this step, leaving one node gene from the green module.
![image](https://github.com/user-attachments/assets/7e17edba-9c0e-4927-ae23-c0cdded67665)

### The selected genes were then used to build the Bayesian Network (BN) using *bnlearn* package
First, the BN was build only with genes without metabolites and metadata. The edges were averaged out after a 100 bootstraps, and only edges with the strength over 0.7 were used to build the final network.  
:rainbow: Here the genes are colored by their attribution to the module:
![image](https://github.com/user-attachments/assets/f0175876-a54c-45a1-acbc-d358c83942b8)

### The genes that comprised the final network were used for GO and KEGG analysis.
While KEGG revealed only one significant pathway, GO discovered an important realtion of the BN genes to actin dynamics.  
:microscope: This is a very important finding as we are working on the role of filamentous actin for a while and have shown that it is involved in epithelial barrier disruption in IBD in mice [Borisova et al., 2020](https://www.nature.com/articles/s41598-020-78141-4) and in humans (data yet unpublished and available per request). And the GO result aligns very well with our data on electron and confocal microscopy regarding actin filamnt dynamics.
![image](https://github.com/user-attachments/assets/003e6b17-9947-4395-979c-4297c3074cc4)

### The BN disocvered here was next tested for relevance and accuracy by building the predictive model.
In order to test the relevance of the BN model, we used a separate dataset from Chinese patients to predict calprotectin level (continious factor in metadata) and disease state (categorial factor from metadata).  
:pencil: Faecal calprotectin has been used to detect intestinal inflammation (colitis or enteritis) and can serve as a biomarker for inflammatory bowel diseases. :pencil:  
Here is the network including calprotectin:
![image](https://github.com/user-attachments/assets/bb6dd785-cf8b-4154-bb7e-7b64678aabad)

### Next, the model was used to predict calprotectin from the same Israeli dataset
:muscle:The most influential gene was *ISX* with bootstrap strength over 0.5.  
There was no calrotectin data for the Chinese test data set, so this prediction vs real measurements correlation for that data is not shown.
<div style="display: flex; gap: 10px;">
  <img src="https://github.com/Elena-Kozhevnikova/Diploma_BI_25_Bayes/blob/main/images/Prediction_correlation.png" alt="Image 1" width="45%" />
  <img src="https://github.com/Elena-Kozhevnikova/Diploma_BI_25_Bayes/blob/main/images/Predictive_genes.png" alt="Image 2" width="45%" />
</div>

### Predict calprotectin levels in the test data form Chinese patients and comapre by disease status
As there was no calrotectin data for the Chinese test data set, we simply predicted it from our model and tested whether it finds any significant difference among the patient groups:  
As was assessed by Wilcoxon test, there was significan difference. As expected, calprotectin is generally higher in patient group.
![Calprot_pred_patient_group](https://github.com/Elena-Kozhevnikova/Diploma_BI_25_Bayes/blob/main/images/Predicted_calprotectin.png)

### Predict disease state in the test data form Chinese patients
Finally, we performed disease predicion using BN. In order to do so, we discretized gene expression data (used "low", "Medium", and "High" levels) and recunstructed the BN using discrete data. After bootstrapping, we found only one wealk parent of the "disease" node, which was also *ISX* gene. The assessment of the model strength was evaluated by predicting the disease state in test data from Chinese patients and shown as confusion matrix.

#### The first matrix shows prediction in the train data, the second - in the test data
<div style="display: flex; gap: 10px;">
  <img src="https://github.com/Elena-Kozhevnikova/Diploma_BI_25_Bayes/blob/main/images/Confusion_matrix_Israel.png" alt="Image 1" width="45%" />
  <img src="https://github.com/Elena-Kozhevnikova/Diploma_BI_25_Bayes/blob/main/images/Confusion_matrix_China.png" alt="Image 2" width="45%" />
</div>

### Metabolites did not incorporate into BN
Finally, we assessed metabolites that were differentially changed based on diagnosis, and found only 6 of those that stringly changes expression upon Crohn's Disease. These metabolites were used to build the BN, but did not form any influantial edges.  
:disappointed_relieved:Thus, metabolites were not further used for analysis.
![image](https://github.com/user-attachments/assets/01c67da6-0142-4806-a2b7-2f700badbc50)
