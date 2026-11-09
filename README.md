# Parkinson's Disease Transcriptomics: Differential Gene Expression Analysis

## 📌 Project Overview
This repository contains an end-to-end RNA-sequencing (RNA-seq) Differential Gene Expression (DGE) pipeline built in **R**. The project analyzes the transcriptomic profile of Parkinson's Disease (PD) subjects compared to healthy controls to identify key molecular dysregulations in dopaminergic neural pathways. 

This project was developed to bridge the gap between traditional wet-lab biological assays and computational data analysis, reflecting my core focus in neuroscience research.

## 🔬 Biological Motivation
My background in biotechnology and my graduation project research on the neuroprotective effects of *Cordyceps militaris* against dopaminergic neuron degeneration highlighted a critical need: to truly understand neurodegeneration and evaluate potential therapies, we must look at the transcriptomic level. This analysis simulates the transcriptomic shift during PD progression, focusing on targeted neurodegenerative and neuroinflammatory markers.

## 🛠️ Methodology & Tools
* **Language:** R
* **Core Packages:** `DESeq2` (Normalization and DGE statistical testing)
* **Visualization:** `ggplot2`, `ggrepel` (Volcano Plots), `pheatmap` (Hierarchical Clustering)
* **Statistical Thresholds:** Adjusted p-value (FDR) < 0.05 and |Log2FoldChange| > 1.

## 📊 Key Biological Insights & Visualizations

### 1. Volcano Plot: Identifying Significant Gene Dysregulation
[Volcano Plot]([Heatmap Plot.png)](https://github.com/mariambana77/Parkinsons-RNAseq-Analysis/blob/37a91aa8909cbe5dffea960fe26d7f9888c31375/Volcano%20plot.png)

**Findings:**
* **Downregulated Pathways:** The significant downregulation of **TH** (Tyrosine Hydroxylase) and **BDNF** reflects the hallmark biological death of dopaminergic neurons and the loss of neurotrophic support in the PD model.
* **Upregulated Pathways:** The sharp upregulation of **SNCA** (Alpha-synuclein) and neuroinflammatory markers like **IL6** aligns with the pathological protein aggregation and reactive gliosis typical in Parkinson's disease.

### 2. Expression Heatmap: Sample Clustering
[Heatmap]([Heatmap Plot.png](https://github.com/mariambana77/Parkinsons-RNAseq-Analysis/blob/37a91aa8909cbe5dffea960fe26d7f9888c31375/Heatmap%20Plot.png)

**Findings:**
The unsupervised hierarchical clustering successfully segregates the PD samples from the Healthy Controls based solely on their genetic expression profiles. The Z-score scaling clearly visualizes the contrast between overactive inflammatory genes and silenced dopaminergic genes across the cohorts.

## 🚀 Future Directions
The next phase of this computational pipeline involves integrating drug-treated transcriptomic data. Specifically, I aim to apply this DGE model to evaluate the efficacy of potential neuroprotective agents—by quantifying whether they can successfully reverse the Log2FoldChange of genes like *TH* and *SNCA* back to healthy baseline levels.

## 👤 About the Author
**Mariam Elbanna** 
Biotechnology Graduate | Aspiring Neuroscientist
Driven by a passion for understanding neural networks and uncovering targeted therapies for neurodegenerative diseases through the integration of biotechnology and advanced data analysis.
