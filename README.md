# oneinall  
A compact and reproducible single-sample GSEA toolkit for R

![R-CMD-check](https://img.shields.io/badge/R%20CMD%20check-passing-brightgreen)
![License](https://img.shields.io/badge/license-MIT-blue)
![Status](https://img.shields.io/badge/status-active-success)

---

## 📦 Overview

**oneinall** is a lightweight R package providing a clean, reproducible, and fully customizable implementation of **single-sample GSEA (Gene Set Enrichment Analysis)**.

It includes:

- A robust single-sample GSEA computation  
- Ready-to-use GSEA-style visualizations  
- Automatic extraction of the **leading-edge subset**  
- Bootstrap-based validation of enrichment  
- Threshold-aware ranking metric visualization  

The goal is to offer a compact, pedagogical, and pipeline-friendly toolkit for transcriptomic analysis.

---

## 🚀 Installation

Install from GitHub:

```r
# install.packages("devtools")
devtools::install_github("cdesterke/oneinall")
```

---

🚀 Main Features

1. gsea_single_vector() — Compute single-sample GSEA

    Enrichment Score (ES)

    Normalized ES (NES)

    Empirical p-value via permutations

    Running ES curve

    Gene set membership tracking

2. plot_gsea_single() — Classical GSEA enrichment plot

    Running ES curve

    Hit positions

    Threshold‑aware gene classification

3. validate_gsea_above_bootstrap() — Bootstrap enrichment validation

    Null distribution of random gene sets

    Empirical p-value

    Annotated histogram

4. plot_rank_metric_thresholded() — Threshold-aware ranking metric plot

(This is the step that was missing in the previous README.)

This visualization shows:

    the ranking metric across the transcriptome

    vertical bars for gene set members

    color‑coding based on whether each gene is above or below a user‑defined threshold

It provides a complementary view to the running ES curve by showing how strongly each gene contributes to the ranking metric.
5. extract_leading_edge() — Leading-edge extraction

    Identifies genes contributing to the maximum absolute ES

    Returns a structured table with ranking and ES information

6. plot_gsea_with_leading_edge() — Full GSEA plot with leading-edge highlighting

    Running ES

    Hit positions

    Leading-edge bars

---

🧪 Full example 

```r
library(oneinall)

# Données simulées
expr <- rnorm(1000)
names(expr) <- paste0("Gene", 1:1000)
geneset <- sample(names(expr), 50)

# 1. Calcul GSEA
res <- gsea_single_vector(expr, geneset)

# 2. Plot GSEA classique
p1 <- plot_gsea_single(res, threshold = 0)
print(p1$plot)

# 3. Validation bootstrap
v <- validate_gsea_above_bootstrap(res, threshold = 0)
print(v$plot)

# 4. Plot du ranking metric (threshold-aware)
p4 <- plot_rank_metric_thresholded(
  result = res,
  threshold = 0,
  title = "Rank metric (threshold-aware)"
)
print(p4)

# 5. Extraction du leading-edge
le <- extract_leading_edge(res)
head(le)

# 6. Plot GSEA complet avec leading-edge
p6 <- plot_gsea_with_leading_edge(res)
print(p6$plot)

```
---

📊 Available Visualizations

    Classical GSEA enrichment curve   
![res](https://github.com/cdesterke/oneinall/tree/main/02_plot.png)

    Ranking gene table
![res](https://github.com/cdesterke/oneinall/tree/main/03_gene_table.png)
    
    Bootstrap validation against null distribution
![res](https://github.com/cdesterke/oneinall/tree/main/04_plot.png)
    
    Threshold-aware ranking metric
![res](https://github.com/cdesterke/oneinall/tree/main/04_plotMetric.png)
    
    Enrichment curve with leading-edge highlighting
![res](https://github.com/cdesterke/oneinall/tree/main/06_plotLE.png)

    


    


---
