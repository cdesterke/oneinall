# oneinall

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

🧠 Main Features
✔ gsea_single_vector()

Computes a single-sample GSEA score:

    Running ES curve

    ES, NES

    Permutation-based p-value

    Sorted ranking

    Hit vector

    Leading-edge implicitly defined

✔ plot_gsea_single()

Generates a classical GSEA plot:

    Running ES curve

    Hit positions

    ES / NES / p-value annotation

✔ plot_rank_metric_thresholded()

Visualizes the ranking metric with:

    Hit bars

    Threshold-aware coloring

    Max ES position

✔ extract_leading_edge()

Extracts the leading-edge subset automatically.
✔ plot_gsea_with_leading_edge()

Full GSEA-style plot including:

    Running ES

    Hit bars

    Highlighted leading-edge genes

    ES / NES / p-value annotation

✔ validate_gsea_above_bootstrap()

Bootstrap validation of enrichment for genes above a threshold.
