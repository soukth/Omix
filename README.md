# 🧬 ProteinFinder: OmiXInsight

**ProteinFinder** is a powerful, interactive GUI-based bioinformatics tool for exploring and analyzing proteomics or transcriptomics datasets. It allows biologists and data scientists to easily visualize, classify, and functionally analyze protein or gene expression data using heatmaps, volcano plots, Venn diagrams, and curated protein groups.

---

## 🚀 Features

- **Interactive GUI (PySimpleGUI)**
- **Upload CSV/XLSX files** for protein or RNA expression data
- **Targeted Analysis** with curated protein groups (e.g., Transporter Proteins, Transcription Factors, etc.)
- **Interactive Heatmaps** using Plotly
- **Volcano Plots** to detect differential expression
- **Venn Diagrams** for comparing datasets
- **Functional Group Classification** (simulated Perplexity AI logic)
- **Automated citation fetching** from PubMed, NCBI, and Scopus (with BeautifulSoup)
- **Excel Report Generation** for categorized results and references

---

## 📂 Input File Format

- First column: Gene or Protein identifiers (e.g., gene symbols like `TP53`, `BRCA1`, etc.)
- Remaining columns: Expression values or statistical measures across conditions
