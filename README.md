# Multi-Omics Integration — Schizophrenia GWAS × RNA-seq

Independent computational genomics project intersecting genetic association signal with transcriptomic differential expression data for schizophrenia, designed as a companion analysis to two prior projects.

## Overview

This project integrates two independent lines of evidence for schizophrenia gene involvement: genes implicated by genome-wide association (GWAS) and genes nominally differentially expressed in postmortem brain tissue (RNA-seq). The goal is to identify genes supported by both genetic and transcriptomic evidence, which would represent stronger candidates for disease relevance than either analysis alone.

## Key Finding

- **195** unique genes identified from GWAS genome-wide significant loci (top 150 independent loci, ±50kb annotation window)
- **200** candidate genes from RNA-seq, ranked by nominal (uncorrected) p-value
- **1** gene overlapping between the two sets: **ZSCAN12**

## Interpretation

This minimal overlap is an expected and informative result, not a failed analysis. The RNA-seq cohort (19 cases, 19 controls) showed no genes surviving genome-wide FDR correction in the companion RNA-seq project, with nominal significance already shown to approximate chance-level expectation (1,361 observed vs. ~1,694 expected genes at p<0.05). Given this limited statistical power, minimal overlap with a well-powered GWAS gene set (130,000+ individuals) is the statistically honest outcome, rather than evidence against shared genetic-transcriptomic mechanisms in schizophrenia.

This finding reinforces the need for substantially larger transcriptomic cohorts (e.g., CommonMind Consortium, PsychENCODE) to meaningfully bridge genetic association and transcriptomic dysregulation in future work.

## Pipeline

| Step | Script | Output |
|---|---|---|
| 1. Load & intersect gene sets | `src/multiomics_integration.py` | `outputs/multiomics_overlap.csv`, `figures/venn_diagram.png` |

## Data Sources

- GWAS gene list: derived from [GWAS-schizophrenia](https://github.com/mhdsuhaanahmed/GWAS-schizophrenia) (PGC3 summary statistics, top 150 independent loci)
- RNA-seq gene list: derived from [RNAseq-schizophrenia](https://github.com/mhdsuhaanahmed/RNAseq-schizophrenia) (GSE87194, DLPFC differential expression)

## Tech Stack

Python · pandas · numpy · matplotlib · matplotlib-venn

## Related Projects

This is the third and final integration step in a three-part schizophrenia genomics investigation:
1. [GWAS-schizophrenia](https://github.com/mhdsuhaanahmed/GWAS-schizophrenia) — genetic association analysis
2. [RNAseq-schizophrenia](https://github.com/mhdsuhaanahmed/RNAseq-schizophrenia) — transcriptomic differential expression
3. **Multiomics-schizophrenia** (this repo) — integration of both evidence layers

## Author

Mohammed Suhaan Ahmed — B.E. Biomedical Engineering, Osmania University
Research interests: Computational Genomics · Psychiatric Genetics · Bioinformatics
