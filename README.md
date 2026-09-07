# PFBA Hormesis Transcriptomics in *Spodoptera frugiperda*

Analysis code associated with the manuscript:

**Early Mitochondrial Transcriptomic Reprogramming Precedes Perfluorobutanoic Acid (PFBA)-Induced Hormetic Growth in *Spodoptera frugiperda***

Eguono W. Omagamre, Tayden Kelly, and Joseph S. Pitula

Environmental Toxicology and Intervention Laboratory (EnviToxIn Lab)  
University of Maryland Eastern Shore

## Overview

This repository contains the R scripts used for transcriptomic, statistical, functional, and figure-supporting analyses in the associated study of early molecular responses to perfluorobutanoic acid (PFBA) in *Spodoptera frugiperda*.

Larvae were exposed to control diet or PFBA at 0.01 or 1 µg/g diet. RNA-seq was performed at 2 and 4 h after exposure with four biological replicates per treatment-time combination (24 libraries total). Phenotypic measurements were evaluated through 6 h.

The analysis examines whether early transcriptomic responses precede the previously observed PFBA-associated hormetic growth phenotype, with particular emphasis on temporal differential expression, mitochondrial respiration and ATP-associated processes, stress/redox responses, and broader developmental regulation.

## Repository organization

```text
scripts/
├── 01_preprocessing/
│   ├── 01_preprocessing_filtering.R
│   └── 02_sample_QC.R
├── 02_atlas/
│   ├── 03_gene_variance_selection.R
│   ├── 04_gene_atlas_construction.R
│   ├── 04B_atlas_topology_zones.R
│   └── 05_functional_gene_atlas.R
├── 03_differential_expression/
├── 04_annotation_enrichment/
├── 05_figures/
├── 06_supplementary/
└── 07_module_validation/
```

## Analysis workflow

### 1. Preprocessing and sample-level QC

The preprocessing scripts construct and filter the gene-level count dataset and perform sample-level quality assessment.

Genes were retained when they had at least 10 raw counts in at least four samples, resulting in 11,661 expressed genes for downstream analysis. Variance-stabilized expression values were used for exploratory sample- and gene-level analyses, whereas differential-expression testing was performed using raw counts in DESeq2.

### 2. Gene-level transcriptomic atlas

Genes were ranked according to variance across the 24-sample dataset. Candidate feature sets of 2,000, 3,000, 4,000, and 6,000 genes were evaluated for stability of the gene-level UMAP representation.

The 3,000-gene atlas was selected as a parsimonious representation at or beyond the variance-curve elbow while retaining strong agreement with the larger reference embedding.

`04_gene_atlas_construction.R` constructs the fixed gene-level UMAP atlas.

`04B_atlas_topology_zones.R` partitions the fixed 3,000-gene atlas into three topology-defined regions using Euclidean distance on the two-dimensional UMAP coordinates followed by single-linkage hierarchical clustering and cutting the dendrogram at k = 3.

This procedure produces Zone 1 (930 genes), Zone 2 (901 genes), and Zone 3 (1,169 genes). The publication reproducibility script was validated against the archived atlas assignments and reproduced all 3,000 gene assignments exactly.

The topology-defined regions are used for descriptive summaries of relative temporal expression behavior and not as independent differential-expression tests.

### 3. Differential-expression analysis

Differential expression was performed with DESeq2 using treatment and time. The workflow includes time-matched treatment-versus-control contrasts, model validation, temporal classification of differentially expressed genes, and cross-dose comparisons.

Genes were classified into temporal response groups including early-only, persistent, and later-emerging responses according to their significance across the 2- and 4-h contrasts. False-discovery-rate correction was performed using the Benjamini-Hochberg procedure.

### 4. Functional annotation and GO enrichment

Functional annotation was integrated with the DESeq2 gene universe.

Gene Ontology Biological Process over-representation analysis is implemented in `11_DEG_GO_enrichment.R`. The analysis uses genes with GO annotation as the testable universe and applies one-sided Fisher's exact tests followed by Benjamini-Hochberg correction.

### 5. Targeted functional modules

The `07_module_validation` scripts examine targeted functional modules used for descriptive interpretation of the transcriptomic response, including stress/redox, mitochondrial respiration, and ATP/bioenergetic genes.

These targeted summaries are descriptive analyses and are not treated as independent pathway-enrichment tests.

### 6. Phenotypic analysis

Phenotypic scripts contain the statistical analyses supporting the early larval weight results and Supplementary Table S1.

Larvae measured at different time points represent independent samples rather than repeated longitudinal measurements of the same individuals.

## Data availability

Raw RNA-sequencing data associated with this study will be deposited in the NCBI Sequence Read Archive (SRA) under a BioProject accession to be added upon deposition.

The RNA-seq dataset comprises 24 paired-end libraries representing three exposure groups, two sampling times, and four biological replicates per treatment-time combination.

Genome and annotation resources required for read alignment and gene-level quantification are not duplicated in this repository. Processed intermediate files are generated by the analysis workflow and are not all included in this source-code repository.

The version of this repository associated with the publication will be permanently archived in Zenodo and assigned a DOI.

## Reproducibility

Scripts are organized approximately in analysis order. Individual scripts document their required input and output paths.

Because the analyses were developed and executed in an HPC environment, some scripts reference project-relative paths to intermediate data objects and external genome/annotation resources. Users reproducing the workflow in another environment should update those paths as appropriate.

Software versions and session information generated during the analysis are retained with the underlying project records.

## Citation

Citation information will be updated following publication.

The archived release corresponding to the submitted manuscript will also be citable through Zenodo.

## Contact

**Eguono W. Omagamre**  
Environmental Toxicology and Intervention Laboratory (EnviToxIn Lab)  
University of Maryland Eastern Shore
