# Overview

This repository contains code snippets for various omics analyses that cover a spectrum of tasks including data preprocessing (i.e., munging), analysis, post-processing, and visualization. To aid in identifying relevant code, scripts and code snippets are loosely organized into the following hierarchical structure:

1. Analysis type: The general bioinformatics analysis type (GWAS, QTL mapping, differential expression testing, etc.).
1. Task type: The stage of data analysis (e.g., munging, analysis, post-processing, or visualization).
1. Analysis software: Software used for the analysis (PLINK, Matrix eQTL, DESeq2, etc.)

```
# Example directory structure
.
└── fwgwas
    ├── analysis
    │   └── lsmm
    ├── munging
    │   └── lsmm
    │       └── bedgraph_overlap.py
    ├── post_processing
    │   └── lsmm
    └── visualization
        └── lsmm
```

