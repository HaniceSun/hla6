<p align="left">
<img src="assets/logo.png" alt="dustle logo" width="150"/>
</p>

# HLA Imputation (hla6)

## Overview

hla6 is designed as an integrated suite of tools for inferring HLA types from either genotyping array data or whole-genome sequencing data. It incorporates the core functionality of SNP2HLA, the CNN-based DEEP*HLA, and the Transformer-based HLARIMNT, providing a unified framework to harmonize their outputs and systematically evaluate their performance, by varying key parameters such as ancestry reference panels and benchmarking against our in-house gold-standard HLA typing results, with the goal of further improving the deep learning models. In addition, hla6 is being applied to the All of Us cohort to investigate type 1 diabetes risk alleles across diverse populations. The "6" in the name reflects the six major HLA loci - HLA-A, HLA-B, and HLA-C from HLA class I, and HLA-DP, HLA-DQ, and HLA-DR from HLA class II - located on human chromosome 6.

## Installation

- using conda

```
git clone git@github.com:HaniceSun/hla6.git
cd hla6
conda env create -f environment.yml
conda activate hla6
```

# Quick Start

```
hla6 -h
hla6 run-snp2hla
```

## Author and License

**Author:** Han Sun

**Email:** hansun@stanford.edu

**License:** [MIT License](LICENSE)
