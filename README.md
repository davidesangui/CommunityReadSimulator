# CommunityReadSimulator
Simple script to generate paired-end shotgun metagenomic sequencing reads from a set of genomes representing a microbial community.
## Installation
To use CommunityReadSimulator, simply download the file CommunityReadSimulator.py from this repository and launch it with Python3. The art_illumina executable is also required, which can be obtained from this repository (linux64 version) or from https://www.niehs.nih.gov/research/resources/software/biostatistics/art.
The script requires numpy, which can be obtained via pip:
```
pip install numpy
```
## How it works
CommunityReadSimulator is essentially a wrapper of ART Illumina read simulator:  
(Weichun Huang, Leping Li, Jason R Myers, and Gabor T Marth. ART: a next-generation sequencing read simulator, Bioinformatics (2012) 28 (4): 593-594).  
Given a set of genomes in fasta format, a relative abundance profile is sampled from a user-defined distribution and paired-end sequencing reads are generated with art_illumina from each genome, according to the relative abundance profile and the desired total number of reads.  
The output consists of two fastq files (forward and reverse reads) and an alignment file (.aln format) which are obtained through the concatenation of the genome-level art_illumina outputs. Moreover, also the ground truth relative abundances are given as output.
