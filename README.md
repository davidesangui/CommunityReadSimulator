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
When multiple cores are available, sequencing reads can be generated from different genomes in parallel.  
The output consists of two fastq files (forward and reverse reads) and an alignment file (.aln format) which are obtained through the concatenation of the genome-level art_illumina outputs. Moreover, also the ground truth relative abundances are given as output.
## Usage
CommunityReadSimulator.py requires four mandatory arguments:
- ```-numReads```: (approximately) the total number of reads that will be generated in the simulated metagenome
- ```-path_to_gens```: path to the directory containing the genomes of the simulated metagenome.
- ```-path_to_ART_illumina```: path to the executable of ART illumina (art_illumina_linux64 in this repository). Only ART Illumina is currently supported, hence no other ART executable should be given.
- ```-output_folder```: the directory to store all output files. If it does not exist, it will be automatically created.
### Example
The data in this repository will be used as example. An example output is stored in the example_output folder. Given the stochasticity in relative abundance sampling and in ART Illumina, your output will be different by the one in this repository. 
```python3 CommunityReadSimulator.py -numReads 10000 -path_to_gens genomes/ -path_to_ART_illumina art_illumina_linux64 -output_folder example_output -p 3```
```-p``` defines the number of processes for parallelization. Setting a number greater than the number of genomes is useless.
