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
The output consists of two fastq files (forward and reverse reads) and two alignment files (.aln format and SAM format) which are obtained through the concatenation of the genome-level art_illumina outputs. Moreover, also the ground truth relative abundances are given as output.
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
## Advanced usage
### Abundance distribution
Relative abundance profiles can be sampled from four different distributions: Lognormal, Exponential, Uniform, Equal.  
If "Equal" distribution is selected, then each genome will have the same relative abundance, i.e., if there are three genomes, then each will have abundance equal to 1/3.  
The distribution can be defined with the option ```-distribution```. The default distribution is Lognormal, which is typically fitted by microbial abundance distirbutions in real communities. 
Mu and Sigma parameter of the Lognormal distribution can be defined with the options ```-mu``` and ```-sigma``` (default values: 1.0 and 2.0 respectively). When the Exponential distribution is selected, then ```-mu``` option defines the lambda parameter of the Exponential distribution. In this case ```-sigma``` option is irrelevant. For all other distributions, the options ```-mu``` and ```-sigma```are irrelevant.  
**Example:**  
```python3 CommunityReadSimulator.py -numReads 10000 -path_to_gens genomes/ -path_to_ART_illumina art_illumina_linux64 -output_folder example_output -p 3 -distribution exponential -mu 1.5```

The user can provide a fixed relative abundance distribution with the option ```-abundances```. The argument of this option is a text file where each line lists a fasta file and, tab-separated, the corresponding **percent** relative abundance. The relative abundance values listed in this file should sum to 100. An example file is given in this repository (predefined_rel_abundances.tsv).  
**Example:**  
```python3 CommunityReadSimulator.py -numReads 10000 -path_to_gens genomes/ -path_to_ART_illumina art_illumina_linux64 -output_folder example_output -p 3 -abundances predefined_rel_abundances.tsv```
### ART Illumina parameters
Via CommunityReadSimulator.py, the user can interact with ART Illumina in a limited way. In particular, only three options can be given to CommunityReadSimulator.py that modify the behaviour of ART Illumina:
- ```-readLength``` Length of simulated sequencing reads generated by ART Illumina. Default: 100. Corresponding art_illumina option: ```-l```
- ```-insertSizeMean``` Mean length of DNA fragment for paired-end reads. Default: 1000. Corresponding art_illumina option: ```-m```
- ```-insertSizeStd``` Standard deviation of length of DNA fragment for paired-end reads. Default: 100. Corresponding art_illumina option: ```-s```
**Example:**  
```python3 CommunityReadSimulator.py -numReads 10000 -path_to_gens genomes/ -path_to_ART_illumina art_illumina_linux64 -output_folder example_output -p 3 -readLength 150 -insertSizeMean 500 -insertSizeStd 50```

However, users may require a more refined ART Illumina command to better mimic their conditions. CommunityReadSimulator.py allows the user to pass in input a text file where the first line corresponds to a user-defined complete art_illumina command, which will be used by CommunityReadSimulator.py for running art_illumina. The flag to pass this input file is ```-ART_command```.  
The command must include the aforementioned ```-l``` ```-m``` ```-s``` options, which must be indicated by the stated flags (-l, -m, -s) and not by possible synonyms allowed by art_illumina (--len, --mflen, --sdev). All other options of art_illumina can be added without restrictions.   
All art_illumina parameters given to CommunityReadSimulator.py will be ignored. Some modifications to the command will be performed: ```-o``` and ```-i``` options will be added or modified according to the arguments given to CommunityReadSimulator.py; ```-f``` option will be added or modified according to the sampled relative abundance; ```-c``` option will be removed; ```-p```, ```-sam```, ```-q``` options will be added if missing.
An example file is given in this repository (art_command.txt).  
**Example:**  
```python3 CommunityReadSimulator.py -numReads 10000 -path_to_gens genomes/ -path_to_ART_illumina art_illumina_linux64 -output_folder example_output -p 3 -ART_command art_command.txt```
