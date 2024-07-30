import os
import argparse
import numpy as np
import multiprocessing
import time
import sys
import subprocess

parser=argparse.ArgumentParser()
parser.add_argument("-numReads",help="total number of reads to be generated. Each read pair corresponds to two reads.",required=True)
parser.add_argument("-path_to_gens",help="path to genome folder",required=True)
group0=parser.add_mutually_exclusive_group(required=True)
group0.add_argument("-path_to_ART_illumina",help="path to the excutable of ART_illumina")
parser.add_argument("-output_folder",help="output folder to store output files. Created if missing.",required=True)
parser.add_argument("-p",help="number of processes for parallelization. Default: 1",type=int,default=1,metavar='')
group=parser.add_mutually_exclusive_group()
group.add_argument("-distribution",help="Distribution to sample relative abundances from. If < equal > then relative abundances will be the same for each genome. Default: lognormal",choices=['lognormal','exponential','uniform','equal'],default='lognormal')
parser.add_argument("-mu",help="mu parameter of lognormal distribution. If <-distribution> is <exponential> this value correspond to the lambda parameter of the exponential distribution. Irrelevant when <-distribution> is <uniform> or <equal>. Default: 1.0",type=float,default=1.0,metavar='')
parser.add_argument("-sigma",help="sigma parameter of lognormal distribution. Irrelevant when <-distribution> is <exponential> , <uniform> or <equal>. Default: 2.0",type=float,default=2.0,metavar='')
group.add_argument("-abundances",help='Path to a text file containing user-defined relative abundances (percent), where each line lists a fasta file in the genome folder and, tab-separated, the corresponding relative abundance (percent). The values should sum to 100.',metavar='')
parser.add_argument("-readLength",help="ART_illumina parameter. Length of simulated reads",default=100,metavar='',type=int)
parser.add_argument("-insertSizeMean",help="ART_illumina parameter. Mean length of DNA fragment for paired-end reads",default=1000,metavar='')
parser.add_argument("-insertSizeStd",help="ART_illumina parameter. Standard deviation of length of DNA fragment for paired-end reads",default=100,metavar='')
group0.add_argument("-ART_command",help="Path to a text file where the first line corresponds to a user-defined complete command to run ART_illumina. This command will be used to run ART_illumina for read generation. All ART_illumina parameters given to "+sys.argv[0]+" will be ignored. Some modifications to the command will be performed: -o and -i options will be added or modified according to the arguments given to "+sys.argv[0]+"; -f option will be added or modified according to the sampled relative abundance; -c option will be removed; -p, -sam, -q options will be added if missing. ",metavar='')
args=parser.parse_args()

def chunks(l, processi):
    lists=[[] for x in range(processi)]
    ind=0
    for ob in l:
        lists[ind].append(ob)
        ind+=1
        if ind==processi:
            ind=0
    end=len(lists)
    for i in range(len(lists)):
        if len(lists[i])==0:
            end=i
            break
    return(lists[:end])

def main(allargs):
    args,genomes,readsfrac_list,genomes_lengths,totreads,artCommand=allargs
    for i in range(len(genomes)):
        genome=genomes[i]
        nreads=round(readsfrac_list[i]*totreads)
        coverage=nreads*(args.readLength if not artCommand else int(artCommand.split(' ')[artCommand.split(' ').index('-l')+1]))/genomes_lengths[i] 
        # launch artbin
        if not artCommand:
            out=subprocess.check_output("{} -p -l {} -m {} -s {} -f {} -o {} -i {} -q -sam".format(args.path_to_ART_illumina,args.readLength,args.insertSizeMean,args.insertSizeStd,coverage,args.output_folder+'/'+genome,args.path_to_gens+'/'+genome), shell=True)
        else:  
            # adjust user-defined ART command
            arguments=artCommand.split(' ')
            if '-o' in arguments: arguments[arguments.index('-o')+1]=args.output_folder+'/'+genome
            elif '--out' in arguments!=-1: arguments[arguments.index('--out')+1]=args.output_folder+'/'+genome
            else: arguments+=['-o',args.output_folder+'/'+genome]
            if '-i' in arguments: arguments[arguments.index('-i')+1]=args.path_to_gens+'/'+genome
            elif '--in' in arguments: arguments[arguments.index('--in')+1]=args.path_to_gens+'/'+genome
            else: arguments+=['-i',args.path_to_gens+'/'+genome]
            if '-f' in arguments: arguments[arguments.index('-f')+1]=str(coverage)
            elif '--fcov' in arguments: arguments[arguments.index('--fcov')+1]=str(coverage)
            else: arguments+=['-f',str(coverage)]
            if '-q' not in arguments and '--quiet' not in arguments: arguments.append('-q')
            if '-p' not in arguments and '--paired' not in arguments: arguments.append('-p')
            if '-sam' not in arguments and '--samout' not in arguments: arguments.append('-sam')
            arguments=[arguments[i] for i in range(len(arguments)) if arguments[i]!='-c' and arguments[i-1]!='-c' and arguments[i]!='--rcount' and arguments[i-1]!='--rcount']
            try:
                subprocess.check_output(' '.join(arguments), shell=True)
            except subprocess.CalledProcessError as e:
                print('----- ART ERROR -----')
                print('The Fasta file causing the error was',genome,'\n')
                

if __name__=='__main__':

    c=1 # counter of log messages
    # check mu and sigma greater than zero
    if args.mu<=0.0 or args.sigma<=0.0: print('Error: -mu and -sigma must be greater than zero.\nAborting...'); exit(1)
    # check output folder/create it
    if not os.path.isdir(args.output_folder):
        if os.path.exists(args.output_folder): print('Error:',args.output_folder,'exists and is not a directory.\nAborting...'); exit(1)
        else: os.mkdir(args.output_folder)

    # abundances distribution
    if not args.abundances:
        gens=os.listdir(args.path_to_gens)
        if args.distribution=='lognormal': rawabs=np.random.lognormal(float(args.mu),float(args.sigma),len(gens))
        elif args.distribution=='exponential': rawabs=np.random.exponential(float(args.mu),len(gens))
        elif args.distribution=='uniform': rawabs=np.random.random(len(gens))
        elif args.distribution=='equal': rawabs=np.ones(len(gens))
        abundances=rawabs/np.sum(rawabs) # relative abundances
    else:
        f=open(args.abundances)
        abundances,gens=[],[]
        for line in f:
            a=line.strip().split('\t')
            abundances.append(float(a[1])/100); gens.append(a[0])
        abundances=np.array(abundances)
    # save abundances to output file
    new=open('{}/relative_abundances.tsv'.format(args.output_folder),'w')
    print('genome\trelative_abundance_%',file=new)
    for gen,abu in zip(gens,abundances):
        print(gen+'\t'+str(abu*100),file=new)
    new.close()

    # get art command if given by user
    if not args.ART_command: artCommand=''
    else: 
        f=open(args.ART_command) 
        artCommand=f.readline().strip(); f.close()
        splitted_command=artCommand.split(' ')
        if not all([x in splitted_command for x in ['-l','-m','-s']]): print('Error: -l, -m and -s options must all be specified in -ART_command.\nAborting...'); exit(1)
        
    # get genome lengths
    print(str(c)+') Calculating genome lengths...'); c+=1
    lengths=[]
    for gen in gens:
        l=0
        f=open(args.path_to_gens+'/'+gen)
        for line in f:
            if line[0]!='>':
                l+=len(line.strip())
        lengths.append(l)
    lenghts=np.array(lengths)
    read_fractions=lenghts*abundances/np.sum(lenghts*abundances)
    lenghts=list(lengths)

    # read generation with parallelization
    print(str(c)+') Running ART_illumina for read generation. This may take a while...'); c+=1
    div_gens=chunks(gens,int(args.p))
    div_readfracts=chunks(read_fractions,int(args.p))
    div_lenghts=chunks(lengths,int(args.p))
    jobs=[]
    with multiprocessing.Pool(args.p) as pool:
        pool.map(main,[(args,div_gens[i],div_readfracts[i],div_lenghts[i],int(args.numReads),artCommand) for i in range(len(div_gens))])

    # remove useless file and cat single-genome outputs into metagenome outputs 
    print(str(c)+') Joining fastq files in single metagenome'); c+=1
    new1=open('{}/metagenome_R1.fq'.format(args.output_folder),'w')
    new2=open('{}/metagenome_R2.fq'.format(args.output_folder),'w')
    newal=open('{}/metagenomeTMP.aln'.format(args.output_folder),'w')
    newsam=open('{}/metagenomeTMP.sam'.format(args.output_folder),'w')
    alnHeader,samHeader=None,None
    alnContigs,samContigs=[],[]
    for gen in gens:
        # cat fastq
        fq=open('{}/{}1.fq'.format(args.output_folder,gen))
        for line in fq: print(line.strip(),file=new1)
        fq.close(); os.remove('{}/{}1.fq'.format(args.output_folder,gen))
        fq=open('{}/{}2.fq'.format(args.output_folder,gen))
        for line in fq: print(line.strip(),file=new2)
        fq.close(); os.remove('{}/{}2.fq'.format(args.output_folder,gen))
        # cat aln
        fal=open('{}/{}1.aln'.format(args.output_folder,gen))
        for line in fal:
            if line.startswith('##') and alnHeader is None: alnHeader=line.strip()
            if line.startswith('@SQ'): alnContigs.append('@SQ'+'\t'+line.strip().split('\t')[1].split(' ')[0]+'\t'+line.strip().split('\t')[2])
            if not (line.startswith('@') or line.startswith('#')): print(line.strip(),file=newal)
        fal.close(); os.remove('{}/{}1.aln'.format(args.output_folder,gen))
        fal=open('{}/{}2.aln'.format(args.output_folder,gen))
        for line in fal: 
            if not (line.startswith('@') or line.startswith('#')): print(line.strip(),file=newal)
        fal.close();  os.remove('{}/{}2.aln'.format(args.output_folder,gen))
        # cat sam
        sam=open('{}/{}.sam'.format(args.output_folder,gen))
        for line in sam:
            if line.startswith('@HD') and samHeader is None: samHeader=line.strip()
            if line.startswith('@SQ'): samContigs.append(('@SQ'+'\t'+line.strip().split('\t')[1].split(' ')[0]+'\t'+line.strip().split('\t')[2]))
            if not line.startswith('@'): print(line.strip(),file=newsam)
        sam.close(); os.remove('{}/{}.sam'.format(args.output_folder,gen)) 
    
    new1.close();new2.close();newal.close();newsam.close()
    # appending header to aln file
    fal=open('{}/metagenomeTMP.aln'.format(args.output_folder))
    newal=open('{}/metagenome.aln'.format(args.output_folder),'w')
    print(alnHeader,file=newal)
    for contig in alnContigs: print(contig,file=newal)
    print('##Header End',file=newal)
    for line in fal: print(line.strip(),file=newal)
    fal.close(); newal.close()
    # appending header to sam file
    sam=open('{}/metagenomeTMP.sam'.format(args.output_folder))
    newsam=open('{}/metagenome.sam'.format(args.output_folder),'w')
    print(samHeader,file=newsam)
    for contig in samContigs: print(contig,file=newsam)
    for line in sam: print(line.strip(),file=newsam)
    sam.close(); newsam.close()
    os.remove('{}/metagenomeTMP.sam'.format(args.output_folder)); os.remove('{}/metagenomeTMP.aln'.format(args.output_folder))
