#! /usr/bin/python

#Analyze read depth in comparison to transcription start

import sys
import argparse
from subprocess import call
import numpy
import scipy
import scipy.stats
import os.path
import multiprocessing


# Calculate mean and confidence intervals ###################################################################################
def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*numpy.array(data)
    n = len(a)
    m, se = numpy.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1+confidence)/2., n-1)
    return m, m-h, m+h

# Calculate mean coverage in 10000bp upstream sequence ###################################################################################
def calcLocalMean(position, chrom, thread_number):
    
    LOG2 = open(args.norm_log2,"r")
    header = LOG2.readline()
    norm = None
    for line in LOG2.readlines():
       info = line.split("\t")
       if chrom == info[0] and position > int(info[1]) and position < int(info[2]):
           norm = numpy.power(2,float(info[3].rstrip()))
           break
    if norm:
        return norm
    else:
        return 1
    LOG2.close()


# Run this for each thread ###################################################################################
def thread_proc(q,thread_number,transcript_list):
    sys.stderr.write("Thread "+str(thread_number)+" started\n")
    sys.stderr.flush()
    coverage_dict=dict()

    #create a list of visited TSS not to count some more than once
    tss_visited = list()

    line_count = 0
    skipped = 0

    #iterate through transcript list from UCSC genome browser
    for transcript in transcript_list:
        coverage_list=dict()
        gene_name = transcript.split("\t")[4].rstrip()
        for i in range(-args.start,0):
            coverage_list[i] = 0
        for i in range(0,args.end+1):
            coverage_list[i] = 0
        line_count += 1
        if (line_count  % 100 == 0):
            sys.stderr.write("\rThread "+str(thread_number)+"\t"+str(line_count)+" genes analyzed")
            sys.stderr.flush()
        #transcription starts are marked at txEnd Field of RefGene Txt file for reverse transcribed genes
        if transcript.split()[1] == '+':
            forward = True
            chrom = transcript.split()[0]
            pos = int(transcript.split()[2])
        else:
            forward = False
            chrom = transcript.split()[0]
            pos = int(transcript.split()[3])
        if chrom+"_"+str(pos) in tss_visited:
            continue
        tss_visited.append(chrom+"_"+str(pos))
        key = gene_name+"\t"+chrom+"_"+str(pos)
        TMP_COVERAGE_BED = open("intermediate/"+os.path.basename(args.bam_file)+str(thread_number)+"tmp_coverage.bed","w")
        call(["samtools","depth","-r",chrom+":"+str(pos-args.start)+"-"+str(pos+args.end),args.bam_file],stdout=TMP_COVERAGE_BED)
        TMP_COVERAGE_BED.close()
    
        TMP_COVERAGE_BED_OUTPUT = open("intermediate/"+os.path.basename(args.bam_file)+str(thread_number)+"tmp_coverage.bed","r")
        content = TMP_COVERAGE_BED_OUTPUT.readlines()
        for i in range(pos-args.start,pos+args.start):
            found = False
            for line in content:
                chrom_found = line.split()[0]
                pos_found = int(line.split()[1])
                if pos_found == i:
                    found = True
                    coverage = int(line.split()[2])
                    if forward:
                        coverage_list[i-pos] = coverage
                    elif not forward:
                        coverage_list[-(i-pos)] = coverage
                    continue
            if not found:
                if forward:
                    coverage_list[i-pos] = 0
                elif not forward:      
                    coverage_list[-(i-pos)] = 0
        TMP_COVERAGE_BED_OUTPUT.close()
        call(["rm","intermediate/"+os.path.basename(args.bam_file)+str(thread_number)+"tmp_coverage.bed"])    
        mean,lower,upper = mean_confidence_interval(coverage_list.values())
        if args.norm:
            control_region_mean = calcLocalMean(pos, chrom, thread_number)
            if control_region_mean == 0:
                control_region_mean = 0.1
            coverage_dict[key] = float(mean)/float(control_region_mean)
        else:
            coverage_dict[key] = mean
   
    
    sys.stderr.write("\rThread "+str(thread_number)+"\t finished\n")
    sys.stderr.flush() 
    q.put(coverage_dict)
#######################################################################################################

# Parse command line arguments ###################################################################################
parser = argparse.ArgumentParser(description='Analyze read depth in comparison to transcription start')
parser.add_argument('-rg','--ref-gene', dest='refgene_file', 
                   help='RefGene file, transcription start should be stored in column 1-3',required=True)
parser.add_argument('-b','--bam', dest='bam_file',
                   help='BAM file',required=True)
parser.add_argument('-s','--start', dest='start',
                   help='Start analyzing coverage at this point before TSS [default:1000]',default=1000,type=int)
parser.add_argument('-e','--end', dest='end',
                   help='Stop analyzing coverage at this point after TSS [default:1000]',default=1000,type=int)
parser.add_argument('-t','--threads', dest='threads',
                   help='Threads to use for computation [default:1]',default=1,type=int)
parser.add_argument('-X','--include-X', dest='x_chrom',
                   help='Include X chromosome genes [default:False]',action="store_true")
parser.add_argument('-gl','--gene-list', dest='gene_list',
                   help='List of gene names',required=True)
parser.add_argument('-norm','--normalize', dest='norm',
                   help='Normalize by local coverage',action="store_true")
parser.add_argument('-norm-file','--normalize-file', dest='norm_log2',
                   help='Normalize by local copynumber from this file')

args = parser.parse_args()
sys.stderr.write("Bam file: "+args.bam_file+"\n")
sys.stderr.write("RefGene file: "+args.refgene_file+"\n")
sys.stderr.write("Genelist: "+str(args.gene_list)+"\n")
sys.stderr.write("Threads: "+str(args.threads)+"\n")
if args.norm and not os.path.isfile(args.norm_log2):
   print "Can't open LOG2-file"
   sys.exit(1)

###############################################################################################
# Analyze data ###################################################################################

try:
    REFGENE = open(args.refgene_file,"r")
except:
    print "Fail to open files specified"
    sys.exit(1)
target_genes = list()
try:
    GENELIST_H = open(args.gene_list,"r")
    for item in GENELIST_H.readlines():
        target_genes.append(item.rstrip())
    GENELIST_H.close()
except:
    print "Failed to open genelist"
    sys.exit(1)
    
#filter genes from genelist if specified
header = REFGENE.readline()
refgene_content = REFGENE.readlines()
target_genes_count = 0
target_content = list()
for i in refgene_content:
    chrom = i.split()[0]
    if chrom.find("_") != -1:
        continue
    if not args.x_chrom:
        if chrom.find("X") != -1:
            continue
    if chrom.find("Y") != -1:
        continue
    if i.split()[4].rstrip() in target_genes:
        target_content.append(i)

#initialize input data
gene_count = 0
sys.stderr.write("\n")
sys.stderr.flush()
thread_input_list = dict()
thread_coverage_list = dict()

for thread_number in range(0,args.threads):
    thread_input_list[thread_number] = list()
    thread_coverage_list[thread_number] = dict()

#dispatch input data
max_gene = len(target_content)
partition_point = max_gene/ args.threads
for thread_number in range(0,args.threads):
    thread_input_list[thread_number] = target_content[(thread_number*partition_point):(thread_number+1)*partition_point]

sys.stderr.write(str(len(thread_input_list[0]))+" genes per thread\n")
sys.stderr.write("--------------------------------------------------\n")
sys.stderr.flush()

#start multiple processes
processes = dict()
queues = dict()
for thread in range(0,args.threads):
    queues[thread] = multiprocessing.Queue()
    processes[thread] = multiprocessing.Process(target=thread_proc,args=(queues[thread],thread,thread_input_list[thread]))
    processes[thread].start()

#wait for processes to finish
for thread in range(0,args.threads):
    thread_coverage_list[thread] = queues[thread].get()
    processes[thread].join()

#collect all data
coverage_dict_all=dict()
for thread in range(0,args.threads):
    for i in  thread_coverage_list[thread].keys():
        coverage_dict_all[i] = thread_coverage_list[thread][i]

print "Gene\tTSS\tCoverage\n"
for gene in coverage_dict_all.keys():
    print gene+"\t"+str(coverage_dict_all[gene])


 
