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
    coverage_list=list()
    #pos = position + 100000
    #step1 get coverage values of TSS-3000 to TSS -1000
    start1 = position - 3000
    end1 = position - 1000
    TMP_COVERAGE_BED = open(args.temp_dir+os.path.basename(args.bam_file)+str(thread_number)+"tmp_coverage1_10000bpupstream.bed","w")
    call(["samtools","depth","-r",chrom+":"+str(start1)+"-"+str(end1),args.bam_file],stdout=TMP_COVERAGE_BED)
    TMP_COVERAGE_BED.close()
    TMP_COVERAGE_BED_OUTPUT = open(args.temp_dir+os.path.basename(args.bam_file)+str(thread_number)+"tmp_coverage1_10000bpupstream.bed","r")
    content = TMP_COVERAGE_BED_OUTPUT.readlines()
    for line in content:
        chrom_found = line.split()[0]
        pos_found = int(line.split()[1])
        coverage = int(line.split()[2])
        coverage_list.append(coverage)
    call(["rm",args.temp_dir+os.path.basename(args.bam_file)+str(thread_number)+"tmp_coverage1_10000bpupstream.bed"])    
    TMP_COVERAGE_BED_OUTPUT.close()

    #step2 get coverage values of TSS+1000 to TSS+3000
    start2 = position + 1000
    end2 = position + 3000
    TMP_COVERAGE_BED = open(args.temp_dir+os.path.basename(args.bam_file)+str(thread_number)+"tmp_coverage2_10000bpupstream.bed","w")
    call(["samtools","depth","-r",chrom+":"+str(start2)+"-"+str(end2),args.bam_file],stdout=TMP_COVERAGE_BED)
    TMP_COVERAGE_BED.close()
    TMP_COVERAGE_BED_OUTPUT = open(args.temp_dir+os.path.basename(args.bam_file)+str(thread_number)+"tmp_coverage2_10000bpupstream.bed","r")
    content = TMP_COVERAGE_BED_OUTPUT.readlines()
    for line in content:
        chrom_found = line.split()[0]
        pos_found = int(line.split()[1])
        coverage = int(line.split()[2])
        coverage_list.append(coverage)
    call(["rm",args.temp_dir+os.path.basename(args.bam_file)+str(thread_number)+"tmp_coverage2_10000bpupstream.bed"])    
    TMP_COVERAGE_BED_OUTPUT.close()
    
    if len(coverage_list) == 0:
        mean = 0.1
    else:
        mean,lower,upper = mean_confidence_interval(coverage_list)
    return mean


# Run this for each thread ###################################################################################
def thread_proc(q,thread_number,transcript_list):
    sys.stderr.write("Thread "+str(thread_number)+" started\n")
    sys.stderr.flush()
    position_lists=dict()

    #create a list of visited TSS not to count some more than once
    tss_visited = list()

    for i in range(-args.start,0):
        position_lists[i] = list()
    for i in range(0,args.end+1):
        position_lists[i] = list()

    line_count = 0
    skipped = 0
    for transcript in transcript_list:
        line_count += 1
        if (line_count  % 100 == 0):
            sys.stderr.write("\rThread "+str(thread_number)+"\t"+str(line_count)+" genes analyzed")
            sys.stderr.flush()
        #use only transcripts on forward strand at that time
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
        if chrom.find("_") != -1:
            continue

        control_region_mean = 1
        if args.norm:
            control_region_mean = calcLocalMean(pos, chrom, thread_number)
        tss_visited.append(chrom+"_"+str(pos))

        TMP_COVERAGE_BED = open(args.temp_dir+os.path.basename(args.bam_file)+str(thread_number)+"tmp_coverage.bed","w")
        call(["samtools","depth","-r",chrom+":"+str(pos-args.start)+"-"+str(pos+args.end),args.bam_file],stdout=TMP_COVERAGE_BED)
        TMP_COVERAGE_BED.close()
    
        TMP_COVERAGE_BED_OUTPUT = open(args.temp_dir+os.path.basename(args.bam_file)+str(thread_number)+"tmp_coverage.bed","r")
        content = TMP_COVERAGE_BED_OUTPUT.readlines()
        for i in range(pos-args.start,pos+args.start):
            found = False
            for line in content:
                chrom_found = line.split()[0]
                pos_found = int(line.split()[1])
                if pos_found == i:
                    found = True
                    #if normalization is switched off control_Region_mean = 1
                    coverage = float(line.split()[2])/ control_region_mean
                    if forward:
                        position_lists[i-pos].append(coverage )
                    elif not forward:
                        position_lists[-(i-pos)].append(coverage)
                    continue
            if not found:
                if forward:
                    position_lists[i-pos].append(0.)
                elif not forward:      
                    position_lists[-(i-pos)].append(0.)
        TMP_COVERAGE_BED_OUTPUT.close()
        call(["rm",args.temp_dir+os.path.basename(args.bam_file)+str(thread_number)+"tmp_coverage.bed"])    

    q.put(position_lists)
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
parser.add_argument('-m','--max-genes', dest='max_genes',
                   help='Maximum numbers of genes to analyze [default:all]',default=0,type=int)
parser.add_argument('-gl','--gene-list', dest='gene_list',
                   help='List of gene names to search for',default="na")
parser.add_argument('-norm','--normalize', dest='norm',
                   help='Normalize by local coverage at 1.000bp upstream',action="store_true")
parser.add_argument('-tmp','--temp-dir', dest='temp_dir',
                   help='Temporary Directory',default="./intermediate/")


args = parser.parse_args()
if args.temp_dir[-1:] != "/":
    args.temp_dir = args.temp_dir+"/"
sys.stderr.write("Bam file: "+args.bam_file+"\n")
sys.stderr.write("RefGene file: "+args.refgene_file+"\n")
sys.stderr.write("Genes: "+str(args.max_genes)+"\n")
sys.stderr.write("Threads: "+str(args.threads)+"\n")
genelist = False
if args.gene_list != "na":
    genelist = True
    sys.stderr.write("Genelist: "+str(args.gene_list)+"\n")
sys.stderr.flush()

###############################################################################################
# Analyze data ###################################################################################

try:
    REFGENE = open(args.refgene_file,"r")
except:
    print "Fail to open files specified"
    sys.exit(1)
target_genes = list()
if genelist:
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
if genelist:
    for i in refgene_content:
        if i.split()[4].rstrip() in target_genes:
            target_content.append(i)
else:
    target_content = refgene_content

#initialize input data
gene_count = 0
sys.stderr.write("\n")
sys.stderr.flush()
thread_input_list = dict()
thread_position_list = dict()
for thread_number in range(0,args.threads):
    thread_input_list[thread_number] = list()
    thread_position_list[thread_number] = list()

#dispatch input data
max_gene = args.max_genes
if args.max_genes > len(target_content) or args.max_genes == 0:
    max_gene = len(target_content)
sys.stderr.write("Searching for "+str(max_gene)+" genes\n")
sys.stderr.flush()
partition_point = max_gene/ args.threads
for thread_number in range(0,args.threads):
    thread_input_list[thread_number] = target_content[(thread_number*partition_point):(thread_number+1)*partition_point]

sys.stderr.write("--------------------------------------------------\n")
sys.stderr.write(str(len(thread_input_list[0]))+" genes per thread\n")
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
    thread_position_list[thread] = queues[thread].get()
    processes[thread].join()

#collect all data
position_lists_all=dict()
for i in range(-args.start,0):
    position_lists_all[i] = list()
    for thread in range(0,args.threads):
        position_lists_all[i].extend(thread_position_list[thread][i])
for i in range(0,args.end+1):
    position_lists_all[i] = list()
    for thread in range(0,args.threads):
        position_lists_all[i].extend(thread_position_list[thread][i])

#output data
sys.stderr.write("--------------------------------------------------\n")
sys.stderr.write("\n"+str(len(position_lists_all[0]))+" TSS analyzed\n")
sys.stderr.flush()
print "Position\tMean Cov\tLowerBound\tUpperBound\tTSS analyzed\n"
for i in range(-args.start,args.end+1):
    if len(position_lists_all[i]) > 3:
        mean,lower_bound, upper_bound = mean_confidence_interval(position_lists_all[i])
        print str(i)+"\t"+str(mean)+"\t"+str(lower_bound)+"\t"+str(upper_bound)+"\t"+str(len(position_lists_all[0]))
    else:
        mean=numpy.mean(position_lists_all[i])
        print str(i)+"\t"+str(mean)+"\t"+str(len(position_lists_all[0]))
