#! /usr/bin/python

#Generate single gene wig files to upload to UCSC Genome browser from BAM file
import sys
import argparse
from subprocess import call
import os

# Parse command line arguments ###################################################################################
parser = argparse.ArgumentParser(description='Analyze read depth in comparison to transcription start')
parser.add_argument('-rg','--ref-gene', dest='refgene_file', 
                   help='RefGene file, transcription start should be stored in column 1-3',required=True)
parser.add_argument('-w','--wig-file', dest='wig_file', 
                   help='Wig file to save results',required=True)
parser.add_argument('-g','--gene', dest='gene', 
                   help='Gene to create WIG file of',required=True)
parser.add_argument('-b','--bam', dest='bam_file',
                   help='BAM file',required=True)
parser.add_argument('-s','--start', dest='start',
                   help='Start analyzing coverage at this point before TSS [default:1000]',default=1000,type=int)
parser.add_argument('-e','--end', dest='end',
                   help='Stop analyzing coverage at this point after TSS [default:1000]',default=1000,type=int)

args = parser.parse_args()
sys.stderr.write("Bam file: "+args.bam_file+"\n")
sys.stderr.write("RefGene file: "+args.refgene_file+"\n")
sys.stderr.write("Gene: "+args.gene+"\n")
sys.stderr.flush()

###############################################################################################
# Read-in refseq data ###################################################################################

try:
    REFGENE = open(args.refgene_file,"r")
except:
    print "Fail to open refGene file specified"
    sys.exit(1) 

gene_starts = dict()
gene_chroms = dict()
gene_dir = dict()
header = REFGENE.readline()
for line in REFGENE.readlines():
    info = line.split("\t")   
    name = info[4].rstrip()
    if name in gene_chroms.keys():
        continue
    gene_chroms[name] = info[0]
    gene_dir[name] = info[1]
    if info[1] == "+":
        gene_starts[name] = int(info[2])
    else:
        gene_starts[name] = int(info[3])

if args.gene not in gene_chroms.keys():
    print "Gene not found in RefGene file"
    sys.exit(1)

chrom = gene_chroms[args.gene]
start = gene_starts[args.gene]-args.start
stop = gene_starts[args.gene]+args.end

call(["./scripts/bam_to_wiggle.py","--chrom="+chrom,"--start="+str(start),"--end="+str(stop),args.bam_file])
intermediate_wig="./intermediate/Wiggle/"+os.path.basename(args.bam_file)+"_"+args.gene
call(["./scripts/bigWigToWig",args.bam_file[:-4]+".bigwig",intermediate_wig])
call(["rm",args.bam_file[:-4]+".bigwig"])
OUTPUT=open(args.wig_file,"w")
OUTPUT.write("track type=wiggle_0 name="+os.path.basename(args.bam_file)+"_"+args.gene+"\n")
INT_WIG = open(intermediate_wig,"r")
for line in INT_WIG.readlines():
    OUTPUT.write(line)

OUTPUT.close()
INT_WIG.close()
