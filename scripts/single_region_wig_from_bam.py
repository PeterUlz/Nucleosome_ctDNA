#! /usr/bin/python

#Generate single gene wig files to upload to UCSC Genome browser from BAM file
import sys
import argparse
from subprocess import call
import os

# Parse command line arguments ###################################################################################
parser = argparse.ArgumentParser(description='Analyze read depth in comparison to transcription start')
parser.add_argument('-chr','--chromosome', dest='chrom', 
                   help='Chromosome',required=True)
parser.add_argument('-s','--start', dest='start',
                   help='Start position',required=True,type=int)
parser.add_argument('-e','--end', dest='end',
                   help='End position',required=True,type=int)
parser.add_argument('-b','--bam', dest='bam_file',
                   help='BAM file',required=True)
parser.add_argument('-w','--wig-file', dest='wig_file', 
                   help='Wig file to save results',required=True)

args = parser.parse_args()
sys.stderr.write("Bam file: "+args.bam_file+"\n")
sys.stderr.write("Region: "+args.chrom+":"+str(args.start)+"-"+str(args.end)+"\n")
sys.stderr.flush()

###############################################################################################
# Read-in refseq data ###################################################################################

chrom = args.chrom
start = str(args.start)
stop = str(args.end)
region=args.chrom+":"+str(args.start)+"-"+str(args.end)

call(["./scripts/bam_to_wiggle.py","--chrom="+chrom,"--start="+str(start),"--end="+str(stop),args.bam_file])
intermediate_wig="./intermediate/Wiggle/"+os.path.basename(args.bam_file)+"_"+region
call(["./scripts/bigWigToWig",args.bam_file[:-4]+".bigwig",intermediate_wig])
call(["rm",args.bam_file[:-4]+".bigwig"])
OUTPUT=open(args.wig_file,"w")
OUTPUT.write("track type=wiggle_0 name="+os.path.basename(args.bam_file)+"_"+region+"\n")
INT_WIG = open(intermediate_wig,"r")
for line in INT_WIG.readlines():
    OUTPUT.write(line)

OUTPUT.close()
INT_WIG.close()
