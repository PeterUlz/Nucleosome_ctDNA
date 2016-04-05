#! /usr/bin/python

# make a text file from wig starting from start and ending at end
#  replace missing entries with 0 entries for further analysis in R
#  remove header and variableStep entries
import argparse
import sys

parser = argparse.ArgumentParser(description='Analyze read depth in comparison to transcription start')
parser.add_argument('-w','--wig', dest='wig_file',
                   help='WIG file',required=True)
parser.add_argument('-o','--output-file', dest='output_file',
                   help='WIG file',required=True)
parser.add_argument('-s','--start', dest='start',
                   help='Start coordinates',required=True,type=int)
parser.add_argument('-e','--end', dest='end',
                   help='End coordinates',required=True,type=int)

args = parser.parse_args()


INFILE = open(args.wig_file,"r")
OUTFILE = open(args.output_file,"w")

content = INFILE.readlines()
last_position = args.start
last_write_position = args.start

for line in content:
    if line[:1] == "t":
        continue
    elif line[:1] == "v":
        continue
    else:
        info = line.rstrip().split("\t")
        position = int(info[0])
        count = info[1]
        if position >= args.start and position <= args.end:
            if (position+1) != last_position:
                for i in range(last_position+1,position):
                    OUTFILE.write(str(i)+"\t0\n")
                    last_write_position = i
            OUTFILE.write(str(position)+"\t"+count+"\n")
            last_write_position = position
            last_position = position

if last_write_position != args.end:
    for i in range(last_write_position+1,args.end):
        OUTFILE.write(str(i)+"\t0\n")

INFILE.close()
OUTFILE.close()




