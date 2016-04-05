#! /usr/bin/python

# Summarize and annotate RNA-Seq output

file1="NonPregnant1"
file2="NonPregnant2"
file3="NonPregnant3"
file4="NonPregnant4"
annot_file = "GPL6244-24073.txt"
output_file= "NonPregnant_annotated.txt"

FILE1_HANDLE = open(file1,"r")
FILE2_HANDLE = open(file2,"r")
FILE3_HANDLE = open(file3,"r")
FILE4_HANDLE = open(file4,"r")
ANNOTATION = open(annot_file,"r")
OUTPUT = open(output_file,"w")

file1_dict = dict()
for line in FILE1_HANDLE.readlines():
    if line[:1] == "#":
        continue
    elif line[:2] == "ID":
        continue
    info = line.rstrip().split("\t")
    file1_dict[info[0]] = float(info[1])

file2_dict = dict()
for line in FILE2_HANDLE.readlines():
    if line[:1] == "#":
        continue
    elif line[:2] == "ID":
        continue
    info = line.rstrip().split("\t")
    file2_dict[info[0]] = float(info[1])

file3_dict = dict()
for line in FILE3_HANDLE.readlines():
    if line[:1] == "#":
        continue
    elif line[:2] == "ID":
        continue
    info = line.rstrip().split("\t")
    file3_dict[info[0]] = float(info[1])

file4_dict = dict()
for line in FILE4_HANDLE.readlines():
    if line[:1] == "#":
        continue
    elif line[:2] == "ID":
        continue
    info = line.rstrip().split()
    file4_dict[info[0]] = float(info[1])

annotate_dict = dict()
for line in ANNOTATION.readlines():
    if line[:1] == "#":
        continue
    elif line[:2] == "ID":
        continue
    info = line.rstrip().split("\t")
    if len(info) < 10 or info[9] == "---":
        gene_name = "NA"
        gene_nm = "NA"
    elif "//" in info[9]:
        gene_def = info[9].split("//")
        gene_name = gene_def[1]
        gene_nm = gene_def[0]
    region = info[2]
    annotate_dict[info[0]] = gene_name+"\t"+gene_nm+"\t"+region

OUTPUT.write("ID\tMean\tNonPregnant1\tNonPregnant2\tNonPregnant3\tNonPregnant4\tGene\tRefSeq\tRegion\n")
for key in file1_dict.keys():
    cum = file1_dict[key]+file2_dict[key]+file3_dict[key]+file4_dict[key]
    mean = cum / 4.
    OUTPUT.write(key+"\t"+str(mean)+"\t"+str(file1_dict[key])+"\t"+str(file2_dict[key])+"\t"+str(file3_dict[key])+"\t"+str(file4_dict[key])+"\t"+annotate_dict[key]+"\n")

FILE1_HANDLE.close()  
FILE2_HANDLE.close()  
FILE3_HANDLE.close()  
FILE4_HANDLE.close()  
ANNOTATION.close()
OUTPUT.close()               
