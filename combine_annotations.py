# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 14:30:32 2021

@author: iss1g18
"""
#import sys


diz = {}

with open("extra_CADDV1.6_hg38_online.tsv",'r') as infile: #input file is your extra annotation
    for line in infile:
        line = line.rstrip().split("\t") #format the line (remove blanks and split by tab)
        query = "\t".join(line[0:-4]) #create your query (e.g. chr-position)
        diz[query] = line[-2] #assign cadd score [raw, phred]

             
with open("master_rawscore.txt",'r') as infile:
    for line in infile:
        line=line.rstrip().split("\t")
        if line[2]=="99999999":
            query="\t".join(line[0:-1])
            try:
                line[2]=str(float(diz[query])) #raw
            except:
                line[2]=str("noanno")
                  
	print("\t".join(line))
