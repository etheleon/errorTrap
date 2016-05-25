#!/usr/bin/env python

from os import listdir, system
from os.path import splitext
import re
import argparse

parser = argparse.ArgumentParser(description='Check if the MDR sequences are similar')
parser.add_argument('MDRpath', metavar='/path/to/mdr/fasta/files', type=str, nargs=1,
                    help='Path to the MDR only fasta files', default="/export2/home/uesu/simulation_fr_the_beginning/howpAssNew/out/mdr/")
args = parser.parse_args()



#Qn: If the MDR only sequences (within the region) are exactly identity to each other

def runCDHIT(f, mypath):
    name = splitext(f)[0]
    cmd = "cd-hit -i %s/%s -c 1 -o %s/%s.cdhitOut > %s/%s.stdOut"%(mypath,f, mypath,name, mypath,name)
    system(cmd)
    outputFile = '%s/%s.stdOut'%(mypath, name)
    with open(outputFile, 'r') as infile:
        for line in infile:
            theMatch = re.search("(\d+)\s+finished\s+(\d+)\s+clusters", line)
            if theMatch is not None:
                status = "same" if theMatch.group(1) == theMatch.group(2) else "diff"
                print "KO: %s, finished: %s clusters: %s status: %s"%(name, theMatch.group(1), theMatch.group(2), status)


[runCDHIT(f, args.MDRpath[0]) for f in listdir(args.MDRpath[0]) if f.endswith(".fna")]
