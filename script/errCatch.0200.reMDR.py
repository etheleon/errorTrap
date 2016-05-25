#!/usr/bin/env python

from collections import deque
from Bio import SeqIO
import re
import argparse
import numpy as np
import pandas as pd
import pprint as pp
import sys

parser = argparse.ArgumentParser(description='Check if the MDR sequences are similar')
parser.add_argument('root',
                    metavar='/path/to/mdr/fasta/files',
                    type=str,
                    nargs=1,
                    help='Path to the MDR only fasta files',
                    default="/export2/home/uesu/simulation_fr_the_beginning/howpAssNew"
                    )
parser.add_argument('ko',
                    metavar='KO',
                    type=str,
                    nargs=1,
                    help='the ko of interest',
                    default="K00001"
                    )
args = parser.parse_args()

root = args.root[0]
ko      = args.ko[0]

class Alignment:
    '''
    KO Specific
    Stores output details from the PADI pipeline for error checking:
    * actual contig sequence
    * MDR contig sequence
    * the start and end locations
    * the reads which accrude into these contig bins
    '''

    start       = 0
    end         = 0
    contigList  = {}
    readInfo    = {}

    def __init__ (self, rootPath, ko):
        self.rootPath = rootPath
        self.ko = ko
        self.getMSALOC()
        self.readContigs()
        self.readMSA()
        self.readStatus()
        self.parseFastQ()

    def getMSALOC(self):
        """
        record the MDR region (not exactly the same)
        """
        counter = 0
        file = self.rootPath + '/out/mdr/' + self.ko + ".fna"
        for record in SeqIO.parse(file, "fasta"):
            if counter > 0:
                break
            theMatch = re.search("msaStart:(\d+) msaEND:(\d+)", record.description)
            counter += 1
            Alignment.start = int(theMatch.group(1))
            Alignment.end = int(theMatch.group(2))
            print "start %s, end %s"%(Alignment.start, Alignment.end)


    def readContigs(self):
        """
        Stores full length contigs
        """
        path = self.rootPath+'/out/newbler/'+ko+"/454AllContigs.fna"
        for record in SeqIO.parse(path, 'fasta'):
            Alignment.contigList[record.id] = {'fullseq': record.seq.upper()}

    def readMSA(self):
        """
        Parses for MDR region based on coordinates outputs the correct strand
        """
        path = self.rootPath + '/out/ntMSA/' + self.ko + ".msa"
        for record in SeqIO.parse(path, "fasta"):
            contig = Alignment.contigList[record.id]['fullseq']
            recseq = str(record.seq)[self.start:self.end]

            shrunk = recseq.replace('-', '').upper()
            substrResult = contig.find(shrunk)
            #its a reverse strand
            if substrResult == -1:
                revshrunk = str(record.seq.reverse_complement())[self.start:self.end].replace('-', '').upper()
                substrResult = contig.find(revshrunk)
                #print ">", record.id,  len(contig), "type: full\tDirection: rev", "\tsubstrResult", substrResult
                #print contig
                backpadding = len(contig) - substrResult - len(revshrunk)
                front = "-" * substrResult
                back = "-" * backpadding
                #print ">", record.id,  len(contig),"type: ntMSA\tDirection: rev", "\tsubstrResult", substrResult
                #print "%s%s%s"%(front,revshrunk,back)
                Alignment.contigList[record.id]['mdr'] = { 'start': substrResult,
                                                            'end' : len(revshrunk) + substrResult,
                                                            'seq' : revshrunk,
                                                            'direction' : 'ntRev'
                }
            #not a reverse strand
            else:
                #print ">", record.id,len(contig),  "type: full\tDirection: ntRev", "\tsubstrResult", substrResult
                #print contig
                backpadding = len(contig) - substrResult - len(shrunk)
                front = "-" * substrResult
                back = "-" * backpadding
                #print ">", record.id,len(contig),  "type: ntMSA\tDirection: ntRev", "\tsubstrResult", substrResult
                #print "%s%s%s"%(front,shrunk,back)
                Alignment.contigList[record.id]['mdr'] = {
                    'start': substrResult,
                    'end' : len(shrunk) + substrResult,
                    'seq' : shrunk,
                    'direction' : 'ntRev'
                }
    def readStatus(self):
        df = pd.read_csv(self.rootPath + "/out/newbler/"+ko+"/454ReadStatus.txt", sep="\t").head(5)
        df.columns=['Accno','ReadStatus','5Contig','5Position','5Strand','3Contig','3Position','3Strand']

        """
        Accno   Read Status     5' Contig       5' Position     5' Strand       3' Contig       3' Position     3' Strand
        simuREAD_62|taxID|191767|loc|7076959-7077060|outpu      Assembled       contig00200     258     -       contig00200     157     +
        simuREAD_332|taxID|201096|loc|478844-478945|output      Assembled       contig02440     104     -       contig02440     8       +
        simuREAD_883|taxID|18|loc|841131-841232|output|s_1      Assembled       contig00107     211     +       contig00107     312     -
        simuREAD_2334|taxID|561|loc|4092117-4092218|output      Assembled       contig00767     304     +       contig00767     404     -
        """
        def splitNstore(row):
            readID = row['Accno'].split("|")[0]
            #print "readID: " + readID
            isAss = row['ReadStatus'] == 'Assembled'
            if isAss:
                isSame = row['5Contig'] == row['3Contig']
                isPos = row['5Position'] == '+'
                if isSame :
                    Alignment.readInfo[readID] = {
                        'parent': row['5Contig'],
                        'startPos':  int(row['5Position']) if isPos else int(row['3Position']),
                        'endPos'      : int(row['3Position']) if isPos else int(row['5Position']),
                        'direction':  'forward' if isPos else 'reverse'
                    }

        df.apply(splitNstore, axis=1)
        #pp.pprint(Alignment.readInfo['simuREAD_62'])

    def parseFastQ(self):
        """
        Parses the fastqfiles and stores the reads into the correct contigs
        """
        for contig in Alignment.contigList:
            Alignment.contigList[contig]['reads'] = deque()

        fq1 = self.rootPath + "/out/newbler/" + ko + "/input/" + ko + ".1.fq"
        fq2 = self.rootPath + "/out/newbler/" + ko + "/input" + ko + ".2.fq"
        for record in SeqIO.parse(fq1, "fastq"):
            readID = record.description.split("|")[0]
            if readID in Alignment.readInfo:
                theParent = Alignment.readInfo[readID]['parent']
            pp.pprint(Alignment.contigList[theParent]['reads'])#.append(str(record))
            Alignment.contigList[theParent]['reads'].append(str(record.seq))
        pp.pprint(Alignment.contigList['contig00001'])

aln = Alignment(root, ko)


