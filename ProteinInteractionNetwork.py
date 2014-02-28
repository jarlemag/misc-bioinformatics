import itertools
import urllib
import xml.etree.ElementTree as ET
import time
import os
import matplotlib.pyplot as plt
import numpy as np
import random

class ProteinInteractionNetwork(object):
    def readInteractions(self,filename,fraction = 1, verbose = True,delimiter = 'pp'):
        with open(filename, 'r') as sourcefile:
            if verbose:
                print 'Reading interactions from',filename
            if fraction < 1:
                pass
                #Select a random subset of the network
            else:
                for line in sourcefile:
                    self.interactions.append(tuple([p.strip() for p in line.split(delimiter)]))
           
    def __init__(self,SIF_file,description = None,size_limit=None,delimiter = 'pp'):
        self.interactions = []
        self.sequencelengths = {}
        self.self_interactions = 0
        self.readInteractions(SIF_file,delimiter  = delimiter)
        self.proteins = set(list(itertools.chain.from_iterable(self.interactions)))
        self.description = description
    def showprogress(self,start_time,target,processed):
        elapsed_time = round(time.time() - start_time,2)
        print 'Processed:',processed,'Elapsed time (s):',elapsed_time
        remaining = target - processed
        print 'Remaining:',remaining
        if processed != 0:
            time_per_entry = elapsed_time / processed
            print 'Estimated time left (min):',round(remaining*time_per_entry/60,1)

            
    def downloadSequenceLengths(self,limit = None, verbose = True,file_in = None,file_out = None):
        if file_in is not None: #If a file (containing previously retrieved sequence lengths) is specified
            if os.path.exists(file_in): #If the output file exists
                with open(file_in,'r') as sourcefile:
                    retrieved = [line.strip().split(' ')[1] for line in sourcefile]
                for protein in self.proteins:
                    if protein not in retrieved:
                        proteins.append(protein) #Append the file with only previously unretrieved sequence lengths              
        else:
            proteins = self.proteins
        self.failedretrieval = []
        self.nosequence = []
        baseurl = "http://www.uniprot.org/uniprot/"
        processed = 0
        start_time = time.time()
        for protein_id in proteins:
            if limit is not None:
                if processed >= limit:
                    break
            xmlurl = baseurl + protein_id + '.xml'
            try:
                xml = urllib.urlopen(xmlurl)
            except:
                self.failedretrieval.append(protein_id)
                if verbose:
                    print 'Failed to access xml for protein',protein_id,'at',xmlurl
                    self.showprogress(start_time,len(self.proteins),processed)
                continue
            try:
                tree = ET.parse(xml)
            except:
                if verbose:
                    print 'Failed to parse xml for protein',protein_id
                    print 'xml URL:',xmlurl
                    self.showprogress(start_time,len(self.proteins),processed)
                continue
            root = tree.getroot()
            try:
                AA_sequence = root.findall('{http://uniprot.org/uniprot}entry/{http://uniprot.org/uniprot}sequence')[0]
            except:
                self.nosequence.append(protein_id)
                continue
            length = AA_sequence.attrib['length']
            self.sequencelengths[protein_id] = int(length)
            if file_out is not None:
                with open(file_out,'a') as target:
                    target.write(protein_id)
                    target.write(' ')
                    target.write(length)
                    target.write('\n')
            processed +=1
    def readSequenceLengths(self,filename,replace = True):
        if replace:
            self.sequencelengths = {}
        with open(filename, 'r') as sourcefile:
            for line in sourcefile:
                key = line.split(' ')[0]
                value= line.split(' ')[1]
                self.sequencelengths[key] = int(value)

    def writeSequenceLengths(self,filename,lengths = None,append = False):
        if lengths is None:
            lengths = self.sequencelengths
        if append:
            target = open(filename, 'a')
        else:
            target = open(filename, 'w')
        for key in self.sequencelengths:
            target.write(key)
            target.write(' ')
            target.write(self.sequencelengths[key])
            target.write('\n')
        target.close()
        return
            
    def findSequenceProportions(self, limit = None, verbose =0,exclude_self_interactions = True):
        proportions = []
        processed = 0
        missing_data = 0
        for interaction in self.interactions:
            if limit is not None:
                if processed >= limit:
                    break
            if exclude_self_interactions:
                if interaction[0] == interaction[1]:
                    if verbose >1:
                        print 'Interaction',interaction,'excluded (self-interaction)'
                    continue
            try:
                lengths = [self.sequencelengths[protein_id] for protein_id in interaction]
            except:
                if verbose > 1:
                    print 'Missing sequence length data for protein(s) in interaction',interaction
                missing_data +=1
                continue
            max_length = max(lengths)
            min_length = min(lengths)
            if verbose > 2:
                print 'Interaction:',interaction
                print 'Longest AA sequence:',max_length
                print 'Shortest AA sequence:',min_length
            proportions.append(round(max_length/float(min_length),3))
            processed +=1
        if verbose >0:
            print 'Sequence length data missing for',missing_data,'interactions','('+str(round((missing_data/float(len(self.interactions))*100),2)),'percent)'
        return proportions
    def countSelfInteractions(self):
        if len(self.interactions) > 0:
            for interaction in self.interactions:
                if interaction[0] == interaction[1]:
                    self.self_interactions +=1

    def randomizeInteractions(self,reject_self_interactions = False,fractionalsize = 1,verbose = True):
        if verbose:
            print 'Randomizing interactions.'
        n = len(self.interactions)
        self.interactions = []
        interaction_number = int(n*fractionalsize)
        start_time = time.time()
        if n > 0:
            processed = 0
            while processed < interaction_number:
                protein_a = random.choice(list(self.proteins))
                protein_b = random.choice(list(self.proteins))
                self.interactions.append((protein_a,protein_b))
                processed +=1
                if verbose:
                    if processed % 1000 == 0:
                      self.showprogress(start_time,interaction_number,processed)

    def randomizeSequenceLengths(self):
        n = len(self.sequencelengths)
        centre = np.median(self.sequencelengths.values())
        randomvalues = list(np.random.normal(centre,size=n))
        randomvalues = [int(value) for value in randomvalues]
        for key in self.sequencelengths:
            self.sequencelengths[key] = randomvalues.pop(random.randrange(len(randomvalues)))
    def writeInteractions(self,filename):
        with open(filename,'w') as target:
            for interaction in self.interactions:
                target.write(interaction[0])
                target.write(' pp ')
                target.write(interaction[1])
                target.write('\n')
        
    def sequenceProportionHistogram(self,percentile_cutoff = 99,normalize = True):
        proportions = self.findSequenceProportions(verbose = 1)
        percentile = np.percentile(proportions,percentile_cutoff)
        fig = plt.figure()
        plt.hist(proportions,bins = 20,range=(1,percentile),normed = normalize)
        plt.xlabel('Relative protein chain length')
        plt.ylabel('Fraction of interacting protein pairs')
        plt.title('Relative chain lengths among interacting protein pairs')
        plt.show()

    def sequenceLengthsHistogram(self,percentile_cutoff = 99,normalize = True,bins=20):
        fig = plt.figure()
        lengths = self.sequencelengths.values()
        percentile = np.percentile(lengths,percentile_cutoff)
        plt.hist(lengths,bins = bins,range=(1,percentile),normed=normalize)
        plt.xlabel('Protein chain length')
        plt.ylabel('Fraction of proteins')
        plt.title('Protein chain length distribution')
        plt.show()
    
if __name__ == "__main__":
    PIN = ProteinInteractionNetwork('Homo sapiens-20121210.sif')
    #PIN.downloadSequenceLengths() #Takes about an hour for 15 000 proteins
    PIN.readSequenceLengths('protlengths.txt')
    PIN.sequenceProportionHistogram()
    #PIN.randomizeInteractions()
    #PIN.sequenceProportionHistogram()
    #yeastPIN = ProteinInteractionNetwork('Saccharomyces cerevisiae-20121210.sif')
    #yeastPIN.downloadSequenceLengths(file_out = 'yeastprotlengths.txt') #Takes about 20 minutes for 6000 proteins
    #yeastPIN.readSequenceLengths('yeastprotlengths.txt')
    #elegansPIN = ProteinInteractionNetwork('Caenorhabditis elegans-20121210.sif')
    #elegansPIN.downloadSequenceLengths(file_out = 'elegansprotlengths.txt')
    #elegansPIN.readSequenceLengths('elegansprotlengths.txt')
    Wu2010PIN = ProteinInteractionNetwork('FIs_043009.txt','Human functional protein interaction network (Wu 2010)',delimiter=' ')
    #Wu2010PIN.downloadSequenceLengths(file_out = 'Wu2010protlengths.txt')
    Wu2010PIN.readSequenceLengths('Wu2010protlengths.txt')
    Wu2010PIN.sequenceLengthsHistogram()
    Wu2010PIN.sequenceProportionHistogram()
