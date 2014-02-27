import itertools
import urllib
import xml.etree.ElementTree as ET
import time
import os
import matplotlib.pyplot as plt
import numpy as np

class ProteinInteractionNetwork(object):
    def __init__(self,SIF_file):
        self.interactions = []
        self.sequencelengths = {}
        with open(SIF_file, 'r') as sourcefile:
            for line in sourcefile:
                self.interactions.append(tuple([p.strip() for p in line.split('pp')]))
            
            self.proteins =  set(list(itertools.chain.from_iterable(self.interactions)))
    def downloadSequenceLengths(self,limit = None, verbose = True,filename = None):
        def showprogress(self):
            elapsed_time = round(time.time() - start_time,2)
            print 'Processed:',processed,'Elapsed time (s):',elapsed_time
            remaining = len(self.proteins) -processed
            print 'Remaining:',remaining
            if processed != 0:
                time_per_entry = elapsed_time / processed
                print 'Estimated time left (min):',round(remaining*time_per_entry/60,1)
        if filename is not None: #If a file (containing previously retrieved sequence lengths) is specified
            if os.path.exists(filename): #If the output file exists
                with open(filename,'r') as sourcefile:
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
                    showprogress(self)
                continue
            try:
                tree = ET.parse(xml)
            except:
                if verbose:
                    print 'Failed to parse xml for protein',protein_id
                    print 'xml URL:',xmlurl
                    showprogress(self)
                continue
            root = tree.getroot()
            try:
                AA_sequence = root.findall('{http://uniprot.org/uniprot}entry/{http://uniprot.org/uniprot}sequence')[0]
            except:
                self.nosequence.append(protein_id)
                continue
            self.sequencelengths[protein_id] = AA_sequence.attrib['length']
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
            
    def findSequenceProportions(self, limit = None, verbose =0):
        proportions = []
        processed = 0
        missing_data = 0
        for interaction in self.interactions:
            if limit is not None:
                if processed >= limit:
                    break
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

if __name__ == "__main__":
    PIN = ProteinInteractionNetwork('Homo sapiens-20121210.sif')
    #PIN.downloadSequenceLengths() #Takes about an hour for 15 000 proteins
    PIN.readSequenceLengths('protlengths.txt')
    prop = PIN.findSequenceProportions(verbose=1)
    percentile = np.percentile(prop,99)
    fig = plt.figure()
    plt.hist(prop,bins = 20,range=(1,percentile),normed = True)
    plt.xlabel('Relative protein size')
    plt.ylabel('Fraction of interacting protein pairs')
    plt.show()
