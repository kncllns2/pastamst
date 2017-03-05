#
# This code takes in the directory that contains subdirectories that 
# contain MAFFT alignments (temporary Pasta files). It does a pairwise
# OPAL alignment between all MAFFT subsets. We then find the minimum 
# spanning tree on this undirected, weighted edge clique where the 
# edge weight corresponds to the quality of the opal alignment.
# It keeps the Opal alignments on the minimum spanning tree and
# does transitivity to merge the remaining alignments. 
#

import sys
from sepp.alignment import MutableAlignment, ExtendedAlignment, _write_fasta
import glob
import numpy as np
from multiprocessing import Pool
import subprocess
import os, json, sys
import itertools

wheres_opal = '/u/sciteam/collins2/software/sate-tools-linux/opal.jar'
nproc = 16

#this was taken from PASTA
def merge_in(self, she):
        '''
        Merges she inside self, assuming we share some common taxa, and the 
        alignment of common taxa is identical across both alignments.
        
        When assumptions are not met, behavior is largely undefined. 
        '''
        #global _T_ID
        #_T_ID += 1
        #ID = _T_ID
        #TIMING_LOG.info("transitivitymerge (%d) started" %ID )
        mykeys = set(self.keys())
        herkeys = set(she.keys())
        #_LOG.debug("Transitive Merge Started. ID:%d - Rows: %d,%d" %(ID,len(mykeys),len(herkeys)))
        shared = mykeys.intersection(herkeys)
        #_LOG.debug("Shared seq: %d" %(len(shared)))
        onlyhers = herkeys - shared
        me_ins = self.get_insertion_columns(shared)
        she_ins = she.get_insertion_columns(shared)
        #_LOG.debug("Insertion Columns: %d,%d" %(len(me_ins),len(she_ins)))

        memap=[]
        shemap=[]

        melen = self.colcount
        shelen =  she.colcount

	me = 0
        ishe = 0
        inew = 0
        while ime < melen or ishe < shelen:
            #print ime,ishe
            if ime in me_ins:
                memap.append(inew)
                ime += 1
            elif ishe in she_ins:
                shemap.append(inew)
                ishe += 1
            else:
                memap.append(inew)
                shemap.append(inew)
                ishe += 1
                ime += 1
            inew += 1

        self.colcount = inew

        for seq in self.itervalues():
            seq.pos = [memap[p] for p in seq.pos]

        for k in onlyhers:
            she[k].pos = [shemap[p] for p in she[k].pos]
            self[k] = she[k]

        #TIMING_LOG.info("transitivitymerge (%d) finished" %ID )
        #_LOG.debug("Transitive Merge Finished. ID:%d; cols after: %d" %(ID,self.colcount))

def MST(opal_dict):
    #assuming min is right
    sort_dict = sorted(opal_dict)
    parent = {}
    rank = {}
    mst = []
    for k in sort_dict:
	for tups in sort_dict[k]:
	    if find(tups[0]) != find(tups[1]):
		union(tups[0],tups[1])		
		mst.append(tups)
     return mst
    
def makeSet(v):
    parent[v] = v
    rank[v] = 0

def find(v):
    if parent.get(v) is None:
	return v
    elif parent[v] is not v:
	parent[v] = find(parent[v])
    return parent[v]

def union(v1, v2):
    if rank[root1] > rank[root2]:
	parent[root2] = root1
    else:
	parent[root1] = root2
    if rank[root1] == rank[root2]:
	rank[root2] += 1

def score(seqs):
    gap = list()
    for i in seqs:
	length = len(seqs[i])
	g = seqs[i].count("-")
	gap.append(g*1.0/length)
    #what should we return? mean/median/average non-gap length/average length
    gap=sort(gap)
    if len(gap)%2==0:
	return gap[len(gap)/2]
    else:
	return (gap[len(gap)/2] + gap[len(gap)/2+1])/2

#MIKE???!!!! I'll need a little help
def opalPairwise(Dict, store_opal):
    mp = Pool(nproc)
    opalDict = dict()
    keyName = keyName + 1
    #NOT SURE HOW MP.MAP works... going to assume it's like & on the command line
    #this k1,k2 is a pairwise look at keys which reference mafft alignments
    for k1, k2 in itertools.combinations(Dict, 2):
	args = ['java', '-Xmx1000m', '-jar', wheres_opal, '--in', Dict[k1], '--in2', Dict[k2], '--out', store_opal] 
	#'--quiet' '--align_method', 'profile'] <-why last 2
    	mp.map(subprocess.call(args))
	seq = MutuableAlignment()
	seq.read_file_object(store_opal)
	key = score(seq)
	#create a dictionary of lists of tuples where the key is the opal alignment score
	#tuple item 0 & 1 are the mafft alignment keys & item 2 is the path to the opal alignment
	#tuples are stored in list in case more than 1 edge has the same score
	t = (k1, k2, opal_path)
	if key in opalDict:
	    opalDict.append(t)
	else:
	opalDict[key] = [t]
    return opalDict

def getMafftAlignment(self, directory):
    #this dict keeps all mafft alignments separate
    dictionary = dict()
    #this dict puts all sequences in mafft alignment into one dictionary
    #this only serves the purpose of making sure we're getting different mafft alignments/sequence names
    #can be removed when certain because it's fairly costly
    bigDictionary = dict()
    keyName = 1
    #this only works for pasta and how it stores the data, it is not generalized
    for file in glob.glob(directory+'/temp*/step2/centroid/r*/d*'):
	seq = MutableAlignment()
	seq.read_file_object(file)
	dictionary[keyName] = seq
	keyName = keyName + 1
	for i in seq:
	    if seq[i] in bigDictionary:
		print ('Error, this sequence already exists or at least this sequence name already exists. Did you repeat mafft alignments?')
	    else:
		bigDictionary[i] = seq[i]
    return dictionary


if __name__ == '__main__':
    #directory of mafft alignments (stored in subdirectories)
    ipdir = str(sys.argv[1])
    #opal output directory, where to store opal alignments
    opdir = str(sys.argv[2])
    mafftDict = getMafftAlignment(ipdir)
    opalDict = opalPairwise(mafftDict, opdir)
    minSpanTree = MST(opalDict)

