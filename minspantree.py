#
# Currently, this code should be placed in the directory that holds
# the pasta output. This code takes in the directory that contains 
# subdirectories that contain MAFFT alignments (temporary Pasta files). 
# It does a pairwise OPAL alignment between all MAFFT subsets. 
# We then find the minimum spanning tree on this undirected, weighted 
# edge clique where the edge weight corresponds to the quality of the 
# opal alignment. It keeps the Opal alignments on the minimum spanning 
# tree and does transitivity to merge the remaining alignments. 
#

import sys, re
import pasta.alignment as all
import glob, operator
import numpy as np
from multiprocessing import Pool
import subprocess
import os, json, sys
import itertools
from random import random
import sys
from dendropy.dataobject.taxon import Taxon
from copy import deepcopy

#wheres_opal = '/u/sciteam/collins2/software/sate-tools-linux/opal.jar'
wheres_opal = '/home/kcollins/Desktop/sate-tools-linux/opal.jar'
nproc = 16
parent = {}
rank = {}

def MST(opal_dict):
    sort_keys = opal_dict.keys()
    sort_keys.sort()
    mst = []
    for k in sort_keys:
	for tups in opal_dict[k]:
	    if find(tups[1]) != find(tups[2]):
		union(tups[1],tups[2])		
		mst.append(tups)
    return mst
    
def makeSet(v):
    parent[v] = v
    rank[v] = 0

def find(v):
    if v not in parent:
	makeSet(v)
    elif parent[v] != v:
    	parent[v] = find(parent[v])
    return parent[v]

def union(v1, v2):
    root1 = find(v1)
    root2 = find(v2)
    if root1 != root2:
    	if rank[root1] > rank[root2]:
	    parent[root2] = root1
    	else:
	    parent[root1] = root2
    	if rank[root1] == rank[root2]:
		rank[root2] = rank[root2] + 1


def score(out):
    seqs = all.Alignment()
    seqs.read_file_object(out)
    gap = list()
    for i in seqs:
	theSeq = list(seqs[i])
	length = 1.0*len(theSeq)
	x = (theSeq.count("-")/length)
	gap.append(x)
    gap.sort()
    if len(gap)%2==0:
	return gap[len(gap)/2]
    else:
	return (gap[len(gap)/2] + gap[len(gap)/2+1])/2.0

def opalPairwise(Dict, storeOpal):
    opalDict = dict()
    argslist = []
    #this k1,k2 is a pairwise look at keys to mafft alignments
    for k1, k2 in itertools.combinations(Dict, 2):
	opalName = '/'+str(k1)+'_'+str(k2)+'.txt'
	argslist.append((Dict[k1],Dict[k2], str(storeOpal+opalName)))
    results = Pool(nproc).map(runOpal, argslist)
    for i in results:
	if i[0] in opalDict:
	    opalDict[i[0]].append((i[1],i[2],i[3]))
    	else:
	    opalDict[i[0]] = [(i[1],i[2],i[3])]
    return opalDict

def runOpal(args):
    mafft1 = str(args[0]) 
    mafft2 = str(args[1])
    output = str(args[2])
    args = ['java', '-Xmx1000m', '-jar', wheres_opal, '--in', mafft1, '--in2', mafft2, '--out', output, '--quiet', '--align_method', 'profile']
    subprocess.call(args)
    x = score(output)
    return [x,output,mafft1,mafft2]


def getMafftAlignment(directory):
    #this dict keeps all mafft alignments separate
    dictionary = dict()
    keyName = 1
    #this only works for pasta and how it stores the data, it is not generalized
    for loc in glob.glob(directory+'/temp*/step2/centroid/r*/d*/*/input.aligned'):	
	dictionary[keyName] = loc
	keyName = keyName + 1
    return dictionary

def trans(minSpanTree):
    print ("We're starting")
    count = 0
    merged = {}
    merged_groups = {}
    alignment = {}
    #I am popping offer the tuple (OPAL, MAFFT1, MAFFT2) as I explore them
    print ('here')
    #I should write this as its own function so I don't have to break a bunch, this is sloppy but works TODO
    while len(minSpanTree) > 0:
	print ('number of alignments: ' + str(len(alignment.keys())))
	#this is to check if I need to make a new alignment ie I have non overlapping OPAL alignments
	print ('Left: ' + str(len(minSpanTree)))
	findHome = False
	#MST is in order, so this is good
	x = minSpanTree[0]
	#now that I have it, I can pop it
        minSpanTree.pop(0)
	#this looks over the various mafft alignments within each alignment made by transitivity
	#I should write this as its own function so I don't have to break a bunch, this is sloppy but works TODO
	for i in merged_groups:
	    #is mafft1 or mafft2 in this specific alignment made by transitivity
	    if x[1] in merged_groups[i] or x[2] in merged_groups[i]:
		#yes? okay, which one? I only need to add the other to merged_groups
		if x[1] in merged_groups[i]:
		    merged_groups[i].append(x[2])
		else:
		    merged_groups[i].append(x[1])
		#now that I have the mafft alignments. I now need to merge_in its opal alignment
		print ('Merging in Opal alignment')
		seq = all.CompactAlignment()
		seq.read_file_object(x[0])
		alignment[i].merge_in(seq)
		#now I add its opal alignment to this list
		merged[i].append(x[0])
		#I've added it (TRUE) and then break because I've added the best opal alignment so far
		findHome = True
		break
	#if there wasn't overlap I need to start a new alignment
	#I should write this as its own function so I don't have to break a bunch, this is sloppy but works TODO
	if findHome is False:
	    print ("Adding new alignment, couldn't merge in")
	    #add its opal alignment
	    merged[count] = [x[0]]
	    #add mafft1
	    merged_groups[count] = [x[1]]
	    #add mafft2
	    merged_groups[count].append(x[2])
	    #add the seq itself to the alignment dictionary
	    seq = all.CompactAlignment()
	    seq.read_file_object(x[0])
	    alignment[count] = seq
	    #now iterate
            count = count + 1
	#now that we've done that, we need to make sure we cannot combine the various alignments
	#this is a pairwise check over the sets of alignments made so far
	#I should write this as its own function so I don't have to break a bunch, this is sloppy but works TODO
	doneLooking = False
	for k1, k2 in itertools.combinations(merged_groups, 2):
	    foundIt = False
	    if doneLooking is True:
		break
	    print ("We're now checking to see overlap")
	    #checks the mafft alignments in the first set
	    for i in merged_groups[k1]:
		if foundIt is True:
		    break
		#checks the mafft alignments in the second set
		for j in merged_groups[k2]:
		    print ('i: ' + i)
		    print ('j: ' + j)
		    #is there a match? move the second set into the first because the first is the "better" one
	    	    if i == j:
			foundIt = True
			doneLooking = True
			print ("There's overlap")
			#remove the name of the mafft alignment that overlaps
	    	    	merged_groups[k2].remove(j)
			print ("1")
			#then extend the first list by the second
		        merged_groups[k1].extend(merged_groups[k2])
			print ("2")
			#now extend the list of the opal alignments
		    	merged[k1].extend(merged[k2])
			print ("3")
			#now merge in the second alignment made by transitivity
		    	alignment[k1].merge_in(alignment[k2])
			print ("4")
			#now remove the second alignment, set of opal and mafft alignments because now its in first one
		    	merged.pop(k2)
			print ("5")
		    	merged_groups.pop(k2)
			print ("6")		
		    	alignment.pop(k2)
			print ("7")
			break
    return alignment[0]

		
if __name__ == '__main__':
    #directory of mafft alignments (stored in subdirectories)
    ipdir = os.getcwd() + '/pastajob'
    #opal output directory, where to store opal alignments
    opdir = os.getcwd() + '/opal'
    mafftDict = getMafftAlignment(ipdir)
    opalDict = opalPairwise(mafftDict, opdir)
    #opalDict = {0.5:[('z','b','h'),('1','a','c'),('2','a','d'),('3','b','c'),('4','b','d')],0.7:[('5','a','f'),('6','c','e'),('7','c','f')],0.9:[('8','b','e'),('9','a','g')],0.6:[('10','b','f'),('11','c','g')],1.1:[('12','c','h'),('13','a','h'),('14','b','g')],1.0:[('15','d','e')],0.8:[('16','e','h')],1.2:[('17','d','f'),('18','f','h')],1.3:[('19','d','g'),('20','g','h'),('21','e','g')],0.1:[('t','a','b'),('u','a','e')],0.3:[('w','c','d')],0.4:[('x','f','e'),('y','f','g')],0.2:[('v','h','d')]}
    minSpanTree = MST(opalDict)
    answer = trans(minSpanTree)
    answer.write_filepath(os.getcwd()+'/alignedMST.fasta')
    
