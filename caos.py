from skbio import TreeNode
from io import StringIO
from Bio import Phylo
from Bio import SeqIO
import re
import sys
import numpy as np
import os
from scipy import stats
import random
import math

def write_fasta(file, accession, sequence):
    file.write(">{}\n".format(accession))
    if len(sequence) <= 70:
        file.write("{}\n".format(sequence))
    else:
        line_num_exact = len(sequence) / 70
        line_num = math.floor(line_num_exact)
        #println(line_num)
        for i in range(0,line_num):
            start = i * 70
            stop = i * 70 + 70
            file.write(sequence[start:stop])
            file.write("\n")


        start = line_num * 70
        stop = len(sequence)
        file.write(sequence[start:stop])
        #file.write("\n")
        file.write("\n")

outputDir = sys.argv[2]
################################################################################
with open('{}/{}_aligned.fasta'.format(outputDir,sys.argv[1])) as f:
    alist = [line.rstrip() for line in f]
seqDict = {}

accession = ''
seq = ''
for line in alist:
    if '>' in line:
        if accession != '':
            seqDict[accession] = seq
        accession = line.replace('>', '')
        accession = accession.replace('_', '')
        accession = accession.replace('.', '')
        seq = ''
    else:
        seq = seq + line
seqDict[accession] = seq
#print(len(seqDict))

homoAccession = sys.argv[1]
homoAccession = homoAccession.replace('_', '')
homoAccession = homoAccession.replace('.', '')
homoSeq = seqDict[homoAccession]
seqLength = len(homoSeq)
#print(homoSeq)

#################################################################################
#Phylogenetic tree
################################################################################
Phylo.convert('{}/{}_trees.nex'.format(outputDir, sys.argv[1]), 'nexus', '{}/{}_tree.tree'.format(outputDir, sys.argv[1]), 'newick')

f = open('{}/{}_tree.tree'.format(outputDir, sys.argv[1]), 'r')
tree = f.readline()
#print(tree)

tree = re.sub(r':\d.\d+', '', tree)
#print(tree)
tree = TreeNode.read(StringIO(tree))
'''
for child in tree.children:
    print(child)
    print(child.is_tip())
exit()
'''
#print(tree)

####################################################################
#find set of characters given a list of taxa l and position number i
####################################################################
def findSet(l, i, seqDict):
    ll = []
    for j in l:
        ll.append(str(seqDict[j][i]))
    ll = list(set(ll))
    return ll

def updateVN(node, child, variation_number, seqDict, length):

    allClades = re.findall(r'[a-zA-Z0-9]+', str(node))
    #print("all")
    #print(allClades)
    for c in child:
        cClades = re.findall(r'[a-zA-Z0-9]+', str(c))
        subClades = list(set(allClades) - set(cClades))
        #print("sub")
        #print(subClades)
        for i in range(length):
            cSet = findSet(cClades, i, seqDict)
            subSet = findSet(subClades, i, seqDict)
            for item in cSet:
                if item not in subSet:
                    variation_number[i] = variation_number[i] + 1
    return variation_number

def generateVN(tree, seqDict, seqLength):
    variation_number = np.zeros((seqLength,), dtype=int)
    queue = []
    queue.append(tree)

    while len(queue) > 0:
        node = queue.pop()
        if node.is_tip():
            continue
        else:
            child = node.children
            variation_number = updateVN(node, child, variation_number, seqDict, seqLength)
            #print('variation')
            #print(variation_number)
            for c in child:
                queue.append(c)
    return variation_number

variation_number = generateVN(tree, seqDict, seqLength)

homoIndexList = []
f_vn = []
for i in range(len(homoSeq)):
    if str(homoSeq[i]) != '-':
        homoIndexList.append(i)
        f_vn.append(variation_number[i])
#print(len(homoIndexList))

outputFile = open("{}/vn_{}.txt".format(outputDir,sys.argv[1]), 'w')

vn_max = max(f_vn)
vn_min = min(f_vn)
#print(vn_min, vn_max)
#exit()
for i in range(len(homoIndexList)):
    j = i + 1
    vn = variation_number[homoIndexList[i]]
    vn = (vn - vn_min) / (vn_max - vn_min)
    outputFile.write('{}\t{}\t{}\n'.format(str(j), homoSeq[homoIndexList[i]], str(vn)))

#outputFile.write("\n")
outputFile.close()
exit()
#cleanSeqDict = {}
outputFile = open("{}/evmutation_{}.txt".format(outputDir,sys.argv[1]), "w")
seq = seqDict[homoAccession]
seq1 = ''
for i in homoIndexList:
    seq1 = seq1 + seq[i]
#cleanSeqDict[key] = seq1
write_fasta(outputFile, homoAccession, seq1)

for key in seqDict.keys():
    #outputFile.write('>{}'.format(key))
    #outputFile.write('\n')
    if key == homoAccession:
        continue
    seq = seqDict[key]

    seq1 = ''
    for i in homoIndexList:
        seq1 = seq1 + seq[i]
    #cleanSeqDict[key] = seq1
    write_fasta(outputFile, key, seq1)
    #outputFile.write(seq1)
    #outputFile.write('\n')
outputFile.close()

#with open("../evmutation/evmutation_fasta.txt", "w") as handle:
#    SeqIO.write(cleanSeqDict.values(), handle, 'fasta')
exit()
