#Created by Claudia
#Jan 28 2015
#Modified Oct 2016 for Beagle output instead of GERMLINE
#MRe-modified to make it simple!! And take bed file of pairwise IBD segments...

import numpy as np
import sys
from sys import argv


#Arguments: path and file names
#datafile: .merged.ibd files
datafilename = sys.argv[1]
print (datafilename)
#Output
mychr = int(sys.argv[2])
print (mychr)
outputfilename = datafilename + '.sharing.by.pos'
print (outputfilename)
bim_file = (sys.argv[3])
print (bim_file)


#load Beagle file
def loadfile (datafile):
    SNP1 = np.genfromtxt(fname = datafile,dtype=int,usecols=5)
    SNP2 = np.genfromtxt(fname = datafile,dtype=int,usecols=6)

    return [SNP1,SNP2]


#load vcf or bim file for map
def loadmapfile (mapfile):

    positions = []
    f=open(mapfile,'r')
    for line in f:
        ichr=int(line.rstrip('\n\r').split()[0])
        if ichr == mychr:
            pos=int(line.rstrip('\n\r').split()[3])
            positions.append(pos)

    return positions


#Write file
def writefile (outputfilename,dict,pos):
    outputfile = open(outputfilename, "w")

    for i in range(len(pos)): #Remove the 1 digit chr no at the beggining of pos
        #if dict[pos[i]] > 0:
        outputfile.write(str(mychr) + "\t" + str(pos[i])+ "\t" + str(dict[pos[i]]) + "\n")       
    outputfile.close()

	

#let's go!
print ('loading data file...')
data = loadfile(datafilename)
print ("let's go!!")

positions=loadmapfile(bim_file)
# print (positions)
print (len(positions))
print (positions[1:10])
print ("Positions done...")
pos_dict = dict(zip(positions, [0]*len(positions)))
# print (pos_dict)
print ("Dict done...looping through IBD segments...")
for i in range(len(data[0])):
    for j in range(positions.index(data[0][i]),(positions.index(data[1][i])+1),1):
        pos_dict[positions[j]]+=1
writefile (outputfilename,pos_dict,positions)

        
print ("Finished!!")
