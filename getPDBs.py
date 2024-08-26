# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 11:21:27 2024

@author: Tigist Tamir, on Machine Name = Luisa
"""

import urllib.request
import re
from multiprocessing import Pool
import glob
from plip.structure.preparation import PDBComplex
from plip.exchange.report import StructureReport
import os.path


f=open('2024_new_PDBs.txt','r')
pdbid=f.read().split(',')


uless=[]	
def getmypdb(pdbid):
	if len(pdbid)>0:
		try:
			urllib.request.urlretrieve('http://files.rcsb.org/download/'+pdbid+'.pdb', pdbid+'.pdb')
		except:
			uless.append(pdbid)

newlist=[]
for i in pdbid:
    if not os.path.exists ('./PDBs/'+i+'.pdb'):
        newlist.append(i)
        getmypdb(i)

for i in uless:
	pdbid.remove(i)

with open("noPDB.txt", "w") as outfile:
    outfile.write("\n".join(str(item) for item in uless))

with open("notdone.txt", "w") as outfile:
    outfile.write("\n".join(str(item) for item in newlist))