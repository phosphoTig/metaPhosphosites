# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 14:53:29 2023

@author: TYT on Luisa
"""

import requests
import re
import glob
import xml.etree.ElementTree as ET
import pandas as pd
import os.path


f=open('NewPDB_IDs.txt','r')
pdb_ids=f.read().split(',')


def search_pdbe_pisa(pdb_ids):
    url ='https://www.ebi.ac.uk/pdbe/pisa/cgi-bin/interfaces.pisa?'+pdb_ids
    response = requests.get(url, stream=True)
    with open('pisafiles/'+pdb_id+'_PISA.xml', 'wb') as f:
        f.write(response.content)
    
        
for pdb_id in pdb_ids:
    if not os.path.exists ('./pisafiles/'+pdb_id+'_PISA.xml'):
        search_pdbe_pisa(pdb_id.lower())
    pdb_ids.remove(pdb_id)


