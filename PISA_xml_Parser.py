import shutil, os, pandas as pd
import os.path
import xml.etree.ElementTree as ET
from lxml import etree
from multiprocessing import Pool, cpu_count
# import codecs
    
def hbond(root):
    # Extract the hydrogen_bond elements
    hbond_elements = root.findall('.//h-bonds')

    # Create a list to store the data
    hbond_data = []
    
    # Iterate over the hydrogen_bond elements and extract the data
    for sub_hbond in hbond_elements:
        if int(sub_hbond.find('n_bonds').text) != 0:
            all_hbond = sub_hbond.findall("bond")
            for hbond in all_hbond:
                data = {
                    "chain1": hbond.find('chain-1').text,
                    "res1": hbond.find('res-1').text,
                    "seqnum1": hbond.find('seqnum-1').text,
                    "atname1": hbond.find('atname-1').text,
                    "chain2": hbond.find('chain-2').text,
                    "res2": hbond.find('res-2').text,
                    "seqnum2": hbond.find('seqnum-2').text,
                    "atname2": hbond.find('atname-2').text,
                    "dist": hbond.find('dist').text
                } 
                hbond_data.append(data)
    
    # Create the DataFrame
    hbondsDf = pd.DataFrame(hbond_data)
    
    return hbondsDf

def saltBridge(root):
    # Extract the hydrogen_bond elements
    bond_elements = root.findall('.//salt-bridges') 

    # Create a list to store the data
    bond_data = []

    # Iterate over the bond elements and extract the data
    for sub_bond in bond_elements:
        if int(sub_bond.find('n_bonds').text) != 0:
            all_bond = sub_bond.findall("bond")
            for bond in all_bond:
                data = {
                    "chain1": bond.find('chain-1').text,
                    "res1": bond.find('res-1').text,
                    "seqnum1": bond.find('seqnum-1').text,
                    "atname1": bond.find('atname-1').text,
                    "chain2": bond.find('chain-2').text,
                    "res2": bond.find('res-2').text,
                    "seqnum2": bond.find('seqnum-2').text,
                    "atname2": bond.find('atname-2').text,
                    "dist": bond.find('dist').text
                }
                bond_data.append(data)

    # Create the DataFrame
    SBbondsDf = pd.DataFrame(bond_data)

    return SBbondsDf

def ssbond(root):
    # Extract the hydrogen_bond elements
    bond_elements = root.findall('.//ss-bonds')

    # Create a list to store the data
    bond_data = []

    # Iterate over the bond elements and extract the data
    for sub_bond in bond_elements:
        if int(sub_bond.find('n_bonds').text) != 0:
            all_bond = sub_bond.findall("bond")
            for bond in all_bond:
                data = {
                    "chain1": bond.find('chain-1').text,
                    "res1": bond.find('res-1').text,
                    "seqnum1": bond.find('seqnum-1').text,
                    "atname1": bond.find('atname-1').text,
                    "chain2": bond.find('chain-2').text,
                    "res2": bond.find('res-2').text,
                    "seqnum2": bond.find('seqnum-2').text,
                    "atname2": bond.find('atname-2').text,
                    "dist": bond.find('dist').text
                }
                bond_data.append(data)

    # Create the DataFrame
    ssbondsDf = pd.DataFrame(bond_data)

    return ssbondsDf

def covalent(root):
    # Extract the hydrogen_bond elements
    bond_elements = root.findall('.//cov-bonds') # root.findall('.//hydrogen_bond')

    # Create a list to store the data
    bond_data = []

    # Iterate over the hydrogen_bond elements and extract the data
    for sub_bond in bond_elements:
        if int(sub_bond.find('n_bonds').text) != 0:
            all_bond = sub_bond.findall("bond")
            for bond in all_bond:
                data = {
                    "chain1": bond.find('chain-1').text,
                    "res1": bond.find('res-1').text,
                    "seqnum1": bond.find('seqnum-1').text,
                    "atname1": bond.find('atname-1').text,
                    "chain2": bond.find('chain-2').text,
                    "res2": bond.find('res-2').text,
                    "seqnum2": bond.find('seqnum-2').text,
                    "atname2": bond.find('atname-2').text,
                    "dist": bond.find('dist').text
                }
                bond_data.append(data)

    # Create the DataFrame
    cbondsDf = pd.DataFrame(bond_data)

    return cbondsDf

def residues(root):
    # Extract the hydrogen_bond elements
    chains = root.findall('.//molecule')
    
    # Create a list to store the data
    residue_data = []
    
    # Iterate over the hydrogen_bond elements and extract the data
    for chain in chains:
        all_residues = chain.findall("residues/residue")
        
        for residues in all_residues:
            if float(residues.find('bsa').text) != 0.0:
                data = {
                    "molecule_id": chain.find('id').text,
                    "chain_id": chain.find('chain_id').text,
                    "class": chain.find('class').text,
                    "tx": chain.find('tx').text,
                    "ty": chain.find('ty').text,
                    "tz": chain.find('tz').text,
                    "int_natoms": chain.find("int_natoms").text,
                    "int_nres": chain.find('int_nres').text,
                    "int_area": chain.find('int_area').text,
                    "int_solv_en": chain.find('int_solv_en').text,
                    "pvalue": chain.find('pvalue').text,
                    "ser_no": residues.find('ser_no').text,
                    "name": residues.find('name').text,
                    "seq_num": residues.find('seq_num').text,
                    "ins_code": residues.find('ins_code').text,
                    "bonds": residues.find('bonds').text,
                    "asa": residues.find('asa').text,
                    "bsa": residues.find('bsa').text,
                    "solv_en": residues.find('solv_en').text}
                residue_data.append(data)

    # Create the DataFrame
    residueDf = pd.DataFrame(residue_data)
    return residueDf


#locate xml files
directory2 = './pisafiles' 

files = []

with open('incompletePISA.txt','r') as f:
    PISA={line.strip() for line in f}

for file in PISA:  
    filename = os.fsdecode(file)
    with open(f"{directory2}/{filename}", 'r', encoding='utf-8', errors='ignore') as f:
        files.append(f"{directory2}/{filename}")
    
    ccat = []
    failedPDB=[]
    # temp={}
    
    for name in files:
        
        resData = []
        # read xml file
        with open(name, 'r', encoding='utf-8', errors='ignore') as f: # with open(name) as f:
            contents = f.read()
            try:
                # Parse the XML data
                xml_data = ET.fromstring(contents)
        
                interface = xml_data.findall('.//interface')
                PDBid=xml_data.find('.//pdb_code').text
        
                for root in interface:
                    # get the interface details
                    interface_id = root.find('id').text
                    interface_type = root.find('type').text
                    nocc = root.find('n_occ').text
                    pval = root.find('pvalue').text
                    area = root.find('int_area').text
                    solv = root.find('int_solv_en').text
                    stabe_en = root.find('stab_en').text
                    css = root.find('css').text
                    overlap = root.find('overlap').text
                    x_rel = root.find('x-rel').text
                    fixed = root.find('fixed').text
                    
        
                    # run functions that extract interface bonds/residues, and add to dictionary by file name + interaction type
                    ss_df = ssbond(root)
                    h_df = hbond(root)
                    co_df = covalent(root)
                    sb_df = saltBridge(root)
                    res_df = residues(root)
        
                    # adding interface properties to the dfs
                    ss_df['PDB'] = PDBid
                    ss_df['bond'] = 'ss'
                    ss_df['interface_id'] = interface_id
                    ss_df['interface_type'] = interface_type
                    ss_df['interface_n_occ'] = nocc
                    ss_df['interface_pval'] = pval
                    ss_df['interface_area'] = area
                    ss_df['interface_solv_en'] = solv
                    ss_df['stab_en'] =stabe_en
                    ss_df['css'] =css
                    ss_df['overlap'] =overlap
                    ss_df['x-rel'] =x_rel
                    ss_df['fixed'] =fixed
        
                    h_df['PDB'] = PDBid
                    h_df['bond'] = 'hBond'
                    h_df['interface_id'] = interface_id
                    h_df['interface_type'] = interface_type
                    h_df['interface_n_occ'] = nocc
                    h_df['interface_pval'] = pval
                    h_df['interface_area'] = area
                    h_df['interface_solv_en'] = solv
                    h_df['stab_en'] =stabe_en
                    h_df['css'] =css
                    h_df['overlap'] =overlap
                    h_df['x-rel'] =x_rel
                    h_df['fixed'] =fixed
        
                    co_df['PDB'] = PDBid
                    co_df['bond'] = 'covalent'
                    co_df['interface_id'] = interface_id
                    co_df['interface_type'] = interface_type
                    co_df['interface_n_occ'] = nocc
                    co_df['interface_pval'] = pval
                    co_df['interface_area'] = area
                    co_df['interface_solv_en'] = solv
                    co_df['stab_en'] =stabe_en
                    co_df['css'] =css
                    co_df['overlap'] =overlap
                    co_df['x-rel'] =x_rel
                    co_df['fixed'] =fixed
        
                    sb_df['PDB'] = PDBid
                    sb_df['bond'] = 'saltBridge'
                    sb_df['interface_id'] = interface_id
                    sb_df['interface_type'] = interface_type
                    sb_df['interface_n_occ'] = nocc
                    sb_df['interface_pval'] = pval
                    sb_df['interface_area'] = area
                    sb_df['interface_solv_en'] = solv
                    sb_df['stab_en'] =stabe_en
                    sb_df['css'] =css
                    sb_df['overlap'] =overlap
                    sb_df['x-rel'] =x_rel
                    sb_df['fixed'] =fixed
                    
                    res_df['PDB'] = PDBid
                    res_df['interface_id'] = interface_id
                    res_df['interface_type'] = interface_type
                    res_df['interface_n_occ'] = nocc
                    res_df['interface_pval'] = pval
                    res_df['interface_area'] = area
                    res_df['interface_solv_en'] = solv
                    res_df['stab_en'] =stabe_en
                    res_df['css'] =css
                    
                    
        
                    # adding the dfs to be concatenated
                    ccat.append(ss_df)
                    ccat.append(h_df)
                    ccat.append(co_df)
                    ccat.append(sb_df)
                    resData.append(res_df)
                    
            
            except:
                failedPDB.append(PDBid)
                with open("20240714_failedPISAs.txt", "w") as outfile:
                    outfile.write("\n".join(str(item) for item in failedPDB))
                return pd.DataFrame()
    output_directory = './PISA resData CSVs'
    df2 = pd.concat(resData)
    df2.to_csv(os.path.join(output_directory,f'{PDBid}_PISA_resData_full.csv'))


df = pd.concat(ccat)

df.to_csv(os.path.join(directory2,'20240604_PISA_interfaceData_full.csv'))

