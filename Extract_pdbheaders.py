#!/usr/bin/env python3
__author__ = "Huan Wang"
__version__ = "1.2"

#adapted by Tigist Tamir on Luisa

from decimal import Decimal, ROUND_HALF_UP
import numpy as np
import pandas as pd
import os, re, sys, time


from Bio.PDB import *
from Bio.PDB import parse_pdb_header
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from itertools import islice

items = ["Resolution", "GoodQ", "Median", "BadQ"]
grade = {"Resolution": [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
                        1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3,
                        2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0,
                        3.1, 3.2, 3.3, 3.4, 3.5, 4.0],
         "GoodQ":  [0.135, 0.145, 0.155, 0.162, 0.185, 0.190, 0.195,
                    0.200, 0.210, 0.215, 0.220, 0.228, 0.232, 0.238,
                    0.242, 0.245, 0.248, 0.250, 0.255, 0.257, 0.260,
                    0.265, 0.268, 0.270, 0.273, 0.275, 0.280],
         "Median": [0.150, 0.162, 0.175, 0.185, 0.200, 0.210, 0.215,
                    0.220, 0.228, 0.232, 0.240, 0.245, 0.250, 0.254,
                    0.258, 0.265, 0.268, 0.270, 0.273, 0.276, 0.280,
                    0.285, 0.290, 0.295, 0.300, 0.305, 0.310],
         "BadQ":   [0.165, 0.185, 0.195, 0.210, 0.220, 0.228, 0.232,
                    0.235, 0.245, 0.250, 0.260, 0.263, 0.266, 0.272,
                    0.275, 0.280, 0.285, 0.290, 0.293, 0.295, 0.297,
                    0.308, 0.310, 0.315, 0.320, 0.330, 0.350]}
rules = pd.DataFrame(grade, columns=items)


def find_PDB_files(path):
    suffix = ".pdb"
    return (f for f in os.listdir(path) if f.endswith(suffix))

    
def deal_round(number, n):
    ''' Rounded the input number (resolution) to n_th digit decimal. (here 0.1)
    For example：
    1.45 is rounded by the function of deal_round("1.45", 0.1),
    it will return 1.5.
    Pay attention, deal_round(1.45, 0.1) would give wrong result, 1.4!
    
    >>> Decimal(1.45)
    Decimal('1.4499999999999999555910790149937383830547332763671875')
    
    >>> Decimal("1.45")
    Decimal('1.45')
    '''
    # VERY IMPORTANT!
    # Here number is a string. It is better to use string type of number!
    # For example, Decimal instance is created from the string '1.45', 
    # and is converted straight to base-10.
    
    val = Decimal(number)  #### here, number is string type
    acc = str(n)  #### n = 0.1 or 0.01 or 0.001. Here, n = 0.1
    return float(Decimal(val.quantize(Decimal(acc), rounding=ROUND_HALF_UP)))


def calc_resolution_grade(resln):
    resln = float(resln)
    if resln < 1.6:
        return "EXCELLENT"
    elif 1.6 <= resln <= 1.79:
        return "EXCELLENT/VERY GOOD"
    elif 1.8 <= resln <= 1.99:
        return "VERY GOOD"
    elif 2.0 <= resln <= 2.29:
        return "VERY GOOD/GOOD"
    elif 2.3 <= resln <= 2.59:
        return "GOOD"
    elif 2.6 <= resln <= 2.89:
        return "GOOD/FAIR"
    elif 2.9 <= resln <= 3.19:
        return "FAIR"
    elif 3.2 <= resln <= 3.49:
        return "FAIR/POOR"
    elif resln >= 3.5:
        return "POOR"
        
        
def calc_R_free_grade(resln, R_free, rules):
    '''Note: here resln and R_free are strings '''
    
    if R_free.upper() == "NULL":
        return "NULL"
    
    elif 1.0 <= float(resln) <= 3.5:
        rounded = deal_round(resln, 0.1)

        array = rules[rules.Resolution.values == rounded]
        GoodQ  = array.GoodQ.values[0]
        Median = array.Median.values[0]
        BadQ   = array.BadQ.values[0]
        
        R_free = float(R_free)    
        if R_free <= (GoodQ - 0.02):
            return "MUCH BETTER THAN AVERAGE at this resolution"
    
        elif (GoodQ - 0.02) < R_free <= ((GoodQ + Median) / 2):
            return "BETTER THAN AVERAGE at this resolution"
    
        elif ((GoodQ + Median) / 2) < R_free <= ((Median + BadQ) / 2):
            return "AVERAGE at this resolution"
    
        elif ((Median + BadQ) / 2) < R_free <= (BadQ + 0.02):
            return "WORSE THAN AVERAGE at this resolution"
    
        elif R_free > (BadQ + 0.02):
            return "UNRELIABLE"
    
        else:
            return "Error!"
        

def parse_info(filename):
    R_value_list = []
    R_free_list  = []
    nonHatoms=[]
    reflec=[]

    str_R_value  = r"3\s+R VALUE\s+\(.+\)\s+:\s*(\d+\.\d+|NULL)"
    str_FREE_R   = r"3\s+FREE R VALUE(\s+|\s+\(.+\)\s+):\s*(\d+\.\d+|NULL)"
    str_B_factor = r"3\s+MEAN B VALUE\s+\(.+:\s+([-+]?\d*\.\d+|NULL)"
    str_bond_len_r = r"BOND LENGTHS REFINED ATOMS\s+\(.+:\s+(\d+)\s*;\s*([\d.]+)\s*;\s*([\d.]+)"
    str_bond_len_o = r"BOND LENGTHS OTHERS\s+\(.+:\s+(\d+)\s*;\s*([\d.]+)\s*;\s*([\d.]+)"
    str_bond_angles_r = r"BOND ANGLES REFINED ATOMS\s+\(.+:\s+(\d+)\s*;\s*([\d.]+)\s*;\s*([\d.]+)"
    str_bond_angles_o = r"BOND ANGLES OTHERS\s+\(.+:\s+(\d+)\s*;\s*([\d.]+)\s*;\s*([\d.]+)"

    pattern_R_value  = re.compile(str_R_value)
    pattern_FREE_R   = re.compile(str_FREE_R)
    pattern_B_factor = re.compile(str_B_factor)
    pattern_bond_len_r = re.compile(str_bond_len_r)
    pattern_bond_len_o = re.compile(str_bond_len_o)
    pattern_bond_angles_r = re.compile(str_bond_angles_r)
    pattern_bond_angles_o = re.compile(str_bond_angles_o)
    
    method = resln = resln_grade = R_value =\
    R_free = R_free_grade = B_value =Bond_len_r =Bond_len_o= Bond_angles_r = Bond_angles_o= "NULL"

    with open(filename, 'r') as fo:
        for line in fo:
            #### extracting method information
            if line.startswith("EXPDTA"):
                method = re.search(r'EXPDTA\s+(.+\w+)', line).group(1)
                print("Expt. Method:\t{:}".format(method))
            
            #### extracting number of reflections used to solve structure
            if 'NUMBER OF REFLECTIONS' in line and not('NULL') in line:
                reflec=re.findall("\d+",line)[1]
            
            #### extracting number of reflections used to solve structure
            if 'NUMBER OF NON-HYDROGEN ATOMS USED IN REFINEMENT.' in line:
                x2=re.split('\n',(''.join(islice(fo,4))))
                datadict={}
                for item in x2:
                    k_v=item.strip().split(':')
                    if len(k_v) == 2:
                        k, v = k_v
                        if 'B VALUE' not in k:
                            if'NULL' not in v:
                                if v!='':
                                    datadict[k]=int(v)
                                else:
                                    datadict[k]=int(0)
                        else:
                            break
                nonHatoms=sum(datadict.values())
                
            #### extracting the highest resolution
            if line.startswith("REMARK   2 RESOLUTION."):
                resln2 = re.search(r'REMARK   2 RESOLUTION.\s+(.+\w+)', line).group(1)
                print(resln2)
                if resln2!="NOT APPLICABLE" and method!="ELECTRON MICROSCOPY" and method!="SOLUTON NMR; THEORETICAL MODEL" and method!="SOLUTION NMR; SOLUTION SCATTERING":
                    resln = re.search(r'[-+]?\d*\.\d+', line).group()
                    print("Resolution:\t{:}".format(resln))
                    
                    resln_grade = calc_resolution_grade(resln)
                    print("Resolution grade:\t{:}".format(resln_grade))
                
            #### dealing with the R_value, R_free value and the average B value
            if line.startswith("REMARK   3"):

                #### extracting the R_value_Working_set
                match_R_work = re.search(pattern_R_value, line)
                if match_R_work:
                    Rvalue = match_R_work.group(1)
                    R_value_list.append(Rvalue)

                #### extracting the R_free value
                match_R_free = re.search(pattern_FREE_R, line)
                if match_R_free:
                    value = match_R_free.group(2)          
                    R_free_list.append(value)
                  
                #### extracting the R_free value
                match_R_free = re.search(pattern_FREE_R, line)
                if match_R_free:
                    value = match_R_free.group(2)
                    R_free_list.append(value)
                
                #### extracting the average B value    
                match_B_val = re.search(pattern_B_factor, line)
                if match_B_val:
                    B_value = match_B_val.group(1)
                    print("Mean B_value:\t{:}".format(B_value))
                
                #### extracting the RMSD Bond lengths refined atoms   
                match_bond_len_r = re.search(pattern_bond_len_r, line)
                if match_bond_len_r:
                    Bond_len_r = match_bond_len_r.group(2)
                    print("RMSD Bond Lengths Refined Atoms:\t{:}".format(Bond_len_r))
                    
                #### extracting the RMSD Bond lengths others    
                match_bond_len_o = re.search(pattern_bond_len_o, line)
                if match_bond_len_o:
                    Bond_len_o = match_bond_len_o.group(2)
                    print("RMSD Bond Lengths Others:\t{:}".format(Bond_len_o))
                
                #### extracting the RMSD Dihedral angles    
                match_bond_angles_r = re.search(pattern_bond_angles_r, line)
                if match_bond_angles_r:
                    Bond_angles_r = match_bond_angles_r.group(2)
                    print("RMSD Bond Angles Refined Atoms:\t{:}".format(Bond_angles_r))
                
                #### extracting the RMSD Dihedral angles    
                match_bond_angles_o = re.search(pattern_bond_angles_o, line)
                if match_bond_angles_o:
                    Bond_angles_o = match_bond_angles_o.group(2)
                    print("RMSD Bond Angles Others:\t{:}".format(Bond_angles_o))
                
                
            if line.startswith("ATOM"):
                break

        if method!= "SOLUTION NMR" and method!= "ELECTRON MICROSCOPY" and method!= "SOLID-STATE NMR" and method!="SOLUTION NMR; THEORETICAL MODEL"and method!="SOLUTION NMR; SOLUTION SCATTERING":
            R_value = min(R_value_list)
            print("R value:\t{:}".format(R_value))

            R_free = min(R_free_list)
            print("R_Free value:\t{:}".format(value))
        
            R_free_grade = calc_R_free_grade(resln, R_free, rules)
            print("R_free_grade:\t{:}".format(R_free_grade))
           
    return (method, resln, resln_grade, R_value, R_free, R_free_grade, B_value, reflec, nonHatoms, Bond_len_r, Bond_len_o, Bond_angles_r, Bond_angles_o)


def main():
    #Add the directory containing PDB files
    path = os.path.normpath('./PDBs')

    initial_time = time.time()
    pdbfiles = find_PDB_files(path)

    data          = {}
    PDB_ids       = []
    Expt_Methods  = []
    Resolutions   = []
    Resln_grades  = []
    R_values      = []
    R_free_values = []
    R_free_grades = []
    Mean_B_values = []
    Reflections_list = []
    ReflecNonHatoms_list = []
    Bond_len_r_list = []
    Bond_len_o_list = []
    Bond_angles_r_list = []
    Bond_angles_o_list = []
    title = ("PDB_id", "Method", "Resolution", "Resolution Grade",
             "R_value", "R_free", "R_free Grade", "Mean B_value", "Reflections", "ReflecNonHatoms", "RMSD_Bond_lengths_r", "RMSD_Bond_lengths_o", "RMSD_Bond_angles_r", "RMSD_Bond_angles_o")

    for i, f in enumerate(pdbfiles):
        start_time = time.time()
        begin = ''.join(("\n", "-" * 50, "\n", "No. {:}, file {:}"))
        print(begin.format(i + 1, f))

        fname = os.path.join(path, f)
        Method, Resolution, Resolution_grade, R_value, R_free, R_free_grade, B_value, Reflections, ReflecNonHatoms,Bond_len_r, Bond_len_o, Bond_angles_r, Bond_angles_o= parse_info(fname)
        
        step_time = time.time() - start_time
        print("\nTime used in this file: {:.3f} Seconds".format(step_time))

        PDB_ids.append(os.path.splitext(f)[0].upper())
        Expt_Methods.append(Method)
        Resolutions.append(Resolution)
        Resln_grades.append(Resolution_grade)
        R_values.append(R_value)
        R_free_values.append(R_free)
        R_free_grades.append(R_free_grade)
        Mean_B_values.append(B_value)
        Reflections_list.append(Reflections)
        ReflecNonHatoms_list.append(ReflecNonHatoms)
        Bond_len_r_list.append(Bond_len_r)
        Bond_len_o_list.append(Bond_len_o)
        Bond_angles_r_list.append(Bond_angles_r)
        Bond_angles_o_list.append(Bond_angles_o)
    print(Reflections_list, ReflecNonHatoms_list)
       

    df = pd.DataFrame(np.column_stack((PDB_ids,
                                       Expt_Methods,
                                       Resolutions,
                                       Resln_grades,
                                       R_values,
                                       R_free_values,
                                       R_free_grades,
                                       Mean_B_values,
                                       Reflections_list,
                                       ReflecNonHatoms_list,
                                       Bond_len_r_list,
                                       Bond_len_o_list,
                                       Bond_angles_r_list,
                                       Bond_angles_o_list)),
                                       columns = title)

    print("The Final Table is \n", df)

    outname = ''.join(("PDB_headers_", str(i + 1), "_2023_new", ".csv"))
    output = os.path.join(path, outname)
    df.to_csv(output, sep=',', index=False)
    
    total_time = time.time() - initial_time
    print("Work Completed. Used Time: {:.3f} Seconds".format(total_time))
    
    
if __name__ == "__main__":
    main()