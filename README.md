# metaPhosphosites
 The collection of scripts in this repository are used to generate structural annotation of phoshphoryaltion sites in protein structures obtained from RCSB PDB, as published in: 
 
 Structural and systems characterization of phosphorylation on metabolic enzymes identifies sex-specific metabolic reprogramming in obesity
 
Tigist Y Tamir, Shreya Chaudhary, Annie X Li, Sonia E Trojan, Cameron T Flower, Paula Vo, Yufei Cui, Jeffrey C Davis, Rachit S Mukkamala, Francesca N Venditti, Alissandra L Hillis, Alex Toker, Matthew G Vander Heiden, Jessica B Spinelli, Norman J Kennedy, Roger J Davis, Forest M White

bioRxiv 2024.08.28.609894; doi: https://doi.org/10.1101/2024.08.28.609894

The following steps include workflow for domain proximity analysis:
 1. Gather PDBs for proteins of interest - getPDBs.py
 2. Evaluate PDB quality - Extract_pdbheaders.py, deltaR_jointPlot.py
 3. Annotate functional domains on proteins - extract_uniprot_Domain_annotations.py
 4. Map functional domains and phosphosites on PDBs -  extract_pdb_seq.py, 
 5. Measure distances from pSTY to domain COM -
 6. Count STY residues within 10 & 25 angstroms of domain COM - count_STY_withDist_DomainCOM.py
 7. Hypergeometric distribution test & multiple hypothesis testing - hypergeomeFDR.py

The following steps include workflow for interface/dimerization domain phosphosite annotaiton:
 1. Gather PISA analysis files in xml format using API - PDBePISA_retriver.py
 2. Parse PISA xml files and store in csv - PISA_xml_Parser.py
 3. Filter PISA csv for phosphosites only - phosFilter.py, PSP_annotated_STY_pdb.csv

Plot results summary - plots_Domains.py
