GLYCO manual 

GLYCO is to calculate number of glycan atoms per surface residue of protein antigen and glycan coverage of epitope residues.

Requirement to run GLYCO
1. Two programs should be installed in your system to run GLYCO, which are FreeSASA and python3.
2. The coordinate section of input PDB file should follow the standard format of PDB.
3. The name of protein residues in PDB file should be as below. Please check if your histidine is defined as one of below histidine names.
   ALA ARG ASN ASP CYS GLN GLU GLY HSD HID HIS HIE ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL


1. If you want to count number of glycan atoms for "all" surface residues of your pathogen,

python3 glycan_coverage.py origin_pdbname cutoff all current_pdb

ex) python3 glycan_coverage.py origin_frame_1.pdb 30 all frame_1.pdb

2. If you want to count number of glycan atoms for "epitope" residues of your pathogen,

python3 glycan_coverage.py origin_pdbname cutoff ep epitope_list

ex) python3 glycan_coverag.py origin_frame_1.pdb 30 ep ep_1.txt

PDB coordinate should follow standard format. Please refer http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html

Protein residue should be in the following set. 
