# GLYCO

GLYCO is to calculate number of glycan atoms per surface residue of protein antigen and glycan coverage of epitope residues.

1.Before you run GLYCO: There are some requirements you may have to check before running program.
 1.1 FreeSASA (https://freesasa.github.io/) and python3 should be installed in your system. 
 1.2 Coordinate section of input PDB file should follow the standard format of PDB (http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html).
 1.3 The name of protein residues in PDB file should be as below. Please check if your histidine is defined as one of below histidine names.
    ALA ARG ASN ASP CYS GLN GLU GLY HSD HID HIS HIE ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL
 1.4 Glycans in the input PDB files should be defined as either ATOM or HETATM

2. Download GLYCO

3. Run GLYCO
There are two modules in GLYCO script and user can choose based on their interest.
 3.1 If you want to count number of glycan atoms for each surface residue on your protein,
    command> python3 glyco.py pdbname.pdb cutoff residue glycan_name
    example> python3 glyco.py 5fyl.pdb 20 residue BMAN AMAN BGLN
 3.2 If you want to calculate glycan coverage (num of glycan atoms/buried surface area) for epitope residues on your protein,
    command> python3 glyco.py pdbname.pdb cutoff epitope epitope_residue_list glycan_name
    example> python3 glyco.py 5fyl.pdb 20 epitope epitope.txt BMAN AMAN BGLN
    
    *epitope.txt should have residue name, chain ID, residue number
    (epitope.txt)
     ARG C 309
     THR A 200
     MET A 196
     HIS C 305
