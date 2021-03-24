# GLYCO

GLYCO is to calculate number of glycan atoms per surface residue of protein antigen and glycan coverage of epitope residues.

1. Before you run GLYCO: There are some requirements you may have to check before running program.
 1.1 FreeSASA (https://freesasa.github.io/) and python3 should be installed in your system.
 1.2 Coordinate section of input PDB file should follow the standard format of PDB (http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html).
 1.3 The name of protein residues in PDB file should be as below. Especially, please check if your histidine is defined as one of below histidine names.
    ALA ARG ASN ASP CYS GLN GLU GLY HSD HID HIS HIE ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL
 1.4 Glycans in the input PDB files should be defined as either ATOM or HETATM

2. Download GLYCO

3. Run GLYCO
  
  3.1 A single frame (pdb)
     If you have a single pdb file and try to run GLYCO, you should follow below. 
     There are two modules in GLYCO script and user can choose based on their interest.
    
    3.1.1 Glycan atoms of each residue:
          (If you want to count number of glycan atoms for each surface residue on your protein,)
    command> python3 glyco.py -pdb pdbname.pdb -cutoff cutoff -module res -glycan glyca names -freesasa path of freesasa executable
    example> python3 glyco.py -pdb 5fyl.pdb -cutoff 20 -moduel res -glycan BMAN AMAN BGLN -freesasa /home/lee/freesasa
    3.1.2 Glycan coverage of epitope regions: 
          (If you want to calculate glycan coverage (num of glycan atoms/buried surface area) for epitope residues on your protein,)
    command> python3 glyco.py -pdb pdbname.pdb -cutoff cutoff -module ep -glycan glyca names -freesasa path of freesasa executable -epitope epitope list
    example> python3 glyco.py -pdb 5fyl.pdb -cutoff 20 -moduel ep -glycan BMAN AMAN BGLN -freesasa /home/lee/freesasa -epitope epitope.txt
    
    *epitope.txt should have following format: residue name, chain ID, residue number
    (epitope.txt)
     ARG C 309
     THR A 200
     MET A 196
     HIS C 305 
 
 3.2 Multiframes
