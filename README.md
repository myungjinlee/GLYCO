# GLYCO

GLYCO (GLYcan COverage) is a program to calculate glycan coverage of glycoproteins (e.g., number of glycan atoms per protein surface residues/epitope residues).

**1. Before you run GLYCO: There are some requirements you may have to check before running the program.<br />**
   - 1.1. FreeSASA(https://freesasa.github.io/) and python3 should be installed.<br />
   - 1.2. Coordinate section of input PDB files should follow the standard ATOM/HETATM record format of PDB (http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html) from column 1 to 54.<br />
   - 1.3. Protein residues in PDB files should be named as below. Especially, please check if your histidine is defined as one of below histidine names.<br />
    ALA ARG ASN ASP CYS GLN GLU GLY HSD HID HIS HIE ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL<br />

**2. Download GLYCO** 

**3. Run GLYCO<br />**
   - GLYCO takes the following arguments. Depending on the module and number of frames you have, the required arguments vary:<br />
    
   | Argument         | Input                                                      | Requirement                  |
   | ---------------- |------------------------------------------------------------| :----------------------------|
   | -pdb             | pdbname.pdb                                                | mandatory when single PDB    |
   | -in_folder       | Input folder that has multiple PDBs                        | mandatory when multiple PDBs |
   | -cutoff          | Cutoff to analyze glycan in Angstrom                       | mandatory                    |
   | -module          | Module name (all or sub)                                   | mandatory                    |
   | -glycan          | Glycan names with comma separators (Do not add space)      | mandatory                    |
   | -out_folder      | Output folder name to save result                          | mandatory                    |
   | -freesasa        | Path of FreeSASA executable                                | mandatory when module "all"  |
   | -residue         | File that has a list of user selected residues             | mandatory when module "sub"  |
   | -probe           | (if module all) Probe radius to define surface             | optional (1.4 &#197; by default)  |
   | -sur_cutoff      | (if module all) Cutoff to define surface residues          | optional (30 A^2 by default) |
   | -num_proc_in     | Number of CPU cores                                        | optional (max. # of cores by default)|
   | -num_parallel    | (if multiple PDBs) Number of frames to submit in parallel  | optional (1 by default)      |
   | -average         | (if multiple PDBs from the same structure) No input        | optional                     |
   
   
   - 3.1. A single frame (PDB)<br />
     - 3.1.1. Glycan coverage of overall protein surface residue -  module: "all":<br />
     
       - Count number of glycan atoms per all surface residues of the protein<br />
       ```
       python3 glyco.py -pdb 5fyl.pdb -cutoff 20 -module all -glycan BMA,AMA,BGL -freesasa /home/lee/freesasa -out_folder output 
       ```
       You can add arguments ```-probe```, ```-sur_cutoff```, ```-num_proc_in``` as needed. <br />
       
       - Output<br /> 
       -- 5fyl_all_glycount.txt: number of glycan atoms per residue<br />
       -- 5fyl_bfactor.pdb: PDB file with glycan atoms as b-factor (You can visualize it with PyMOL.) <br />
       
     - 3.1.2. Glycan coverage of user selected residues - module: "sub":<br />
       
       - Calculate number of glycan atoms of selected residues<br />
       ```
       python3 glyco.py -pdb 5fyl.pdb -cutoff 20 -module sub -glycan BMA,AMA,BGL -residue residuelist.txt -out_folder output
       ```
       You can add argument ```-num_proc_in``` as needed. <br /><br />
         &nbsp; &nbsp; *residuelist.txt should be in the following format: residue name, chain ID, residue number<br />
            &nbsp; &nbsp; (residuelist.txt)<br />
            &nbsp; &nbsp;  ARG&nbsp; A&nbsp; 309<br />
            &nbsp; &nbsp;  THR&nbsp; A&nbsp; 200<br />
            &nbsp; &nbsp;   MET&nbsp; B&nbsp; 196<br />
            &nbsp; &nbsp;   ASP&nbsp; C&nbsp; 305<br />
       
        - Output<br /> 
        -- sub_glysum.txt: summation of number of glycan atoms of the input residues. This value excludes the overlapped, redundant glycan atoms shared among a residue list. <br /><br />
 
   - 3.2. Multiframes: If you have multiple frames of pdb files, you can submit multiple jobs in parallel.<br />
     - 3.1.1. Glycan coverage of overall protein surface residue - module: "all":<br />
       - Count number of glycan atoms per all surface residues of the proteinbr />
        *Input PDBs should be named as PREFIX_INDEX.pdb (e.g., frame_1.pdb, frame_2.pdb) and placed in a folder.<br />
        *(num_proc_in x num_parallel) should not exceed the total/available number of CPUs in your system.<br />
       ```
       python3 glyco.py -in_folder input -cutoff 20 -module all -glycan BMA,AMA,BGL -freesasa /data/leem/freesasa -num_proc_in 22 -num_parallel 2 -out_folder results -average
       ```
       
       You can add arguments ```-probe```, ```-sur_cutoff``` as needed. <br />
       - Output<br /> 
        -- PREFIX_INDEX_all_glycount.txt: number of glycan atoms per residue <br />
        -- PREFIX_INDEX_bfactor.pdb: PDB file with glycan atoms as b-factor (You can visualize it with PyMOL.) <br />
        -- ave_all_glycount.txt: averaged number of glycan atoms per residue <br /> 
     
     - 3.1.2. Glycan coverage of user selected residues - module: "sub":<br />
       - Count number of glycan atoms per input residue
         
       ```
       python3 glyco.py -in_folder input -cutoff 20 -module sub -glycan BMA,AMA,BGL -num_proc_in 28 -num_parallel 2 -residue residuelist.txt -out_folder results -average
       ```
       - Output<br /> 
        -- ave_sub_glycount.txt: averaged number of glycan atoms per a set of input residue <br />  
       
 Please report any bugs or questions to Myungjin Lee, Ph.D. (myungjin.lee@nih.gov)
      
