# GLYCO (version 1)

GLYCO (GLYcan COverage) is a program to calculate glycan coverage of glycoproteins (e.g., number of glycan atoms per protein surface residues/epitope residues).

**1. Before you run GLYCO: There are some requirements you may have to check before running the program.<br />**
   - 1.1. FreeSASA(https://freesasa.github.io/) 
   - 1.2. Python3 with modules (sys, os, pathlib, re, time, numpy, argparse, datetime, glob, fileinput, collections, multiprocessing, multiprocessing.pool, concurrent.futures)<br />
   - 1.3. Input: PDB of glycosylated protein
      - 1.3.1 Coordinate section of input PDB files should follow the standard ATOM/HETATM record format of PDB (http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html) from column 1 to 54.<br />
      - 1.3.2 Protein residues in PDB files should be named as below. Especially, please check if your histidine is defined as one of below histidine names.<br />
    ALA ARG ASN ASP CYS GLN GLU GLY HSD HIS HIE HID HSE HSP ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL<br />
      - Example input files use mannose-5 glycosylated protein with glycan name BGL, BMA, AMA.
        (Command, Input, Output described below can be found in the folder "Examples".)<br /><br />
        <img src="https://github.com/myungjinlee/GLYCO/blob/main/images/man5.png" width="350" height="230"><br />
        
**2. Download GLYCO** <br />
   - Download glyco.py<br /><br />

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
   | -sur_cutoff      | (if module all) Cutoff to define surface residues          | optional (30 &#197;<sup>2</sup> by default) |
   | -num_proc_in     | Number of CPU cores                                        | optional (max. # of cores by default)|
   | -num_parallel    | (if multiple PDBs) Number of frames to submit in parallel  | optional (1 by default)      |
   | -average         | (if multiple PDBs from the same structure) No input        | optional                     |
   
     
   - 3.1. A single frame (PDB)<br />
     - 3.1.1. Glycan coverage of overall protein surface residue -  module: "all":<br />
     
       - Count number of glycan atoms per all surface residues of the protein<br />
       ```
       python3 glyco.py -pdb frame_1.pdb -cutoff 20 -module all -glycan BMA,AMA,BGL -freesasa /home/lee/freesasa -out_folder output 
       ```
       You can add arguments ```-probe```, ```-sur_cutoff```, ```-num_proc_in``` as needed. <br />

       - Output<br /> 
         - frame_1_all_glycount.txt: number of glycan atoms per residue<br />
           <img src="https://github.com/myungjinlee/GLYCO/blob/main/images/F1.png" width="300" height="70"> 
         - frame_1_all_bfactor.pdb: number of glycans from “frame_1_all_glycount.txt” are embedded in the PDB as B-factors (You can visualize it with PyMOL.)<br /> 
           <img src="https://github.com/myungjinlee/GLYCO/blob/main/images/F2.png" width="680" height="70">                         
         - frame_1_all_glysum.txt: summation of number of glycans in the PDB<br />
           <img src="https://github.com/myungjinlee/GLYCO/blob/main/images/F3.png" width="200" height="70">
         - frame_1_all_outer.rsa: output of running Freesasa <br />
        GLYCO utilizes the 5th column (e.g., All-atoms, ABS) to evaluate if the residue is a surface residue or not <br />
          <img src="https://github.com/myungjinlee/GLYCO/blob/main/images/F4.png" width="550" height="70"> <br />
          
     - 3.1.2. Glycan coverage of user selected residues - module: "sub":<br />
       
       - Calculate number of glycan atoms of selected residues<br />
       ```
       python3 glyco.py -pdb frame_1.pdb -cutoff 20 -module sub -glycan BMA,AMA,BGL -residue residuelist.txt -out_folder output
       ```
       You can add argument ```-num_proc_in``` as needed. <br /><br />
         &nbsp; &nbsp; *residuelist.txt should be in the following format: residue name, chain ID, residue number<br />
            &nbsp; &nbsp; (residuelist.txt)<br />
            &nbsp; &nbsp;  LEU&nbsp; A&nbsp;  95<br />
            &nbsp; &nbsp;  VAL&nbsp; A&nbsp;  97<br />
            &nbsp; &nbsp;  THR&nbsp; A&nbsp;  98<br />

        - Output<br /> 
          - frame_1_sub_glysum.txt: summation of number of glycan atoms of the input residue list. This value excludes the overlapped, redundant glycan atoms shared among a residue list. <br />
            <img src="https://github.com/myungjinlee/GLYCO/blob/main/images/F5.png" width="200" height="70">
          
   - 3.2. Multiframes: If you have multiple frames of pdb files, you can submit multiple jobs in parallel.<br />
     - 3.2.1. Glycan coverage of overall protein surface residue - module: "all":<br />
       - Count number of glycan atoms per all surface residues of the protein <br />
        *Input PDBs should be named as PREFIX_INDEX.pdb (e.g., frame_1.pdb, frame_2.pdb) and placed in a folder.<br />
        *(num_proc_in x num_parallel) should not exceed the total/available number of CPUs in your system.<br />
       ```
       python3 glyco.py -in_folder input -cutoff 20 -module all -glycan BMA,AMA,BGL -freesasa /data/leem/freesasa -num_proc_in 22 -num_parallel 2 -out_folder results -average
       ```
       
       You can add arguments ```-probe```, ```-sur_cutoff``` as needed. <br />
       - Output<br /> 
         - frame_1..5_all_glycount.txt: number of glycan atoms per residue<br />
         - frame_1..5_bfactor.pdb: number of glycans from “frame_1_all_glycount.txt” are embedded in the PDB as B-factors (You can visualize it with PyMOL.)<br /> 
         - ave_all_glycount.txt: averaged number of glycan atoms per surface residue over number of PDBs<br />
           <img src="https://github.com/myungjinlee/GLYCO/blob/main/images/F6.png" width="300" height="70">
         - ave_all_glysum.txt: averaged summation of number of glycans for all surface residues over number of PDBs <br /> 
           <img src="https://github.com/myungjinlee/GLYCO/blob/main/images/F7.png" width="200"  height="70">
        
     - 3.2.2. Glycan coverage of user selected residues - module: "sub":<br />
       - Count number of glycan atoms per input residue
         
       ```
       python3 glyco.py -in_folder input -cutoff 20 -module sub -glycan BMA,AMA,BGL -num_proc_in 28 -num_parallel 2 -residue residuelist.txt -out_folder results -average
       ```
       - Output<br /> 
         - frame_1_sub_glysum.txt: summation of number of glycan atoms of the input residue list.<br />
         - ave_sub_glysum.txt: averaged number of glycan atoms for input residue list <br />
           <img src="https://github.com/myungjinlee/GLYCO/blob/main/images/F8.png" width="180" height="70"> <br /> 
 
 Please report any bugs or questions to Myungjin Lee, Ph.D. (myungjin.lee@nih.gov)
      
