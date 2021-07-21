# GLYCO

GLYCO (GLYcan COverage) is a program to calculate glycan coverage of glycoproteins (e.g., number of glycan atoms per protein surface residues/epitope residues).

**1. Before you run GLYCO: There are some requirements you may have to check before running the program.<br />**
   - 1.1. FreeSASA (https://freesasa.github.io/) and python3 should be installed.<br />
   - 1.2. Coordinate section of input PDB files should follow the standard ATOM/HETATM record format of PDB (http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html) from column 1 to 54.<br />
   - 1.3. Protein residues in PDB files should be named as below. Especially, please check if your histidine is defined as one of below histidine names.<br />
    ALA ARG ASN ASP CYS GLN GLU GLY HSD HID HIS HIE ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL<br />

**2. Download GLYCO** 

**3. Run GLYCO<br />**
   - GLYCO takes the following arguments. Depending on the module and number of frames you have, the required arguments vary:<br />
 Sample code<br />
&nbsp;&nbsp;&nbsp;&nbsp;5th position in an really ugly code  
    5th position in a clear an readable code  
    Again using non-breaking spaces :)
    
   &nbsp;&nbsp;&nbsp;-----------------------------------------------------------------------------------------------------<br />
       &nbsp; &nbsp; &nbsp; -pdb&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; pdbname.pdb<br />
       &nbsp; &nbsp; &nbsp; -cutoff&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; cutoff in $\mathring A$   $\text{\r \textit{g}}$ \AA Angstrom<br />
       &nbsp; &nbsp; &nbsp; -module &emsp; &emsp; &emsp; &emsp;res or ep<br />
       &nbsp; &nbsp; &nbsp; -glycan &emsp; &emsp; &emsp; &emsp;list glycan names with comma separators<br />
       &nbsp; &nbsp; &nbsp; -freesasa&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp;path of FreeSASA executable<br />
       &nbsp; &nbsp; &nbsp; -probe &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; (optional)probe radius (1.4 A by default)<br />
       &nbsp; &nbsp; &nbsp; -sur_cutoff&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; &nbsp; (optional)cutoff to define surface (30 A^2 by default)<br />
       &nbsp; &nbsp; &nbsp; -epitope&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;(if module is "ep")File with epitope residues<br />
       &nbsp; &nbsp; &nbsp; -in_folder &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(if multiple PDBs)Input folder with multiple PDBs<br />
       &nbsp; &nbsp; &nbsp; -out_folder&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;Output folder to save results<br />
       &nbsp; &nbsp; &nbsp; -num_parallel &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; (if multiple PDBs)number of frames to submit in parallel (1 by default)<br />
       &nbsp; &nbsp; &nbsp; -num_proc_in&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; (optional)number of CPU cores to allocate (maximum number of cores by default)<br />
       &nbsp; &nbsp; &nbsp; -average &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; (optional)(if multiple PDBs from the same structure)to average glycan numbers over multiple PDBs<br />
   &nbsp;&nbsp;&nbsp;-----------------------------------------------------------------------------------------------------<br />
 
     
   If you want to ignore silence warnings from Python, you could use ``` python3 -W glyco.py ``` and add argments after that.<br />
   
   - 3.1. A single frame (PDB)<br />
     - 3.1.1. Glycan atoms of each residue -  module: "res":<br />
     
       - Count number of glycan atoms for each surface residue on your protein<br />
       ```
       python3 glyco.py -pdb 5fyl.pdb -cutoff 20 -module res -glycan BMA,AMA,BGL -freesasa /home/lee/freesasa -out_folder output 
       ```
       You can add arguments ```-probe -sur_cutoff -num_proc_in  ```as needed. <br />
       
       - Output<br /> 
       -- 5fyl_res_glycount.txt: number of glycan atoms per residue<br />
       -- 5fyl_bfactor.pdb: PDB file with glycan atoms as b-factor (You can visualize it with PyMOL.) <br />
       
     - 3.1.2. Glycan atoms of epitope residues - module: "ep":<br />
       
       - Calculate number of glycan atoms of epitope residues<br />
       ```
       python3 glyco.py -pdb 5fyl.pdb -cutoff 20 -module ep -glycan BMA,AMA,BGL -epitope epitope.txt -out_folder output
       ```
       You can add argument ```-num_proc_in  ```as needed. <br />
       *epitope.txt should be in the following format: residue name, chain ID, residue number<br />
         (epitope.txt)<br />
          ARG&nbsp; A&nbsp; 309<br />
          THR&nbsp; A&nbsp; 200<br />
          MET&nbsp; B&nbsp; 196<br />
          ASP&nbsp; C&nbsp; 305<br />
       
        - Output<br /> 
        -- ep_glysum.txt: summation of number of glycan atoms of the input epitope residues. This value excludes the overlapped, redundant glycan atoms shared among an epitope. <br /><br />
       Once you calculate buried surface area of your epitope and divide ep_glysum by the buried surface area, you can get epitope-glycan coverage. Calculating buried surface area is not provided by GLYCO, but there are many ways you can estimate the epitope-buried surface area such as Pisa(https://www.ebi.ac.uk/msd-srv/prot_int/cgi-bin/piserver) or making your own script with FreeSASA output. 
 
   - 3.2. Multiframes: If you have multiple frames of pdb files, you can submit multiple jobs in parallel.<br />
     - 3.1.1. Glycan atoms of each residues - module: "res":<br />
       - Count number of glycan atoms<br />
        *Input PDBs should be named as PREFIX_INDEX.pdb (e.g., frame_1.pdb, frame_2.pdb).<br />
        *(num_proc_in x num_parallel) should not exceed the total/available number of CPUs in your system.<br />
       ```
       python3 glyco.py -in_folder input -cutoff 20 -module res -glycan BMA,AMA,BGL -freesasa /data/leem/freesasa -num_proc_in 22 -num_parallel 2 -out_folder results -average
       ```
       - Output<br /> 
        -- PREFIX_INDEX_res_glycount.txt: number of glycan atoms per residue <br />
        -- PREFIX_INDEX_bfactor.pdb: PDB file with glycan atoms as b-factor (You can visualize it with PyMOL.) <br />
        -- ave_res_glycount.txt: averaged number of glycan atoms per residue <br /> 
     
     - 3.1.2. Glycan atoms of epitope regions - module: "ep":<br />
       - Count number of glycan atoms per epitope residue
         
       ```
       python3 glyco.py -in_folder input -cutoff 20 -module ep -glycan BMA,AMA,BGL -num_proc_in 28 -num_parallel 1 -epitope /home/leem/glyco/multiframes/epitope/epitope.txt -out_folder results -average
       ```
       - Output<br /> 
        -- ave_ep_glycount.txt: averaged number of glycan atoms per epitope <br />  
       
 Please report any bugs or questions to Myungjin Lee, Ph.D. (myungjin.lee@nih.gov)
      
