# GLYCO

GLYCO is a program to calculate number of glycan atoms per surface residue of proteins ("res" module) and number of glycan atoms per epitope residues ("ep" module).

**1. Before you run GLYCO: There are some requirements you may have to check before running the program.<br />**
   - 1.1. FreeSASA (https://freesasa.github.io/) and python3 should be installed.<br />
   - 1.2. Coordinate section of input PDB files should follow the standard ATOM/HETATM record format of PDB (http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html) from column 1 to 54.<br />
   - 1.3. Protein residues in PDB files should be named as below. Especially, please check if your histidine is defined as one of below histidine names.<br />
    ALA ARG ASN ASP CYS GLN GLU GLY HSD HID HIS HIE ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL<br />
   - 1.4. Glycans in input PDB files should be defined as either ATOM or HETATM.<br />

**2. Download GLYCO** 
   - Download glyco.py, template folder in your current working directory.<br />

**3. Run GLYCO<br />**
   - GLYCO takes the following arguments. Depending on the module and number of frames you have, the required arguments vary:<br />
 
   &nbsp;&nbsp;&nbsp;------------------------------------------------------------------<br />
       &nbsp; &nbsp; &nbsp; -pdb&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; pdbname.pdb<br />
       &nbsp; &nbsp; &nbsp; -cutoff&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; cutoff in Angstrom<br />
       &nbsp; &nbsp; &nbsp; -module&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; res or ep<br />
       &nbsp; &nbsp; &nbsp; -glycan&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; list glycan names with comma separators<br />
       &nbsp; &nbsp; &nbsp; -freesasa&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp;path of FreeSASA executable<br />
       &nbsp; &nbsp; &nbsp; -path&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; path of current working directory<br />
       &nbsp; &nbsp; &nbsp; -frame_start&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp;&nbsp;index of first frame<br />
       &nbsp; &nbsp; &nbsp; -frame_end&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; &nbsp; index of last frame<br />
       &nbsp; &nbsp; &nbsp; -frame_gap&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; &nbsp; number of frames to bundle<br />
       &nbsp; &nbsp; &nbsp; -num_proc&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;number of CPU cores to allocate(optional) (maximum number of cores by default)<br />
       &nbsp; &nbsp; &nbsp; -probe&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; probe radius(optional) (1.4 A by default)<br />
       &nbsp; &nbsp; &nbsp; -sur_cutoff&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;cutoff to define surface (optional) (30 A^2 by default)<br />
   &nbsp;&nbsp;&nbsp;------------------------------------------------------------------<br />
   
   - 3.1. A single frame (PDB): If you have a single PDB file, you should follow below.<br />
     - 3.1.1. Glycan atoms of each residue -  module: "res":<br />
     
       - Count number of glycan atoms for each surface residue on your protein<br />
       ```
       python3 glyco.py -pdb 5fyl.pdb -cutoff 20 -module res -glycan BMA,AMA,BGL -num_proc 28 -freesasa /home/lee/freesasa -probe 1.4 -surf_cutoff 30
       ```
       There are a bunch of output files, but you want to focus on "res_count.txt" that has number of glycan atoms per residue.<br />
       
       - You can visualize it with bfactor script as shown below.<br /> 
       ```
       python3 bfactor.py res_count.txt frame_1.pdb
       ```
       You can open the output "frame_1_bfactor.pdb" with PyMOL. 
       
     - 3.1.2. Glycan atoms of epitope residues - module: "ep":<br />
       
       - Calculate number of glycan atoms of epitope residues<br />
       ```
       python3 glyco.py -pdb 5fyl.pdb -cutoff 20 -module ep -glycan BMA,AMA,BGL -num_proc 28 -epitope epitope.txt 
       ```
       *epitope.txt should be in the following format: residue name, chain ID, residue number<br />
         (epitope.txt)<br />
          ARG&nbsp; A&nbsp; 309<br />
          THR&nbsp; A&nbsp; 200<br />
          MET&nbsp; B&nbsp; 196<br />
          ASP&nbsp; C&nbsp; 305<br />
       
       The output is "ep_glysum.txt" which has summation of number of glycan atoms of the input epitope residues. This value excludes the overlapped, redundant glycan atoms shared among an epitope. 
       Once you calculate buried surface area of your epitope and divide ep_glysum by the buried surface area, you can get epitope-glycan coverage. Calculating buried surface area is not provided by GLYCO, but there are many ways you can estimate the epitope-buried surface area such as Pisa(https://www.ebi.ac.uk/msd-srv/prot_int/cgi-bin/piserver) or making your own script with FreeSASA output. 
 
   - 3.2. Multiframes: If you have multiple frames of pdb files, you can submit multiple jobs in your HPC system. <br />
     - 3.1.1. Glycan atoms of each residues - module: "res":<br />
       - Count number of glycan atoms
         1) Input PDBs should be named as frame_INDEX.pdb (e.g., frame_1.pdb, frame_2.pdb) and located in folder "input"
         2) The "input" folder, "template" folder, and glyco.py should be in your current working directory. 
         3) The current working directory in 2) should be entered in the argument "-path". 
         4) Change line ## in the code that works for your HPC system.
       ```
       bash multi_glyco.sh -cutoff 20 -frame_start 1 -frame_end 50 -path /home/leem/glyco/multiframes -glycan BMA,AMA,BGL -num_proc 28 -freesasa /data/leem/freesasa -module res
       ```
       - Average number of glycan atoms over multiple frames: Once you finish calculating number of glycan atoms per each frame, you can average "res_count.txt" over the multiple frames. You have to run it in where all directories, (e.g., "frames") are located. ($WORKING_DIR/$CUTOFF/res/)<br /> 
       ```
       python3 ave_res.py -frame_start 1 -frame_end 50 
       ```
       The output is "ave_res_count.txt"
     
       - You can visualize it with bfactor script as shown below.<br /> 
       ```
       python3 bfactor.py ave_res_count.txt frame_1.pdb
       ```
       You can open the output "frame_1_bfactor.pdb" with PyMOL to display the glycan coverage per protein surface residue. 
     - 3.1.2. Glycan atoms of epitope regions - module: "ep":<br />
       - Count number of glycan atoms per epitope residue
         1) Should follow the same instruction in Section 3.1.1. 1)-4).
         2) The script groups several jobs in one bundle and submit each bundled job in one node. You should specify the number of jobs for a bundle in argument "-frame_gap". (Each epitope-glycan coverage normally finishes in a short time, which may reduce the efficiency of calculation in HPC system if each job was distributed per node. In this case, bundling can solve the problem.)
       ```
       bash multi_glyco.sh -cutoff 20 -frame_start 1 -frame_end 50 -frame_gap 10 -module ep -path /home/leem/glyco/multiframes -glycan BMA,AMA,BGL -num_proc 28 -epitope /home/leem/glyco/multiframes/epitope/epitope.txt
       ```
       - Average number of glycan atoms over multiple frames: Once you finish calculating number of glycan atoms per each frame, you can average "ep_glysum.txt" over the multiple frames. You have to run it in where all directories, (e.g., "frames") are located. ($WORKING_DIR/$CUTOFF/ep/)<br /> 
       ```
       python3 ave_ep.py -frame_start 1 -frame_end 50 
       ```
       The output is "ave_ep_gly.txt"     
       
 Please report any bugs or questions to Myungjin Lee, Ph.D. (myungjin.lee@nih.gov)
      
