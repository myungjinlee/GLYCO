# GLYCO

GLYCO is to calculate number of glycan atoms per surface residue of protein ("res" module) and glycan coverage of epitope residues ("ep" module).

**1. Before you run GLYCO: There are some requirements you may have to check before running program.<br />**
   - 1.1. FreeSASA (https://freesasa.github.io/) and python3 should be installed in your system.<br />
   - 1.2. Coordinate section of input PDB file should follow the standard format of PDB (http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html).<br />
   - 1.3. The name of protein residues in PDB file should be as below. Especially, please check if your histidine is defined as one of below histidine names.<br />
    ALA ARG ASN ASP CYS GLN GLU GLY HSD HID HIS HIE ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL<br />
   - 1.4. Glycans in the input PDB files should be defined as either ATOM or HETATM.<br />
   - 1.5. (Only if you want to run multiple frames) name of the pdb files should be frame_1.pdb, frame_2.pdb, frame_3.pdb ... etc.

**2. Download GLYCO**

**3. Run GLYCO<br />**
   - GLYCO takes the following arguments. Depending on the module and number of frames, you need to enter all of these or some:<br />
 
   ------------------------------------------------------------------<br />
       &nbsp; &nbsp; &nbsp; -pdb&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; pdbname.pdb;<br />
       &nbsp; &nbsp; &nbsp; -cutoff&nbsp; &nbsp; &nbsp; &nbsp; cutoff in Angstrom<br />
       &nbsp; &nbsp; &nbsp; -module&nbsp; &nbsp; &nbsp; res or ep<br />
       &nbsp; &nbsp; &nbsp; -glycan&nbsp; &nbsp; &nbsp; list glycan names with comma separator<br />
       &nbsp; &nbsp; &nbsp; -freesasa&nbsp; &nbsp; &nbsp;path of freesasa executable<br />
       &nbsp; &nbsp; &nbsp; -path&nbsp; &nbsp; &nbsp; &nbsp; path of current working directory<br />
   ------------------------------------------------------------------<br />
   
   - 3.1. A single frame (pdb): If you have a single pdb file, you should follow below.<br />
     - 3.1.1. Glycan atoms of each residue:<br />
     
       - Count number of glycan atoms for each surface residue on your protein<br />
       ```
       python3 glyco.py -pdb name.pdb -cutoff cutoff -module res -glycan glycan names -freesasa path of freesasa executable
       ```
       example)
       ```
       python3 glyco.py -pdb 5fyl.pdb -cutoff 20 -module res -glycan BMAN,AMAN,BGLN -freesasa /home/lee/freesasa
       ```
       There are a bunch of output files, but you want to focus on "res_count.txt" that has number of glycan atoms per residue.<br />
       
       - You can visualize it with bfactor script as shown below.<br /> 
       ```
       python3 bfactor.py res_count.txt pdbname.pdb
       ```
       example)
       ```
       python3 bfactor.py res_count.txt frame_1.pdb
       ```
       
     - 3.1.2. Glycan coverage of epitope regions:<br />
       
       - Calculate glycan coverage (num of glycan atoms/buried surface area) for epitope residues of your protein)<br />
       ```
       python3 glyco.py -pdb name.pdb -cutoff cutoff -module ep -glycan glycan names -epitope epitope list
       ```
       example)
       ```
       python3 glyco.py -pdb 5fyl.pdb -cutoff 20 -module ep -glycan BMAN,AMAN,BGLN -epitope epitope.txt
       ```
       *epitope.txt should have following format: residue name, chain ID, residue number<br />
         (epitope.txt)<br />
          ARG C 309<br />
          THR A 200<br />
          MET A 196<br />
          HIS C 305<br />
 
   - 3.2. Multiframes: If you have multiple frames of pdb files, you can submit multiple jobs in your HPC system. <br />
     - 3.1.1. Glycan atoms of each residues:<br />
       - Count number of glycan atoms
         1) Input pdbs should be named as frame_INDEX.pdb such as frame_1.pdb, frame_2.pdb etc and deposit them in folder name "input"
         2) Folders input, template, and script glyco.py should be all in your current working directory
       ```
       bash multi_res_run.sh -frame_start index of first frame -path path of current working directory -glycan glycan names (comma separated) -cutoff cutoff -freesasa path of freesasa executable
       ```
       example)
       ```
       bash multi_res_run.sh -frame_start 1 -frame_end 50 -path /home/leem/glyco/multiframes -glycan BMA,AMA -cutoff 20 -freesasa /data/leem/freesasa
       ```
       - Average number of glycan atoms over multiple frames: Once you finish calculating number of glycan atoms per each frame, you can average "res_count.txt" over the frames in this step. You have to run it in where all directories, "frames" are located. ($WORKING_DIR/$CUTOFF/res/)<br /> 
       ```
       python3 ave_mult.py -frame_start index of first frame -frame_end index of last frame
       ```
       example) 
       ```
       python3 ave_mult.py -frame_start 1 -frame_end 50 
       ```
       The output is "ave_res_count.txt"
     
       - You can visualize it with bfactor script as shown below.<br /> 
       ```
       python3 bfactor.py ave_res_count.txt pdbname.pdb
       ```
       example)
       ```
       python3 bfactor.py ave_res_count.txt frame_1.pdb
       ```
     - 3.1.2. Glycan coverage of epitope regions:<br />
       ```
       bash multi_res_run.sh -frame_start index of first -frame_end index of last frame -module ep -path path of working directory -glycan glycan names
       ```
       example)
       ```
       bash multi_res_run.sh -frame_start 1 -frame_end 50 -module ep -path /home/leem/glyco/multiframes -glycan BMA,AMA
       ```
       
     
