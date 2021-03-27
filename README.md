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
   - 3.1. A single frame (pdb): If you have a single pdb file, you should follow below.<br />
     - 3.1.1. Glycan atoms of each residue:<br />
          (If you want to count number of glycan atoms for each surface residue on your protein,)<br />
       - command> python3 glyco.py -pdb pdbname.pdb -cutoff cutoff -module res -glycan glycan names -freesasa path of freesasa executable<br />
       - example> python3 glyco.py -pdb 5fyl.pdb -cutoff 20 -module res -glycan BMAN AMAN BGLN -freesasa /home/lee/freesasa<br />
       - There are a bunch of output files, but you want to focus on "res_count.txt" that has number of glycan atoms per residue. 
       - If you execute below command, you can visualize 
       - command> bfactor.py res_count.txt pdbname.pdb
       - example> bfactor.py res_count.txt frame_1.pdb
     - 3.1.2. Glycan coverage of epitope regions:<br />
          (If you want to calculate glycan coverage (num of glycan atoms/buried surface area) for epitope residues on your protein,)<br />
       - command> python3 glyco.py -pdb pdbname.pdb -cutoff cutoff -module ep -glycan glycan names -freesasa path of freesasa executable -epitope epitope list <br />
       - example> python3 glyco.py -pdb 5fyl.pdb -cutoff 20 -module ep -glycan BMAN AMAN BGLN -freesasa /home/lee/freesasa -epitope epitope.txt<br />
         *epitope.txt should have following format: residue name, chain ID, residue number<br />
         (epitope.txt)<br />
         ARG C 309<br />
         THR A 200<br />
         MET A 196<br />
         HIS C 305<br />
 
   - 3.2. Multiframes: If you have multiple frames of pdb files, you can submit multiple jobs in your HPC system. <br />
     - 3.1.1. Glycan atoms of each residues:<br />
       - command> bash multi_res_sub.sh 20 /data/SBIS/leem25/glycan_density/glyco/multiframes 3 5 
       - 
