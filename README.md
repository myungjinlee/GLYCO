This script is to calculate glycan atoms per surface residue of protein antigen.

Before you run the script, you need to set up your VMD with "move_to_origin.tcl"

1. If you want to count number of glycan atoms for "all" surface residues of your pathogen,

python3 glycan_coverage.py origin_pdbname cutoff all current_pdb

ex) python3 glycan_coverage.py origin_frame_1.pdb 30 all frame_1.pdb

2. If you want to count number of glycan atoms for "epitope" residues of your pathogen,

python3 glycan_coverage.py origin_pdbname cutoff ep epitope_list

ex) python3 glycan_coverag.py origin_frame_1.pdb 30 ep ep_1.txt
