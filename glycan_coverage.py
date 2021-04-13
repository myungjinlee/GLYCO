#!/usr/bin/env python
# coding: utf-8

import sys
import os
from pathlib import Path
import re
import time
import numpy
import argparse
import datetime
from collections import defaultdict
from collections import OrderedDict
import numpy as np

# For multiprocessing
import multiprocessing 
from multiprocessing.pool import Pool

# For progress bars
from tqdm import tqdm


# ---------------------------------------------------------------------------------#
#                                 Helper methods                                  #
# ---------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------#
# Generate a surface residue dictionary 
# ---------------------------------------------------------------------------------#
def outermost(inpdb, infloat, freesasa_path):
    print('PROGRESS: Excluding buried residues of the input protein', flush=True)
     
    outrsa = Path(inpdb).name.replace(".pdb", "") + "_outer.rsa"
    
    os.system(freesasa_path + ' ' + str(inpdb) + " -p 1.4 --rsa > " + str(outrsa))
    keys = []
    outdict = {}
    surf_dict = {}
    with open(outrsa, "r") as f:
        lines = f.readlines()
        for line in lines:
            atype = line[0:3].strip()
            resname = line[3:8].strip()
            chain = line[8].strip()
            resnum = line[9:14].strip()
            surf_area = line[14:22].strip()
            # RES GLU A 185A  119.57
            if atype == 'RES' and resname in ('ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU',
                                              'GLY', 'HSD', 'HIS', 'HIE', 'HID', 'ILE', 'LEU',
                                              'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
                                              'TYR', 'VAL'):
                surf_key = (resname, chain, resnum)
                surf_dict[surf_key] = float(surf_area)
                if float(surf_area) >= infloat:
                    keys.append([resname, chain, resnum])
                    # 'RES     ALA      A       1    108.06'
                    
    keys = numpy.array(keys)
    with open(inpdb, "r") as f:
        lines = f.readlines()
        for line in tqdm(lines, desc="Excluding buried residues of the input protein"):
            if line[17:20] in ('ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HSD',
                               'HIS', 'HIE', 'HID', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                               'SER', 'THR', 'TRP', 'TYR', 'VAL'):
                for k in keys:
                    atomnum = line[6:11].strip()
                    atomtype = line[12:16].strip()
                    resname = line[17:20].strip()
                    chain = line[21:22].strip()
                    resid = line[22:28].strip()
                    # ATOM  14823  N   ARG A1000      -8.710  -0.790   7.690
                    if k[0] == resname and k[2] == resid: 
                        key = (atomnum, atomtype, resname, chain, resid)
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        if atomtype.find('H') == -1:
                            outdict[key] = [x, y, z]
                                           
    return outdict, surf_dict


# --------------------------------------------------------------------------------- #
# Epitope file to dictionary
# --------------------------------------------------------------------------------- #
def input_to_dict(infile, indict):
    outdict = {}
    for line in open(infile):
        line = line.split()
        if len(line) == 3:
            # LYS A  121
            chain = line[1]
            resid = line[2]
            for key, val in indict.items():
                if chain == key[3] and resid == key[4]:
                    outdict[key] = val
    return outdict


# --------------------------------------------------------------------------------- #
# Generate a protein dictionary that has key and coordinate as values
# and, generate another protein dictionary that has x as key and the rest as values
# the second dictionary is for checking crossing with protein region
# input: pdb
# output1: dict[(atomnum, atomtype, resname, chain, resid)] = [x,y,z] 
# output2: dict[x] = [y,z,(atomnum, atomtype, resname, chain, resid)]
# --------------------------------------------------------------------------------- #
def gen_pro_dict(inpdb):
    print('PROGRESS: Generating a protein dictionary')
    outdict1 = {}
    outdict2 = defaultdict(list) 
    for line in open(inpdb):
        if line[0:4] == 'ATOM':
            if line[17:20] in ('ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HSD',
                               'HIS', 'HIE', 'HID', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                               'SER', 'THR', 'TRP', 'TYR', 'VAL'):
                resname = line[17:21].strip()
                chain = line[21].strip()
                resid = line[22:28].strip()
                atomnum = line[6:11].strip()
                atomtype = line[12:16].strip()
                key = (atomnum, atomtype, resname, chain, resid)
                xf = float(line[30:38].strip())
                yf = float(line[38:46].strip())
                zf = float(line[46:54].strip())
                x = int(float(line[30:38].strip()))
                y = int(float(line[38:46].strip()))
                z = int(float(line[46:54].strip()))
                if atomtype.find('H') == -1:
                    outdict1[key] = [xf, yf, zf]
                    outdict2[x].append([y, z, key])  # to check crossing
    return outdict1, outdict2


# ---------------------------------------------------------------------------------#
# Generate a glycan dictionary that has glycan keys and values as glycan coordinate
# input: pdb
# output: dict[(atomnum, atomtype, resname, chain, resid)] = [x,y,z]
# ---------------------------------------------------------------------------------#
def gen_gly_dict(inpdb, input_glycan):
    print('PROGRESS: Generating a glycan dictionary') 
    outdict = {}
    for line in open(inpdb):
        if line[0:4] == 'ATOM' or line[0:4] == 'HETA':
            if line[17:21].strip() in input_glycan:
                resname = line[17:21].strip()
                chain = line[21].strip()
                resid = line[22:28].strip()
                atomnum = line[6:11].strip()
                atomtype = line[12:16].strip()
                key = (atomnum, atomtype, resname, chain, resid)
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                if atomtype.find('H') == -1:
                    outdict[key] = [x, y, z]
    return outdict


# ---------------------------------------------------------------------------------#
# Extract epitope residues from a dictionary and generate an epitope dictionary
# input: epitope array, total dictionary
# output: dict[(atomnum, atomtype, resname, chain, resid)] = [x,y,z]
# ---------------------------------------------------------------------------------#
def ep_to_dict(inarr, indict):
    outdict = {}
    for ep in inarr:
        for key, val in indict.items():
            if ep[0] == key[2] and ep[1] == key[4]:
                outdict[key] = val
    return outdict


# ---------------------------------------------------------------------------------#
# Calculate distances between two coordinates
# input: list (coordinates)
# output: float (distance)
# ---------------------------------------------------------------------------------#
def cal_dist(inlist1, inlist2):
    inarr1 = numpy.array(inlist1)
    inarr2 = numpy.array(inlist2)
    outfloat = numpy.linalg.norm(inarr1-inarr2)
    return outfloat


# ---------------------------------------------------------------------------------#
# Generate a vector between protein residue and glycan atom
# input: 
# output:  
# ---------------------------------------------------------------------------------#
def vector(inlist1, inlist2):
    outlist = []
    inarr1 = numpy.array(inlist1)
    inarr2 = numpy.array(inlist2)
    leng = numpy.linalg.norm(inarr2-inarr1)
    leng = int(round(leng))
    for i in range(1, int(round(leng))): 
        t = i/float(leng)
        point = inarr1 + t*(inarr2-inarr1) 
        outlist.append(point)
    return outlist


def check_protein(inlist, intuple, pro_dict, ori_pro_dict):

    flag = True
    tmparr = numpy.array(inlist)
    tmparr = numpy.around(tmparr, decimals=0)
    for tmp in tmparr:
        # These are the coordinates of the parametrized line
        x = int(tmp[0])
        y = int(tmp[1])
        z = int(tmp[2]) 

        # To speed up things, only look at atoms in slices of x-axis
        filter_x = pro_dict.get(x, []) + pro_dict.get(x - 1, []) + pro_dict.get(x + 1, [])  
        
        if filter_x:
            for elem in filter_x:
                # We used the integer to filter by x slices, now retrieve the real float coordinates
                elem_float_x, elem_float_y, elem_float_z = ori_pro_dict[elem[2]]
            
                if abs(y - elem_float_y) <= 1 and abs(z - elem_float_z) <= 1 and abs(x - elem_float_x) <= 1:  
                        if not (elem[2][2] == intuple[2] and elem[2][3] == intuple[3] and elem[2][4] == intuple[4]):
                            flag = False

    return flag


def merge_dict_lists(added_dict=None, inplace_dict=None):    
    for k, v in added_dict.items(): 
        if k in inplace_dict.keys(): 
            inplace_dict[k] += v 
        else: 
            inplace_dict[k] = v 
            
            
# ---------------------------------------------------------------------------------#
#                                  Main methods                                    #
# ---------------------------------------------------------------------------------#


def process_cmd_args():
    # --------------------------------------------------------------------------------- #
    # Taking parameters
    # --------------------------------------------------------------------------------- #
    
    args_ = {}    

    parser = argparse.ArgumentParser()
    parser.add_argument(' ')


    parser = argparse.ArgumentParser(usage = '\n'+ '\n' +
             'module1: To calculate number of glyan atoms of residues on the protein surface' + '\n'
             'python3 glyco.py -pdb pdbname.pdb -cutoff 30 -glycan BMA AMA -module res -freesasa /home/lee/freesasa' + '\n' + '\n' +
             'module2: To calculate glyan coverage of epitope residues' + '\n'
             'python3 glyco.py -pdb pdbname.pdb -cutoff 30 -glycan BMA AMA -module ep -freesasa /home/lee/freesasa -epitope epitope_list.txt')
    parser.add_argument('-pdb', type=str, help='Enter your name of pdb. ex) input.pdb', required=True)
    parser.add_argument('-cutoff', type=float, help='Enter your distance cutoff in Angstrom. ex) 30', required=True)
    parser.add_argument('-glycan', type=str, help='Enter your name of glycans with space for separators. ex) BMA AMA BGLN', required=True)
    parser.add_argument('-module', type=str, help='Enter your module, either res or ep. ex) res', required=True)
    parser.add_argument('-epitope', type=str, help='Enter your epitope list. ex) epitope_list.txt')
    parser.add_argument('-freesasa', type=str, help='Enter your path of Freesasa, ex) /home/lee/freesasa')
    parser.add_argument('-nproc', type=int, help='Enter the number of workers to use', default=multiprocessing.cpu_count())

    args = parser.parse_args()

    if args.module == 'res' and args.freesasa == None:
        print('Please add freesasa path')
        syse.exit(1) 
    elif args.module == 'ep' and args.epitope == None:
        print('Please add epitope list')
        sys.exit(1)

    args_["sur_cutoff"] = 30
    args_["origin_pdb"] = args.pdb
    args_["dist_cutoff"] = int(args.cutoff)
    args_["input_glycan"] = args.glycan.replace(',', ' ')
    args_["res_or_ep"] = args.module
    args_["input_epitope"] = args.epitope
    args_["freesasa_path"] = args.freesasa 
    args_["nproc"] = args.nproc
    args = args_
    
    return args


def preprocess_dicts(origin_pdb, sur_cutoff, input_glycan, input_epitope, res_or_ep, freesasa_path):
    # ---------------------------------------------------------------------------------#
    # Generate dictionaries with input files
    # ---------------------------------------------------------------------------------#
    
    ori_pro_dict, pro_dict = gen_pro_dict(origin_pdb)
    # pro_dict = gen_pro_dict(origin_pdb)[1]
    gly_dict = gen_gly_dict(origin_pdb, input_glycan)
    # surf_dict = outermost(origin_pdb, sur_cutoff)[1]

    res_dict = {}
    if res_or_ep == 'res':
        res_dict = outermost(origin_pdb, sur_cutoff, freesasa_path)[0]

    elif res_or_ep == 'ep':
        res_dict = input_to_dict(input_epitope, ori_pro_dict)
        
    return ori_pro_dict, pro_dict, gly_dict, res_dict


def launch_workers(gly_dict, res_dict, pro_dict, nproc, dist_cutoff, ori_pro_dict):
    
    glycan_data_list = list(gly_dict.items())
    protein_data_list = list(res_dict.items())

    # Set the longest of {protein, glycan} lists to be the outter loop when generating pairs, and partition outter into
    # sublists to send to multipel CPUs with full inner list
    if len(glycan_data_list) > len(protein_data_list):
        outter_data = glycan_data_list
        inner_data = protein_data_list
        glycan_out = True
    else:
        outter_data = protein_data_list 
        inner_data = glycan_data_list
        glycan_out = False

    # Split calculations into multiple bundles
    outter_data = np.array(outter_data, dtype=object)
    bundle_num = nproc * 2
    bundle_outter_data = np.array_split(outter_data, bundle_num, axis=0) 

    # Store results 
    all_results = {}

    # Define a pool: nproc CPUs will compute simultaneously at all times, and taking care of all jobs in a queue
    tp = Pool(processes=nproc) 

    # Loop through all our jobs and add them to the queue 
    for job_id, outter_data in enumerate(bundle_outter_data): 

        if glycan_out:
            glycan_mini = outter_data
            protein_mini = inner_data
        else:
            glycan_mini = inner_data
            protein_mini = outter_data

        # Add worker to queue
        res = tp.apply_async(worker, (job_id, protein_mini, glycan_mini, pro_dict, dist_cutoff, ori_pro_dict)) 

        # Save AsyncResult object
        all_results[job_id] = res 

    # Close and join pool, waits for all jobs to be added to queue on main thread   
    print("Added all jobs to queue.\n") 
    tp.close()
    tp.join() 
    
    return all_results


# Performs the main work
def worker(job_id, res_dict, gly_dict, pro_dict, dist_cutoff, ori_pro_dict):
    # ---------------------------------------------------------------------------------#
    # Main loop
    # First  for-loop: Go through surface protein
    # Second for-loop: Go through glycan atoms
    # Inside for-loop: 1.Select and count glycans within cutoff based on pair distances
    #                  2.Check if glycans are covered by protein and exclude them
    # ---------------------------------------------------------------------------------#  
    print("Starting worker {}".format(job_id), flush=True)
    
    count_dict = {}
    num_dict = {}
    tmp_dict = defaultdict(list)
    for res_dict_item in res_dict:
        key1, val1 = res_dict_item
        res_key = (key1[2], key1[3], key1[4])
        # if res_key not in tmp_dict.keys():
        #    tmp_dict[res_key] = []
        count = 0 
        gly_list = []
        for gly_dict_item in gly_dict:
            key2, val2 = gly_dict_item
            dist = cal_dist(val1, val2)
            if dist < dist_cutoff:

                line = vector(val1, val2)
                flag = check_protein(line, key1, pro_dict, ori_pro_dict)
                
                if flag:
                    count += 1
                    gly_list.append(key2[0])
              #  else:
              #      pass
                
        tmp_dict[res_key].extend(gly_list)
        count_dict[key1] = gly_list
        num_dict[key1] = count

    # to remove redundant glycan atoms
    final_dict = {}
    for key, val in tmp_dict.items():
        final_dict[key] = list(OrderedDict.fromkeys(val))
        
    print("Finishing worker {}".format(job_id), flush=True)
        
    return count_dict, num_dict, final_dict


def aggregate_worker_res(all_results):
    
    # Retrieve results from calls
    print("\nCollecting results...")
    all_results = {k: v.get() for k, v in all_results.items()}
    
    # Aggregate calls from each worker
    count_dict = {}
    num_dict = {}
    final_dict = {}

    for all_results_i in all_results.values():
        count_dict_i, num_dict_i, final_dict_i = all_results_i

        merge_dict_lists(added_dict=count_dict_i, inplace_dict=count_dict)
        merge_dict_lists(added_dict=num_dict_i, inplace_dict=num_dict)
        merge_dict_lists(added_dict=final_dict_i, inplace_dict=final_dict)
    
    return count_dict, num_dict, final_dict


def write_outputs(res_or_ep, count_dict, num_dict, final_dict):
    
    if res_or_ep == 'res':
    
        # ---------------------------------------------------------------------------------#
        # Write an output 'atom_gly.txt'
        # protein atom, [list of glycan atom number per protein atom]  
        # ---------------------------------------------------------------------------------#
    #    with open('atom_gly.txt', 'w') as f_obj_1:
    #        for key, val in count_dict.items():
    #            line = str(key) + ' ' + str(val) + '\n'
    #            f_obj_1.write(line)

        # ---------------------------------------------------------------------------------#
        # Write an output 'atom_count.txt'
        # protein atom, number of glycan atoms per protein atom 
        # ---------------------------------------------------------------------------------#
    #    with open('atom_count.txt', 'w') as f_obj_2:
    #        for key, val in num_dict.items():
    #            line = str(key) + ' ' + str(val) + '\n' 
    #            f_obj_2.write(line)

        # ---------------------------------------------------------------------------------#
        # Write outputs 'res_glycount.txt'
        # Protein residue, [list of glycan atom number per protein residue] 
        # Protein residue, number of glycan atoms 
        # ---------------------------------------------------------------------------------#
    #    with open('res_gly.txt', 'w') as f_obj_3, open('res_count.txt', 'w') as f_obj_4:
        with open('res_glycount.txt', 'w') as f_obj_res:
            for key, val in final_dict.items():
    #            line1 = str(key) + ' ' + str(val) + '\n'
                line = str('{:22s}'.format(str(key))) + ' ' + str('{:4d}'.format(len(val))) + '\n'
    #            f_obj_3.write(line1)
                f_obj_res.write(line)
        
    elif res_or_ep == 'ep':
        # --------------------------------------------------------------------------------- #
        # Write an output 'ep_glysum.txt'
        # Protein residue, number of glycan atoms 
        # --------------------------------------------------------------------------------- #
        with open('ep_glysum.txt', 'w') as f_obj_ep:        
            total_glycan_list = []
            for key, val in final_dict.items():
                total_glycan_list.extend(val)
            total_glycan = list(OrderedDict.fromkeys(total_glycan_list))
            # f_obj_5.write(str(input_epitope) + ' ' + str(len(total_glycan)))
            f_obj_ep.write(str(len(total_glycan)))


# --------------------------------------------------------------------------------- #
#                                      MAIN                                         #
# --------------------------------------------------------------------------------- #


def main():
    print("Starting")

    # Extract arguments from command line
    args = process_cmd_args()
    sur_cutoff = args["sur_cutoff"]
    origin_pdb = args["origin_pdb"]
    dist_cutoff = args["dist_cutoff"]
    input_glycan = args["input_glycan"]
    res_or_ep = args["res_or_ep"]
    input_epitope = args["input_epitope"]
    freesasa_path = args["freesasa_path"]
    nproc = args["nproc"]

    """
    # Debug inputs
    sur_cutoff = 30
    origin_pdb = "input/origin_frame_100_renumber-HXB2.pdb"
    dist_cutoff = 26
    input_glycan = ["BMA", "AMA", "BGLN"]
    res_or_ep = "res"
    input_epitope = None 
    freesasa_path = "/home/reveizm2/glyco_2/dependencies/freesasa"
    nproc = 28
    """

    # python3 glyco_1.py -pdb input/origin_frame_100_renumber-HXB2.pdb -cutoff 26 -module res -glycan BMA,AMA,BGL -freesasa /home/reveizm2/glyco_2/dependencies/freesasa -nproc 28

    # Generate dictionaries with input files
    ori_pro_dict, pro_dict, gly_dict, res_dict = preprocess_dicts(origin_pdb, sur_cutoff, input_glycan, input_epitope,
                                                                  res_or_ep, freesasa_path)

    # Launch workers
    all_results = launch_workers(gly_dict, res_dict, pro_dict, nproc, dist_cutoff, ori_pro_dict)

    # Retrieve and aggregate data from all workers
    count_dict, num_dict, final_dict = aggregate_worker_res(all_results)

    # Write results
    write_outputs(res_or_ep, count_dict, num_dict, final_dict)


if __name__ == "__main__":
    start = time.time()

    main()

    end = time.time()
    print('TIME SPENT: ', str(datetime.timedelta(seconds=end-start)))
