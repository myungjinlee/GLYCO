#!/usr/bin/env python
# coding: utf-8

#-------------------------------------------------------------#
# GLYCO
#
# 07/27/2021
# Developed by Myungjin Lee, Ph.D, Mateo Reveiz
# Please contact myungjin.lee@nih.gov for any question or bug.
#-------------------------------------------------------------#

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
import glob  
import fileinput

# For multiprocessing
import multiprocessing 
from multiprocessing.pool import Pool
from concurrent.futures import ProcessPoolExecutor as Pool_out


# ---------------------------------------------------------------------------------#
#                                 Helper methods                                   #
# ---------------------------------------------------------------------------------#

# Generate a surface residue dictionary 
def outermost(inpdb, infloat, freesasa_path, inprobe, out_folder):
    print('PROGRESS: Excluding buried residues of the input protein', flush=True)
     
    outrsa = out_folder + "/" + Path(inpdb).name.replace(".pdb", "") + "_outer.rsa"
    
    os.system(freesasa_path + ' ' + str(inpdb) + " -p " + str(inprobe) + " --rsa > " + str(outrsa))
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
            if atype == 'RES' and resname in ('ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU',
                                              'GLY', 'HSD', 'HIS', 'HIE', 'HID', 'ILE', 'LEU',
                                              'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
                                              'TYR', 'VAL'):
                surf_key = (resname, chain, resnum)
                surf_dict[surf_key] = float(surf_area)
                if float(surf_area) >= infloat:
                    keys.append([resname, chain, resnum])
    keys = numpy.array(keys)
    with open(inpdb, "r") as f:
        lines = f.readlines()
        for line in open(inpdb):
            if line[17:20] in ('ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HSD',
                               'HIS', 'HIE', 'HID', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                               'SER', 'THR', 'TRP', 'TYR', 'VAL'):
                atomnum = line[6:11].strip()
                atomtype = line[12:16].strip()
                resname = line[17:20].strip()
                chain = line[21:22].strip()
                resid = line[22:28].strip()
            
                for k in keys:                    
                    if k[0] == resname and k[1] == chain and k[2] == resid: 
                        key = (atomnum, atomtype, resname, chain, resid)
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        if atomtype.find('H') == -1:
                            outdict[key] = [x, y, z]
    
    return outdict, surf_dict


# Convert an subset residue file to a dictionary
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


# Generate a protein dictionary that has key and coordinate as values
# and, generate another protein dictionary that has x as key and the rest as values
# the second dictionary is for checking crossing with protein region
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


# Generate a glycan dictionary that has glycan keys and values as glycan coordinate
def gen_gly_dict(inpdb, input_glycan):
    print('PROGRESS: Generating a glycan dictionary') 
    outdict = {}
    for line in open(inpdb):
        if line[0:4] == 'ATOM' or line[0:4] == 'HETA':
            if line[17:21].strip() in input_glycan.split():
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


# Extract residues from a dictionary and generate an subset residue dictionary
def ep_to_dict(inarr, indict):
    outdict = {}
    for ep in inarr:
        for key, val in indict.items():
            if ep[0] == key[2] and ep[1] == key[4]:
                outdict[key] = val
    return outdict


# Calculate distances between two coordinates
def cal_dist(inlist1, inlist2):
    inarr1 = numpy.array(inlist1)
    inarr2 = numpy.array(inlist2)
    outfloat = numpy.linalg.norm(inarr1-inarr2)
    return outfloat


def check_protein(inar_diff, inarr1, dist, intuple, pro_dict, ori_pro_dict):

    leng = 1 / dist
    
    for t in numpy.arange(0.0, 1.0, leng): 
        point = inarr1 + t * (inar_diff) 
        
        # These are the coordinates of the parametrized line
        x = int(point[0])
        y = int(point[1])
        z = int(point[2]) 

        # To speed up things, only look at atoms in slices of x-axis
        filter_x = pro_dict.get(x, []) + pro_dict.get(x - 1, []) + pro_dict.get(x + 1, [])  
        
        for elem in filter_x:
            # We used the integer to filter by x slices, now retrieve the real float coordinates
            elem_float_x, elem_float_y, elem_float_z = ori_pro_dict[elem[2]]
        
            if abs(y - elem_float_y) <= 1 and abs(z - elem_float_z) <= 1 and abs(x - elem_float_x) <= 1:  
                    if not (elem[2][2] == intuple[2] and elem[2][3] == intuple[3] and elem[2][4] == intuple[4]):
                        return False

    return True


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


    parser = argparse.ArgumentParser(usage= '\n'+ '\n' +
             'Single PDB------------------------------------------------------------------------------------------------------------------------------------------' + '\n'
             'module "all": Glycan coverage of the entire protein surface' + '\n'
             'python3 glyco.py -pdb frame.pdb -cutoff 26 -sur_cutoff 30 -probe 1.5 -module all -glycan BMA,AMA,BGL -freesasa /home/bin/freesasa -num_proc_in 22 -out_folder results' + '\n' + '\n'
             'module "sub": Glycan coverage of the subset of protein residues' + '\n'
             'python3 glyco.py -pdb frame.pdb -cutoff 26 -module sub -glycan BMA,AMA,BGL -residue residuelist.txt -num_proc_in 22 -num_parallel 1 -out_folder results' + '\n' 
             '----------------------------------------------------------------------------------------------------------------------------------------------------' + '\n' + '\n'
             'Multiple PDBs---------------------------------------------------------------------------------------------------------------------------------------' + '\n'
             'module "all": ' + '\n'
             'python3 glyco.py -in_folder input -cutoff 26 -sur_cutoff 30 -probe 1.5 -module all -glycan BMA,AMA,BGL -freesasa /home/bin/freesasa -num_proc_in 22 -num_parallel 2 -out_folder results -average'  + '\n'  + '\n'
             'module "sub": ' + '\n'
             'python3 glyco.py -in_folder input -cutoff 26  -module ep -glycan BMA,AMA,BGL -residue subset_residue.txt -num_proc_in 22 -num_parallel 2 -out_folder results -average' + '\n'
             '----------------------------------------------------------------------------------------------------------------------------------------------------' + '\n' + '\n')
    parser.add_argument('-pdb', type=str, help='Enter your name of pdb. ex) input.pdb')
    parser.add_argument('-cutoff', type=float, help='Enter your distance cutoff for glycans in Angstrom. ex) 30', required=True)
    parser.add_argument('-glycan', type=str, help='Enter your name of glycans with space for separators. ex) BMA,AMA,BGLN', required=True)
    parser.add_argument('-module', type=str, help='Enter your module name, either all or sub. ex) all', required=True)
    parser.add_argument('-residue', type=str, help='Enter your protein residue list. ex) residuelist.txt')
    parser.add_argument('-freesasa', type=str, help='Enter your path of Freesasa executable, ex) /home/lee/freesasa')
    parser.add_argument('-num_proc_in', type=int, help='Enter the number of workers to use per frame. ex) 9', default=multiprocessing.cpu_count())
    parser.add_argument('-num_parallel', type=int, help='Enter the number of frames to do in parallel. ex) 5', default=1)
    parser.add_argument('-probe', type=float, help='Enter probe radius of FreeSASA in Angstrom. ex) 1.5', default=1.4)
    parser.add_argument('-sur_cutoff', type=float, help='Enter surface area cutoff of FreeSASA in Angstrom^2. ex) 40', default=30)
    parser.add_argument('-average', action='store_true', help='Add if using input folder with multiple pdbs of the same protein.')
    parser.add_argument('-in_folder', type=str, help='Input folder where input pdb files are located. ex) input')
    parser.add_argument('-out_folder', type=str, help='Output folder where results will be saved. ex) output', required=True)

    args = parser.parse_args()

    if args.module == 'all' and args.freesasa == None:
        print('Please add freesasa path')
        sys.exit(1) 
    elif args.module == 'sub' and args.residue == None:
        print('Please add residue list')
        sys.exit(1)
        
        
    if args.pdb is not None:
        if args.in_folder is not None:
            print("Please use either -pdb OR -in_folder. Both cannot be used at the same time.")
            sys.exit(1)
        else:
            pass
    else:
        if args.in_folder is not None:
            input_files = list(glob.glob(os.path.join(args.in_folder, "*.pdb")))
            
            if args.num_parallel > len(input_files):
                print("The number of parallel jobs (-num_parallel {}) exceeds the number of pdb files found in the input folder ({}).".format(args.num_parallel, len(input_files)))
                sys.exit(1) 
            
        else:
            print("Please use either -pdb OR -in_folder.")
            sys.exit(1)
        
         
    if args.average:
        if args.pdb is not None:
            print("Please do not use -pdb when -average is true. All input pdbs should be in the /input folder.")
            sys.exit(1)
        if args.in_folder is None:
            print("-average requires -in_folder and out_folder to be provided.")
            
        
    if Path(args.out_folder).exists():
        print("Output folder '{}' already exists. Please prodive another folder.".format(args.out_folder))
        sys.exit(1)
        
    if args.num_proc_in * args.num_parallel > multiprocessing.cpu_count():
        print("The number of CPUs requested ({}) exceeds the available CPUs detected ({}).".format(args.num_proc_in * args.num_parallel, multiprocessing.cpu_count()))
        sys.exit(1)
        

    args_["sur_cutoff"] = args.sur_cutoff
    args_["origin_pdb"] = args.pdb
    args_["dist_cutoff"] = int(args.cutoff)
    args_["input_glycan"] = args.glycan.replace(',', ' ')
    args_["all_or_sub"] = args.module
    args_["input_sub"] = args.residue
    args_["freesasa_path"] = args.freesasa 
    args_["nproc_in"] = args.num_proc_in
    args_["nproc_out"] = args.num_parallel
    args_["probe"] = args.probe
    args_["average"] = args.average
    args_["in_folder"] = args.in_folder
    args_["out_folder"] = args.out_folder
    args = args_
    
    return args

# Generate dictionaries with input files
def preprocess_dicts(origin_pdb, sur_cutoff, input_glycan, input_sub, all_or_sub, freesasa_path, inprobe, out_folder):
    
    ori_pro_dict, pro_dict = gen_pro_dict(origin_pdb)
    gly_dict = gen_gly_dict(origin_pdb, input_glycan)

    res_dict = {}
    if all_or_sub == 'all':
        res_dict = outermost(origin_pdb, sur_cutoff, freesasa_path, inprobe, out_folder)[0]

    elif all_or_sub == 'sub':
        res_dict = input_to_dict(input_sub, ori_pro_dict)
        
    return ori_pro_dict, pro_dict, gly_dict, res_dict


def launch_workers(gly_dict, res_dict, pro_dict, nproc, dist_cutoff, ori_pro_dict):
    
    glycan_data_list = list(gly_dict.items())
    protein_data_list = list(res_dict.items())

    # Set the longest of {protein, glycan} lists to be the outter loop when generating pairs, and partition outter into
    # sublists to send to multiple CPUs with full inner list
    if len(glycan_data_list) > len(protein_data_list):
        outter_data = glycan_data_list
        inner_data = protein_data_list
        glycan_out = True
    else:
        outter_data = protein_data_list 
        inner_data = glycan_data_list
        glycan_out = False

    # Split calculations into multiple bundles
    outter_data = numpy.array(outter_data, dtype=object)
    bundle_num = nproc * 2
    bundle_outter_data = numpy.array_split(outter_data, bundle_num, axis=0) 

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
    print("Starting worker {}".format(job_id), flush=True)
    
    count_dict = {}
    num_dict = {}
    tmp_dict = defaultdict(list)
    for res_dict_item in res_dict:
        key1, val1 = res_dict_item
        res_key = (key1[2], key1[3], key1[4])

        count = 0 
        gly_list = []
        for gly_dict_item in gly_dict:
            key2, val2 = gly_dict_item
            
            inarr1 = numpy.array(val1)
            inarr2 = numpy.array(val2)
            inar_diff = inarr2 - inarr1
            
            dist = numpy.linalg.norm(inar_diff)
            
            if dist < dist_cutoff:

                flag = check_protein(inar_diff, inarr1, dist, key1, pro_dict, ori_pro_dict)
                
                if flag:
                    count += 1
                    gly_list.append(key2[0])
                else:
                    pass
                    
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


def write_outputs(all_or_sub, count_dict, num_dict, final_dict, out_folder, prefix, origin_pdb):
    
    if all_or_sub == 'all':
        
        print(out_folder + "/" + prefix + '_all_glycount.txt')
        with open(out_folder + "/" + prefix + '_all_glycount.txt', 'w') as f_obj_res:
            for key, val in final_dict.items():
                line = str('{:22s}'.format(str(key))) + ' ' + str('{:4d}'.format(len(val))) + '\n'
                f_obj_res.write(line)

        print(origin_pdb)
        print(out_folder + "/{}_bfactor.pdb".format(prefix))
          
        os.system("cp " + origin_pdb + ' ' + out_folder + "/" + prefix + '_bfactor.pdb')
        for line in fileinput.input(out_folder + "/" + prefix + '_bfactor.pdb', inplace=1):
            line = line.strip()
            if line.startswith('ATOM'):
                origin_line = line
                key = (str(line[17:21].strip()),str(line[21:22].strip()),str(line[22:26].strip()))
                if key in final_dict.keys():
                    count = len(final_dict[key])
                    val = str('{:>6.2f}'.format(count))
                    new_line = origin_line[:60] + val + origin_line[66:]
                    print(new_line)
                else:
                    val = str('{:>6.2f}'.format(0))
                    new_line = origin_line[:60] + val + origin_line[66:]
                    print(new_line)
            else:
                print(line)

    elif all_or_sub == 'sub':
        with open(out_folder + "/" + prefix + '_sub_glysum.txt', 'w') as f_obj_ep:        
            total_glycan_list = []
            for key, val in final_dict.items():
                total_glycan_list.extend(val)
            total_glycan = list(OrderedDict.fromkeys(total_glycan_list))
            f_obj_ep.write(str(len(total_glycan)))


def process_request(args, origin_pdb):

    sur_cutoff = args["sur_cutoff"]
    dist_cutoff = args["dist_cutoff"]
    input_glycan = args["input_glycan"]
    all_or_sub = args["all_or_sub"]
    input_sub = args["input_sub"]
    freesasa_path = args["freesasa_path"]
    nproc_in = args["nproc_in"]
    probe = args["probe"]
    average = args["average"]
    out_folder = args["out_folder"]

    # Generate dictionaries with input files
    ori_pro_dict, pro_dict, gly_dict, res_dict = preprocess_dicts(origin_pdb, sur_cutoff, input_glycan, input_sub,
                                                                  all_or_sub, freesasa_path, probe, out_folder)
    # Launch workers
    all_results = launch_workers(gly_dict, res_dict, pro_dict, nproc_in, dist_cutoff, ori_pro_dict)

    # Retrieve and aggregate data from all workers
    count_dict, num_dict, final_dict = aggregate_worker_res(all_results)

    # Write results
    prefix = Path(origin_pdb).name.replace(".pdb", "")
    write_outputs(all_or_sub, count_dict, num_dict, final_dict, out_folder, prefix, origin_pdb)
   
   
def average_frames_res(out_folder):
    
    total_dict = {}
    output_files = list(glob.glob(os.path.join(out_folder, "*all_glycount.txt")))
    
    for output_file in output_files:
        for line in open(output_file):
            line = line.split()
            resname = " ".join(re.findall("[a-zA-Z]+", line[0]))
            chain = " ".join(re.findall("[a-zA-Z]+", line[1]))
            resid = line[2][1:-2]
            resid = str(resid)
            count = line[3]
            key = (resname, chain, resid)
            if key not in total_dict.keys():
                total_dict[key]=[]
            total_dict[key].append(float(count.strip()))
    
    with open(out_folder + '/ave_all_glycount.txt','w') as file_obj:
        for key in sorted(total_dict, key=lambda x: (x[1])):
            if len(total_dict[key]) > 0:
                ave_count = sum(total_dict[key])/len(total_dict[key])
                file_obj.write(str("{:20s}".format(str(key))) + str("{:>9.2f}".format(ave_count))+'\n')

def average_frames_ep(out_folder):
    
    output_files = list(glob.glob(os.path.join(out_folder, "*sub_glysum.txt")))
    glysum = 0
    total_average = 0
    
    for output_file in output_files:
        for line in open(output_file):
            line = line.strip()
            glysum += float(line)
    total_average = glysum / len(output_files)
    
    with open(out_folder + '/ave_sub_gly.txt', 'w') as file_obj:
        file_obj.write(str("{:>9.2f}".format(total_average)) + '\n')



# --------------------------------------------------------------------------------- #
#                                      MAIN                                         #
# --------------------------------------------------------------------------------- #

 

def main():
    print("Starting")

    # Extract arguments from command line
    args = process_cmd_args()
    
    # Create output folder
    out_folder = args["out_folder"]
    os.mkdir(out_folder) 
    
    if args["in_folder"] is None:
        # Single pdb to calculate
        process_request(args, args["origin_pdb"])
    else:
        # Multiple pdbs found in the in_folder
    
        # Get files from input folder
        in_folder = args["in_folder"]
        input_files = list(glob.glob(os.path.join(in_folder, "*.pdb")))
    
        # Create pool
        nproc_out = args["nproc_out"]
        args["origin_pdb"] = None
        OUTER_POOL = Pool_out(nproc_out)
        
        # Submit all jobs
        for input_file in input_files:
            OUTER_POOL.submit(process_request, args, input_file)
    
        # Shutdown and wait for children to finish
        OUTER_POOL.shutdown(wait=True)
        
        # Average frames if requested
        if args["average"]:
            print("Averaging frames...")
            
            if args["all_or_sub"] == "all":
                average_frames_res(out_folder)
            elif args["all_or_sub"] == "sub":
                average_frames_ep(out_folder)
            else:
                print("Unknown module!")



if __name__ == "__main__":
    start = time.time()

    main()

    end = time.time()
    print('TIME SPENT: ', str(datetime.timedelta(seconds=end-start)))

