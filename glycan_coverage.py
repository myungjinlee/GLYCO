#!/home/leem25/software/miniconda3/bin/python

import sys
import os
import re
import time
import numpy
import datetime
from collections import defaultdict
from collections import OrderedDict

start = time.time()

f_obj_1 = open('all_gly.txt', 'w')
f_obj_2 = open('all_count.txt', 'w')
f_obj_3 = open('res_gly.txt', 'w')
f_obj_4 = open('res_count.txt', 'w')
f_obj_5 = open('epitope_total.txt', 'w')


def outermost(inpdb, infloat):
    outrsa = 'outer.rsa'
    os.system("freesasa " + str(inpdb) + " -p 1.4 --rsa > " + str(outrsa))
    keys = []
    outdict = {}
    surf_dict = {}
    out_surf_acc_dict = {}
    for line in open(outrsa):
        line=line.split()
        if line[0] == 'RES' and line[1] in ('ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HSD','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'):
            surf_key = (line[1],line[2],line[3])
            surf_dict[surf_key] = float(line[4])
            if float(line[4]) >= infloat:
                keys.append([line[1],line[2],line[3]])
                    # 'RES     ALA      A       1    108.06'
    keys=numpy.array(keys)
    for line in open(inpdb):
        if line[17:20] in ('ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HSD','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'):
            for k in keys:
                resname = line[17:20].strip()
                chain = line[21:22].strip()
                resid = line[22:28].strip()
                if k[0] == resname and k[2] == resid: #k[1] == chain and k[2] == resid:
                    atomnum = line[6:11].strip()
                    atomtype = line[12:16].strip() 
                    key = (atomnum, atomtype, resname, chain, resid)
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    if atomtype.find('H') == -1:
                        outdict[key]=[x,y,z]
    return outdict, surf_dict

def read_epit(instr):
    outlist = []
    inlist=instr.split(', ')
    r=re.compile("([a-zA-Z]+)([0-9]+)")
    inlist=[r.match(string).groups() for string in inlist]
    for ep in inlist:
        if ep[0] == 'V':
            index = inlist.index(ep)
            inlist.remove(ep)
            ep = ep[:0]+ep[1:]
            ep = ep[:0]+('VAL',)+ep[0:]
            inlist.insert(index, ep)
        if ep[0] == 'Y':
            index = inlist.index(ep)
            inlist.remove(ep)
            ep = ep[:0]+ep[1:]
            ep = ep[:0]+('TYR',)+ep[0:]
            inlist.insert(index, ep)
        if ep[0] == 'W':
            index = inlist.index(ep)
            inlist.remove(ep)
            ep = ep[:0]+ep[1:]
            ep = ep[:0]+('TRP',)+ep[0:]
            inlist.insert(index, ep)
        if ep[0] == 'T':
            index = inlist.index(ep)
            inlist.remove(ep)
            ep = ep[:0]+ep[1:]
            ep = ep[:0]+('THR',)+ep[0:]
            inlist.insert(index, ep)
        if ep[0] == 'S':
            index = inlist.index(ep)
            inlist.remove(ep)
            ep = ep[:0]+ep[1:]
            ep = ep[:0]+('SER',)+ep[0:]
            inlist.insert(index, ep)
        if ep[0] == 'P':
            index = inlist.index(ep)
            inlist.remove(ep)
            ep = ep[:0]+ep[1:]
            ep = ep[:0]+('PRO',)+ep[0:]
            inlist.insert(index, ep)
        if ep[0] == 'F':
            index = inlist.index(ep)
            inlist.remove(ep)
            ep = ep[:0]+ep[1:]
            ep = ep[:0]+('PHE',)+ep[0:]
            inlist.insert(index, ep)
        if ep[0] == 'M':
            index = inlist.index(ep)
            inlist.remove(ep)
            ep = ep[:0]+ep[1:]
            ep = ep[:0]+('MET',)+ep[0:]
            inlist.insert(index, ep)
        if ep[0] == 'K':
            index = inlist.index(ep)
            inlist.remove(ep)
            ep = ep[:0]+ep[1:]
            ep = ep[:0]+('LYS',)+ep[0:]
            inlist.insert(index, ep)
        if ep[0] == 'L':
            index = inlist.index(ep)
            inlist.remove(ep)
            ep = ep[:0]+ep[1:]
            ep = ep[:0]+('LEU',)+ep[0:]
            inlist.insert(index, ep)
        if ep[0] == 'I':
            index = inlist.index(ep)
            inlist.remove(ep)
            ep = ep[:0]+ep[1:]
            ep = ep[:0]+('ILE',)+ep[0:]
            inlist.insert(index, ep)
        if ep[0] == 'H':
            index = inlist.index(ep)
            inlist.remove(ep)
            ep = ep[:0]+ep[1:]
            ep = ep[:0]+('HSD',)+ep[0:]
            inlist.insert(index, ep)
        if ep[0] == 'G':
            index = inlist.index(ep)
            inlist.remove(ep)
            ep = ep[:0]+ep[1:]
            ep = ep[:0]+('GLY',)+ep[0:]
            inlist.insert(index, ep)
        if ep[0] == 'Q':
            index = inlist.index(ep)
            inlist.remove(ep)
            ep = ep[:0]+ep[1:]
            ep = ep[:0]+('GLN',)+ep[0:]
            inlist.insert(index, ep)
        if ep[0] == 'E':
            index = inlist.index(ep)
            inlist.remove(ep)
            ep = ep[:0]+ep[1:]
            ep = ep[:0]+('GLU',)+ep[0:]
            inlist.insert(index, ep)
        if ep[0] == 'C':
            index = inlist.index(ep)
            inlist.remove(ep)
            ep = ep[:0]+ep[1:]
            ep = ep[:0]+('CYS',)+ep[0:]
            inlist.insert(index, ep)
        if ep[0] == 'D':
            index = inlist.index(ep)
            inlist.remove(ep)
            ep = ep[:0]+ep[1:]
            ep = ep[:0]+('ASP',)+ep[0:]
            inlist.insert(index, ep)
        if ep[0] == 'N':
            index = inlist.index(ep)
            inlist.remove(ep)
            ep = ep[:0]+ep[1:]
            ep = ep[:0]+('ASN',)+ep[0:]
            inlist.insert(index, ep)
        if ep[0] == 'R':
            index = inlist.index(ep)
            inlist.remove(ep)
            ep = ep[:0]+ep[1:]
            ep = ep[:0]+('ARG',)+ep[0:]
            inlist.insert(index, ep)
        if ep[0] == 'A':
            index = inlist.index(ep)
            inlist.remove(ep)
            ep = ep[:0]+ep[1:]
            ep = ep[:0]+('ALA',)+ep[0:]
            inlist.insert(index, ep)
    outarr = numpy.array(inlist)
    return outarr


def gen_pro_dict(inpdb):
    outdict1 = {}
    outdict2 = defaultdict(list) 
    for line in open(inpdb):
        if line[0:4] == 'ATOM':
            if line[17:20] in ('ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HSD','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'):
                resname = line[17:21].strip()
                chain = line[21].strip()
                resid = line[22:28].strip()
                atomnum = line[6:11].strip()
                atomtype = line[12:16].strip()
                key = (atomnum, atomtype, resname, chain, resid)
                x = int(float(line[30:38].strip()))
                y = int(float(line[38:46].strip()))
                z = int(float(line[46:54].strip()))
                if atomtype.find('H') == -1:
                    outdict1[key]=[x,y,z]
                    outdict2[x].append([y,z,key]) #check crossing
    return outdict1, outdict2

def gen_gly_dict(inpdb):
    outdict = {}
    for line in open(inpdb):
        if line[0:4] == 'ATOM' and line[17:21].strip() in 'BGLC, BMAM, AMAM, BGAL, AFUC, ANE5':
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
                    outdict[key]=[x,y,z]
    return outdict

def ep_to_dict(inarr, indict):
    outdict={}
    for ep in inarr:
        for key, val in indict.items():
            if ep[0]==key[2] and ep[1]==key[4]:
                outdict[key]=val
    print(outdict)
    return outdict

def get_ca(indict):
    outdict={}
    for key, val in indict.items():
        if key[1]=='CA':
            outdict[key]=val
    return outdict

def cal_dist(inlist1,inlist2):
    inarr1=numpy.array(inlist1)
    inarr2=numpy.array(inlist2)
    outstr = numpy.linalg.norm(inarr1-inarr2)
    return outstr

def vector(inlist1, inlist2):
    outlist = []
    inarr1=numpy.array(inlist1)
    inarr2=numpy.array(inlist2)
    leng = numpy.linalg.norm(inarr2-inarr1)
    leng = int(round(leng))
    for i in range(1, int(round(leng))): 
        t = i/float(leng)
        point = inarr1 + t*(inarr2-inarr1) 
        outlist.append(point)
    return outlist

def check_protein(inlist, indict, intuple):
    flag = True
    tmparr = numpy.array(inlist)
    tmparr = numpy.around(tmparr, decimals=0)
    for tmp in tmparr:
        x = int(tmp[0])
        y = int(tmp[1])
        z = int(tmp[2])
        filter_x = pro_dict.get(x)
        if filter_x:
            for elem in filter_x:
                if y == elem[0]:
                    if (elem[1]-1) <= z and z <= (elem[1]+1):
                            if elem[2][2] == intuple[2] and elem[2][3] == intuple[3] and elem[2][4] == intuple[4]:
                                flag = True
                            else:
                                flag = False
    return flag

sur_cutoff = 30
origin_pdb = sys.argv[1]
dist_cutoff = int(sys.argv[2])
all_or_ep = sys.argv[3]

ori_pro_dict = gen_pro_dict(origin_pdb)[0]
pro_dict = gen_pro_dict(origin_pdb)[1]
gly_dict = gen_gly_dict(origin_pdb)
surf_dict = outermost(origin_pdb, sur_cutoff)[1]

res_dict = {}
if all_or_ep == 'all':
    res_dict = outermost(origin_pdb, sur_cutoff)[0]
elif all_or_ep == 'ep':
    input_epitope = sys.argv[4]
    file_obj2=open(input_epitope, 'r')
    epitope = ''
    for line in file_obj2:
        epitope = line.strip()
    ep_site = read_epit(epitope)
    res_dict = ep_to_dict(ep_site, ori_pro_dict)

count_dict = {}
num_dict = {}
tmp_dict = {}
total_list = []
for key1, val1 in res_dict.items():
    res_key = (key1[2],key1[3],key1[4])
    if res_key not in tmp_dict.keys():
        tmp_dict[res_key]=[]
    count = 0 
    gly_list = []
    for key2, val2 in gly_dict.items():
        dist = cal_dist(val1, val2)
        if dist < dist_cutoff:
            flag = True
            line = vector(val1, val2)
            flag = check_protein(line, pro_dict, key1)
            if flag == False:
               print(flag, ' PROTEIN', key1, val1)
               print(' gly', key2, val1)
            if flag == True:
                count += 1
                gly_list.append(key2[0])
    tmp_dict[res_key].extend(gly_list)
    total_list.extend(gly_list)
    count_dict[key1] = gly_list
    num_dict[key1] = count

final_dict = {}
for key, val in tmp_dict.items():
    area = float(surf_dict[key])
    if area == 0.00:
        area = 0.000000000000001
    final_dict[key]=list(OrderedDict.fromkeys(val))

for key, val in count_dict.items():
    line = str(key) + ' ' + str(val) + '\n'
    f_obj_1.write(line)

for key, val in num_dict.items():
    line = str(key) + ' ' + str(val) + '\n' 
    f_obj_2.write(line)

for key, val in final_dict.items():
    area = float(surf_dict[key])
    if area == 0.00:
        area = 0.000000000000001
    line1 = str(key) + ' ' + str(val) + '\n'
    line2 = str(key) + ' '+ str('{:4d}'.format(len(val))) +'  ' + str('{:8.3f}'.format(len(val)/area)) + '\n'
    f_obj_3.write(line1)
    f_obj_4.write(line2)

if all_or_ep == 'ep':
    total_glycan_list = []
    for key, val in final_dict.items():
        total_glycan_list.extend(val)
    total_glycan = list(OrderedDict.fromkeys(total_glycan_list))
    f_obj_5.write(str(input_epitope) + ' ' + str(len(total_glycan)))

end = time.time()
print('TIME SPENT: ', str(datetime.timedelta(seconds=end-start)))
print ('SUCCESSFULLY DONE!')
