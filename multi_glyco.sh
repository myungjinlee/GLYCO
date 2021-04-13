#!/bin/bash

# ------------------------------------------------------- #
# Myungjin Lee, Ph.D. National Institutes of Health
# Released date: 04/07/2021
# Description: multi_glyco.sh is a bash script to
# ------------------------------------------------------- #

#ls /data/SBIS/leem25/glycan_density/hiv_man9_wt/epitope > 'ep_list.txt'
#ls $2/epitope > 'ep_list.txt'

# cmd> bash multi_ep_sub.sh 10 /data/SBIS/leem25/glycan_density/glyco/multiframes 3 5
# path: where all inputs are located (GLYCO, template folder, input pdb folder)

let arg=$#

while [ $# -gt 0 ]; 
    do
           case "$1" in
                -cutoff)
                    shift
                    cutoff=$1
                    shift
                    ;;
                -module)
                    shift
                    module=$1
                    shift
                    ;;
                -glycan) 
                    shift
                    glycan=$1 
                    shift 
                    ;;
                -frame_start)
                    shift
                    frame_start=$1
                    shift
                    ;;
                -frame_end)
                    shift
                    frame_end=$1
                    shift
                    ;;
                -frame_gap)
                    shift
                    frame_gap=$1
                    shift
                    ;;
                -path)
                    shift
                    path=$1
                    shift
                    ;;
                -freesasa)
                    shift
                    freesasa=$1
                    shift
                    ;;
                -num_proc)
                    shift
                    freesasa=$1
                    shift
                    ;;
                -epitope)
                    shift
                    epitope=$1
                    shift
                    ;;
                *)
                   echo "$1 is not a recognized flag! please see the following example"
                   echo "bash multi_res_run.sh -cutoff 10 -frame_start 1 -frame_end 3 -path /home/leem/multiframes -glycan BMA,AMA,BGL"
                   exit 1;
                   ;;
          esac
    done

# ------------------------------------------------------- #
# Removing slash (/) at the working directory path
# ------------------------------------------------------- #

if [ ${path: -1} == '/' ]; then
path=${path:0:-1}
fi

# ------------------------------------------------------- #
# -------------------------MAIN-------------------------- #
# ------------------------------------------------------- #

# ------------------------------------------------------- #
# When input module is "res"
# ------------------------------------------------------- #
if [ $module == 'res' ]; then
    if [ $arg -ne 16 ]; then
        echo "---------------------------------------ERROR: Wrong number of arguments, please see the following example-------------------------------------------------------"
        echo "bash multi_glyco.sh -cutoff 10 -frame_start 1 -frame_end 10 -module res -path /home/leem/multiframes -freesasa /home/leem/freesasa -glycan BMA,AMA,BGL -num_proc 28"
        exit 1
    elif [ ! $cutoff ] ||[ ! $frame_start ] || [ ! $frame_end ] || [ ! $module ] || [ ! $path ] || [ ! $freesasa ] || [ ! $glycan ]; then
        echo "---------------------------------------ERROR: Missing arguments, please see the following example---------------------------------------------------------------"
        echo "bash multi_glyco.sh -cutoff 10 -frame_start 1 -frame_end 10 -module res -path /home/leem/multiframes -freesasa /home/leem/freesasa -glycan BMA,AMA,BGL -num_proc 28"
    exit 1
    fi

    echo module: $module, GLYCO will calculate glycan atoms of all residues on the protein
    mkdir $cutoff
    cd $cutoff
        mkdir res
        cd res
            cp $path/template/res/*sh .
            sed -i "s/CUTOFF/$cutoff/g" gen_sub.sh
            sed -i "s/CUTOFF/$cutoff/g" run.sh
            sed -i "s/GLYCAN/$glycan/g" run.sh
            sed -i "s#FREESASA_PATH#$freesasa#g" run.sh
            sed -i "s#WORKING_DIR#$path#g" run.sh
            sed -i "s/NUM_PROC/$num_proc/g" run.sh
            for indx in `eval echo {$frame_start..$frame_end}`; do
                mkdir frame_$indx
                cd frame_$indx
                    cp ../run.sh .
                    cp $path/input/frame_${indx}.pdb . 
                    sed -i "s/INDEX/$indx/g" run.sh
                cd ../
            done
            wait
            bash gen_sub.sh $path $frame_start $frame_end > swarm.sub
            # Please edit below line as to your system
            swarm -g 56 --partition quick -t ${num_proc} -f swarm.sub --time 03:59:00 --job-name=${cutoff}_res
        cd ../
    cd ../

# ------------------------------------------------------- #
# When input module is "ep"
# ------------------------------------------------------- #
elif [[ $module == 'ep' ]]; then
    
    if [ $arg -ne 18 ]
    then
        echo "---------------------------------------ERROR: Wrong number of arguments, please see the following example---------------------------------------------"
        echo "bash multi_glyco.sh -cutoff 10 -frame_start 1 -frame_end 10 -frame_gap 3 -module ep -path /home/leem/multiframes -glycan BMA,AMA,BGL -epitope /home/leem/multiframe/epitope/ep.txt -num_proc 28"
        exit 1
    elif [ ! $cutoff ] || [ ! $frame_start ] || [ ! $frame_end ] || [ ! $frame_gap ] || [ ! $module ] || [ ! $path ] || [ ! $glycan ] || [ ! $epitope ]; then
        echo "---------------------------------------ERROR: Missing arguments, please see the following example-----------------------------------------------------"
        echo "bash multi_glyco.sh -cutoff 10 -frame_start 1 -frame_end 10 -frame_gap 3 -module ep -path /home/leem/multiframes -glycan BMA,AMA,BGL -epitope /home/leem/multiframe/epitope/ep.txt -num_proc 28"
        exit 1
    fi

    mkdir $cutoff
    cd $cutoff
        mkdir ep
        cd ep/
    #while IFS= read -r line
    #    do
    #    echo "$line"
    #    mkdir ${line}
    #    cd ${line}
            cp $path/template/ep/* .
            sed -i "s#WORKING_DIR#$path#g" run.sh 
            sed -i "s/CUTOFF/$cutoff/g" run.sh
            sed -i "s/GLYCAN/$glycan/g" run.sh
            sed -i "s#EPITOPE#$epitope#g" run.sh
            sed -i "s/NUM_PROC/$num_proc/g" run.sh
   #        sed -i "s/ZZZ/${line}/g" run.sh
            for indx in `eval echo {$frame_start..$frame_end}`; do
                mkdir frame_$indx
                cd frame_$indx
                    cp ../run.sh .
                    cp $path/input/frame_${indx}.pdb . 
                    sed -i "s/INDEX/$indx/g" run.sh
                    cd ../
            done
            for (( index=$frame_start; index<=$frame_end; index+=$frame_gap )); do
                if [[ $index != $frame_start ]];then
                    let f=$((index-$frame_gap))
                    let l=$index
                    for ind in `eval echo {$f..$l}`; do cat $PWD/frame_$ind/run.sh ; done > $f.sh
                    cp temp.sub swarm.sub
                    sed -i "s/BBB/$f/g" swarm.sub 
                    swarm -g 32 --partition quick -t ${num_proc} -f swarm.sub --time 3:59:00 --job-name=$f
                fi
            done
            if [ $((index-$frame_gap)) != $frame_end ]; then
                let f=$((index-$frame_gap))
                let l=$index
                for ind in `eval echo {$f..$frame_end}`; do cat $PWD/frame_$ind/run.sh ; done > $f.sh
                cp temp.sub swarm.sub
                sed -i "s/BBB/$f/g" swarm.sub
                swarm -g 32 --partition quick -t ${num_proc} -f swarm.sub --time 3:59:00 --job-name=$f
            fi
        cd ../
        # done < /data/SBIS/leem25/glycan_density/hiv_man9_wt/ep_list.txt
    cd ../../
fi
