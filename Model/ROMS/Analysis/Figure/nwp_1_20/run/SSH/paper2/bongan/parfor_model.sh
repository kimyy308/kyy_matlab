#!/bin/bash
 scr_name=${1%%.m}
 models=()
 models+=(IPSL-CM5A-LR IPSL-CM5A-MR MPI-ESM-LR NorESM1-M)
 nmodels=${#models[@]}

 for ((ii=0; ii<nmodels; ii++)); do
   model=${models[$ii]}
   tmp_scr=${scr_name}_$ii
   echo $scr_name > /tmp/log.$tmp_scr
   awk '/^Models =/ {printf("Models = [\"%s\"];\n",model); next} //' model=$model  ${scr_name}.m > ${tmp_scr}.m
   runmatlab.sh ${tmp_scr} --clear &
 done

 while [ 1 ]; do
   sleep 10
   procs=`ls -1 ${scr_name}_[0-9]* 2>/dev/null | wc -l`
   if [[ $procs -le 0 ]]; then break; fi
   echo "$procs matlab scripts are running.."
 done
 echo "All complete!"

