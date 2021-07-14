#!/bin/bash
 scr_name=${1%%.m}
 Experiments=()
 Experiments+=(historical rcp85 rcp45)
 nExperiments=${#Experiments[@]}
 Tables=(Omon interp)
 nTables=${#Tables[@]}

 for ((ii=0; ii<nExperiments; ii++)); do
   Experiment=${Experiments[$ii]}
 for ((jj=0; jj<nTables; jj++)); do
   Table=${Tables[$jj]}
   scr_no=$((ii*nTables+jj))
   tmp_scr=${scr_name}_${scr_no}
   echo $scr_name > /tmp/log.$tmp_scr
   awk '/^Experiments =/ {printf("Experiments = [\"%s\"];\n",Experiment); next} /^Tables =/ {printf("Tables = [\"%s\"];\n",Table); next} //' Experiment=$Experiment Table=$Table  ${scr_name}.m > ${tmp_scr}.m
   runmatlab.sh ${tmp_scr} --clear &
 done
 done

 while [ 1 ]; do
   sleep 10
   procs=`ls -1 ${scr_name}_[0-9]* 2>/dev/null | wc -l`
   if [[ $procs -le 0 ]]; then break; fi
   echo "$procs matlab scripts are running.."
 done
 echo "All complete!"

