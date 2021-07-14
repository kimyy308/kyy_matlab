#!/bin/bash

src_name=${1%%.m}
start_yr=1985
end_yr=1999
num_yr=$((end_yr-start_yr+1))
echo $src_name $start_yr $end_yr $num_yr
printf "" > tmp/pid_file
for ((i_yr=1; i_yr<=num_yr; i_yr++)); do
  yr=$((i_yr-1+start_yr))
  echo $yr > tmp/log_${i_yr}
  awk '/for year =/ {printf("for year = %d\n",int(yr)); next} //' yr=$yr  ${src_name}.m > ${src_name}_${i_yr}.m
  ./runmatlab.sh $src_name _$i_yr &
  sleep 1
done

while [ 1 ]; do
  sleep 10
  procs=`ls ${src_name}_* 2>/dev/null|wc -l`
  if [[ $procs -le 0 ]]; then break; fi
  echo "$procs matlab scripts are running.."
done
echo "All complete!"
rm tmp/pid_file

