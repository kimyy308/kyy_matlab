#!/bin/bash

src_name=${1%%.m}
start_yr=1995
end_yr=2004
num_yr=$((end_yr-start_yr+1))
echo $src_name $start_yr $end_yr $num_yr
printf "" > tmp/pid_file
for ((i_yr=1; i_yr<=num_yr; i_yr++)); do
  yr=$((i_yr-1+start_yr))
  echo $yr > tmp/log_${i_yr}
#  awk '/for year =/ {printf("for year = %d\n",int(yr)); next} / for c_year =/ {printf("for year = 1:0\n"); next} //' yr=$yr  ${src_name}.m > ${src_name}_${i_yr}.m
  awk '/for year =/ {printf("for year = %d\n",int(yr)); next} // ' yr=$yr  ${src_name}.m > ${src_name}_${i_yr}.m
#  ./runmatlab.sh $src_name _$i_yr &
done

#while [ 1 ]; do
#  sleep 10
#  procs=`ls ${src_name}_* 2>/dev/null|wc -l`
#  if [[ $procs -le 0 ]]; then break; fi
#  echo "$procs matlab scripts are running.."
#done
#echo "Climate script starts.."
#awk '/for yr =/ {printf("for yr = 1:0\n"); next} //' ${src_name}.m > ${src_name}_0.m
#./runmatlab.sh $src_name _0
echo "All complete!"
mkdir -p tmp
rm tmp/pid_file

