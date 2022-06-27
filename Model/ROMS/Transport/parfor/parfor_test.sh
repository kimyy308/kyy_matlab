#!/bin/bash

src_name=${1%%.m}
start_i=1
end_i=2
num_i=$((end_i-start_i+1))
echo $src_name $start_i $end_i $num_i
printf "" > tmp/pid_file
for ((i_i=1; i_i<=num_i; i_i++)); do
  i=$((i_i-1+start_i))
  echo $i > tmp/log_${i_i}
  awk '/for testnameind2=/ {printf("for testnameind2 = %d\n",int(i)); next} / for c_testnameind2 =/ {printf("for testnameind2 = 1:0\n"); next} //' i=$i  ${src_name}.m > ${src_name}_${i_i}.m
  ./runmatlab.sh $src_name _$i_i &
done
#while [ 1 ]; do
#  sleep 10
#  procs=`ls ${src_name}_* 2>/dev/null|wc -l`
#  if [[ $procs -le 0 ]]; then break; fi
#  echo "$procs matlab scripts are running.."
#done
#echo "Climate script starts.."
#awk '/for i =/ {printf("for i = 1:0\n"); next} //' ${src_name}.m > ${src_name}_0.m
#./runmatlab.sh $src_name _0
echo "All complete!"
rm tmp/pid_file

