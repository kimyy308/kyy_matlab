#!/bin/bash
src_name=$1
func=$2

matlab_cmd="matlab -nosplash -nodesktop -nodisplay -nojvm -r \"${src_name}${func}; exit;\""
echo $$ >> tmp/pid_file
date >> tmp/log${func}
#echo $matlab_cmd >> tmp/log${func}
eval $matlab_cmd >> tmp/log${func}
rm -f ${src_name}${func}.m
exit
