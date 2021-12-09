#!/bin/bash       
                  
#SBATCH -p lgpu   
#SBATCH --gres=gpu         
                  
source /home/eduard/miniconda3/bin/activate openmm
#source /home/eduard/miniconda3/condabin/activate cuda
job_name=`cat myjob.def | sed -n '1p'  | awk '{print $2}'`  
temp=`cat myjob.def | sed -n '2p'  | awk '{print $2}'`  

last=25          
echo $job_name
echo $temp 
python3 npt.py $1   #${job_name} ${temp} 
ret_val=$?                                                                                  
           
if [ $ret_val -eq 0 -a $1 -lt $last ] ; then
    let next=$1+1          
    sbatch -J ${job_name::2}_$temp_$next -o out/npt_$next.out run_npt.sh $next
elif [ ! $ret_val -eq 0 ] ; then
    echo "Error in run $1 .." >> out/error.log
fi                         
                  
                  




