#!/bin/bash -l
nreaz=$1
ncore=$2
pflotran_exe=$3

name="1dthermal"
dir="pflotran/"

cd $dir

srun -n $ncore $pflotran_exe -pflotranin $name".in" -stochastic -num_realizations $nreaz -num_groups $nreaz -screen_output off ;
#srun -n 50 ~/pflotran-cori -pflotranin 1dthermal.in -stochastic -num_realizations 50 -num_groups 50 -screen_output off


#for ireaz in $(seq $(($ncore+1)) $nreaz)
#do
#    while [ $(ls -ltr *.txt | wc -l) -lt $(($ireaz-$ncore)) ]
#    do
#	sleep 0.01
#    done
#
#    (cd $ireaz ;
#	$pflotran_exe -pflotranin $mc_name".in"  -screen_output off ;
#	cd ../ ;
#	echo $ireaz > $ireaz.txt) &
#
#done

#while  [ $(ls -ltr *.txt | wc -l) -lt $nreaz ]
#do
#    sleep 1
#done

wait

cd ../
