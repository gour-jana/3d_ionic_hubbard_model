#!/bin/bash 

datafile=input.dat
pathfile=path.txt
sub=job.sge

mkdir -p 3d_ionic_DOS
cd 3d_ionic_DOS

for U1 in  1.0 2.0 3.0 
do
mkdir -p U_$U1
cd U_$U1

for tp in 0.0 0.4 0.6
do 
mkdir -p tp_$tp
cd tp_$tp

for seed in 87959 #98352 79869 69857 96387
do 
mkdir -p seed_$seed
cd seed_$seed

for strnth in 0.0 0.5 1.0
do
mkdir -p strnth_$strnth
cd strnth_$strnth


cp ../../../../3d_tca_ttp_ionic_Hubbard_model.f90 .
cp ../../../../mtfort90.f .


#!******************************************************************************
echo   "8                         !d "             > $datafile
echo   "4                          !dc "            >> $datafile
echo   "-1.0d0                     !t1     "        >> $datafile
echo   $tp"0d0                      !t2     "        >> $datafile
echo   "29                         !Tgrid_max "      >> $datafile
echo   "2000                       !MCSW "          >> $datafile
echo   "10                          !intrvl   "      >> $datafile
echo   $U1"d0                      !U1 "             >> $datafile  
echo   "1.0d0                      !filling "       >> $datafile
echo   "0.02                       !gama_m "          >> $datafile 
echo   $strnth"0d0                       !strnth "        >> $datafile 
echo   "$seed                      !seed "          >> $datafile
#!******************************************************************************
echo "../../../../../main_data/U_$U1/tp_$tp/seed_$seed/strnth_$strnth"    >> $pathfile
#!******************************************************************************

gfortran -o run.x  mtfort90.f  3d_tca_ttp_ionic_Hubbard_model.f90  -llapack -lblas

#********************************************************************************
echo "#PBS  -N   3DI_"$U1"_tp_"$tp                            > $sub
echo "#PBS -l nodes=1:ppn=1"                        >>$sub
echo "#PBS -j oe "                                  >> $sub   
echo "#PBS -o out.log"                                  >> $sub 
echo "#PBS -e err.log"                                  >> $sub           
echo "cd"    \$PBS_O_WORKDIR                         >> $sub
echo "date"                                         >> $sub
echo "./run.x <input.dat> main_out "                    >> $sub
echo "date"                                         >> $sub
#******************************************************************************
chmod 777 job.sge
qsub  job.sge
cd ..
done

cd ..
done 

cd ..
done

cd ..
done
