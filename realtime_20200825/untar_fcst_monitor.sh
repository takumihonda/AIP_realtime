#/bin/bash


TOP=/data_ballantine02/miyoshi-t/amemiya/SCALE-LETKF-rt-archive/result/ope/d4_500m
#flist=`ls $TOP/2020083*/dafcst/nc*.gz`
flist=`ls $TOP/202009*/dafcst/grads_ref3d.tar`
flist=`ls $TOP/2020090[3-6]*/dafcst/grads_ref3d.tar`
flist=`ls $TOP/202008[2,3]*/dafcst/grads_ref3d.tar`

wk=`pwd`

fcst_tmp=$wk/fcst
mkdir -p $fcst_tmp

echo $wk

for fname in $flist
do
  echo $fname
  tar -xvf $fname -C $fcst_tmp

done

exit


