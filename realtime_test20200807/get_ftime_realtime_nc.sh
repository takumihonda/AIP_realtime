#/bin/bash


TOP=/data_ballantine02/miyoshi-t/amemiya/SCALE-LETKF-rt-archive/result/ope/d4_500m
flist=`ls $TOP/2*/dafcst/nc*.gz`

#TOP=/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_test20200807/dafcst_nc
#flist=$TOP


ofile=/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_test20200807/amemiya.txt

tmp=/home/honda/SCALE-LETKF/AIP_REALTIME/python/tmp

echo $ofile

rm -f $ofile
touch $ofile

for fname in $flist
do

  rm -rf $tmp
  mkdir -p $tmp
  cp $fname $tmp/
  cd $tmp

  tar zxvf $tmp/*.gz 


  nclist=`ls *.nc`
  for nc in $nclist
  do
    list=`ls -l --time-style="+%Y-%m-%d %H:%M:%S" $nc | awk '{print $6, $7,$8}'`
    echo $list
    echo $list >> $ofile
  done

  cd -

#  list=`ls -l --time-style="+%Y-%m-%d %H:%M:%S" *.nc | awk '{print $6, $7,$8}'`
#  echo $list

#  list=`zcat ${fname} | tar tv | awk '{print $5, $6}'` 


done

echo $ofile

