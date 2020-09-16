#/bin/bash


TOP=/data_ballantine02/miyoshi-t/amemiya/SCALE-LETKF-rt-archive/result/ope/d4_500m
#flist=`ls $TOP/2020083*/dafcst/nc*.gz`
flist=`ls $TOP/202009*/dafcst/nc*.gz`
flist=`ls $TOP/202008[2,3]*/dafcst/nc*.gz`

#echo $flist
#exit
#flist=$TOP


ofile=/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_20200825/result.txt

tmp=/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_20200825/tmp

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
    echo $list >> $ofile
    echo $list
exit
  done

  cd -

#  list=`ls -l --time-style="+%Y-%m-%d %H:%M:%S" *.nc | awk '{print $6, $7,$8}'`
#  echo $list

#  list=`zcat ${fname} | tar tv | awk '{print $5, $6}'` 


done

echo $ofile

