#/bin/bash


ofile=/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_test20200807/honda.txt


echo $ofile

rm -f $ofile
touch $ofile


nclist=`ls *.nc`
for nc in $nclist
do
  list=`ls -l --time-style="+%Y-%m-%d %H:%M:%S" $nc | awk '{print $6, $7,$8}'`
  echo $list
  echo $list >> $ofile
done


