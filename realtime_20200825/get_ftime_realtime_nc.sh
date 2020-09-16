#/bin/bash


TOP=/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_20200825/as-of-pm16-0827/dafcst_nc
TOP=/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_SAFE/realtime_20200825/dafcst_nc



ofile=$TOP/../result.txt

echo $ofile

rm -f $ofile
touch $ofile

cd $TOP

#nclist=`ls *.nc`
nclist=`ls 2020090[6,7]*.nc`


for nc in $nclist
do
  list=`ls -l --time-style="+%Y-%m-%d %H:%M:%S" $nc | awk '{print $6, $7,$8}'`
  echo $list
  echo $list >> $ofile
done

cd -


echo $ofile

