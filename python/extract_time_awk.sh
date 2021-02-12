#/bin/bash


#DIR=/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY/20201117/D4_500m_CTRL_NP1024/exp/3660515_cycle_20190824150000
DIR=/data_ballantine02/miyoshi-t/honda/SCALE-LETKF/AIP_D4_VERIFY/20201117/D4_500m_CTRL_NP1024/exp/3662073_cycle_20190824150000

IFILE=$DIR/job.o
OFILE=$DIR/DA.txt

rm -f $OFILE
awk '{if($1=="[Info:DA]")print $7 }' $IFILE > $OFILE

OFILE=$DIR/FCST.txt
awk '{if($1=="[Info:fcst]" && $2=="End")print $8 }' $IFILE > $OFILE
