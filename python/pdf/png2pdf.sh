#!/bin/bash

mkdir -p png

cd png

for fn in `ls *.png`; do
  fn=$( basename $fn '.png')
  echo $fn
  convert ${fn}.png ${fn}.pdf
#  mv ${fn}.png png/
done

cd -
