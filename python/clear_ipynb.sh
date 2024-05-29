#!/bin/bash

in_dir=$1
echo "Clearing output cells for all Jupyter Notebooks in $in_dir"

if [ ! -d $in_dir ]
then
  echo "$in_dir doesn't exist. Abort!"
  exit
fi

if [ ! -d $in_dir/python ]
then
  echo "$in_dir/python doesn't exist. Abort!"
  exit
fi

flist=$(ls $in_dir/python/*.ipynb)
for fsc in $flist 
do
  echo "  Processing $fsc"
  jupyter nbconvert --clear-output --inplace $fsc
done

echo "Done"
