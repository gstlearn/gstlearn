#!/bin/bash

in_dir=$1
out_dir=$2
runner=$3
echo "Processing $in_dir to $out_dir"

if [ -d $out_dir ]
then
  if [ -z "$(ls -A $out_dir)" ]
  then
    echo "$out_dir exists and is empty."
  else
    echo "$out_dir exists and full of stuff."
  fi
  while true; do
    read -p "Do you want to continue (Y/N) ?" yn

    case $yn in 
	y ) break;;
	Y ) break;;
	n ) exit;;
	N ) exit;;
	* ) ;;
    esac
  done
fi

if [ ! -d $in_dir/python ]
then
  echo "$in_dir/r doesn't exist. Abort!"
  exit
fi

mkdir -p $out_dir

flist=$(ls $in_dir/python/*.ipynb)
for fsc in $flist 
do
  echo "  Processing $fsc"
  python3 $runner $fsc $out_dir html
done

echo "Done"
