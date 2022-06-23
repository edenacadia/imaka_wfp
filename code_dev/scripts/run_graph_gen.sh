#!/bin/bash
# bash script for reading in a run directory
# Eden McEwen 
# July 24, 2020

# take in a dir
dirname=$1
dirtemp="${dirname}tmp/"

#check to see if input dir exists
if [ ! -d "$dirname" ]; then
    echo "Input dir doesn't exist: ${dirname}"
    exit -1
fi

#making temp file directory
if [ ! -d "${dirtemp}" ]; then
    echo "Creating tmp dir: ${dirtemp}"
    mkdir $dirtemp
fi

#split into tmp files for running
# if the file is too large, split rows
for filename in ${dirname}RUN*.txt; do
    echo "creating temps for ${filename}"
    file1="${filename%%.*}"
    name="${file1##*/}"
    tmpprefix="${dirtemp}${file1##*/}_tmp"
    split -l 2 --numeric-suffixes $filename $tmpprefix --additional-suffix=.txt
    line=$(head -n 1 "$filename")
    sed -i "s;^;${line}\n;" $tmpprefix*.txt
done


# background runs on all files
for filename in ${dirtemp}RUN*.txt; do
    echo "processing ${filename}"
    file1="${filename%%.*}"
    screen="${file1##*/}"
    echo "  Starting screen:  log to ${screen}"
    #screen -X -S $screen quit
    screen -S $screen -d -m 
    screen -S $screen -X stuff $"python3 wp_pipeline.py -d $filename\n"
done