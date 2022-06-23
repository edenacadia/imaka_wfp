#!/bin/bash
filename=$1
outdir=$2
outfile="wspd250_"$filename

# generateurls from a run#.txt file
#newfile="wspd250_"+$filename
#while read line; do
# reading each line
#chrlen=${#line}
#if [$chrlen]
#then
#  $prefix
#echo $line
#done < $filename

# wget chmod u+x scriptfiles
prefix='http://mkwc.ifa.hawaii.edu/forecast/mko/models/gfs/'
wspd_suffix=".wspd.txt"
wdir_suffix=".wdir.txt"


while read line; do
echo $line
wget -N -P "$outdir"wind_250/ "$prefix""$line"00/gfs."$line"00.wspd.txt
wget -N -P "$outdir"wind_250/ "$prefix""$line"00/gfs."$line"00.wdir.txt
wget -N -P "$outdir"wind_250/ "$prefix""$line"06/gfs."$line"06.wspd.txt
wget -N -P "$outdir"wind_250/ "$prefix""$line"06/gfs."$line"06.wdir.txt
wget -N -P "$outdir"wind_250/ "$prefix""$line"12/gfs."$line"12.wspd.txt
wget -N -P "$outdir"wind_250/ "$prefix""$line"12/gfs."$line"12.wdir.txt
wget -N -P "$outdir"wind_250/ "$prefix""$line"18/gfs."$line"18.wspd.txt
wget -N -P "$outdir"wind_250/ "$prefix""$line"18/gfs."$line"18.wdir.txt
done < $filename

