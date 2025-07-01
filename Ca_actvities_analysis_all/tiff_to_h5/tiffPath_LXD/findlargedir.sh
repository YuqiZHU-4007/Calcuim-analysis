#!/bin/bash

# example: bash findlargedir.sh 1000

pf=$(find . -type f -regex '.*\.tif$' | sed -r 's|/[^/]+$||' |uniq -c ) # count file numbers and get folder paths
num=$1 #CHANGE FISRT PARAMETER the amont tiff files should be greater than or equal to $num
ORIIFS=$IFS
IFS=$'\n'
if [[ -e largetiff.csv ]] ; then # remove existed file paths
    rm largetiff.csv
fi
for s in $pf; do
    s="${s#"${s%%[![:space:]]*}"}" # trim header spaces
    n=$(echo "$s" | sed -r 's/[^0-9]*([0-9]+).*$/\1/') # file numbers
    if [[ $n -ge $num ]] ; then # file numbers greater than or equal to given number
        echo $s | sed -r 's/^[0-9]+\ //1'>> largetiff.csv # put a new folder path into largetiff.csv
    fi
done
IFS=$ORIIFS