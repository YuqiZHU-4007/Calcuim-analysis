#!/bin/bash

dirs=(
    "./fear_conditioning/huc/20190514/fish3/rec25106988/raw tiff"
    "./fear_conditioning/huc/20190514/fish3/rec25106988/rawtiff"
    "./fear_conditioning/huc/20190514/raw tiff"
    "./fear_conditioning/RAW tiff"
    "./fear_conditioning/huc/20190520/fish/rec25106988/RAWTIFF"
    "./fear_conditioning/huc/20190520/fish/rec25106920/"
    "./fear_conditioning/huc/20190514/fish/rec25106988/raw tiff/rawtiff"
)
for dir in "${dirs[@]}"; do
    mkdir -p "$dir"
    for n in {1..20}; do
        foldernum=$( printf 'z%02d' $n )
        mkdir "$dir/$foldernum"
        for m in {1..5}; do
            filenum=$( printf 'p%02d.tiff' $m )
            touch "$dir/$foldernum/$filenum"
        done
    done
done

dirs=(
    "./fear_conditioning/contain"
    );
for dir in $dirs; do
    mkdir -p "$dir"
    for m in {1..38000}; do
        filenum=$( printf 'p%05d.tiff' $m )
        touch "$dir/$filenum"
    done
done