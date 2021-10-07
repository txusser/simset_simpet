#!/bin/bash
# runs the programs mat, calc_ad, and calc_comp_ad
# for each of the materials  given as arguments.
# 
# sample_usage:  run_mat BGO tungsten
# This command line would take the input files
# BGO.dat and tungsten.dat and create the outputs
#

for filename in $*
do
        echo generating "$filename" data files - may take up to an hour
        nice mat < "$filename".dat
        nice calc_ad < "$filename".dat > "$filename".ad
        nice calc_comp_ad < "$filename".dat > "$filename".comp_ad
        echo
done

echo done done done

