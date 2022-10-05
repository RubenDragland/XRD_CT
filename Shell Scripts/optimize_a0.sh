#! /usr/bin/bash

while getopts m:i:d: flag
do
    case "${flag}" in
        m) mode=${OPTARG};;
        i) interpolation=${OPTARG};;
        f) div=${OPTARG};;
    esac
done

echo "Running optimisation of a0 Spherical Harmonics coefficient"
matlab -nojvm -nodisplay -nodesktop -r "./SASTT_optimization_SH_a0(mode, interpolation, div)" ; exit # Make file so that it comes with arguments