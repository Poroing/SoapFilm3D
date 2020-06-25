#!/bin/bash

fmmtl_expansion_order=5
fmmtl_theta=0.5
fmmtl_minimum_cell_size=1
winding_expansion_order=2
winding_beta=2
number_sources=1000
deltas="0.1 0.01 0.001 0.0001 0.00001"

run_evaluation_flag=1
retry_on_fail_flag=
while getopts nrt options; do
    case $options in
        n) run_evaluation_flag=;;
        r) run_evaluation_flag=1;;
        t) retry_on_fail_flag=1;;
    esac
done

echo $run_evaluation_flag
if [[ -n $run_evaluation_flag ]]; then
    for delta in $deltas; do
        echo "../build/Apps/SoapFilm3D/Tests \
            $delta \
            $fmmtl_expansion_order $fmmtl_theta $fmmtl_minimum_cell_size \
            $winding_expansion_order $winding_beta \
            $number_sources"

        while \
                ! ../build/Apps/SoapFilm3D/Tests \
                    $delta \
                    $fmmtl_expansion_order $fmmtl_theta $fmmtl_minimum_cell_size \
                    $winding_expansion_order $winding_beta \
                    $number_sources \
                > delta_evaluation_$delta.csv \
            && \
            [[ -n $retry_on_fail_flag ]]
        do
            echo "Retrying"
        done
    done
fi

methods="BarnesHut Fmmtl"
deltas=$(echo $deltas | tr " " ",")
echo $deltas
for method in $methods; do
    pdflatex "\def\method{$method} \def\deltaparameters{$deltas} \input{delta_evaluation.tex}";
    mv delta_evaluation.pdf delta_evaluation_${method}_no_outliers.pdf;
    pdflatex "\def\showoutliers{} \def\method{$method} \def\deltaparameters{$deltas} \input{delta_evaluation.tex}";
    mv delta_evaluation.pdf delta_evaluation_${method}_outliers.pdf;
done
