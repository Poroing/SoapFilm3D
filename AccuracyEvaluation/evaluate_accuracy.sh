#!/bin/bash

#${deltas} stores the different regularization coefficients for which the
#relative error of the Barnes-Hut method and fast multipole method relative
#error will be computed. Each regularization coefficient must be separated
#by a space (This is to be able to iterate over them with bash for loop).
deltas="0.01 0.001 0.0001 0.00001 0"

#${regularization_coefficients} stores the different regularization coeffcients
#for which the relative error of the fast multipole method during the
#simulation will be computation. It must have the same format as ${deltas}
regularization_coefficients="0.5,0.25,0.125"

LATEX_FILENAME="delta_evaluation.tex"

USAGE="Usage: $0 [-e | -E] [-t] [-s | -S] [-p | -P] [-h] [-c]"

HELP=$(cat << EOF
${USAGE}

   Compute the relative error of some methods that compute the Biot-Savart law and produces
   boxplots for comparison. The methods evaluated is the Barnes-Hut method (implemented by
   modifying the libigl code for fast winding number) and the Fast Multipole Method (implemented
   through the fmmtl library. The relative error is computed by comparing the method with the
   naive method for the computation of the Biot-Savart law. The evaluation is made with different
   regularization parameters for the Biot-Savart law and on multiple randomly generated set of
   sources.

   This also computes and produces boxplot of the relative error of the Fast Multipole Method
   compared to the naive method for each computation of the Biot-Savart law in the
   simulation of 330 time steps of a bubble lattice made of 20 bubbles with umbrella laplacian
   smoothing, smoothing coefficient 0.01 and target edge length 0.3

   The plots produced can be changed by modifying this script or the file ${LATEX_FILENAME}. See
   inside of those file for documentation.

ARGUMENTS

   -e (defaut) Run the different methods to obtain the data necessary to make the plots.
      Those files have the form delta_evaluation_\${delta}.csv where \${delta} is the
      regularization coefficient used.

   -E (not -e). This should be used if the data is already present"

   -t Rerun a method if it fails for any reason.

   -s (default) Run the simulation to obtain the data necessary to make the plots. Those files
      have the form delta_evaluation_\${delta}.csv where \${delta} is the regularization
      coefficient used

   -S (not -s). This should be used if the data is already data present

   -p (default) Produces the boxplots as pdf files from the produced data. Two boxplots are
      produced per methods, one with the ouliers and one without.

      The pdf files names have the
      form delta_evaluation_\${method}_\${outliers}.pdf, where \${method} is "BarnesHut" for
      the relative error of the "Barnes-Hut" method, "Fmmtl" for the Fast Multipole Method
      and Error for the simulation, \${outliers} is "outliers" for the boxplots with the
      outliers and "no_outliers" for those without.

      The boxplots are made by compiling the LaTeX file ${LATEX_FILENAME} (See the documentation
      in the specified file).

   -P (not -p). This should be used to debug this script. This allows quicker debuging of the
      the data production.

   -c Remove every file produced by the script (implies -E -S -P).

   -h Show this help, then quit.
EOF
)


produce_boxplots() {
    pdflatex "\def\method{$1} \def\deltaparameters{$2} \input{$LATEX_FILENAME}";
    mv delta_evaluation.pdf delta_evaluation_${1}_no_outliers.pdf;
    pdflatex "\def\showoutliers{} \def\method{$1} \def\deltaparameters{$2} \input{$LATEX_FILENAME}";
    mv delta_evaluation.pdf delta_evaluation_${1}_outliers.pdf;
}

## Method Parameters

fmmtl_expansion_order=5
fmmtl_theta=0.5
fmmtl_minimum_cell_size=128
winding_expansion_order=2
winding_beta=2
number_sources=1000

run_evaluation_flag=1
retry_on_fail_flag=
run_simulation_flag=1
produce_pdf_flag=1
help_flag=
clean_flag=
while getopts eEtsSpPhc options; do
    case $options in
        E) run_evaluation_flag=;;
        e) run_evaluation_flag=1;;
        t) retry_on_fail_flag=1;;
        s) run_simulation_flag=1;;
        S) run_simulation_flag=;;
        p) produce_pdf_flag=1;;
        P) produce_pdf_flag=;;
        h) help_flag=1;;
        c) clean_flag=1;;
    esac
done

if [[ -n $help_flag ]]; then
    echo -e "${HELP}"
    exit
fi

if [[ -n $clean_flag ]]; then
    rm -f *.csv
    rm -f *.aux
    rm -f *.pdf
    rm -f *.log
    rm -rf MergedBubbleLatticeDelta
    exit
fi

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

if [[ -n $produce_pdf_flag ]]; then
    methods="BarnesHut Fmmtl"
    deltas=$(echo $deltas | tr " " ",")
    echo $deltas
    for method in $methods; do
        produce_boxplots $method "$deltas"
    done
fi

if [[ -n $run_simulation_flag ]]; then
    ../run_simulations.py \
        -o scene=mergedbubblelattice \
        -o mesh-size-m=20 \
        -o remeshing-resolution=0.3 \
        -o regularization-coefficient=$regularization_coefficients \
        -o fast-summation=fmmtl \
        -o smoothing-coef=0.01 \
        -o smoothing-type=umbrella \
        -o fmmtl-minimum-cell-size=$fmmtl_minimum_cell_size \
        -o fmmtl-expansion-order=$fmmtl_expansion_order \
        -o print-biot-savart-error=1 #Enable the computation and output of the relative error. \
        -t 3.3 \
        MergedBubbleLatticeDelta/
fi

../process_results.py \
    -o regularization-coefficient=$regularization_coefficients \
    MergedBubbleLatticeDelta \
    csv \
    --mode execution_time \
    --filename-stem biot-savart

for regularization_coefficient in $(echo $regularization_coefficients | tr "," " "); do
    mv \
        MergedBubbleLatticeDelta/Regularization-coefficient$regularization_coefficient/biot-savart.csv \
        delta_evaluation_${regularization_coefficient}.csv
done

if [[ -n $produce_pdf_flag ]]; then
    produce_boxplots Error "$regularization_coefficients"
fi

