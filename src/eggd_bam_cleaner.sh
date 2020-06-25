#!/bin/bash
# eggd_bam_cleaner 1.0.0
set -e -x -o pipefail

main() {

    echo "Value of bam: '$bam'"
    echo "Value of multi_bool: '$multi_bool'"
    echo "Value of MAPQ threshold: '$mapq"

    # The following line(s) use the dx command-line tool to download your file
    # inputs to the local file system using variable names for the filenames. To
    # recover the original filenames, you can use the output of "dx describe
    # "$variable" --name".

    mkdir files && cd files

    dx download "$bam" -o input_bam

    # sort bam by name, req. for fixmate
    echo "sorting bam by name"
    samtools sort -n -o name_sort.bam input_bam

    # fixing bam with fixmate
    echo "fixing mate coords and insert size fields"
    samtools fixmate -m name_sort.bam fixmate.bam

    # sorting bam by pos, req. for markdup
    echo "sorting bam by coordinates"
    samtools sort -o pos_sort.bam fixmate.input_bam

    # removing duplicates with markdup
    echo "removing duplicates"
    samtools markdup -r pos_sort.bam rm_dup.bam

    if "$multi_bool"; then
        # multi bool true => remove multi-mapped reads
        # default threshold 20 unless specified
        echo "removing multi mapped reads"
        output_bam="${bam}_rm_dup_multimap.bam"
        
        samtools view -q $mapq -b rm_dup.bam > "$output_bam"
        samtools index $output_bam
    else
        # multi_bool false => not removing mult-mapped reads
        output_bam="${bam}_rm_dup.bam"
        samtools index $output_bam


    bam=$(dx upload output_bam --brief)
    index=$(dx upload "{output_bam}.bai" --brief)

    # The following line(s) use the utility dx-jobutil-add-output to format and
    # add output variables to your job's output as appropriate for the output
    # class.  Run "dx-jobutil-add-output -h" for more information on what it
    # does.

    dx-jobutil-add-output bam "$bam" --class=file
    dx-jobutil-add-output index "$index" --class=file
}
