#!/bin/bash
# eggd_bam_cleaner 1.0.0
set -e -x -o pipefail

main() {

    # input variables
    echo "Value of bam: '$bam'"
    echo "Value of multi_bool: '$multi_bool'"
    echo "Value of MAPQ threshold: '$mapq"

    # remove _markdup if present in name
    if [[ $bam_name =~ "_markdup" ]]
    then
        echo "markdup in name"
        bam_name="${bam_name/_markdup/}"
    else
        echo "not in name"
    fi

    # unpack and build samtools
    cd /packages
    tar -jxvf samtools-1.10.tar.bz2
    cd samtools-1.10
    ./configure --prefix=/packages
    make
    make install
    export PATH=/packages/bin:$PATH

    echo "finished building samtools, beginning analysis"

    # working dir to download bams to
    mkdir ~/files && cd ~/files

    dx download "$bam" -o input_bam

    # sort bam by name, req. for fixmate
    echo "sorting bam by name"
    samtools sort -n -o name_sort.bam input_bam

    # fixing bam with fixmate
    echo "fixing mate coords and insert size fields"
    samtools fixmate -m name_sort.bam fixmate.bam

    # sorting bam by pos, req. for markdup
    echo "sorting bam by coordinates"
    samtools sort -o pos_sort.bam fixmate.bam

    # removing duplicates with markdup
    echo "removing duplicates"
    samtools markdup -r pos_sort.bam rm_dup.bam

    if [[ "$multi_bool" ]]; then
        # multi bool true => remove multi-mapped reads
        # default threshold ($mapq) 20 unless specified 
        echo "removing multi mapped reads"
        output_bam="${bam_name/.bam/_rm_dup_rm_multimap.bam}"
        
        samtools view -q $mapq -b rm_dup.bam > "$output_bam"
        samtools index $output_bam
    else
        # multi_bool false => not removing mult-mapped reads
        output_bam="${bam_name/.bam/_rm_dup.bam}"
        samtools index $output_bam
    fi

    bam=$(dx upload $output_bam --brief)
    index=$(dx upload "${output_bam}.bai" --brief)

    # The following line(s) use the utility dx-jobutil-add-output to format and
    # add output variables to your job's output as appropriate for the output
    # class.  Run "dx-jobutil-add-output -h" for more information on what it
    # does.

    dx-jobutil-add-output bam "$bam" --class=file
    dx-jobutil-add-output index "$index" --class=file
}
