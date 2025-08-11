#Snakemake workflow for cleaning raw metagenomes 

import os

#Preparing files 

configfile: "config/config.yml"

INPUT = config["input"]
with open(INPUT) as f:
    SAMPLES = [line.strip() for line in f if line.strip()]

READS = config["reads"]
ASSEMBLY = config["assembly"]
ALIGNMENT = config["alignment"]
METABAT2 = config["metabat2"]
CONCOCT = config["concoct"]
MAXBIN2 = config["maxbin2"]
SEMIBIN2 = config["semibin2"]
CONTIGS_10K = config["10k-contigs"]
CONCOCT_FINAL = config["concoct-final"]
CHECKM = config["checkM"]
OUTPUT = config["output"]

###### Protocol ######

rule all:
    input:
        expand("{sample}-all_done.txt", sample = SAMPLES)
        
rule minimap2_align:
    input:
        assembly = ASSEMBLY + "{sample}.fa",
        fwd = READS + "{sample}_1.fastq.gz",
        rev = READS + "{sample}_2.fastq.gz"
    output:
        bam = ALIGNMENT + "{sample}-sorted.bam", 
        flag = "{sample}-align_done.txt"
    conda:
        "envs/metabinning.yaml"
    shell:
        """
        minimap2 -ax sr -t 6 {input.assembly} {input.fwd} {input.rev} | \
        samtools view -bS -@ 6 | \
        samtools sort -@ 6 -o {output.bam}
        touch {output.flag}
        """

rule cut_concoct_10K:
    input:
        ASSEMBLY + "{sample}.fa"
    output:
        assembly_10k = CONTIGS_10K + "{sample}_10K.fasta",
        bed = ALIGNMENT + "{sample}_10K.bed",
        flag = "{sample}_10K-cut_done.txt"
    conda:
        "envs/metabinning.yaml"
    shell:
        """
        cut_up_fasta.py \
        {input} \
        -c 10000 \
        -o 0 -m \
        --merge_last \
        -b {output.bed} > {output.assembly_10k}
        touch {output.flag}
        """

rule depth_calc_concoct:
    input:
        cut = "{sample}_10K-cut_done.txt",
        align = "{sample}-align_done.txt"
    params:
        bed = ALIGNMENT + "{sample}_10K.bed",
        bam = ALIGNMENT + "{sample}-sorted.bam"
    output:
        depth = ALIGNMENT + "{sample}-coverage_table.tsv",
        flag = "{sample}-depth_10K-done.txt"
    conda:
        "envs/metabinning.yaml"
    shell:
        """
        scripts/concoct_coverage_table.py {params.bed} {params.bam} > {output.depth}
        touch {output.flag}
        """

rule index_and_depth_calc:
    input:
        ALIGNMENT + "{sample}-sorted.bam"
    output:
        index = ALIGNMENT + "{sample}-sorted.bam.bai",
        depth = ALIGNMENT + "{sample}-depth.txt", 
        flag = "{sample}-index-depth-done.txt"
    conda:
        "envs/metabinning.yaml"
    shell:
        """
        samtools index {input}
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input}
        touch {output.flag}
        """
        
rule maxbin2:
    input:
        assembly = ASSEMBLY + "{sample}.fa",
        depth = ALIGNMENT + "{sample}-depth.txt",
        bam = ALIGNMENT + "{sample}-sorted.bam"
    params:
        maxbin2 = directory(MAXBIN2 + "{sample}"),
        final_output = directory(OUTPUT + "{sample}"),
        maxbin2_dir = directory(MAXBIN2),
        sample = "{sample}" 
    output:
        flag = "{sample}-maxbin2-done.txt"
    conda:
        "envs/metabinning.yaml"
    shell:
        """
        mkdir -p {params.maxbin2}
	run_MaxBin.pl \
        -contig {input.assembly} \
        -abund {input.depth} \
        -out {params.maxbin2} \
        -thread 8
        
        # Ensure final output directory exists
        mkdir -p {params.final_output}

        #Rename and move the files
        for values in $(ls {params.maxbin2_dir} | grep {params.sample} | grep fasta | awk -F'.' '{{print $2}}')
        do

        mv {params.maxbin2_dir}{params.sample}.$values.fasta {params.final_output}/{params.sample}-maxbin2-$values.fa

        done

        #Eliminate the intermediate files
        rm -r {params.maxbin2_dir}{params.sample}*
        
        touch {output.flag}
        """

rule metabat2:
    input:
        assembly = ASSEMBLY + "{sample}.fa",
        depth = ALIGNMENT + "{sample}-depth.txt",
        bam = ALIGNMENT + "{sample}-sorted.bam"
    params:
        metabat2 = directory(METABAT2 + "{sample}"),
        final_output = directory(OUTPUT + "{sample}"),
        metabat2_dir = directory(METABAT2),
        sample = "{sample}"
    output:
        flag = "{sample}-metabat2-done.txt"
    conda:
        "envs/metabinning.yaml"
    shell:
        """
        metabat2 \
        -i {input.assembly} \
        -a {input.depth} \
        -o {params.metabat2} \
        -t 8
        
        # Ensure final output directory exists
        mkdir -p {params.final_output}

        #Rename and move the files
        for values in $(ls {params.metabat2_dir} | grep {params.sample} | awk -F'.' '{{print $2}}' | sed -e 's|BinInfo||')
        do
            mv {params.metabat2_dir}{params.sample}.$values.fa {params.final_output}/{params.sample}-metabat2-$values.fa
        done

        #Eliminate intermediate files
        rm {params.metabat2_dir}{params.sample}.BinInfo.txt
        
        touch {output.flag}
        """
        
rule concoct:
    input:
        assembly_10k = CONTIGS_10K + "{sample}_10K.fasta",
        depth = ALIGNMENT + "{sample}-coverage_table.tsv",
        assembly = ASSEMBLY + "{sample}.fa"
    params:
        bins = directory(CONCOCT + "{sample}"),
        bins_final = directory(CONCOCT_FINAL + "{sample}"),
        final_output = directory(OUTPUT + "{sample}"),
        sample = "{sample}"
    output:
        flag = "{sample}-concoct-done.txt"
    conda:
        "envs/metabinning.yaml"
    shell:
        """
        #Create output directories
        mkdir -p {params.bins}
	mkdir -p {params.bins_final}
	
	#Run concoct
	concoct \
        --composition_file {input.assembly_10k} \
        --coverage_file {input.depth} \
        -t 16 \
        -b {params.bins} 
        
        #Merge subcontigs
        scripts/merge_cutup_clustering.py \
        {params.bins}/clustering_gt1000.csv \
        > {params.bins}/clustering_merged.csv
        
        #Extract bins
        scripts/extract_fasta_bins.py \
        {input.assembly} \
        {params.bins}/clustering_merged.csv \
        --output_path {params.bins_final}

        # Ensure final output directory exists
        mkdir -p {params.final_output}
        
        #Rename and move the files
        for values in $(ls {params.bins_final} | awk -F'.' '{{print $1}}')
        do

        mv {params.bins_final}/$values.fa {params.final_output}/{params.sample}-concoct-$values.fa

        done

        #Eliminate intermediate files
        rm -r {params.bins_final}
        rm -r {params.bins}
        
        touch {output.flag}
        """
        
rule semibin2:
    input:
        assembly = ASSEMBLY + "{sample}.fa",
        bam = ALIGNMENT + "{sample}-sorted.bam"
    params:
        final_output = directory(OUTPUT + "{sample}"),
        semibin2 = directory(SEMIBIN2 + "{sample}")
    output:
        flag = "{sample}-semibin2-done.txt"
    conda:
        "envs/metabinning.yaml"
    shell:
        """
        SemiBin2 single_easy_bin \
        --self-supervised \
        -i {input.assembly} \
        -b {input.bam} \
        --threads 16 \
        -o {params.semibin2}

        # Ensure final output directory exists
        mkdir -p {params.final_output}
        
        #Rename and move the files
        gunzip {params.semibin2}/output_bins/*.gz 

        for f in {params.semibin2}/output_bins/SemiBin_*.fa; do
            id=$(basename "$f" | cut -d'_' -f2 | cut -d'.' -f1)
            mv "$f" {params.final_output}/{wildcards.sample}-semibin2-${{id}}.fa
        done

        #Eliminate intermediate files
        rm -r {params.semibin2}
        
        touch {output.flag}
        """

rule checkM: 
    input:
        semibin2 = "{sample}-semibin2-done.txt",
        concoct = "{sample}-concoct-done.txt",
        metabat2 = "{sample}-metabat2-done.txt",
        maxbin2 = "{sample}-maxbin2-done.txt"
    params:
        final_output = directory(OUTPUT + "{sample}/"),
        alignment = directory(ALIGNMENT),
        sample = "{sample}", 
	assembly_10k = CONTIGS_10K + "{sample}_10K.fasta"
    output:
        checkM = directory(CHECKM + "{sample}"), 
        flag = "{sample}-checkM-done.txt"
    shell:
        """
        #Eliminate the alignment files that we will not need anymore
        rm {params.alignment}{params.sample}*
	rm {params.assembly_10k}

	#Create the output directory
	mkdir -p {output.checkM}

        checkm lineage_wf \
	-x fa \
	-t 24 \
	--tab_table \
	--quiet \
	{params.final_output} \
	{output.checkM} > {output.checkM}/checkm_quality_summary.tsv
        touch {output.flag}
        """   

rule end:
    input:
        checkM = "{sample}-checkM-done.txt"
    output:
        done = "{sample}-all_done.txt"
    shell:
        """
        touch {output.done}
        """
