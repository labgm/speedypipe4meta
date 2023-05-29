# /home/victor/miniconda3/bin:/home/victor/miniconda3/condabin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/usr/lib/jvm/java-8-oracle/bin:/usr/lib/jvm/java-8-oracle/db/bin:/usr/lib/jvm/java-8-oracle/jre/bin:/home/victor/local/bin
# /home/victor/mambaforge/bin/conda:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/usr/lib/jvm/java-8-oracle/bin:/usr/lib/jvm/java-8-oracle/db/bin:/usr/lib/jvm/java-8-oracle/jre/bin
configfile: "config.yaml"

workdir: config["workdir"]

timezone: 'America/Sao_Paulo'

rule all:
    input:
        fastqc_forward = ["results/" + sample + "/fastqc/" + sample + "_1_fastqc.html" for sample in config["samples"]],
        fastqc_revers = ["results/" + sample + "/fastqc/" + sample + "_2_fastqc.html" for sample in config["samples"]],
        classification = ["results/" + sample + "/kraken/" + sample + "_report.html" for sample in config["samples"]],
        anotation = ["results/" + sample + "/prokka/" + sample + ".gbk" for sample in config["samples"]],
        anotation_bin = ["results/" + sample + "/prokka_bin/" for sample in config["samples"]],
        # megares = ["results/" + sample + "/MEGARes/" for sample in config["samples"]],
        metaquast = ["results/" + sample + "/metaquast/combined_reference/report.tsv" for sample in config["samples"]],
        best_assembly = ["results/" + sample + "/assembly/best_assembly" for sample in config["samples"]]

# TODO: remember to remove files extracted at the end of the pipeline

# Step 1: Quality Analysis
rule fastqc:
    input:
        forward = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["forward"]),
        revers = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["revers"])
    params:
        outdir = "results/{sample}/fastqc"
    output:
        forward = "results/{sample}/fastqc/{sample}_1_fastqc.html",
        revers = "results/{sample}/fastqc/{sample}_2_fastqc.html"
    log:
        stdout = "results/{sample}/fastqc/log-stdout.txt",
        stderr = "results/{sample}/fastqc/log-stderr.txt"
    conda:
        "envs/fastqc.yaml"
    benchmark:
        "results/{sample}/fastqc/benchmark.txt"
    threads:
        config["threads"]
    shell:
        "fastqc --threads {threads} --outdir {params.outdir} {input.forward} {input.revers} > {log.stdout} 2> {log.stderr}"

# Step 2:  Trim and Crop
rule trimmomatic:
    input:
        forward = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["forward"]),
        revers = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["revers"])
    output:
        "results/{sample}/trimmomatic/{sample}_forward_paired.fq.gz", # for paired forward reads
        "results/{sample}/trimmomatic/{sample}_forward_unpaired.fq.gz", # for unpaired forward reads
        "results/{sample}/trimmomatic/{sample}_reverse_paired.fq.gz", # for paired reverse reads
        "results/{sample}/trimmomatic/{sample}_reverse_unpaired.fq.gz"  # for unpaired reverse reads 
    threads:
        config["threads"]
    shell:
        """
        trimmomatic PE \
            -threads {threads} \
            {input[0]} {input[1]} \
            {output[0]} {output[1]} \
            {output[2]} {output[3]} \
            ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:70 \
        """

# Step 3: Assembly Megahit
rule megahit:
    input:
        forward = "results/{sample}/trimmomatic/{sample}_forward_paired.fq.gz",
        revers = "results/{sample}/trimmomatic/{sample}_reverse_paired.fq.gz"
    params:
        klist = config["megahit"]['klist'],
        output = "results/{sample}/assembly/megahit"
    output:
        "results/{sample}/assembly/megahit/final.contigs.fa"
    log:
        stdout = "results/{sample}/assembly/megahit/log-stdout.txt",
        stderr = "results/{sample}/assembly/megahit/log-stderr.txt"
    conda:
        "envs/megahit.yaml"
    # benchmark:
    #     "results/{sample}/spades/benchmark.txt"
    resources:
        mem_gb = 100
    threads:
        config["threads"]
    shell:
        """
            megahit -f -1 {input[0]} -2 {input[1]} -t {threads} --presets meta-large -o {params.output} --min-contig-len 300
        """
        # megahit -1 {input[0]}  -2 {input[1]}  -t {threads} -m {resources.mem_gb} -o {params.output} (Ver com o professor)

# Step 3: Assembly MetaSPAdes
rule metaspades:
    input:
        forward = "results/{sample}/trimmomatic/{sample}_forward_paired.fq.gz",
        revers = "results/{sample}/trimmomatic/{sample}_reverse_paired.fq.gz"
    output:
        "results/{sample}/assembly/metaspades/scaffolds.fasta"
    params:
        klist = "27,37,47,57,67,77,87,97,107,117,127",
        output = "results/{sample}/assembly/metaspades"
    log:
        stdout = "results/{sample}/assembly/metaspades/log-stdout.txt",
        stderr = "results/{sample}/assembly/metaspades/log-stderr.txt"
    # conda:
    #     "envs/metaspades.yaml"
    resources:
        mem_gb = 100
    threads:
        config["threads"]
    shell:
        """
            metaspades.py -o {params.output} -1 {input.forward} -2 {input.revers} -k {params.klist} --threads {threads} --memory {resources.mem_gb}
        """
        # -k {config["metaspades"]["kmer_sizes"]}
# Step 3: Assembly IDBA-UD
rule pre_idba:
    input:
        forward = "results/{sample}/trimmomatic/{sample}_forward_paired.fq.gz",
        revers = "results/{sample}/trimmomatic/{sample}_reverse_paired.fq.gz"
    output:
        forward = "results/{sample}/trimmomatic/{sample}_forward_paired.fq",
        revers = "results/{sample}/trimmomatic/{sample}_reverse_paired.fq"
    shell:
        """
            gzip -dk {input.forward} {input.revers}
        """

rule fq2fa:
    input:
        forward = "results/{sample}/trimmomatic/{sample}_forward_paired.fq",
        revers = "results/{sample}/trimmomatic/{sample}_reverse_paired.fq"
    output:
        "results/{sample}/assembly/idba/{sample}_1_2.fa"
    shell: 
        """
            fq2fa --merge {input.forward} {input.revers} {output}
        """

rule idba:
    input:
        "results/{sample}/assembly/idba/{sample}_1_2.fa"
    output:
        "results/{sample}/assembly/idba/scaffold.fa"
    params:
        outdir = "results/{sample}/assembly/idba",
    log:
        stdout = "results/{sample}/assembly/idba/log-stdout.txt",
        stderr = "results/{sample}/assembly/idba/log-stderr.txt"
    # conda:
        # "envs/idba.yaml"
    resources:
        mem_gb = 100
    threads:
        config["threads"]
    shell:
        """
            idba_ud -r {input} -o {params.outdir} --mink 27 --step 10 --maxk 127 --num_threads {threads}
        """
            #idba_ud -r {input.forward} -l {input.revers} -o {output} --num_threads {threads} --mink {config["idba"]["min_kmer_size"]} --maxk {config["idba"]["max_kmer_size"]} --step {config["idba"]["kmer_step"]} --pre_correction --no_correct --min_contig {config["idba"]["min_contig_length"]}
            # fq2fa --merge F1_forward_paired.fq F1_reverse_paired.fq F1_1_2.fa
# Step 3 Quast Analysis
rule metaquast:
    input: 
        megahit = "results/{sample}/assembly/megahit/final.contigs.fa",
        metaspades = "results/{sample}/assembly/metaspades/scaffolds.fasta",
        idba = "results/{sample}/assembly/idba/scaffold.fa"
    output:
        "results/{sample}/metaquast/combined_reference/report.tsv"
    params:
        outdir = "results/{sample}/metaquast"
    shell:
        """
            metaquast.py -o {params.outdir} {input.megahit} {input.metaspades} {input.idba} -l Megahit,SPades,Idba
        """

rule best_assembly:
    input:
        "results/{sample}/metaquast/combined_reference/report.tsv"
    output:
        directory("results/{sample}/assembly/best_assembly")
    params:
        "results/{sample}/assembly/final.contigs.fa"
    shell:
        """
            python scripts/best_assembly.py --input {input}
            
        """

# Step 4 (Classificação)
rule kraken:
    input:
        forward = "results/{sample}/trimmomatic/{sample}_forward_paired.fq.gz",
        revers = "results/{sample}/trimmomatic/{sample}_reverse_paired.fq.gz"
    output:
        report_kraken = "results/{sample}/kraken/{sample}_report.txt",
        report_krona = "results/{sample}/kraken/{sample}_report.krona",
        report_html = "results/{sample}/kraken/{sample}_report.html",
        
    log:
        stdout = "results/{sample}/kraken/log-stdout.txt",
        stderr = "results/{sample}/kraken/log-stderr.txt"
    params:
        db = config["kraken"]['db']
    threads:
        config["threads"]
    shell:
        """
            kraken2 --db {params.db} --threads {threads} --paired {input.forward} {input.revers} --report {output.report_kraken} > {log.stdout} 2> {log.stderr}
            kreport2krona.py -r {output.report_kraken} -o {output.report_krona}
            ktImportText {output.report_krona} -o {output.report_html}
        """

rule prokka:
    input:
        # prokka = "results/bin/prokka/binaries/linux/tbl2asn",
        fasta = "results/{sample}/assembly/best_assembly/final.contigs.fa"
    params:
        outdir = "results/{sample}/prokka",
        prefix = "{sample}",
        # prokka = "results/bin/prokka/binaries/linux"
    output:
        "results/{sample}/prokka/{sample}.gbk"
    log:
        stdout = "results/{sample}/prokka/log-stdout.txt",
        stderr = "results/{sample}/prokka/log-stderr.txt"
    conda:
        "envs/prokka.yaml"
    benchmark:
        "results/{sample}/prokka/benchmark.txt"
    threads:
        config["threads"]
    shell:
        # """
        # export PATH={params.prokka}:$PATH
        # prokka --force --cpus {threads} --outdir {params.outdir} --prefix {params.prefix} {input.fasta} --centre X --compliant > {log.stdout} 2> {log.stderr}
        # """
        """
        prokka --force --cpus {threads} --outdir {params.outdir} --prefix {params.prefix} {input.fasta} --centre X --compliant > {log.stdout} 2> {log.stderr}
        """
# Step 6 (Anotação Contings)
rule bowtie2:
    input:
        contig = "results/{sample}/assembly/best_assembly/final.contigs.fa",
        forward = "results/{sample}/trimmomatic/{sample}_forward_paired.fq.gz",
        revers = "results/{sample}/trimmomatic/{sample}_reverse_paired.fq.gz"
    output:
        sam = "results/{sample}/bowtie2/{sample}.sam",
        db = directory("results/{sample}/bowtie2/")
    params:
        db = "{sample}db"
    threads:
        config["threads"]
    shell:
        """
            bowtie2-build -f {input.contig} {output.db}/{params.db}
            bowtie2 -x  {output.db}/{params.db} -1 {input.forward} -2 {input.revers} --no-unal -p {threads} -S {output.sam} 
        """
# https://sourceforge.net/projects/bbmap/files/latest/download
rule pileup_install:
    output:
        "results/bin/bbmap/pileup.sh"
    shell:
        """
            wget https://sinalbr.dl.sourceforge.net/project/bbmap/BBMap_39.01.tar.gz
            tar -xvzf BBMap_39.01.tar.gz
            if [ ! -d ~/results/bin/bbmap ]; then
                echo 'Movendo diretório...'
                mv bbmap/ results/bin/
            else
                echo 'Diretório já existe'
            fi
        """

rule pileup:
    input:
        pileup = "results/bin/bbmap/pileup.sh",
        sam = "results/{sample}/bowtie2/{sample}.sam"
    output:
        coverage = "results/{sample}/pileup/{sample}_cov.txt",
        abundance = "results/{sample}/pileup/{sample}_abundance.txt"
    params:
        var = "{print $1\"\t\"$5}"
    shell:
        """
            {input.pileup} -Xmx5G in={input.sam} out={output.coverage}
            awk '{params.var}' {output.coverage} | grep -v '^#'> {output.abundance}
        """

import os
import fnmatch

rule maxbin2:
    input:
        contig = "results/{sample}/assembly/best_assembly/final.contigs.fa",
        abundance = "results/{sample}/pileup/{sample}_abundance.txt"
    output:
        folder = directory("results/{sample}/maxbin/"),
        # files = "results/{sample}/maxbin/maxbin.{id_bins}.fasta",
        listBins = "results/{sample}/maxbin/list_bins.txt"
    conda:
        "envs/maxbin.yaml"
    # wildcard_constraints:
    #     id_bins = "[0-9]",
    shell:
        """
            /home/victor/storage/MaxBin-2.2.7/run_MaxBin.pl -min_contig_length 300 -thread {threads} -contig {input.contig} -out {output.folder}/maxbin -abund {input.abundance}
            ls {output.folder}/*.fasta > {output.listBins}
        """
# ls {output.folder}/*.fasta > {output.listBins}
# PREFIXES = [int(f.split(".")[1]) for f in glob_wildcards("results/M1/maxbin.{prefix}.fasta")]
rule prokka_bin:
    input:
        "results/{sample}/maxbin/"
        # "results/{sample}/maxbin/list_bins.txt"
    output:
        directory("results/{sample}/prokka_bin/")
    log:
        stdout = "results/{sample}/prokka_bin/log-stdout.txt",
        stderr = "results/{sample}/prokka_bin/log-stderr.txt"
    params:
        outdir = "results/{sample}/prokka_bin/",
        prefix = "{sample}",
        prokka = "results/bin/prokka/binaries/linux"
    threads:
        config["threads"]
    shell:
        """
            export PATH={params.prokka}:$PATH
            for file in {input}/*.fasta; do 
                prefix=$(basename $file | cut -d '.' -f 2) 
                prokka --force --cpus {threads} --outdir {params.outdir}/$prefix --prefix $prefix $file --centre X --compliant > {log.stdout} 2> {log.stderr}; 
            done
        """

rule megares:
    input:
        'amrplusplus_v2/main_AmrPlusPlus_v2_withRGI.nf'
    output:
        directory("results/{sample}/MEGARes/")
    shell:
        """
            nextflow run {input}  -profile singularity --output {output} -resume
        """