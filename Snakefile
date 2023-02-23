configfile: "config.yaml"

workdir: config["workdir"]

rule all:
    input:
        #reativar fastqc_forward = ["results/" + sample + "/fastqc/" + sample + "_1_fastqc.html" for sample in config["samples"]],
        #reativar fastqc_revers = ["results/" + sample + "/fastqc/" + sample + "_2_fastqc.html" for sample in config["samples"]],
        # trim_forward_paired = ["results/" + sample + "/trimmomatic/" + sample + "_forward_paired.fq.gz" for sample in config["samples"]],
        # trim_forward_unpaired = ["results/" + sample + "/trimmomatic/" + sample + "_forward_unpaired.fq.gz" for sample in config["samples"]],
        # trim_reverse_paired = ["results/" + sample + "/trimmomatic/" + sample + "_reverse_paired.fq.gz" for sample in config["samples"]],
        # trim_reverse_unpaired = ["results/" + sample + "/trimmomatic/" + sample + "_reverse_unpaired.fq.gz" for sample in config["samples"]]
        # scaffolds = ["results/" + sample + "/megahit/final.contigs.fa" for sample in config["samples"]]
        classification = ["results/" + sample + "/kraken/resultado.txt" for sample in config["samples"]],
        anotation = ["results/" + sample + "/prokka/" + sample + ".gbk" for sample in config["samples"]],
        # bowtie2 = ["results/" + sample + "/bowtie2/" + sample + ".sam"  for sample in config["samples"]],
        abundance = ["results/" + sample + "/pileup/" + sample + "_abundance.txt"  for sample in config["samples"]],
        binning = ["results/" + sample + "/maxbin/" for sample in config["samples"]]
        # filter = ["results/" + sample + "/mob_recon/chromosome_filtered.fasta" for sample in config["samples"]],
        # prokka = ["results/" + sample + "/prokka/" + sample + ".gbk" for sample in config["samples"]],
        # quast = ["results/" + sample + "/quast/report.tsv" for sample in config["samples"]]

# TODO: remember to remove files extracted at the end of the pipeline

# (CARD 01)------------------------------------------------------------------------------------------        
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
        TrimmomaticPE \
            -threads {threads} \
            {input[0]} {input[1]} \
            {output[0]} {output[1]} \
            {output[2]} {output[3]} \
            ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:70 \
        """
# Comando em Linha Única
# nTrimmomaticPE -threads 4 /home/victorf/pipeline-metagenomic/data/M1_1.fq.gz /home/victorf/pipeline-metagenomic/data/M1_2.fq.gz results/M1/trimmomatic/M1_forward_paired.fq.gz results/M1/trimmomatic/M1_forward_unpaired.fq.gz results/M1/trimmomatic/M1_reverse_paired.fq.gz results/M1/trimmomatic/M1_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:70

# Step 3: Pre-Assembly
rule kmerstream:
    input:
        forward = "results/{sample}/adapterremoval/{sample}_1.fastq",
        revers = "results/{sample}/adapterremoval/{sample}_2.fastq",
        singleton = "results/{sample}/adapterremoval/{sample}_singleton.fastq"
    params:
        collapsed = "results/{sample}/adapterremoval/{sample}_collapsed.fastq",
        collapsed_truncated = "results/{sample}/adapterremoval/{sample}_collapsed_truncated.fastq"
    output:
        "results/{sample}/kmerstream/ar-{sample}.tsv"
    log:
        stdout = "results/{sample}/kmerstream/log-stdout.txt",
        stderr = "results/{sample}/kmerstream/log-stderr.txt"
    conda:
        "envs/kmerstream.yaml"
    benchmark:
        "results/{sample}/kmerstream/benchmark.txt"
    threads:
        config["threads"]
    shell:
        """
        params=()
        if [[ -f {params.collapsed} && -f {params.collapsed_truncated} ]]; then
            params+=({params.collapsed} {params.collapsed_truncated})
        fi
        KmerStream --kmer-size=7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81,83,85,87,89,91,93,95,97,99,101,103,105,107,109,111,113,115,117,119,121,123,125,127 --output={output} --threads={threads} --tsv {input.forward} {input.revers} {input.singleton} "${{params[@]}}"
        """

# (CARD 02)------------------------------------------------------------------------------------------        
# Step 4: Assembly Megahit
rule megahit:
    input:
        forward = "results/{sample}/trimmomatic/{sample}_forward_paired.fq.gz",
        revers = "results/{sample}/trimmomatic/{sample}_reverse_paired.fq.gz"
    params:
        klist = config["megahit"]['klist'],
        output = "results/{sample}/megahit"
    #     prefix = "results/{sample}/spades",
    #     collapsed = "results/{sample}/adapterremoval/{sample}_collapsed.fastq",
    #     collapsed_truncated = "results/{sample}/adapterremoval/{sample}_collapsed_truncated.fastq"
    output:
        "results/{sample}/megahit/final.contigs.fa"
    log:
        stdout = "results/{sample}/megahit/log-stdout.txt",
        stderr = "results/{sample}/megahit/log-stderr.txt"
    conda:
        "envs/megahit.yaml"
    # benchmark:
    #     "results/{sample}/spades/benchmark.txt"
    # threads:
    #     config["threads"]
    resources:
        mem_gb = 100
    threads:
        config["threads"]
    shell:
        """
            megahit -f -1 {input[0]} -2 {input[1]} -t {threads} --k-list {params.klist} -o {params.output}
        """
        # megahit -1 {input[0]}  -2 {input[1]}  -t {threads} -m {resources.mem_gb} -o {params.output}
# Step 5 (Classificação)

rule kraken:
    input:
        "results/{sample}/megahit/final.contigs.fa"
    output:
        "results/{sample}/kraken/resultado.txt"
    params:
        db = config["kraken"]['db']
    threads:
        config["threads"]
    shell:
        "kraken2 --db {params.db} --threads {threads} --report {output} {input}"

# Cloning the prokka repository because tbl2asn in bioconda is old and throws errors
rule install_prokka:
    output:
        "results/bin/prokka/binaries/linux/tbl2asn"
    conda:
        "envs/prokka.yaml"
    shell:
        """
        rm -rf results/bin/prokka
        git clone https://github.com/tseemann/prokka.git results/bin/prokka > /dev/null 2> /dev/null
        """

rule prokka:
    input:
        # prokka = "results/bin/prokka",
        prokka = "results/bin/prokka/binaries/linux/tbl2asn",
        # chromosome = "results/{sample}/mob_recon/chromosome.fasta"
        chromosome = "results/{sample}/megahit/final.contigs.fa"
    params:
        outdir = "results/{sample}/prokka",
        prefix = "{sample}",
        prokka = "results/bin/prokka/binaries/linux",
        # prokka2 = "results/bin/prokka/bin",
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
        # prokka --force --cpus {threads} --outdir {params.outdir} --prefix {params.prefix} {input.chromosome} --centre X --compliant > {log.stdout} 2> {log.stderr}
        # """
        """
        export PATH={params.prokka}:$PATH
        prokka --force --cpus {threads} --outdir {params.outdir} --prefix {params.prefix} {input.chromosome} --centre X --compliant > {log.stdout} 2> {log.stderr}
        """

rule bowtie2:
    input:
        contig = "results/{sample}/megahit/final.contigs.fa",
        forward = "results/{sample}/trimmomatic/{sample}_forward_paired.fq.gz",
        revers = "results/{sample}/trimmomatic/{sample}_reverse_paired.fq.gz"
    output:
        sam = "results/{sample}/bowtie2/{sample}.sam",
        # db = "results/{sample}/bowtie2/"
    params:
        db = "{sample}db"
    threads:
        config["threads"]
    shell:
        """
            bowtie2-build -f {input.contig} {params.db}
            bowtie2 -x {params.db} -1 {input.forward} -2 {input.revers} --no-unal -p {threads} -S {output.sam} 
        """
# https://sourceforge.net/projects/bbmap/files/latest/download
rule pileup_install:
    output:
        "results/bin/bbmap/pileup.sh"
    shell:
        """
            wget https://sinalbr.dl.sourceforge.net/project/bbmap/BBMap_39.01.tar.gz
            tar -xvzf BBMap_39.01.tar.gz
            mv bbmap/ results/bin/
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
            {input.pileup} -Xmx1G in={input.sam} out={output.coverage}
            awk '{params.var}' {output.coverage} | grep -v '^#'> {output.abundance}
        """

# def BINS(wildcards):
#     return [str(wildcards.bin_id)]

rule maxbin2:
    input:
        contig = "results/{sample}/megahit/final.contigs.fa",
        abundance = "results/{sample}/pileup/{sample}_abundance.txt",
        # bins = expand("{bin}.fasta", bin=BINS)
    output:
        directory("results/{sample}/maxbin/")
    params:
        wildcards = "{sample}"
    wildcard_constraints:
        bin_id="[0-9]+"
    conda:
        "envs/maxbin.yaml"
    shell:
        """
            mkdir results/{params.wildcards}/maxbin
            /home/victor/storage/MaxBin-2.2.7/run_MaxBin.pl -min_contig_length 500 -thread {threads} -contig {input.contig} -out {output}/maxbin -abund {input.abundance}
        """

rule prokka_bin:
    input:
        # forward = lambda wildcards: fnmatch.filter(os.listdir("results/{sample}/maxbin/"), '*.fasta')
        fna_files = expand("results/{sample}/maxbin.{bin_id}.fasta", bin_id=BINS)
    output:
        gff_files = expand("results/{sample}/{bin_id}.gff", bin_id=BINS),
        fna_files = expand("results/{sample}/{bin_id}.fna", bin_id=BINS),
        faa_files = expand("results/{sample}/{bin_id}.faa", bin_id=BINS)
    wildcard_constraints:
        bin_id="[0-9]+"
    threads:
        config["threads"]
    shell:
        "prokka --cpus {threads} {params.prokka_args} {input.fna_files} && "


# rule prokka_bin:
#     input:
#         fna_files = expand("results/{sample}/maxbin.{bin_id}.fasta", bin_id=BINS)
#     output:
#         gff_files = expand("results/{sample}/{bin_id}.gff", bin_id=BINS),
#         fna_files = expand("results/{sample}/{bin_id}.fna", bin_id=BINS),
#         faa_files = expand("results/{sample}/{bin_id}.faa", bin_id=BINS)
#     params:
#         prokka_args = "--kingdom Bacteria --outdir {PROKKA_OUTDIR}"
#     threads: 4
#     shell:
#         "prokka --cpus {threads} {params.prokka_args} {input.fna_files} && "
#         "mv {input.fna_files}/*.gff {output.gff_files} && "
#         "mv {input.fna_files}/*.fna {output.fna_files} && "
#         "mv {input.faa_files}/*.faa {output.faa_files}"


# rule prokkabin:
#     input:
#         # prokka = "results/bin/prokka",
#         prokka = "results/bin/prokka/binaries/linux/tbl2asn",
#         # chromosome = "results/{sample}/mob_recon/chromosome.fasta"
#         chromosome = "results/{sample}/megahit/final.contigs.fa"
#     params:
#         outdir = "results/{sample}/prokka",
#         prefix = "{sample}",
#         prokka = "results/bin/prokka/binaries/linux",
#         # prokka2 = "results/bin/prokka/bin",
#     output:
#         "results/{sample}/prokka/{sample}.gbk"
#     log:
#         stdout = "results/{sample}/prokka/log-stdout.txt",
#         stderr = "results/{sample}/prokka/log-stderr.txt"
#     conda:
#         "envs/prokka.yaml"
#     benchmark:
#         "results/{sample}/prokka/benchmark.txt"
#     threads:
#         config["threads"]
#     shell:
#         # """
#         # prokka --force --cpus {threads} --outdir {params.outdir} --prefix {params.prefix} {input.chromosome} --centre X --compliant > {log.stdout} 2> {log.stderr}
#         # """
#         """
#         export PATH={params.prokka}:$PATH
#         prokka --force --cpus {threads} --outdir {params.outdir} --prefix {params.prefix} {input.chromosome} --centre X --compliant > {log.stdout} 2> {log.stderr}
#         """
# rule metabat:
#     input:
#         contigs="contigs.fa",
#         bam="leituras.bam"
#     output:
#         "binarios.fa"
#     shell:
#         "metabat2 -i contigs.fa -a leituras.bam -o binarios.fa"


# Step 
# rule extract:
#     input:
#         forward = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["forward"]),
#         revers = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["revers"])
#     output:
#         forward = "results/{sample}/extract-file/{sample}_1.fastq",
#         revers = "results/{sample}/extract-file/{sample}_2.fastq"
#     conda:
#         "envs/extract-file.yaml"
#     benchmark:
#         "results/{sample}/extract-file/benchmark.txt"
#     shell:
#         """
#         ./scripts/extract-file.sh {input.forward} {output.forward}
#         ./scripts/extract-file.sh {input.revers} {output.revers}
#         """

# rule mob_recon:
#     input:
#         "results/{sample}/cdhit/contigs.fasta"
#     params:
#         output = "results/{sample}/mob_recon"
#     output:
#         "results/{sample}/mob_recon/chromosome.fasta"
#     log:
#         stdout = "results/{sample}/mob_recon/log-stdout.txt",
#         stderr = "results/{sample}/mob_recon/log-stderr.txt"
#     conda:
#         "envs/mobsuite.yaml"
#     benchmark:
#         "results/{sample}/mob_recon/benchmark.txt"
#     threads:
#         config["threads"]
#     shell:
#         """
#         mob_recon --force -u -c -t -n {threads} -i {input} -o {params.output} > {log.stdout} 2> {log.stderr}
#         """

# rule filter_contigs:
#     input:
#         "results/{sample}/mob_recon/chromosome.fasta"
#     output:
#         "results/{sample}/mob_recon/chromosome_filtered.fasta"
#     conda:
#         "envs/filter_contigs.yaml"
#     params:
#         length = config['contigs']['minlength']
#     shell:
#         """
#         ./scripts/filter_contigs.py {params.length} {input} {output}
#         """

# Cloning the prokka repository because tbl2asn in bioconda is old and throws errors
