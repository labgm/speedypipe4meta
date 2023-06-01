#!/usr/bin/env python3

from os import path
import argparse
import glob
import os

parser = argparse.ArgumentParser(description='Takes all fastq files from a folder and creates a snakemake config file for processing the datasets.')
parser.add_argument('input_folder', type=str, help='Input folder containing the FASTQ files')
parser.add_argument('--config_file', type=str, default='config.yaml', help='File to write config file [default: config.yaml]')
parser.add_argument('--threads', type=int, default=40, help='Maximum number of threads for each task in the pipeline [default: 4]')
parser.add_argument('--mem_mb', type=int, default=256000, help='Maximum amount of RAM in MB for each task in the pipeline [default: 256000]')
args = parser.parse_args()

if path.isdir(args.input_folder):
    config = open(args.config_file, "w")
    config.write("workdir: %s\n\n" % (os.getcwd()))
    config.write("threads: %d\n" % (args.threads))
    config.write("mem_mb: %d\n\n" % (args.mem_mb))
    config.write("# MEGAHIT parameters\n")
    config.write("megahit:\n")
    config.write("    klist: '65'\n\n")
    config.write("# KRAKEN parameters\n")
    config.write("kraken:\n")
    config.write("    db: '/usr/local/kraken2/db_standard'\n\n")
    config.write("samples:\n")
    forward = glob.iglob(args.input_folder + "/**/*_1.*", recursive=True)
    for f in forward:
        partition = f.partition("_1.")
        sample = partition[0].split("/")[-1]
        r = partition[0] + "_2." + partition[2]
        config.write("    '%s':\n" % (sample))
        config.write("        forward: '%s'\n" % (f))
        config.write("        revers: '%s'\n" % (r))
    config.write("\n# Mininum length of contigs\n")
    config.write("contigs:\n")
    config.write("    minlength: 200\n")
    config.write("DB:\n")
    config.write("    MEGARES_ANN: 'https://www.meglab.org/downloads/megares_v3.00/megares_annotations_v3.00.csv'\n")
    config.write("    MEGARES_FASTA: 'https://www.meglab.org/downloads/megares_v3.00/megares_database_v3.00.fasta'\n")
    config.close()
else:
    print("Please provide a folder containing FASTQ files")
