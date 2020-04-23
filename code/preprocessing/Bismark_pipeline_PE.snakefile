# Parallelisation options
import multiprocessing
import sys
import itertools
import os
import collections
import json
import glob
cpuCount = (multiprocessing.cpu_count() - 2)

# Import config and make report
configfile: "config_demultiplex.yaml"
report: "report/workflow.rst"

if config["sequencer"] == "NextSeq":
    OpticalDupsPixelDistance = 100
elif config["sequencer"] == "HiSeq":
    OpticalDupsPixelDistance = 2500
elif config["sequencer"] == "NovaSeq":
    OpticalDupsPixelDistance = 12000
else:
    sys.exit("Specify sequencer (NextSeq, HiSeq, NovaSeq)")

if config["dualindex"] == False:
    ReadRegex = '[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+)_[0-9]+:[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+'
else:
    ReadRegex = '[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+)_[0-9]+:[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+\+[a-zA-Z0-9]+'

if config["demultiplex"] == True:
    # Sample list
    
    f = open("SampleSheet.csv", 'r')
    NAMES ={}
    SAMPLES =[]
    TARGET = []
    TARGET_BASE = []
    LIBRARIES = []
    lineno = 1

    def has_duplicates(my_list):
        for index in range(len(my_list)):
            if index!=0: # check only from 2nd element onwards
                item = my_list[index] # get the current item
                if item in my_list[:index]: # check if current item is in the beginning of the list uptill the current element
                    return True # duplicate exist
        return False # duplicate does not exist

    for line in f:
            if 'Sample_ID' in line:
                column = line.split(",")
                sampleID    =column.index('Sample_ID')
                sampleName  =column.index('Sample_Name')
                libraryID   =column.index('Sample_Project')
                for line in f:
                    column = line.split(",")
                    SAMPLES += [column[sampleID]]
                    NAMES[column[sampleID]] = column[sampleName]
                    LIBRARIES += [column[libraryID]]
                    LIBRARIES = list(filter(None, LIBRARIES))
                    if (has_duplicates(SAMPLES) == False) and (len(LIBRARIES) > 0):
                        # placeholder for demultiplex with library folder
                        TARGET += [column[sampleID]+"_S"+str(lineno)+"_R1_001"]
                        TARGET += [column[sampleID]+"_S"+str(lineno)+"_R2_001"]
                        TARGET_BASE += [column[sampleID]+"_S"+str(lineno)]
                        lineno = lineno + 1
                    elif (has_duplicates(SAMPLES) == False) and (len(LIBRARIES) == 0):
                        TARGET += [column[sampleID]+"_S"+str(lineno)+"_R1_001"]
                        TARGET += [column[sampleID]+"_S"+str(lineno)+"_R2_001"]
                        TARGET_BASE += [column[sampleID]+"_S"+str(lineno)]
                        lineno = lineno + 1

    IDS = tuple(TARGET)
    R1IDS = tuple(TARGET_BASE)
    # LIBRARY_IDS = tuple(list(set(LIBRARIES)))

    rule all:
        input:
            "multiqc_report.html"

    rule demultiplex:
        input:
            samplesheet = "SampleSheet.csv"
        output:
            flag = "demultiplexed_reads/bcl2fastq.SUCCESS",
            all = expand("demultiplexed_reads/{sample}.fastq.gz", sample = IDS),
        threads: 8
        params:
            runfolder = config["runfolder"],
            dir = "demultiplexed_reads",
        shell:
            "ml purge && ml bcl2fastq2/2.20.0-intel-2019a;"
            "bcl2fastq "
            "    --runfolder-dir {params.runfolder}"
            "    --output-dir {params.dir}"
            "    --interop-dir {params.dir}/InterOp "
            "    --sample-sheet {input} "
            "    --loading-threads {threads}"
            "    --processing-threads {threads} "
            "    --writing-threads {threads} "
            "    --barcode-mismatches 0 "
            "    --no-lane-splitting "
            "    --fastq-compression-level 6 "
            "    2>&1 | tee bcl2fastq.log && touch {output.flag} && "
            "    find {params.dir} -iname *gz -exec mv {{}} {params.dir} \;"

    rule adaptortrimming:
        input:
            flag = "demultiplexed_reads/bcl2fastq.SUCCESS",
            fq1 = "demultiplexed_reads/{sample}_R1_001.fastq.gz",
            fq2 = "demultiplexed_reads/{sample}_R2_001.fastq.gz"
        output:
            temp("trimmed_reads/{sample}_R1_001_val_1.fq.gz"),
            temp("trimmed_reads/{sample}_R2_001_val_2.fq.gz")
        run:
            if config["darkcycles"] == True:
                shell("""
                ml purge && ml Trim_Galore/0.6.0-foss-2018b-Python-3.6.6; ml FastQC/0.11.8-Java-1.8;
                trim_galore --paired {input.fq1} {input.fq2} --fastqc --gzip {config[rrbs]} -o ./trimmed_reads
                """)
            elif config["darkcycles"] == False:
                shell("""
                ml purge && ml Trim_Galore/0.6.0-foss-2018b-Python-3.6.6; ml FastQC/0.11.8-Java-1.8;
                trim_galore --paired {input.fq1} {input.fq2} --fastqc --clip_R1 3 --clip_R2 4 --three_prime_clip_R2 1 --gzip {config[rrbs]} -o ./trimmed_reads
                """)
            elif not config["darkcycles"]:
                sys.exit("Specify if dark cycles (True or False)")

else:
    IDS, = glob_wildcards("demultiplexed_reads/{sample}.fastq.gz")
    R1IDS, = glob_wildcards("demultiplexed_reads/{sample}_R1_001.fastq.gz")

    rule all:
        input:
            "multiqc_report.html"

    rule adaptortrimming:
        input:
            fq1 = "demultiplexed_reads/{sample}_R1_001.fastq.gz",
            fq2 = "demultiplexed_reads/{sample}_R2_001.fastq.gz"
        output:
            temp("trimmed_reads/{sample}_R1_001_val_1.fq.gz"),
            temp("trimmed_reads/{sample}_R2_001_val_2.fq.gz")
        run:
            if config["darkcycles"] == True:
                shell("""
                ml purge && ml Trim_Galore/0.6.0-foss-2018b-Python-3.6.6; ml FastQC/0.11.8-Java-1.8;
                trim_galore --paired {input.fq1} {input.fq2} --fastqc --three_prime_clip_R1 1 --three_prime_clip_R2 1 --gzip {config[rrbs]} -o ./trimmed_reads
                """)
            elif config["darkcycles"] == False:
                shell("""
                ml purge && ml Trim_Galore/0.6.0-foss-2018b-Python-3.6.6; ml FastQC/0.11.8-Java-1.8;
                trim_galore_6.0.3 --paired {input.fq1} {input.fq2} --fastqc  --three_prime_clip_R1 1 --three_prime_clip_R2 1 --clip_R1 3 --clip_R2 3 --gzip {config[rrbs]} -o ./trimmed_reads
                """)
            elif not config["darkcycles"]:
                sys.exit("Specify if dark cycles (True or False)")


rule mapping_bismark:
    input:
        genome = config["humangenome"],
        fq1 = "trimmed_reads/{sample}_R1_001_val_1.fq.gz",
        fq2 = "trimmed_reads/{sample}_R2_001_val_2.fq.gz"
    output:
        temp("mapped_reads/{sample}_R1_001_val_1_bismark_bt2_pe.bam"),
        temp(directory("bismarkTemp_{sample}"))
    threads: 4
    params:
        tempdir = "bismarkTemp_{sample}"
    shell:
        "ml purge && ml Bismark/0.20.1-intel-2018b; "
        "bismark --nucleotide_coverage --temp_dir {params.tempdir} --multicore {threads} {input.genome} -1 {input.fq1} -2 {input.fq2} -o ./mapped_reads/"


rule query_sort_samtools:
    input:
        bam = "mapped_reads/{sample}_R1_001_val_1_bismark_bt2_pe.bam"
    output:
        bam = temp("query_sorted_reads/{sample}_R1_001_val_1_bismark_bt2_pe.bam"),
    threads: 4
    shell:
        "ml purge && ml SAMtools/1.9-intel-2018b; "
        "samtools sort -n -@ {threads} {input.bam} -o {output.bam}; "

rule rm_optical_dups:
    input: "query_sorted_reads/{sample}_R1_001_val_1_bismark_bt2_pe.bam"
    output: temp("rm_optical_dups/{sample}_R1_001_val_1_bismark_bt2_pe.bam")
    params:
        OpticalDupsPixelDistance = OpticalDupsPixelDistance,
        ReadRegex = ReadRegex
    shell:
        "ml purge && ml picard/2.18.27-Java-1.8; "
        "java -jar $EBROOTPICARD/picard.jar MarkDuplicates "
        "I={input} "
        "O={output} "
        "M=rm_optical_dups/{wildcards.sample}_MarkDuplicates.txt "
        "OPTICAL_DUPLICATE_PIXEL_DISTANCE={params.OpticalDupsPixelDistance} "
        "REMOVE_SEQUENCING_DUPLICATES=true "
        "TAGGING_POLICY=OpticalOnly "
        "READ_NAME_REGEX='{params.ReadRegex}' "
        "ASSUME_SORT_ORDER=queryname; "

rule sort_samtools:
    input:
        bam = "rm_optical_dups/{sample}_R1_001_val_1_bismark_bt2_pe.bam"
    output:
        bam = "sorted_reads/{sample}_R1_001_val_1_bismark_bt2_pe.bam",
        bai = "sorted_reads/{sample}_R1_001_val_1_bismark_bt2_pe.bam.bai"
    threads: 4
    shell:
        "ml purge && ml SAMtools/1.9-intel-2018b; "
        "samtools sort -@ {threads} {input.bam} -o {output.bam}; "
        "samtools index {output.bam} {output.bai} "

if config["genomebuild"] == "GRCh37" or config["genomebuild"] == "GRCh38" or config["genomebuild"] == "GRCm38":
    rule hsmetrics:
        input:
            bam = "sorted_reads/{sample}_R1_001_val_1_bismark_bt2_pe.bam",
            bai = "sorted_reads/{sample}_R1_001_val_1_bismark_bt2_pe.bam.bai"
        output:
            txt = "{sample}_hs_metrics.txt"
        threads: 1
        params:
            genomepath = config["humangenomefa"],
            targetpath = config["target"],
        shell:
            "ml purge && ml picard/2.18.27-Java-1.8; "
            "java -jar $EBROOTPICARD/picard.jar CollectHsMetrics "
            "I={input.bam} "
            "O={output.txt} "
            "R={params.genomepath} "
            "TARGET_INTERVALS={params.targetpath} "
            "BAIT_INTERVALS={params.targetpath} "
elif not config["genomebuild"]:
    sys.exit("Specify genomebuild (GRCh37 or GRCh38)")

rule mapping_to_lambda_bismark:
    input:
        genome = config["lambdagenome"],
        fq1 = "trimmed_reads/{sample}_R1_001_val_1.fq.gz",
        fq2 = "trimmed_reads/{sample}_R2_001_val_2.fq.gz"
    output:
        "mapped_reads/lambda.{sample}_R1_001_val_1_bismark_bt2_pe.bam",
        temp(directory("bismarkTempLambda_{sample}"))
    threads: 2
    params:
        tempdir = "bismarkTempLambda_{sample}"
    shell:
        "ml purge && ml Bismark/0.20.1-intel-2018b;"
        "bismark --nucleotide_coverage --prefix lambda --temp_dir {params.tempdir} --multicore {threads} {input.genome} -1 {input.fq1} -2 {input.fq2} -o ./mapped_reads/"

rule extract_methylation_calls:
    input:
        "rm_optical_dups/{sample}_R1_001_val_1_bismark_bt2_pe.bam"
    output:
        "methylation_extractor_output/{sample}_R1_001_val_1_bismark_bt2_pe.bismark.cov.gz"
    threads: 6
    shell:
        "ml purge && ml Bismark/0.20.1-intel-2018b; "
        "bismark_methylation_extractor -p --bedGraph --gzip --multicore {threads} {input} -o ./methylation_extractor_output/"

rule deduplicate:
    input:
        "rm_optical_dups/{sample}_R1_001_val_1_bismark_bt2_pe.bam"
    output:
        temp("dedupl_reads/{sample}_R1_001_val_1_bismark_bt2_pe.deduplicated.bam")
    shell:
        "ml purge && ml Bismark/0.20.1-intel-2018b;"
        "deduplicate_bismark -p --bam {input} --output_dir ./dedupl_reads"

rule extract_methylation_calls_deduplicated:
    input:
        "dedupl_reads/{sample}_R1_001_val_1_bismark_bt2_pe.deduplicated.bam"
    output:
        "methylation_extractor_output/{sample}_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
    threads: 5
    shell:
        "ml purge && ml Bismark/0.20.1-intel-2018b;"
        "bismark_methylation_extractor --bedGraph --gzip --multicore {threads} {input} -o ./methylation_extractor_output/"

rule rename:
    input:
        expand("methylation_extractor_output/lambda.{sample}_R1_001_val_1_bismark_bt2_pe.bismark.cov.gz", sample=R1IDS),
        expand("methylation_extractor_output/{sample}_R1_001_val_1_bismark_bt2_pe.{dupl}bismark.cov.gz", sample=R1IDS, dupl = ['deduplicated.', ''])
    output:
        expand("methylation_extractor_output/{sample}{dupl}.cov.gz", sample=R1IDS, dupl = ['dedupl', '']),
        expand("methylation_extractor_output/lambda.{sample}.cov.gz", sample=R1IDS)
    shell:
        "cd methylation_extractor_output/ ;"
        "rename _R1_001_val_1_bismark_bt2_pe.bismark.cov.gz .cov.gz *_R1_001_val_1_bismark_bt2_pe.bismark.cov.gz; "
        "rename _R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz dedupl.cov.gz *_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz;"
        "rm *txt.gz"

if config["genomebuild"] == "GRCh37" or config["genomebuild"] == "GRCh38":
    rule multiqc:
        input:
            expand("{sample}_hs_metrics.txt", sample = R1IDS),
            expand("methylation_extractor_output/{sample}{dupl}.cov.gz", sample=R1IDS, dupl = ['dedupl', '']),
            expand("methylation_extractor_output/lambda.{sample}.cov.gz", sample=R1IDS),
        output:
            report("multiqc_report.html")
        shell:
            "ml purge && ml MultiQC/1.7-intel-2018b-Python-3.6.6 && ml swap matplotlib/2.2.3-intel-2018b-Python-3.6.6; "
            "multiqc -f ."
elif not config["genomebuild"]:
    sys.exit("Specify genomebuild (GRCh37 or GRCh38)")
