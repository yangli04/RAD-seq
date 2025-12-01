configfile: "workflow/config.yaml"

from pathlib import Path

SAMPLE=["AKD1C","AKD1N","AKD1xC","AKD1xN","AKD2C","AKD2N","AKD2xC","AKD2xN","ANC1C","ANC1N","ANC1xC","ANC1xN","ANC2C","ANC2N","ANC2xC","ANC2xN","C1C","C1N","C1xC","C1xN","C2C","C2N","C2xC","C2xN","Ch1C","Ch1N","Ch1xC","Ch1xN","Ch2C","Ch2N","Ch2xC","Ch2xN","Chh1C","Chh1N","Chh1xC","Chh1xN","Chh2C","Chh2N","Chh2xC","Chh2xN","CO1C","CO1N","CO1xC","CO1xN","CO2C","CO2N","CO2xC","CO2xN","COR1C","COR1N","COR1xC","COR1xN","COR2C","COR2N","COR2xC","COR2xN"]

INTERDIR = Path("preprocess/interdir")
OUTDIR = Path("preprocess/outdir")

wildcard_constraints:
    sample="|".join(SAMPLE)

rule all:
    input:
        "report_qc/report_falco_before.html",
        "report_qc/report_falco_after.html",
        "report_qc/mapping_rRNA.html",
        "report_qc/mapping_genome.html",
        "report_qc/cutadapt_PE_qc.html",
        OUTDIR / "feature_counts/counts_hg_genome.txt",
        expand(OUTDIR / "calculated_table/{reftype}/{sample}.pileup.tsv",sample=SAMPLE, reftype = ['genome', 'rRNA']),


rule falco_before_PE:
    input:
        r1="raw_data/{sample}_S1_R1_001.fastq.gz",
        r2="raw_data/{sample}_S1_R2_001.fastq.gz",
    output:
        INTERDIR / "quality_control/falco_before/{sample}_r1/fastqc_report.html",
        INTERDIR / "quality_control/falco_before/{sample}_r2/fastqc_report.html",
        INTERDIR / "quality_control/falco_before/{sample}_r1/fastqc_data.txt",
        INTERDIR / "quality_control/falco_before/{sample}_r2/fastqc_data.txt",
        INTERDIR / "quality_control/falco_before/{sample}_r1/summary.txt",
        INTERDIR / "quality_control/falco_before/{sample}_r2/summary.txt",
    params:
        r1=lambda wildcards: INTERDIR / f"quality_control/falco_before/{wildcards.sample}_r1",
        r2=lambda wildcards: INTERDIR / f"quality_control/falco_before/{wildcards.sample}_r2",
        falco=config['falco'],
    shell:
        """
        {params.falco} -o {params.r1} {input.r1}
        {params.falco} -o {params.r2} {input.r2}
        """


rule report_falco_before:
    input:
        expand(INTERDIR / "quality_control/falco_before/{sample}_r1/fastqc_data.txt", sample=SAMPLE),
        expand(INTERDIR / "quality_control/falco_before/{sample}_r2/fastqc_data.txt", sample=SAMPLE),
    output:
        "report_qc/report_falco_before.html",
    threads: 2
    resources:
        mem_mb=8000,
    shell:
        "multiqc -f -m fastqc -n {output} {input}"

rule cutadapt_PE:
    input:
        fq_1="raw_data/{sample}_S1_R1_001.fastq.gz",
        fq_2="raw_data/{sample}_S1_R2_001.fastq.gz",
    output:
        fq_1=INTERDIR / "cutadapt_PE/temp/temp_{sample}_S1_R1.fastq.gz",
        fq_2=INTERDIR / "cutadapt_PE/temp/temp_{sample}_S1_R2.fastq.gz",

        inter_1=INTERDIR / "cutadapt_PE/temp/inter_{sample}_S1_R1.fastq.gz",
        inter_2=INTERDIR / "cutadapt_PE/temp/inter_{sample}_S1_R2.fastq.gz",

        trimmed_1=INTERDIR / "cutadapt_PE/trimmed_{sample}_S1_R1.fastq.gz",
        trimmed_2=INTERDIR / "cutadapt_PE/trimmed_{sample}_S1_R2.fastq.gz",

        report1=INTERDIR / "cutadapt_PE/log/temp_{sample}.cutadapt.log",
        report2=INTERDIR / "cutadapt_PE/log/inter_{sample}.cutadapt.log",
        report3=INTERDIR / "cutadapt_PE/log/trimmed_{sample}.cutadapt.log",

    threads: 10
    resources:
        mem_mb=4000,
    params:
        cutadapt=config['cutadapt']
    shell:
        """
        {params.cutadapt} -j {threads} -U 11 --rename='{{id}}_{{r2.cut_prefix}} {{comment}}' -m 15\
                --max-n=0 -e 0.15 -q 20 --nextseq-trim=20 -O 1 --pair-filter=any \
                -a AGATCGGAAGAGCACACGTCTG -A AGATCGGAAGAGCGTCGTGT \
                -o {output.fq_1} -p {output.fq_2} \
                {input.fq_1} {input.fq_2} > {output.report1}

        {params.cutadapt} -j {threads} -m 15 -u -11 -n 5 -O 1 \
                -g ACACGACGCTCTTCCGATCT -a AGATCGGAAGAGCGTCGTGT --pair-filter=any \
                -G ACACGACGCTCTTCCGATCT -A AGATCGGAAGAGCGTCGTGT \
                -o {output.inter_1} -p {output.inter_2} \
                {output.fq_1} {output.fq_2} > {output.report2}

        {params.cutadapt} -j {threads} -m 15:15 -o {output.trimmed_1} -p {output.trimmed_2} \
                {output.inter_1} {output.inter_2} > {output.report3}
        """

rule report_cutadapt_PE:
    input:
        expand(INTERDIR / "cutadapt_PE/log/temp_{sample}.cutadapt.log", sample=SAMPLE),
        expand(INTERDIR / "cutadapt_PE/log/inter_{sample}.cutadapt.log", sample=SAMPLE),
        expand(INTERDIR / "cutadapt_PE/log/trimmed_{sample}.cutadapt.log", sample=SAMPLE),
    output:
        "report_qc/cutadapt_PE_qc.html",
    shell:
        "multiqc -f -m cutadapt -n {output} {input}"

rule falco_after_PE:
    input:
        r1=INTERDIR /"cutadapt_PE/trimmed_{sample}_S1_R1.fastq.gz",
        r2=INTERDIR / "cutadapt_PE/trimmed_{sample}_S1_R2.fastq.gz",
    output:
        INTERDIR / "quality_control/falco_after/{sample}_PE_r1/fastqc_report.html",
        INTERDIR / "quality_control/falco_after/{sample}_PE_r2/fastqc_report.html",
        INTERDIR / "quality_control/falco_after/{sample}_PE_r1/fastqc_data.txt",
        INTERDIR / "quality_control/falco_after/{sample}_PE_r2/fastqc_data.txt",
        INTERDIR / "quality_control/falco_after/{sample}_PE_r1/summary.txt",
        INTERDIR / "quality_control/falco_after/{sample}_PE_r2/summary.txt",
    params:
        r1=lambda wildcards: INTERDIR / f"quality_control/falco_after/{wildcards.sample}_PE_r1",
        r2=lambda wildcards: INTERDIR / f"quality_control/falco_after/{wildcards.sample}_PE_r2",
        falco=config['falco'],
    shell:
        """
        {params.falco} -o {params.r1} {input.r1}
        {params.falco} -o {params.r2} {input.r2}
        """

rule report_falco_after:
    input:
        expand(INTERDIR / "quality_control/falco_after/{sample}_PE_r1/fastqc_data.txt",sample=SAMPLE),
        expand(INTERDIR / "quality_control/falco_after/{sample}_PE_r2/fastqc_data.txt",sample=SAMPLE),
    output:
        "report_qc/report_falco_after.html",
    threads: 2
    resources:
        mem_mb=8000,
    shell:
        "multiqc -f -m fastqc -n {output} {input}"

rule map_rRNA_PE:
    input:
        r1=INTERDIR / "cutadapt_PE/trimmed_{sample}_S1_R1.fastq.gz",
        r2=INTERDIR / "cutadapt_PE/trimmed_{sample}_S1_R2.fastq.gz",
    output:
        sam=INTERDIR / "map_rRNA/raw_sam/{sample}.rRNA.sam",
        fq1=INTERDIR / "rRNA_depleted/{sample}_R1.fastq.gz",
        fq2=INTERDIR / "rRNA_depleted/{sample}_R2.fastq.gz",
        summary=INTERDIR / "map_rRNA/hisat_summary/{sample}.rRNA.summary",
    params:
        un=lambda wildcards: INTERDIR / f"rRNA_depleted/{wildcards.sample}_R%.fastq.gz",
        ref_rRNA=config["ref_rRNA"],
        tmp=INTERDIR / "map_rRNA/tmp",
    threads: 12
    resources:
        mem_mb=5000,
    shell:
        """
        hisat-3n -q -x {params.ref_rRNA} \
        --summary-file {output.summary} --new-summary -1 {input.r1} -2 {input.r2} \
        --base-change A,G -p {threads} \
        --un-conc-gz {params.un} --no-unal \
        | samtools sort -o {output.sam} -O sam
        """

rule flag_sort_index_depth_rRNA:
    input:
        sam=INTERDIR / "map_rRNA/raw_sam/{sample}.rRNA.sam",
    output:
        temp_bam=INTERDIR / temp("map_rRNA/temp_sam/{sample}.rRNA.0.bam"),
        bam=INTERDIR / "map_rRNA/raw_bam/{sample}.rRNA.bam", 
        flagstat=INTERDIR / "map_rRNA/flagstat/{sample}.rRNA.raw.flagstat",
        index=INTERDIR / "map_rRNA/raw_bam/{sample}.rRNA.bam.bai",
    threads: 2
    resources:
        mem_mb=10000,
    shell:
        """
        samtools view -Shb {input} > {output.temp_bam}
        samtools sort {output.temp_bam} -o {output.bam}
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.flagstat}
        """

rule dedup_rRNA:
    input:
        bam=INTERDIR / "map_rRNA/raw_bam/{sample}.rRNA.bam",
        index=INTERDIR / "map_rRNA/raw_bam/{sample}.rRNA.bam.bai",
    output:
        bam=INTERDIR / temp("map_rRNA/dedup_bam/{sample}.rRNA.dedup.tmp.bam"),
        log=INTERDIR /"map_rRNA/dedup_logs/{sample}.rRNA.dedup.log",
    resources:
        mem_mb=60000,
    params:
        umi_tools=config['umi_tools']
    shell:
        "{params.umi_tools} dedup --paired --stdin={input.bam} --log={output.log} --stdout={output.bam}"

rule flag_sort_index_depth_rRNA_dedup:
    input:
        bam=INTERDIR /"map_rRNA/dedup_bam/{sample}.rRNA.dedup.tmp.bam",
    output:
        bam=INTERDIR /"map_rRNA/dedup_bam/{sample}.rRNA.dedup.bam",
        index=INTERDIR /"map_rRNA/dedup_bam/{sample}.rRNA.dedup.bam.bai",
        flagstat=INTERDIR /"map_rRNA/flagstat_dedup/{sample}.rRNA.dedup.flagstat",
        sam=INTERDIR /"map_rRNA/dedup_sam/{sample}.rRNA.dedup.sam",
    threads: 2
    resources:
        mem_mb=10000,
    shell:
        """
        samtools sort {input} -o {output.bam}
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.flagstat}
        samtools view -h {output.bam} > {output.sam}
        """

rule mapping_rRNA_qc:
    input:
        expand(INTERDIR /"map_rRNA/hisat_summary/{sample}.rRNA.summary", sample=SAMPLE),
        expand(INTERDIR /"map_rRNA/flagstat_dedup/{sample}.rRNA.dedup.flagstat", sample=SAMPLE),
        expand(INTERDIR /"map_rRNA/dedup_logs/{sample}.rRNA.dedup.log", sample=SAMPLE),
        expand(INTERDIR /"map_rRNA/flagstat/{sample}.rRNA.raw.flagstat", sample=SAMPLE),
    output:
        "report_qc/mapping_rRNA.html",
    shell:
        """
        multiqc -f -n {output} {input}
        """

rule hisat_table_rRNA:
    input:
        INTERDIR /"map_rRNA/dedup_sam/{sample}.rRNA.dedup.sam",
    output:
        INTERDIR /"hisat3n_table/rRNA/{sample}.hisat3n_table.tsv",
    threads: 16
    resources:
        mem_mb=3000,
    params:
        ref_rRNA_fa=config["ref_rRNA_fa"],
    shell:
        """
        hisat-3n-table -p {threads} --alignments {input} --ref {params.ref_rRNA_fa} --output-name {output} --base-change A,G
        """


rule calculate_table_rRNA:
    input:
        INTERDIR /"hisat3n_table/rRNA/{sample}.hisat3n_table.tsv",
    output:
        temp=temp(INTERDIR /"calculated_table/rRNA/{sample}.temp.hisat3n_table.bed6"),
        table=OUTDIR / "calculated_table/rRNA/{sample}.pileup.tsv",
    resources:
        mem_mb=30000,
    params:
        ref_rRNA_fa=config["ref_rRNA_fa"],
    shell:
        """
        awk 'BEGIN{{FS="\\t"; OFS="\\t"}} NR > 1 && $7+$5 > 2  {{ print $1,$2,$3,$5/($7+$5),$7+$5 }}' {input} > {output.temp}
        variant motif -i {output.temp} -o {output.table} -c 1,2,3 -f {params.ref_rRNA_fa}
        """

rule map_genome_PE:
    input:
        r1=INTERDIR / "rRNA_depleted/{sample}_R1.fastq.gz",
        r2=INTERDIR / "rRNA_depleted/{sample}_R2.fastq.gz",
    output:
        sam=INTERDIR / "map_genome/raw_sam/{sample}.genome.sam",
        fq1=INTERDIR / "genome_depleted/{sample}_R1.fastq.gz",
        fq2=INTERDIR / "genome_depleted/{sample}_R2.fastq.gz",
        summary=INTERDIR / "map_genome/hisat_summary/{sample}.genome.summary",
    params:
        un=lambda wildcards: INTERDIR / f"genome_depleted/{wildcards.sample}_R%.fastq.gz",
        ref_genome=config["ref_genome"],
        tmp=INTERDIR / "map_genome/tmp",
    threads: 12
    resources:
        mem_mb=5000,
    shell:
        """
        hisat-3n -q -x {params.ref_genome} \
        --summary-file {output.summary} --new-summary -1 {input.r1} -2 {input.r2} \
        --base-change A,G -p {threads} \
        --un-conc-gz {params.un} --no-unal \
        | samtools sort -o {output.sam} -O sam
        """

rule flag_sort_index_depth_genome:
    input:
        sam=INTERDIR / "map_genome/raw_sam/{sample}.genome.sam",
    output:
        temp_bam=INTERDIR / temp("map_genome/temp_sam/{sample}.genome.0.bam"),
        bam=INTERDIR / "map_genome/raw_bam/{sample}.genome.bam", 
        flagstat=INTERDIR / "map_genome/flagstat/{sample}.genome.raw.flagstat",
        index=INTERDIR / "map_genome/raw_bam/{sample}.genome.bam.bai",
    threads: 2
    resources:
        mem_mb=10000,
    shell:
        """
        samtools view -Shb {input} > {output.temp_bam}
        samtools sort {output.temp_bam} -o {output.bam}
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.flagstat}
        """

rule dedup_genome:
    input:
        bam=INTERDIR / "map_genome/raw_bam/{sample}.genome.bam",
        index=INTERDIR / "map_genome/raw_bam/{sample}.genome.bam.bai",
    output:
        bam=INTERDIR / temp("map_genome/dedup_bam/{sample}.genome.dedup.tmp.bam"),
        log=INTERDIR /"map_genome/dedup_logs/{sample}.genome.dedup.log",
    resources:
        mem_mb=60000,
    params:
        umi_tools=config['umi_tools']
    shell:
        "{params.umi_tools} dedup --paired --stdin={input.bam} --log={output.log} --stdout={output.bam}"

rule flag_sort_index_depth_genome_dedup:
    input:
        bam=INTERDIR /"map_genome/dedup_bam/{sample}.genome.dedup.tmp.bam",
    output:
        bam=INTERDIR /"map_genome/dedup_bam/{sample}.genome.dedup.bam",
        index=INTERDIR /"map_genome/dedup_bam/{sample}.genome.dedup.bam.bai",
        flagstat=INTERDIR /"map_genome/flagstat_dedup/{sample}.genome.dedup.flagstat",
        sam=INTERDIR /"map_genome/dedup_sam/{sample}.genome.dedup.sam",
    threads: 2
    resources:
        mem_mb=10000,
    shell:
        """
        samtools sort {input} -o {output.bam}
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.flagstat}
        samtools view -h {output.bam} > {output.sam}
        """

rule featureCounts_hg_genome:
    input:
        dedup_bam=expand(INTERDIR /"map_genome/dedup_bam/{sample}.genome.dedup.bam", sample=SAMPLE),
    output:
        counts=OUTDIR / "feature_counts/counts_hg_genome.txt",
        summary=OUTDIR / "feature_counts/counts_hg_genome.txt.summary",
    params:
        ref_hg_genome_gtf=config["ref_hg_genome_gtf"],
        featureCounts=config["featureCounts"],
    threads: 12
    resources:
        mem_mb=4000,
    shell:
        """
        {params.featureCounts} -T {threads} --countReadPairs -p -t exon -g gene_id -a {params.ref_hg_genome_gtf} -o {output.counts} {input}
        """

rule mapping_genome_qc:
    input:
        expand(INTERDIR /"map_genome/hisat_summary/{sample}.genome.summary", sample=SAMPLE),
        expand(INTERDIR /"map_genome/flagstat_dedup/{sample}.genome.dedup.flagstat", sample=SAMPLE),
        expand(INTERDIR /"map_genome/dedup_logs/{sample}.genome.dedup.log", sample=SAMPLE),
        expand(INTERDIR /"map_genome/flagstat/{sample}.genome.raw.flagstat", sample=SAMPLE),
    output:
        "report_qc/mapping_genome.html",
    shell:
        """
        multiqc -f -n {output} {input}
        """

rule hisat_table_genome:
    input:
        INTERDIR /"map_genome/dedup_sam/{sample}.genome.dedup.sam",
    output:
        INTERDIR /"hisat3n_table/genome/{sample}.hisat3n_table.tsv",
    threads: 16
    resources:
        mem_mb=3000,
    params:
        ref_genome_fa=config["ref_genome_fa"],
    shell:
        """
        hisat-3n-table -p {threads} --alignments {input} --ref {params.ref_genome_fa} --output-name {output} --base-change A,G
        """


rule calculate_table_genome:
    input:
        INTERDIR /"hisat3n_table/genome/{sample}.hisat3n_table.tsv",
    output:
        temp=temp(INTERDIR /"calculated_table/genome/{sample}.temp.hisat3n_table.bed6"),
        table=OUTDIR / "calculated_table/genome/{sample}.pileup.tsv",
    resources:
        mem_mb=30000,
    params:
        ref_genome_fa=config["ref_genome_fa"],
    shell:
        """
        awk 'BEGIN{{FS="\\t"; OFS="\\t"}} NR > 1 && $7+$5 > 2  {{ print $1,$2,$3,$5/($7+$5),$7+$5 }}' {input} > {output.temp}
        variant motif -i {output.temp} -o {output.table} -c 1,2,3 -f {params.ref_genome_fa}
        """
