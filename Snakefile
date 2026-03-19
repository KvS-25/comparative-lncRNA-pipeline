# ============================================================
# comparative-lncRNA-pipeline Snakefile
# Snakemake workflow for lncRNA stress response analysis
#
# Usage:
#   snakemake --cores 4                    # local run
#   snakemake --profile profiles/slurm     # SLURM cluster run
# ============================================================

import os
import yaml

# ---- Load config ----
configfile: "config/config.yaml"

SAMPLES     = config["samples"]["names"]
OUT_BASE    = config["output"]["base"]
OUT_PAF     = config["output"]["paf"]
OUT_BED     = config["output"]["bed"]
OUT_FASTA   = config["output"]["fasta"]
OUT_GO      = config["output"]["go"]
OUT_KEGG    = config["output"]["kegg"]
OUT_PLOTS   = config["output"]["plots"]
SAMPLES_DIR = config["samples"]["directory"]
GENOME      = config["genome"]["reference"]
THREADS     = config["slurm"]["threads"]
INDEX_SIZE  = config["minimap2"]["index_size"]
PRESET      = config["minimap2"]["preset"]

# ============================================================
# Rule: all — defines final targets
# ============================================================
rule all:
    input:
        # Alignment outputs
        expand("{out}/{sample}.paf", out=OUT_PAF, sample=SAMPLES),
        expand("{out}/{sample}.bed", out=OUT_BED, sample=SAMPLES),
        # Region comparison
        os.path.join(OUT_BASE, "multiinter_output.bed"),
        os.path.join(OUT_BASE, "conserved.bed"),
        os.path.join(OUT_BASE, "needle_specific.bed"),
        os.path.join(OUT_BASE, "root_specific.bed"),
        os.path.join(OUT_BASE, "cold_specific.bed"),
        os.path.join(OUT_BASE, "drought_specific.bed"),
        # GO enrichment
        expand("{out}/{cat}_BP_GO_enrichment.txt", out=OUT_GO,
               cat=["needle_specific","root_specific","cold_specific","drought_specific"]),
        # KEGG enrichment
        expand("{out}/{cat}_KEGG_enrichment.txt", out=OUT_KEGG,
               cat=["needle_specific","root_specific","cold_specific","drought_specific"]),
        # Plots
        os.path.join(OUT_PLOTS, "upset_plot.png"),
        os.path.join(OUT_PLOTS, "region_counts_bar.png"),

# ============================================================
# Rule: align — minimap2 alignment per sample
# ============================================================
rule align:
    input:
        genome = GENOME,
        fasta  = os.path.join(SAMPLES_DIR, "candidate_transcript_{sample}.fasta")
    output:
        paf = os.path.join(OUT_PAF, "{sample}.paf"),
        bed = os.path.join(OUT_BED, "{sample}.bed")
    log:
        "logs/align_{sample}.log"
    threads: THREADS
    resources:
        runtime  = "2880",
        mem_mb   = 40000,
        slurm_partition = config["slurm"]["partition"],
        slurm_account   = config["slurm"]["account"]
    conda:
        "envs/alignment.yaml"
    shell:
        """
        echo "[$(date)] Starting alignment for {wildcards.sample}" >> {log}

        minimap2 \
            -t {threads} \
            -I {INDEX_SIZE} \
            -x {PRESET} \
            {input.genome} \
            {input.fasta} \
            > {output.paf} 2>> {log}

        if [ ! -s {output.paf} ]; then
            echo "ERROR: PAF file empty for {wildcards.sample}" >> {log}
            exit 1
        fi

        echo "[$(date)] Converting PAF to BED" >> {log}

        awk 'BEGIN{{OFS="\\t"}} {{print $6, $8, $9, $1, $12}}' {output.paf} \
            | sort -k1,1 -k2,2n \
            | bedtools merge -c 4,5 -o distinct,mean \
            > {output.bed} 2>> {log}

        echo "[$(date)] Finished {wildcards.sample}" >> {log}
        """

# ============================================================
# Rule: multiinter — bedtools multi-sample comparison
# ============================================================
rule multiinter:
    input:
        beds = expand("{out}/{sample}.bed", out=OUT_BED, sample=SAMPLES)
    output:
        multiinter     = os.path.join(OUT_BASE, "multiinter_output.bed"),
        conserved      = os.path.join(OUT_BASE, "conserved.bed"),
        needle_specific= os.path.join(OUT_BASE, "needle_specific.bed"),
        root_specific  = os.path.join(OUT_BASE, "root_specific.bed"),
        cold_specific  = os.path.join(OUT_BASE, "cold_specific.bed"),
        drought_specific=os.path.join(OUT_BASE, "drought_specific.bed")
    log:
        "logs/multiinter.log"
    conda:
        "envs/alignment.yaml"
    shell:
        """
        echo "[$(date)] Running bedtools multiinter" >> {log}

        bedtools multiinter \
            -i {input.beds} \
            -names {SAMPLES} \
            > {output.multiinter} 2>> {log}

        echo "[$(date)] Filtering regions" >> {log}

        awk '$4 == 4' {output.multiinter} > {output.conserved}
        awk '$6==1 && $7==0 && $8==1 && $9==0' {output.multiinter} > {output.needle_specific}
        awk '$6==0 && $7==1 && $8==0 && $9==1' {output.multiinter} > {output.root_specific}
        awk '$6==1 && $7==1 && $8==0 && $9==0' {output.multiinter} > {output.cold_specific}
        awk '$6==0 && $7==0 && $8==1 && $9==1' {output.multiinter} > {output.drought_specific}

        echo "[$(date)] Done" >> {log}
        """

# ============================================================
# Rule: go_analysis — topGO enrichment
# ============================================================
rule go_analysis:
    input:
        needle = os.path.join(OUT_BASE, "needle_specific.bed"),
        root   = os.path.join(OUT_BASE, "root_specific.bed"),
        cold   = os.path.join(OUT_BASE, "cold_specific.bed"),
        drought= os.path.join(OUT_BASE, "drought_specific.bed")
    output:
        enrichment = expand("{out}/{cat}_BP_GO_enrichment.txt", out=OUT_GO,
               cat=["needle_specific","root_specific","cold_specific","drought_specific"]),
        refgenes = expand("{out}/{cat}_refgenes.txt", out=OUT_GO,
               cat=["needle_specific","root_specific","cold_specific","drought_specific"])
    log:
        "logs/go_analysis.log"
    conda:
        "envs/goanalysis.yaml"
    shell:
        "Rscript scripts/03_go_analysis.R > {log} 2>&1"

# ============================================================
# Rule: kegg_analysis — KEGG pathway enrichment
# ============================================================
rule kegg_analysis:
    input:
        expand("{out}/{cat}_refgenes.txt", out=OUT_GO,
               cat=["needle_specific","root_specific","cold_specific","drought_specific"])
    output:
        expand("{out}/{cat}_KEGG_enrichment.txt", out=OUT_KEGG,
               cat=["needle_specific","root_specific","cold_specific","drought_specific"])
    log:
        "logs/kegg_analysis.log"
    conda:
        "envs/goanalysis.yaml"
    shell:
        "Rscript scripts/04_kegg_analysis.R > {log} 2>&1"

# ============================================================
# Rule: plots — generate all visualisation figures
# ============================================================
rule plots:
    input:
        multiinter = os.path.join(OUT_BASE, "multiinter_output.bed"),
        go_results = expand("{out}/{cat}_BP_GO_enrichment.txt", out=OUT_GO,
                            cat=["needle_specific","root_specific","cold_specific","drought_specific"]),
        kegg_results = expand("{out}/{cat}_KEGG_enrichment.txt", out=OUT_KEGG,
                              cat=["needle_specific","root_specific","cold_specific","drought_specific"])
    output:
        upset  = os.path.join(OUT_PLOTS, "upset_plot.png"),
        counts = os.path.join(OUT_PLOTS, "region_counts_bar.png")
    log:
        "logs/plots.log"
    conda:
        "envs/goanalysis.yaml"
    shell:
        "Rscript scripts/05_plots.R > {log} 2>&1"
