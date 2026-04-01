# ============================================================
# comparative-lncRNA-pipeline Snakefile
#
# Select species at runtime:
#   snakemake --config species=spruce --profile profiles/slurm --use-conda
#   snakemake --config species=pine   --profile profiles/slurm --use-conda
#
# Dry run:
#   snakemake --config species=spruce --dry-run --cores 4
# ============================================================

import os

configfile: "config/config.yaml"

# ---- Resolve active species block ----
SPECIES = config.get("species", config.get("default_species", "pine"))
if SPECIES not in config:
    raise ValueError(
        f"Species '{SPECIES}' not found in config.yaml. "
        f"Available: {[k for k in config if isinstance(config[k], dict)]}"
    )
SP = config[SPECIES]

# ---- Pull values from active species block ----
SAMPLES     = SP["samples"]["names"]
SAMPLES_DIR = SP["samples"]["directory"]
GENOME      = SP["genome"]["reference"]
ANNOTATION  = SP["genome"]["annotation"]
THREADS     = SP["slurm"]["threads"]
OUT_BASE    = SP["output"]["base"]
OUT_PAF     = SP["output"]["paf"]
OUT_BED     = SP["output"]["bed"]
OUT_FASTA   = SP["output"]["fasta"]
OUT_GO      = SP["output"]["go"]
OUT_KEGG    = SP["output"]["kegg"]
OUT_PLOTS   = SP["output"]["plots"]
LOGS_DIR    = SP["slurm"]["logs"]
ERRORS_DIR  = SP["slurm"]["errors"]

# ---- Shared parameters ----
INDEX_SIZE = config["minimap2"]["index_size"]
PRESET     = config["minimap2"]["preset"]

# ---- Detect which comparison categories are possible --------
stress_samples  = [s for s in SAMPLES if s not in ("SSE", "SZE")]
HAS_NEEDLE  = any(s.endswith("N") for s in stress_samples)
HAS_ROOT    = any(s.endswith("R") for s in stress_samples)
HAS_COLD    = any(s[1] == "C"    for s in stress_samples)
HAS_DROUGHT = any(s[1] == "D"    for s in stress_samples)
HAS_EMBRYO  = "SSE" in SAMPLES and "SZE" in SAMPLES

CATEGORIES = []
if HAS_NEEDLE and HAS_ROOT:
    CATEGORIES += ["needle_specific", "root_specific"]
if HAS_COLD and HAS_DROUGHT:
    CATEGORIES += ["cold_specific", "drought_specific"]
if HAS_EMBRYO:
    CATEGORIES += ["somatic_embryo_specific", "zygotic_embryo_specific"]

# ============================================================
# Rule: all
# ============================================================
rule all:
    input:
        expand(os.path.join(OUT_PAF,   "{sample}.paf"),   sample=SAMPLES),
        expand(os.path.join(OUT_BED,   "{sample}.bed"),   sample=SAMPLES),
        expand(os.path.join(OUT_FASTA, "{sample}.fasta"), sample=SAMPLES),
        os.path.join(OUT_BASE, "multiinter_output.bed"),
        os.path.join(OUT_BASE, "conserved.bed"),
        expand(os.path.join(OUT_BASE, "{cat}.bed"),                        cat=CATEGORIES),
        expand(os.path.join(OUT_GO,   "{cat}_BP_GO_enrichment.txt"),       cat=CATEGORIES),
        expand(os.path.join(OUT_KEGG, "{cat}_KEGG_enrichment.txt"),        cat=CATEGORIES),
        os.path.join(OUT_PLOTS, "upset_plot.png"),

# ============================================================
# Rule: align — minimap2 per sample (SLURM job)
# ============================================================
rule align:
    input:
        fasta  = os.path.join(SAMPLES_DIR, "candidate_transcript_{sample}.fasta"),
        genome = GENOME,
    output:
        paf = os.path.join(OUT_PAF, "{sample}.paf"),
        bed = os.path.join(OUT_BED, "{sample}.bed"),
    params:
        preset     = PRESET,
        index_size = INDEX_SIZE,
        threads    = THREADS,
        logs_dir   = LOGS_DIR,
        errors_dir = ERRORS_DIR,
    resources:
        slurm_account   = SP["slurm"]["account"],
        slurm_partition = SP["slurm"]["partition"],
        runtime         = SP["slurm"]["runtime"],
        mem_mb          = 40000,
        cpus_per_task   = THREADS,
    log:
        os.path.join(LOGS_DIR, "align_{sample}.log"),
    shell:
        """
        mkdir -p {params.logs_dir} {params.errors_dir} \
                 $(dirname {output.paf}) $(dirname {output.bed})

        minimap2 \
            -x {params.preset} \
            -I {params.index_size} \
            -t {params.threads} \
            {input.genome} {input.fasta} \
            > {output.paf} 2> {log}

        awk 'BEGIN{{OFS="\\t"}} $5 != "*" {{
            score = int($12 * 1000 / 60);
            if (score > 1000) score = 1000;
            print $6, $8, $9, $1, score, $5
        }}' {output.paf} | sort -k1,1 -k2,2n > {output.bed}
        """

# ============================================================
# Rule: extract_fasta — bedtools getfasta (login node)
# ============================================================
rule extract_fasta:
    input:
        bed    = os.path.join(OUT_BED,   "{sample}.bed"),
        genome = GENOME,
    output:
        fasta  = os.path.join(OUT_FASTA, "{sample}.fasta"),
    resources:
        runtime = 60,
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        bedtools getfasta \
            -fi {input.genome} \
            -bed {input.bed} \
            -fo {output.fasta}
        """

# ============================================================
# Rule: multiinter — bedtools region comparison (login node)
# ============================================================
rule multiinter:
    input:
        beds = expand(os.path.join(OUT_BED, "{sample}.bed"), sample=SAMPLES),
    output:
        multiinter              = os.path.join(OUT_BASE, "multiinter_output.bed"),
        conserved               = os.path.join(OUT_BASE, "conserved.bed"),
        needle_specific         = os.path.join(OUT_BASE, "needle_specific.bed"),
        root_specific           = os.path.join(OUT_BASE, "root_specific.bed"),
        cold_specific           = os.path.join(OUT_BASE, "cold_specific.bed"),
        drought_specific        = os.path.join(OUT_BASE, "drought_specific.bed"),
        somatic_embryo_specific = os.path.join(OUT_BASE, "somatic_embryo_specific.bed"),
        zygotic_embryo_specific = os.path.join(OUT_BASE, "zygotic_embryo_specific.bed"),
    params:
        names       = " ".join(SAMPLES),
        n           = len(SAMPLES),
        has_needle  = HAS_NEEDLE,
        has_root    = HAS_ROOT,
        has_cold    = HAS_COLD,
        has_drought = HAS_DROUGHT,
        has_embryo  = HAS_EMBRYO,
        col         = {s: (i + 6) for i, s in enumerate(SAMPLES)},
        out_base    = OUT_BASE,
    resources:
        runtime = 60,
    run:
        import subprocess, os
        os.makedirs(params.out_base, exist_ok=True)

        bed_files = " ".join(input.beds)

        subprocess.run(
            f"bedtools multiinter -i {bed_files} -names {params.names} "
            f"> {output.multiinter}",
            shell=True, check=True
        )

        subprocess.run(
            f"awk '$4=={params.n}' {output.multiinter} > {output.conserved}",
            shell=True, check=True
        )

        def present(s): return f"${params.col[s]}==1"
        def absent(s):  return f"${params.col[s]}==0"

        def write_specific(pos, neg, outfile):
            p = " && ".join(present(s) for s in pos)
            n = " && ".join(absent(s)  for s in neg)
            expr = f"{p} && {n}" if n else p
            subprocess.run(
                f"awk '{expr}' {output.multiinter} > {outfile}",
                shell=True, check=True
            )

        needle_s  = [s for s in SAMPLES if s.endswith("N") and s not in ("SSE","SZE")]
        root_s    = [s for s in SAMPLES if s.endswith("R") and s not in ("SSE","SZE")]
        cold_s    = [s for s in SAMPLES if len(s) > 1 and s[1]=="C" and s not in ("SSE","SZE")]
        drought_s = [s for s in SAMPLES if len(s) > 1 and s[1]=="D" and s not in ("SSE","SZE")]

        if params.has_needle and params.has_root:
            write_specific(needle_s,  root_s,    output.needle_specific)
            write_specific(root_s,    needle_s,  output.root_specific)
        else:
            open(output.needle_specific, 'w').close()
            open(output.root_specific,   'w').close()

        if params.has_cold and params.has_drought:
            write_specific(cold_s,    drought_s, output.cold_specific)
            write_specific(drought_s, cold_s,    output.drought_specific)
        else:
            open(output.cold_specific,    'w').close()
            open(output.drought_specific, 'w').close()

        if params.has_embryo:
            write_specific(["SSE"], ["SZE"], output.somatic_embryo_specific)
            write_specific(["SZE"], ["SSE"], output.zygotic_embryo_specific)
        else:
            open(output.somatic_embryo_specific, 'w').close()
            open(output.zygotic_embryo_specific, 'w').close()

# ============================================================
# Rule: go_analysis (login node)
# ============================================================
rule go_analysis:
    input:
        bed        = os.path.join(OUT_BASE, "{cat}.bed"),
        annotation = ANNOTATION,
    output:
        gene_go  = os.path.join(OUT_GO, "gene_to_GO.txt"),
        mstrg    = os.path.join(OUT_GO, "mstrg_to_refgene.txt"),
        refgenes = os.path.join(OUT_GO, "{cat}_refgenes.txt"),
        result   = os.path.join(OUT_GO, "{cat}_BP_GO_enrichment.txt"),
    conda: "envs/goanalysis.yaml"
    resources:
        runtime = 120,
    log:
        os.path.join(LOGS_DIR, "go_{cat}.log"),
    shell:
        "Rscript scripts/03_go_analysis.R "
        "{input.bed} {input.annotation} "
        "{output.gene_go} {output.mstrg} {output.refgenes} {output.result} "
        "> {log} 2>&1"

# ============================================================
# Rule: kegg_analysis (login node)
# ============================================================
rule kegg_analysis:
    input:
        refgenes   = os.path.join(OUT_GO,   "{cat}_refgenes.txt"),
        annotation = ANNOTATION,
    output:
        gene_kegg     = os.path.join(OUT_KEGG, "gene_to_KEGG.txt"),
        pathway_names = os.path.join(OUT_KEGG, "kegg_pathway_names.txt"),
        result        = os.path.join(OUT_KEGG, "{cat}_KEGG_enrichment.txt"),
    conda: "envs/goanalysis.yaml"
    resources:
        runtime = 120,
    log:
        os.path.join(LOGS_DIR, "kegg_{cat}.log"),
    shell:
        "Rscript scripts/04_kegg_analysis.R "
        "{input.refgenes} {input.annotation} "
        "{output.gene_kegg} {output.pathway_names} {output.result} "
        "> {log} 2>&1"

# ============================================================
# Rule: plots (login node)
# ============================================================
rule plots:
    input:
        multiinter   = os.path.join(OUT_BASE, "multiinter_output.bed"),
        go_results   = expand(os.path.join(OUT_GO,   "{cat}_BP_GO_enrichment.txt"), cat=CATEGORIES),
        kegg_results = expand(os.path.join(OUT_KEGG, "{cat}_KEGG_enrichment.txt"),  cat=CATEGORIES),
    output:
        upset = os.path.join(OUT_PLOTS, "upset_plot.png"),
    params:
        plots_dir = OUT_PLOTS,
        names     = " ".join(SAMPLES),
    conda: "envs/goanalysis.yaml"
    resources:
        runtime = 60,
    log:
        os.path.join(LOGS_DIR, "plots.log"),
    shell:
        "Rscript scripts/05_plots.R "
        "{input.multiinter} {params.plots_dir} {params.names} "
        "> {log} 2>&1"
