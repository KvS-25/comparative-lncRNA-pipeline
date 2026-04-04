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
        expand(os.path.join(OUT_BED,  "{sample}.bed"),  sample=SAMPLES),
        expand(os.path.join(OUT_PAF,  "{sample}.paf"),  sample=SAMPLES),
        os.path.join(OUT_BASE, "multiinter_output.bed"),
        os.path.join(OUT_BASE, "conserved.bed"),
        expand(os.path.join(OUT_BASE, "{cat}.bed"),                        cat=CATEGORIES),
        os.path.join(OUT_GO,   "gene_to_GO.txt"),
        os.path.join(OUT_GO,   "mstrg_to_refgene.txt"),
        os.path.join(OUT_KEGG, "gene_to_KEGG.txt"),
        os.path.join(OUT_KEGG, "kegg_pathway_names.txt"),
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
    conda: "envs/alignment.yaml"
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
    conda: "envs/alignment.yaml"
    params:
        names    = " ".join(SAMPLES),
        n        = len(SAMPLES),
        out_base = OUT_BASE,
    resources:
        runtime = 60,
    shell:
        """
        bedtools multiinter \
            -i {input.beds} \
            -names {params.names} \
            > {output.multiinter}

        awk '$4=={params.n}' {output.multiinter} > {output.conserved}

        python3 - "{output.multiinter}" "{params.out_base}" {params.names} << 'PYEOF'
import sys, subprocess, os

bed_file = sys.argv[1]
out_base = sys.argv[2]
samples  = sys.argv[3:]

col = dict((s, i + 6) for i, s in enumerate(samples))

def awk_filter(present, absent, outfile):
    p = " && ".join("$" + str(col[s]) + "==1" for s in present)
    a = " && ".join("$" + str(col[s]) + "==0" for s in absent)
    expr = (p + " && " + a) if a else p
    cmd = "awk '" + expr + "' " + bed_file + " > " + outfile
    subprocess.run(cmd, shell=True, check=True)
    open(outfile, 'a').close()

stress  = [s for s in samples if s not in ("SSE", "SZE")]
needle  = [s for s in stress if s.endswith("N")]
root    = [s for s in stress if s.endswith("R")]
cold    = [s for s in stress if len(s) > 1 and s[1] == "C"]
drought = [s for s in stress if len(s) > 1 and s[1] == "D"]

if needle and root:
    awk_filter(needle, root,   out_base + "/needle_specific.bed")
    awk_filter(root,   needle, out_base + "/root_specific.bed")
else:
    open(out_base + "/needle_specific.bed", "w").close()
    open(out_base + "/root_specific.bed",   "w").close()

if cold and drought:
    awk_filter(cold,    drought, out_base + "/cold_specific.bed")
    awk_filter(drought, cold,    out_base + "/drought_specific.bed")
else:
    open(out_base + "/cold_specific.bed",    "w").close()
    open(out_base + "/drought_specific.bed", "w").close()

if "SSE" in samples and "SZE" in samples:
    awk_filter(["SSE"], ["SZE"], out_base + "/somatic_embryo_specific.bed")
    awk_filter(["SZE"], ["SSE"], out_base + "/zygotic_embryo_specific.bed")
else:
    open(out_base + "/somatic_embryo_specific.bed", "w").close()
    open(out_base + "/zygotic_embryo_specific.bed", "w").close()
PYEOF
        """

# ============================================================
# Rule: build_go_maps — build gene→GO and MSTRG maps once
#       from the full annotation (no wildcard, runs once)
# ============================================================
rule build_go_maps:
    input:
        annotation = ANNOTATION,
    output:
        gene_go = os.path.join(OUT_GO, "gene_to_GO.txt"),
        mstrg   = os.path.join(OUT_GO, "mstrg_to_refgene.txt"),
    conda: "envs/goanalysis.yaml"
    resources:
        runtime = 60,
    log:
        os.path.join(LOGS_DIR, "build_go_maps.log"),
    shell:
        "Rscript scripts/03_go_analysis.R "
        "build_maps {input.annotation} "
        "{output.gene_go} {output.mstrg} "
        "> {log} 2>&1"

# ============================================================
# Rule: go_analysis — per-category GO enrichment
# ============================================================
rule go_analysis:
    input:
        bed      = os.path.join(OUT_BASE, "{cat}.bed"),
        gene_go  = os.path.join(OUT_GO, "gene_to_GO.txt"),
        mstrg    = os.path.join(OUT_GO, "mstrg_to_refgene.txt"),
    output:
        refgenes = os.path.join(OUT_GO, "{cat}_refgenes.txt"),
        result   = os.path.join(OUT_GO, "{cat}_BP_GO_enrichment.txt"),
    conda: "envs/goanalysis.yaml"
    resources:
        runtime = 120,
    log:
        os.path.join(LOGS_DIR, "go_{cat}.log"),
    shell:
        "Rscript scripts/03_go_analysis.R "
        "enrich {input.bed} {input.gene_go} {input.mstrg} "
        "{output.refgenes} {output.result} "
        "> {log} 2>&1"

# ============================================================
# Rule: build_kegg_maps — build gene→KEGG map once
# ============================================================
rule build_kegg_maps:
    input:
        annotation = ANNOTATION,
    output:
        gene_kegg     = os.path.join(OUT_KEGG, "gene_to_KEGG.txt"),
        pathway_names = os.path.join(OUT_KEGG, "kegg_pathway_names.txt"),
    conda: "envs/goanalysis.yaml"
    resources:
        runtime = 60,
    log:
        os.path.join(LOGS_DIR, "build_kegg_maps.log"),
    shell:
        "Rscript scripts/04_kegg_analysis.R "
        "build_maps {input.annotation} "
        "{output.gene_kegg} {output.pathway_names} "
        "> {log} 2>&1"

# ============================================================
# Rule: kegg_analysis — per-category KEGG enrichment
# ============================================================
rule kegg_analysis:
    input:
        refgenes      = os.path.join(OUT_GO,   "{cat}_refgenes.txt"),
        gene_kegg     = os.path.join(OUT_KEGG, "gene_to_KEGG.txt"),
        pathway_names = os.path.join(OUT_KEGG, "kegg_pathway_names.txt"),
    output:
        result = os.path.join(OUT_KEGG, "{cat}_KEGG_enrichment.txt"),
    conda: "envs/goanalysis.yaml"
    resources:
        runtime = 120,
    log:
        os.path.join(LOGS_DIR, "kegg_{cat}.log"),
    shell:
        "Rscript scripts/04_kegg_analysis.R "
        "enrich {input.refgenes} {input.gene_kegg} {input.pathway_names} "
        "{output.result} "
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
    conda: "envs/goanalysis.yaml"
    resources:
        runtime = 60,
    log:
        os.path.join(LOGS_DIR, "plots.log"),
    shell:
        "Rscript scripts/05_plots.R > {log} 2>&1"