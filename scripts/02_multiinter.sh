#!/bin/bash
# ============================================================
# Script 02: bedtools multiinter region comparison
#
# Auto-detects which categories to produce based on sample
# names in config/config.yaml:
#   - needle/root specific   (samples ending in N vs R)
#   - cold/drought specific  (2nd char C vs D, stress samples)
#   - embryo specific        (SSE vs SZE, if both present)
#
# Usage:
#   bash scripts/02_multiinter.sh
# Requires:
#   config/config.yaml, BED files from 01_align.sh
#   bedtools in PATH (activate alignment env first)
# ============================================================

set -euo pipefail

CONFIG="config/config.yaml"

# ---- Parse config -------------------------------------------
OUT_BED=$(python3  -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['output']['bed'])")
OUT_BASE=$(python3 -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['output']['base'])")
LOGS=$(python3     -c "import yaml; c=yaml.safe_load(open('$CONFIG')); print(c['slurm']['logs'])")

mapfile -t SAMPLES < <(python3 -c "
import yaml
c = yaml.safe_load(open('$CONFIG'))
for s in c['samples']['names']:
    print(s)
")

mkdir -p "$OUT_BASE" "$LOGS"
exec > >(tee -a "$LOGS/multiinter.log") 2>&1

echo "[$(date)] Samples: ${SAMPLES[*]}"

# ---- Validate BED files -------------------------------------
BED_FILES=()
for s in "${SAMPLES[@]}"; do
    f="${OUT_BED}/${s}.bed"
    [ -f "$f" ] || { echo "ERROR: Missing BED: $f" >&2; exit 1; }
    BED_FILES+=("$f")
done

# ---- Run bedtools multiinter --------------------------------
echo "[$(date)] Running bedtools multiinter"
bedtools multiinter \
    -i "${BED_FILES[@]}" \
    -names "${SAMPLES[@]}" \
    > "${OUT_BASE}/multiinter_output.bed"

echo "[$(date)] multiinter done: $(wc -l < "${OUT_BASE}/multiinter_output.bed") regions"

N=${#SAMPLES[@]}

# ---- Conserved (present in all samples) ---------------------
awk -v n="$N" '$4 == n' "${OUT_BASE}/multiinter_output.bed" \
    > "${OUT_BASE}/conserved.bed"
echo "  conserved:          $(wc -l < "${OUT_BASE}/conserved.bed")"

# ---- Helper: build awk filter expression --------------------
# multiinter cols: chr start end num_samples name_list | sample_cols from col 6
# Build column-index map via python so order always matches config
python3 - "${OUT_BASE}/multiinter_output.bed" "${OUT_BASE}" "${SAMPLES[@]}" << 'PYEOF'
import sys, subprocess, os

bed_file   = sys.argv[1]
out_base   = sys.argv[2]
samples    = sys.argv[3:]

# column index in awk (1-based): col 6 = first sample
col = {s: (i + 6) for i, s in enumerate(samples)}

def awk_filter(present, absent, infile, outfile):
    parts = [f"${col[s]}==1" for s in present] + [f"${col[s]}==0" for s in absent]
    expr  = " && ".join(parts)
    subprocess.run(f"awk '{expr}' {infile} > {outfile}", shell=True, check=True)
    count = sum(1 for _ in open(outfile))
    print(f"  {os.path.basename(outfile):<40} {count}")

# Stress samples (exclude embryo)
stress = [s for s in samples if s not in ("SSE", "SZE")]
needle = [s for s in stress if s.endswith("N")]
root   = [s for s in stress if s.endswith("R")]
cold   = [s for s in stress if len(s) >= 2 and s[1] == "C"]
drought= [s for s in stress if len(s) >= 2 and s[1] == "D"]

# Needle vs Root
if needle and root:
    awk_filter(needle,  root,   bed_file, f"{out_base}/needle_specific.bed")
    awk_filter(root,    needle, bed_file, f"{out_base}/root_specific.bed")
else:
    open(f"{out_base}/needle_specific.bed", "w").close()
    open(f"{out_base}/root_specific.bed",   "w").close()
    print("  needle/root: skipped (not enough tissue groups)")

# Cold vs Drought
if cold and drought:
    awk_filter(cold,    drought, bed_file, f"{out_base}/cold_specific.bed")
    awk_filter(drought, cold,    bed_file, f"{out_base}/drought_specific.bed")
else:
    open(f"{out_base}/cold_specific.bed",    "w").close()
    open(f"{out_base}/drought_specific.bed", "w").close()
    print("  cold/drought: skipped (not enough condition groups)")

# Embryo
if "SSE" in samples and "SZE" in samples:
    awk_filter(["SSE"], ["SZE"], bed_file, f"{out_base}/somatic_embryo_specific.bed")
    awk_filter(["SZE"], ["SSE"], bed_file, f"{out_base}/zygotic_embryo_specific.bed")
else:
    open(f"{out_base}/somatic_embryo_specific.bed", "w").close()
    open(f"{out_base}/zygotic_embryo_specific.bed", "w").close()
    print("  embryo: skipped (SSE and/or SZE not in sample list)")

PYEOF

echo "[$(date)] Done"