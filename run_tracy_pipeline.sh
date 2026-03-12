#!/usr/bin/env bash
# Tracy Sanger Variant Pipeline (HFE H63D / C282Y)
# Author: Adrian Ilich (refactored for CLI)
# Requires: bash >= 4
set -euo pipefail

if [[ -z "${BASH_VERSION-}" ]]; then
  echo "Please run this script with bash (e.g., 'bash run_tracy_pipeline.sh'), not 'source'." 1>&2
  return 1 2>/dev/null || exit 1
fi

########################################
# Defaults (override via CLI)
########################################
AB1_DIR="samples"
OUTDIR="tracy_out"

# Reference FASTAs (one contig each; see README for headers)
REF_H63D="refs/HFE_H63D.clean.fa"
REF_C282Y="refs/HFE_C282Y.clean.fa"

# Loci to process (comma-separated -> array)
LOCI=("C282Y")

# Behavior
SANITIZE_MODE="copy"   # copy|move
DO_BASECALL=true
USE_BCFTOOLS=true

########################################
# Parse CLI args
########################################
print_help() {
  cat <<'EOF'
Usage: run_tracy_pipeline.sh [options]

Core:
  --ab1-dir PATH         Input AB1 directory (default: samples)
  --outdir PATH          Output directory (default: tracy_out)
  --loci LIST            Comma-separated loci (default: C282Y). Example: C282Y,H63D
  --ref-H63D PATH        FASTA path for H63D (default: refs/HFE_H63D.clean.fa)
  --ref-C282Y PATH       FASTA path for C282Y (default: refs/HFE_C282Y.clean.fa)

Behavior:
  --move                 Move (not copy) raw AB1 files into samples/safe/
  --copy                 Copy raw AB1 files into samples/safe/ (default)
  --no-basecall          Skip 'tracy basecall'
  --no-bcftools          Do not convert BCF -> VCF (keep BCF only)
  --help, -h             Show this help and exit
EOF
}

if [[ $# -gt 0 ]]; then
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --ab1-dir)   AB1_DIR="${2:?}"; shift 2 ;;
      --outdir)    OUTDIR="${2:?}"; shift 2 ;;
      --loci)      IFS=',' read -r -a LOCI <<< "${2:?}"; shift 2 ;;
      --ref-H63D)  REF_H63D="${2:?}"; shift 2 ;;
      --ref-C282Y) REF_C282Y="${2:?}"; shift 2 ;;
      --move)      SANITIZE_MODE="move"; shift ;;
      --copy)      SANITIZE_MODE="copy"; shift ;;
      --no-basecall) DO_BASECALL=false; shift ;;
      --no-bcftools) USE_BCFTOOLS=false; shift ;;
      --help|-h)   print_help; exit 0 ;;
      *) echo "Unknown option: $1"; print_help; exit 1 ;;
    esac
  done
fi

########################################
# Setup
########################################
AB1_SAFE_DIR="${AB1_DIR}/safe"
mkdir -p "${OUTDIR}/fastq" "${OUTDIR}/vcf" "${OUTDIR}/logs" "${AB1_SAFE_DIR}"

# bcftools detection (optional)
BCFTOOLS_BIN=""
if $USE_BCFTOOLS && command -v bcftools >/dev/null 2>&1; then
  BCFTOOLS_BIN="$(command -v bcftools)"
fi

# tracy required
if ! command -v tracy >/dev/null 2>&1; then
  echo "ERROR: tracy not found. Activate conda/mamba env with tracy installed."
  exit 1
fi

echo "✔ Tracy pipeline starting..."
echo "Using AB1 directory : ${AB1_DIR}"
echo "AB1 safe directory  : ${AB1_SAFE_DIR}"
echo "Output directory    : ${OUTDIR}"
echo "bcftools            : ${BCFTOOLS_BIN:-not available (BCF only)}"
echo "Selected loci       : ${LOCI[*]}"

# Check references for selected loci
for locus in "${LOCI[@]}"; do
  fasta_var="REF_${locus}"
  fasta="${!fasta_var:-}"
  if [[ -z "${fasta}" ]]; then
    echo "ERROR: No reference variable found for locus ${locus} (expected ${fasta_var})."
    exit 1
  fi
  if [[ ! -s "${fasta}" ]]; then
    echo "ERROR: FASTA for ${locus} not found or empty: ${fasta}"
    exit 1
  fi
  echo " - ${locus} : ${fasta}"
done
echo

########################################
# PREP: sanitize filenames into ${AB1_DIR}/safe/
########################################
shopt -s nullglob nocaseglob
copied=0
for src in "${AB1_DIR}"/*.ab1 "${AB1_DIR}"/*.AB1; do
  base="$(basename "$src")"
  safe_base="$(echo "$base" | sed 's/[^A-Za-z0-9._-]/_/g')"
  safe_base="${safe_base%.*}.ab1"    # normalize extension to lowercase .ab1
  dst="${AB1_SAFE_DIR}/${safe_base}"

  if [[ -e "$dst" ]]; then
    # Collision handling
    i=1
    while [[ -e "${AB1_SAFE_DIR}/${safe_base%.ab1}_dup${i}.ab1" ]]; do
      i=$((i+1))
    done
    echo "WARNING: name collision '${base}' -> '${safe_base}'. Using '${safe_base%.ab1}_dup${i}.ab1'."
    dst="${AB1_SAFE_DIR}/${safe_base%.ab1}_dup${i}.ab1"
  fi

  if [[ "${SANITIZE_MODE}" == "move" ]]; then
    mv -f "$src" "$dst"
  else
    cp -p "$src" "$dst"
  fi
  copied=$((copied+1))
done
shopt -u nocaseglob
echo "[PREP] Placed ${copied} sanitized AB1 file(s) in ${AB1_SAFE_DIR}"

# From here on, only read from AB1_SAFE_DIR
AB1_DIR="${AB1_SAFE_DIR}"

########################################
# MAIN PIPELINE
########################################
processed=0
skipped=0
SKIP_LIST="${OUTDIR}/skipped_samples.txt"
: > "${SKIP_LIST}"

# Helper function: call one locus with retry for trim error
call_one_locus () {
  local ab1="$1"       # path to input AB1
  local safe="$2"      # sanitized sample base
  local locus="$3"     # e.g., H63D
  local fasta="$4"     # path to single-contig FASTA

  local outprefix="${OUTDIR}/vcf/${safe}_${locus}"
  local log="${OUTDIR}/logs/${safe}_${locus}_decompose.log"

  if ! tracy decompose -v -a homo_sapiens -r "${fasta}" -o "${outprefix}" "${ab1}" \
        > "${log}" 2>&1 ; then
    if grep -qi "sum of the left and right trim size is larger than the trace" "${log}" 2>/dev/null; then
      echo "Retrying ${safe} (${locus}) with --trim 0 ..."
      if ! tracy decompose -v -a homo_sapiens --trim 0 -r "${fasta}" -o "${outprefix}" "${ab1}" \
            >> "${log}" 2>&1 ; then
        echo "SKIP: decompose failed for ${safe} (${locus}) even with --trim 0 (see ${log})."
        echo "${safe},${locus}" >> "${SKIP_LIST}"
        return 1
      fi
    else
      echo "SKIP: decompose failed for ${safe} (${locus}) (see ${log})."
      echo "${safe},${locus}" >> "${SKIP_LIST}"
      return 1
    fi
  fi

  # Optional BCF -> VCF
  if [[ -n "${BCFTOOLS_BIN}" ]]; then
    "${BCFTOOLS_BIN}" view "${outprefix}.bcf" -Ov -o "${outprefix}.vcf" || \
      echo "WARNING: bcftools failed on ${outprefix}.bcf; leaving BCF as-is."
  fi
  return 0
}

shopt -s nullglob nocaseglob
found_any=false
for f in "${AB1_DIR}"/*.ab1 "${AB1_DIR}"/*.AB1; do
  [[ -e "$f" ]] || continue
  found_any=true

  fname="${f##*/}"
  fname="${fname%.*}"  # strip extension case-insensitively
  safe="$(echo "$fname" | sed 's/[^A-Za-z0-9._-]/_/g')"

  echo "--------------------------------------------------"
  echo "Processing: $fname  ->  ${safe}"
  echo "--------------------------------------------------"

  # 1) Optional basecall to FASTQ for QC
  if $DO_BASECALL; then
    tracy basecall -f fastq -o "${OUTDIR}/fastq/${safe}.fastq" "$f" \
      > "${OUTDIR}/logs/${safe}_basecall.log" 2>&1 || true
  fi

  any_success=false
  for locus in "${LOCI[@]}"; do
    fasta_var="REF_${locus}"
    fasta="${!fasta_var}"
    if call_one_locus "${f}" "${safe}" "${locus}" "${fasta}"; then
      any_success=true
    fi
  done

  if [[ "${any_success}" != true ]]; then
    echo "SKIP: ${safe} failed all selected loci."
    skipped=$((skipped+1))
    continue
  fi

  echo "✔ Completed: ${safe}"
  processed=$((processed+1))
done
shopt -u nocaseglob
shopt -u nullglob

if [[ "${found_any}" != true ]]; then
  echo "WARNING: No AB1 files found in ${AB1_DIR}."
fi

echo "--------------------------------------------------"
echo "DONE! Processed samples: ${processed}; Samples failed all loci: ${skipped}; Loci skipped (total): $(wc -l < \"${SKIP_LIST}\" | tr -d ' ')"
echo "Outputs in: ${OUTDIR}"
if [[ -s "${SKIP_LIST}" ]]; then
  echo "⚠️  Loci skipped after retry (no usable trace):"
  cat "${SKIP_LIST}"
fi
echo "--------------------------------------------------"

########################################
# SUMMARY CSV (BCF-only)
########################################
# Prereqs: pandas, cyvcf2 (provided by environment.yml)

export PY_AB1_DIR="${AB1_DIR}"
export PY_OUTDIR="${OUTDIR}"
export PY_LOCI="$(IFS=, ; echo "${LOCI[*]}")"

python << 'PYCODE'
import os, re, glob
import pandas as pd
from cyvcf2 import VCF

OUTDIR = os.environ.get("PY_OUTDIR", "tracy_out")
AB1_DIR = os.environ.get("PY_AB1_DIR", "samples/safe")
SEL = [s for s in os.environ.get("PY_LOCI", "H63D").split(",") if s]

# Locus metadata (contig name must match the single-contig FASTA used by tracy decompose)
METADATA = {
    "H63D":  {"contig":"HFE_H63D_hg38_chr6_26090750_26091252",  "local_pos":202, "rsid":"rs1799945"},
    "C282Y": {"contig":"HFE_C282Y_hg38_chr6_26092607_26093122", "local_pos":307, "rsid":"rs1800562"},
}

EXPECTED_POS = {locus: None for locus in SEL}

def sanitize(name: str) -> str:
    return re.sub(r'[^A-Za-z0-9._-]', '_', name)

def parse_from_basename(base: str):
    m = re.match(r'^([A-H]\d{2})_([^_]+)_(.+)$', base)
    if m:
        return m.group(1), m.group(2)
    parts = base.split('_', 2)
    if len(parts) >= 2 and re.match(r'^[A-H]\d{2}$', parts[0]):
        return parts[0], parts[1]
    return None, None

def detect_control(basename: str):
    u = basename.upper()
    if "NEGATIVE" in u or "_NTC_" in u or u.endswith("_NTC"): return "negative"
    if "POSITIVE" in u: return "positive"
    return "none"

def gt_to_zyg(gt: str):
    if gt in {"0/0","0|0"}: return "hom-ref"
    if gt in {"0/1","1/0","0|1","1|0"}: return "het"
    if gt in {"1/1","1|1"}: return "hom-alt"
    return "ambiguous"

# Build sample list from AB1 filenames
samples = []
for apath in sorted(glob.glob(os.path.join(AB1_DIR, "*.ab1"))):
    base = os.path.basename(apath).rsplit(".",1)[0]
    well, sample_id = parse_from_basename(base)
    ctrl = detect_control(base)
    samples.append((base, well, sample_id, ctrl))

records = []
vcf_dir = os.path.join(OUTDIR, "vcf")

for base, well, sample_id, ctrl in samples:
    safe = sanitize(base)
    row = {
        "well": well or "",
        "sample_id": sample_id or "",
        "control_type": ctrl,
        "qc_status": "NA",
        "qc_notes": ""
    }

    for locus in SEL:
        row[f"{locus}_genotype"] = "."
        row[f"{locus}_zygosity"] = "ambiguous"
        row[f"{locus}_mutation_present"] = "UNKNOWN"

    for locus in SEL:
        meta = METADATA.get(locus)
        if not meta:
            row["qc_notes"] += f" missing metadata for {locus};"
            continue

        bcf_path = os.path.join(vcf_dir, f"{safe}_{locus}.bcf")
        if not os.path.exists(bcf_path):
            continue

        try:
            v = VCF(bcf_path)
        except Exception as e:
            row["qc_notes"] += f" failed to open {locus} BCF: {e};"
            continue

        target_records = []
        if meta["contig"] in v.seqnames:
            for rec in v(f"{meta['contig']}:{meta['local_pos']}-{meta['local_pos']}"):
                target_records.append(rec)

        if not target_records and meta.get("rsid"):
            for rec in v:
                if (rec.ID or "").lower() == meta["rsid"]:
                    target_records.append(rec); break

        if target_records:
            rec = target_records[0]
            gt_arr = rec.genotypes[0] if rec.genotypes else None
            gt_str = f"{gt_arr[0]}/{gt_arr[1]}" if gt_arr and None not in gt_arr[:2] else "."
            zyg = gt_to_zyg(gt_str)
            mut = "YES" if zyg in ("het","hom-alt") else ("NO" if zyg=="hom-ref" else "UNKNOWN")
            row[f"{locus}_genotype"] = gt_str
            row[f"{locus}_zygosity"] = zyg
            row[f"{locus}_mutation_present"] = mut
        else:
            row[f"{locus}_genotype"] = "0/0"
            row[f"{locus}_zygosity"] = "hom-ref"
            row[f"{locus}_mutation_present"] = "NO"

    if ctrl == "negative":
        any_mut = any(row.get(f"{l}_mutation_present")=="YES" for l in SEL)
        row["qc_status"] = "FAIL" if any_mut else "PASS"
        row["qc_notes"]  = "NTC shows variant signal" if any_mut else "NTC OK"
    elif ctrl == "positive":
        status="NA"; notes=[]
        for locus in SEL:
            exp = EXPECTED_POS.get(locus)
            if exp:
                got = row[f"{locus}_genotype"].replace("|","/")
                hit = (got == exp)
                status = "PASS" if (status in {"NA","PASS"} and hit) else ("FAIL" if not hit else status)
                if not hit: notes.append(f"{locus} exp {exp}, got {got}")
        if status!="NA":
            row["qc_status"]=status
            row["qc_notes"]="; ".join(notes) if notes else "POS control matched"
        else:
            row["qc_notes"]="POS control present; no expected genotypes configured"

    records.append(row)


def well_key(w):
    if isinstance(w, str) and re.match(r"^[A-H]\d{2}$", w): return (w[0], int(w[1:]))
    return ("Z", 999)

import pandas as pd

df = pd.DataFrame(records)
if not df.empty:
    df = df.sort_values(by=["well","sample_id"], key=lambda c: c.map(well_key) if c.name=="well" else c)
else:
    df = pd.DataFrame(columns=["well","sample_id","control_type","qc_status","qc_notes"])  # empty scaffold
out_path = os.path.join(OUTDIR, "hfe_results.csv")
df.to_csv(out_path, index=False)
print(f"Wrote {out_path} with {len(df)} rows; loci: {', '.join(SEL)}")
PYCODE
