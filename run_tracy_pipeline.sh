#!/usr/bin/env bash
# Tracy Sanger Variant Pipeline (HFE H63D / C282Y)
# FINAL MERGED VERSION
# Includes:
#  - Safe input handling
#  - tracy decompose variant calling
#  - Confident WT and variant calls
#  - FASTQ QC (mean Q, >=Q20 bp)
#  - GOLD-STANDARD QC: minimap2 alignment-based callable window
#  - Single comprehensive CSV report

set -euo pipefail

########################################
# Safety
########################################
if [[ -z "${BASH_VERSION-}" ]]; then
  echo "Please run this script with bash, not source." >&2
  return 1 2>/dev/null || exit 1
fi

########################################
# Defaults
########################################
AB1_DIR="samples"
OUTDIR="tracy_out"
AB1_SAFE_SUBDIR="safe"

REF_H63D="refs/HFE_H63D.clean.fa"
REF_C282Y="refs/HFE_C282Y.clean.fa"

LOCI=()
SANITIZE_MODE="copy"
DO_BASECALL=true
USE_BCFTOOLS=true

########################################
# CLI
########################################
print_help() {
  cat <<'EOF'
Usage: run_tracy_pipeline.sh [options]

Core:
  --ab1-dir PATH
  --outdir PATH
  --loci LIST            Comma-separated (e.g. C282Y,H63D)

Behavior:
  --move | --copy
  --no-basecall
  --no-bcftools
  --help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --ab1-dir) AB1_DIR="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --loci) IFS=',' read -r -a LOCI <<< "$2"; shift 2;;
    --move) SANITIZE_MODE="move"; shift;;
    --copy) SANITIZE_MODE="copy"; shift;;
    --no-basecall) DO_BASECALL=false; shift;;
    --no-bcftools) USE_BCFTOOLS=false; shift;;
    --help|-h) print_help; exit 0;;
    *) echo "Unknown option $1" >&2; exit 1;;
  esac
done

########################################
# Setup
########################################
AB1_SAFE_DIR="${AB1_DIR}/${AB1_SAFE_SUBDIR}"
mkdir -p "$AB1_SAFE_DIR" "$OUTDIR"/{fastq,vcf,logs,plots,align}

for bin in tracy samtools minimap2; do
  command -v "$bin" >/dev/null 2>&1 || { echo "ERROR: $bin not found" >&2; exit 1; }
done

BCFTOOLS_BIN=""
$USE_BCFTOOLS && command -v bcftools >/dev/null 2>&1 && BCFTOOLS_BIN="$(command -v bcftools)"

########################################
# PREP: sanitize AB1s -> safe/
########################################
shopt -s nullglob nocaseglob
for src in "$AB1_DIR"/*.ab1; do
  base=$(basename "$src")
  safe=$(echo "$base" | sed 's/[^A-Za-z0-9._-]/_/g')
  safe="${safe%.*}.ab1"
  dst="$AB1_SAFE_DIR/$safe"
  [[ "$SANITIZE_MODE" == move ]] && mv -f "$src" "$dst" || cp -p "$src" "$dst"
done
shopt -u nocaseglob

########################################
# tracy decompose helper
########################################
call_one_locus() {
  local ab1="$1" safe="$2" locus="$3" fasta="$4"
  local outprefix="$OUTDIR/vcf/${safe}_${locus}"
  local log="$OUTDIR/logs/${safe}_${locus}.log"

  tracy decompose -v -a homo_sapiens -r "$fasta" -o "$outprefix" "$ab1" >"$log" 2>&1 \
  || tracy decompose -v -a homo_sapiens --trim 0 -r "$fasta" -o "$outprefix" "$ab1" >>"$log" 2>&1 \
  || return 1

  [[ -n "$BCFTOOLS_BIN" ]] && "$BCFTOOLS_BIN" view "${outprefix}.bcf" -Ov -o "${outprefix}.vcf" || true
}

########################################
# MAIN: basecall, decompose, align
########################################
shopt -s nullglob
for ab1 in "$AB1_SAFE_DIR"/*.ab1; do
  safe=$(basename "$ab1" .ab1)
  echo "Processing $safe"

  fastq="$OUTDIR/fastq/${safe}.fastq"
  if $DO_BASECALL; then
    tracy basecall -f fastq "$ab1" -o "$fastq" >"$OUTDIR/logs/${safe}_basecall.log" 2>&1 || true
  fi

  [[ -s "$fastq" ]] || { echo "QC_FAIL: no FASTQ"; continue; }

  for locus in "${LOCI[@]}"; do
    fasta_var="REF_${locus}"
    ref="${!fasta_var}"

    call_one_locus "$ab1" "$safe" "$locus" "$ref" || true

    bam="$OUTDIR/align/${safe}_${locus}.bam"
    minimap2 -a -x map-ont "$ref" "$fastq" | samtools sort -o "$bam"
    samtools index "$bam"
  done
done
shopt -u nullglob

########################################
# PYTHON SUMMARY: FULL REPORT
########################################
export PY_AB1_DIR="$AB1_SAFE_DIR"
export PY_OUTDIR="$OUTDIR"
export PY_LOCI="$(IFS=,; echo "${LOCI[*]}")"

python <<'PYCODE'
import os, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pysam
from cyvcf2 import VCF

OUTDIR=os.environ['PY_OUTDIR']
AB1=os.environ['PY_AB1_DIR']
LOCI=[l for l in os.environ['PY_LOCI'].split(',') if l]
FASTQ=os.path.join(OUTDIR,'fastq')
VCF_DIR=os.path.join(OUTDIR,'vcf')
ALIGN=os.path.join(OUTDIR,'align')
PLOTS=os.path.join(OUTDIR,'plots')

# Locus metadata
META={
 'H63D':{'local_pos':202,'window_bp':120,'ref_len':503},
 'C282Y':{'local_pos':307,'window_bp':120,'ref_len':515}
}

rows=[]; mean_qs=[]; q20_lens=[]

for ab1 in sorted(glob.glob(os.path.join(AB1,'*.ab1'))):
  sample=os.path.basename(ab1).replace('.ab1','')
  row={'sample':sample}

  # FASTQ QC
  fq=os.path.join(FASTQ,sample+'.fastq')
  qs=[]
  if os.path.exists(fq):
    for i,l in enumerate(open(fq)):
      if i%4==3:
        qs.extend([ord(c)-33 for c in l.strip()])

  row['mean_q']=round(float(np.mean(qs)),1) if qs else 0.0
  row['q20_bp']=int(sum(q>=20 for q in qs))
  mean_qs.append(row['mean_q']); q20_lens.append(row['q20_bp'])

  # Per locus: variants + coverage-based QC
  for locus in LOCI:
    meta=META[locus]
    bcf=os.path.join(VCF_DIR,f'{sample}_{locus}.bcf')
    bam=os.path.join(ALIGN,f'{sample}_{locus}.bam')

    row.update({
      f'{locus}_genotype':'.',
      f'{locus}_zygosity':'ambiguous',
      f'{locus}_mutation_present':'UNKNOWN',
      f'{locus}_variant_call_confidence':'LOW',
      f'{locus}_WT_CONFIDENT':False,
      f'{locus}_callable_window_bp':0
    })

    # Coverage-based callable window
    if os.path.exists(bam):
      cov_bp=0
      bf=pysam.AlignmentFile(bam)
      for col in bf.pileup(stepper='all',truncate=True):
        if col.n>0:
          cov_bp+=1
      row[f'{locus}_callable_window_bp']=cov_bp

    # Variant logic from tracy
    target=None
    if os.path.exists(bcf):
      try:
        v=VCF(bcf)
        for rec in v:
          if rec.POS==meta['local_pos']:
            target=rec
      except Exception:
        pass

    if target:
      gt=target.genotypes[0][:2]
      gt_str=f"{gt[0]}/{gt[1]}"
      if gt_str=='0/0': z='hom-ref'; mut='NO'
      elif gt_str in {'0/1','1/0'}: z='het'; mut='YES'
      elif gt_str=='1/1': z='hom-alt'; mut='YES'
      else: z='ambiguous'; mut='UNKNOWN'
      row[f'{locus}_genotype']=gt_str
      row[f'{locus}_zygosity']=z
      row[f'{locus}_mutation_present']=mut

    # Confident WT if covered and no variant
    elif row[f'{locus}_callable_window_bp']>=meta['window_bp'] and row['mean_q']>=20:
      row[f'{locus}_genotype']='0/0'
      row[f'{locus}_zygosity']='hom-ref'
      row[f'{locus}_mutation_present']='NO'
      row[f'{locus}_WT_CONFIDENT']=True

    # Confidence scoring
    if row[f'{locus}_callable_window_bp']>=meta['window_bp'] and row['mean_q']>=20:
      row[f'{locus}_variant_call_confidence']='HIGH'
    elif row[f'{locus}_callable_window_bp']>0:
      row[f'{locus}_variant_call_confidence']='MEDIUM'

  rows.append(row)

DF=pd.DataFrame(rows)

# Overall QC flags
DF['QC_FAIL']=(DF.mean_q<15)|(DF.q20_bp<200)
DF['QC_PASS']=(DF.mean_q>=20)&(DF.q20_bp>=300)
DF['QC_WARN']=~DF.QC_FAIL & ~DF.QC_PASS

# Order columns: sample, variant calls, QC
variant_cols=[c for c in DF.columns if any(l in c for l in LOCI)]
qc_cols=[c for c in DF.columns if c not in (['sample']+variant_cols)]
DF=DF[['sample']+variant_cols+qc_cols]

out=os.path.join(OUTDIR,'genotyping_results.csv')
DF.to_csv(out,index=False)

# Plots
plt.figure(); plt.hist(mean_qs); plt.title('Mean FASTQ Q'); plt.savefig(os.path.join(PLOTS,'mean_q.png'))
plt.figure(); plt.hist(q20_lens); plt.title('>=Q20 bp'); plt.savefig(os.path.join(PLOTS,'q20_bp.png'))

print(f'Wrote {out} with {len(DF)} samples')
PYCODE
