#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------------------------
# Download hg38 (either whole genome or chr6 only), index with samtools,
# and extract HFE H63D and C282Y windows into single-contig FASTAs whose
# headers match the Tracy pipeline metadata.
#
# Usage:
#   bash get_hg38_and_make_hfe_refs.sh \
#     --out refs \
#     [--source ucsc|ensembl] \
#     [--whole-genome] \
#     [--skip-download]
#
# Defaults:
#   --out refs
#   --source ucsc
#   chr6-only download unless --whole-genome is set
#
# Outputs:
#   refs/HFE_H63D.clean.fa
#   refs/HFE_C282Y.clean.fa
#   (plus downloaded FASTA + .fai in refs/genome/)
#
# Coordinates (hg38):
#   H63D  : chr6:26090750-26091252
#   C282Y : chr6:26092607-26093122
# ------------------------------------------------------------------------------

OUT="refs"
SRC="ucsc"              # ucsc|ensembl  (we implement ucsc chr6-only and ucsc whole)
WHOLE=false
SKIP_DL=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out) OUT="${2:?}"; shift 2 ;;
    --source) SRC="${2:?}"; shift 2 ;;
    --whole-genome) WHOLE=true; shift ;;
    --skip-download) SKIP_DL=true; shift ;;
    -h|--help)
      grep '^# ' "$0" | sed 's/^# \{0,1\}//'
      exit 0
      ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

mkdir -p "${OUT}/genome"

# --------------------------
# Decide download targets
# --------------------------
FASTA_PATH=""
FASTA_URL=""
MD5_URL=""

if [[ "${SRC}" == "ucsc" ]]; then
  if $WHOLE; then
    # Whole-genome FASTA (UCSC hg38.fa.gz)
    # https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/  (hg38.fa.gz)
    FASTA_URL="https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz"
    MD5_URL="https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/md5sum.txt"
    FASTA_PATH="${OUT}/genome/hg38.fa.gz"
  else
    # Minimal: chr6-only FASTA (sufficient for HFE windows)
    # https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/  (chr6.fa.gz)
    FASTA_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr6.fa.gz"
    MD5_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/md5sum.txt"
    FASTA_PATH="${OUT}/genome/chr6.fa.gz"
  fi
else
  echo "Currently implemented source is --source ucsc. If you want Ensembl primary assembly instead, I can add it."
  exit 1
fi

# --------------------------
# Download (unless skipped)
# --------------------------
if ! $SKIP_DL; then
  echo "Downloading FASTA from: ${FASTA_URL}"
  wget -q -O "${FASTA_PATH}" "${FASTA_URL}"

  # Fetch md5 list and verify (best-effort)
  if [[ -n "${MD5_URL}" ]]; then
    echo "Fetching md5 list for verification (best-effort): ${MD5_URL}"
    wget -q -O "${OUT}/genome/md5sum.txt" "${MD5_URL}" || true
    # Try to verify if present
    if [[ -s "${OUT}/genome/md5sum.txt" ]]; then
      echo "Verifying md5..."
      fname="$(basename "${FASTA_PATH}")"
      # md5sum.txt contains lines like: <md5>  <filename>
      expected_md5="$(grep -E "[[:space:]]${fname}$" "${OUT}/genome/md5sum.txt" | awk '{print $1}' || true)"
      if [[ -n "${expected_md5}" ]]; then
        # Compute md5 of the downloaded file
        if command -v md5sum >/dev/null 2>&1; then
          actual_md5="$(md5sum "${FASTA_PATH}" | awk '{print $1}')"
        elif command -v md5 >/dev/null 2>&1; then
          actual_md5="$(md5 -q "${FASTA_PATH}")"
        else
          actual_md5=""
        fi
        if [[ -n "${actual_md5}" && "${actual_md5}" != "${expected_md5}" ]]; then
          echo "WARNING: md5 mismatch for ${fname}. Expected ${expected_md5}, got ${actual_md5}."
        else
          echo "md5 OK for ${fname}"
        fi
      else
        echo "md5 entry for ${fname} not found in md5sum.txt (continuing)."
      fi
    fi
  fi
else
  echo "Skipping download as requested (--skip-download). Expecting ${FASTA_PATH} to exist."
fi

# --------------------------
# Prepare FASTA & index
# --------------------------
if [[ ! -s "${FASTA_PATH}" ]]; then
  echo "ERROR: FASTA not found at ${FASTA_PATH}"
  exit 1
fi

# Uncompress if needed
if [[ "${FASTA_PATH}" == *.gz ]]; then
  echo "Uncompressing FASTA..."
  gunzip -c "${FASTA_PATH}" > "${OUT}/genome/$(basename "${FASTA_PATH%.gz}")"
  FASTA_PATH="${OUT}/genome/$(basename "${FASTA_PATH%.gz}")"
fi

# Index with samtools
if ! command -v samtools >/dev/null 2>&1; then
  echo "ERROR: samtools not found on PATH. Please install samtools."
  exit 1
fi

echo "Indexing FASTA with samtools faidx..."
samtools faidx "${FASTA_PATH}"

# --------------------------
# Extract HFE windows
# --------------------------
# Coordinates (hg38):
#   H63D  : chr6:26090750-26091252
#   C282Y : chr6:26092607-26093122
#
# Output headers must match the pipeline's metadata keys:
#   >HFE_H63D_hg38_chr6_26090750_26091252
#   >HFE_C282Y_hg38_chr6_26092607_26093122

mkdir -p "${OUT}"

# Helper to extract and reheader a window into a single-contig FASTA
extract_window () {
  local region="$1"   # e.g., chr6:26090750-26091252
  local name="$2"     # e.g., HFE_H63D_hg38_chr6_26090750_26091252
  local outfa="${OUT}/${name}.clean.fa"

  echo "Extracting ${region} -> ${outfa}"
  # samtools faidx prints header like >chr6:26090750-26091252
  samtools faidx "${FASTA_PATH}" "${region}" \
    | sed "1s|^>.*$|>${name}|" > "${outfa}"

  # Quick sanity
  if [[ ! -s "${outfa}" ]]; then
    echo "ERROR: Failed to create ${outfa}"
    exit 1
  fi
}

extract_window "chr6:26090750-26091252" "HFE_H63D_hg38_chr6_26090750_26091252"
extract_window "chr6:26092607-26093122" "HFE_C282Y_hg38_chr6_26092607_26093122"

echo
echo "Done. Created:"
echo "  ${OUT}/HFE_H63D.clean.fa"
echo "  ${OUT}/HFE_C282Y.clean.fa"
echo
echo "Genome FASTA & index at:"
echo "  ${FASTA_PATH}"
echo "  ${FASTA_PATH}.fai"
