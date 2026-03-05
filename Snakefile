#############################################
# CONFIG / PATHS
#############################################
import os, re, glob, yaml, pathlib, hashlib, shutil, csv
from snakemake.io import directory

# Load config: precedence → snakemake.config → config.yaml → env
cfg = {}
if "snakemake" in globals() and hasattr(snakemake, "config") and snakemake.config:
    cfg = dict(snakemake.config)
elif os.path.exists("config.yaml"):
    cfg = yaml.safe_load(open("config.yaml")) or {}

AB1_DIR = os.environ.get("AB1_DIR", cfg.get("ab1_dir", "ab1"))  # raw root
OUTDIR  = os.environ.get("OUTDIR",  cfg.get("outdir",  "tracy_out"))
PANEL   = os.environ.get("PANEL",   cfg.get("panel",   "panel.yaml"))
SAMPLE_SHEET = os.environ.get("SAMPLE_SHEET", cfg.get("sample_sheet", ""))
LOCI_SUBSET  = [x for x in os.environ.get("LOCI_SUBSET", ",".join(cfg.get("loci_subset", []))).split(",") if x]

# Absolute paths
AB1_DIR = os.path.abspath(AB1_DIR)
AB1_SAMPLES_DIR = os.path.join(AB1_DIR, "samples")
AB1_QC_DIR      = os.path.join(AB1_DIR, "qc")

SAFE_AB1_DIR     = os.path.join(AB1_DIR, "safe")
SAFE_SAMPLES_DIR = os.path.join(SAFE_AB1_DIR, "samples")
SAFE_QC_DIR      = os.path.join(SAFE_AB1_DIR, "qc")

OUTDIR  = os.path.abspath(OUTDIR)
PANEL   = os.path.abspath(PANEL)
MANIFEST = os.path.join(OUTDIR, "ab1_manifest.csv")

# Panel / species
panel = yaml.safe_load(open(PANEL))
species = panel.get("species", "homo_sapiens")

# Loci to process
LOCI = [l["locus_id"] for l in panel["loci"]]
if LOCI_SUBSET:
    LOCI = [l for l in LOCI if l in set(LOCI_SUBSET)]

# ------------ helpers ------------
def sanitize_name(base: str) -> str:
    # Minimal sanitization: replace whitespace with underscores
    return re.sub(r"\s+", "_", base)

def detect_control_from_base(base: str) -> str:
    u = base.upper()
    if "NEGATIVE" in u or "_NTC_" in u or u.endswith("_NTC"): return "negative"
    if "POSITIVE" in u or "_POS_" in u: return "positive"
    return "unknown"  # for files under qc/ with unclear names

# Pre-create out dirs
pathlib.Path(os.path.join(OUTDIR, "vcf")).mkdir(parents=True, exist_ok=True)
pathlib.Path(os.path.join(OUTDIR, "logs")).mkdir(parents=True, exist_ok=True)
pathlib.Path(os.path.join(OUTDIR, "fastq")).mkdir(parents=True, exist_ok=True)

#############################################
# DAG / TARGETS
#############################################
rule all:
    input:
        os.path.join(OUTDIR, "hfe_results.csv"),
        os.path.join(OUTDIR, "software_versions.txt"),
        MANIFEST

#############################################
# CHECKPOINT: sanitize_ab1  → copies + classifies + manifest + dirs
#############################################
checkpoint sanitize_ab1:
    """
    Copy raw AB1s from:
      - ab1/samples/ (unsafe names)
      - ab1/qc/      (unsafe names)
    into:
      - ab1/safe/samples/ (sanitized)
      - ab1/safe/qc/      (sanitized)
    and write OUTDIR/ab1_manifest.csv.
    """
    output:
        manifest     = MANIFEST,
        safe_samples = directory(SAFE_SAMPLES_DIR),
        safe_qc      = directory(SAFE_QC_DIR)
    run:
        raw_sample_ab1 = sorted(glob.glob(os.path.join(AB1_SAMPLES_DIR, "*.ab1")))
        raw_qc_ab1     = sorted(glob.glob(os.path.join(AB1_QC_DIR, "*.ab1")))
        raw_list = [(p, "samples") for p in raw_sample_ab1] + [(p, "qc") for p in raw_qc_ab1]

        os.makedirs(SAFE_SAMPLES_DIR, exist_ok=True)
        os.makedirs(SAFE_QC_DIR, exist_ok=True)
        rows, seen = [], set()

        for raw_path, group in raw_list:
            base  = os.path.basename(raw_path)[:-4]
            sbase = sanitize_name(base)
            ctrl  = detect_control_from_base(base) if group == "qc" else "none"
            dst_dir = SAFE_QC_DIR if group == "qc" else SAFE_SAMPLES_DIR
            dst = os.path.join(dst_dir, f"{sbase}.ab1")

            # collision guard
            if dst in seen or (os.path.exists(dst) and not os.path.samefile(raw_path, dst)):
                h = hashlib.sha1(base.encode("utf-8")).hexdigest()[:8]
                sbase = f"{sbase}_{h}"
                dst = os.path.join(dst_dir, f"{sbase}.ab1")

            shutil.copy2(raw_path, dst)
            rows.append([base, sbase, ctrl, group, dst, raw_path])
            seen.add(dst)

        with open(output.manifest, "w", newline="") as fo:
            w = csv.writer(fo)
            w.writerow(["original_base","sanitized","control_type","group","safe_path","raw_path"])
            w.writerows(rows)

#############################################
# Functions that depend on checkpoint's manifest
#############################################
def samples_from_manifest(wildcards):
    manifest = checkpoints.sanitize_ab1.get(**wildcards).output.manifest
    with open(manifest) as f:
        r = csv.DictReader(f)
        return [row["sanitized"] for row in r]

def safe_path_from_manifest(wc):
    manifest = checkpoints.sanitize_ab1.get(**wc).output.manifest
    with open(manifest) as f:
        r = csv.DictReader(f)
        for row in r:
            if row["sanitized"] == wc.sample:
                return row["safe_path"]
    raise ValueError(f"Sample {wc.sample} not found in manifest")


#############################################
# RULES: basecall_fastq / decompose / summarize / software_versions
#############################################
rule basecall_fastq:
    input:
        ab1 = safe_path_from_manifest
    output:
        fq  = os.path.join(OUTDIR, "fastq", "{sample}.fastq"),
        log = os.path.join(OUTDIR, "logs",  "{sample}_basecall.log")
    conda:
        "envs/base.yml"
    threads: 1
    shell:
        r"""
        set -euo pipefail
        tracy basecall -f fastq -o "{output.fq}" "{input.ab1}" > "{output.log}" 2>&1 || true
        """

rule decompose:
    input:
        ab1   = safe_path_from_manifest,
        fasta = lambda w: next(l["reference_fasta"] for l in panel["loci"] if l["locus_id"]==w.locus)
    output:
        bcf = os.path.join(OUTDIR, "vcf", "{sample}_{locus}.bcf")
    log:
        L = os.path.join(OUTDIR, "logs", "{sample}_{locus}_decompose.log")
    conda:
        "envs/base.yml"
    threads: 1
    shell:
        r"""
        set -euo pipefail
        if ! tracy decompose -v -a {species} -r "{input.fasta}" \
              -o "{OUTDIR}/vcf/{wildcards.sample}_{wildcards.locus}" "{input.ab1}" \
              > "{log.L}" 2>&1 ; then
          if grep -qi "sum of the left and right trim size" "{log.L}"; then
            echo "Retrying with --trim 0" >> "{log.L}"
            if ! tracy decompose -v -a {species} --trim 0 -r "{input.fasta}" \
                  -o "{OUTDIR}/vcf/{wildcards.sample}_{wildcards.locus}" "{input.ab1}" \
                  >> "{log.L}" 2>&1 ; then
              echo "[WARN] Tracy failed after --trim 0; writing placeholder BCF." >> "{log.L}"
              : > "{output.bcf}"
            fi
          else
            echo "[WARN] Tracy failed (non-trim); writing placeholder BCF." >> "{log.L}"
            : > "{output.bcf}"
          fi
        fi
        """

rule summarize:
    input:
        manifest = MANIFEST,
        bcfs = lambda wc: expand(
            os.path.join(OUTDIR, "vcf", "{sample}_{locus}.bcf"),
            sample=samples_from_manifest(wc), locus=LOCI
        ),
        panel = PANEL
    output:
        csv = os.path.join(OUTDIR, "hfe_results.csv")
    conda:
        "envs/base.yml"
    threads: 1
    params:
        ab1_dir     = SAFE_AB1_DIR,
        outdir      = OUTDIR,
        loci        = ",".join(LOCI),
        sample_sheet= SAMPLE_SHEET,
        manifest    = MANIFEST
    script:
        "scripts/summarize.py"

rule software_versions:
    output:
        os.path.join(OUTDIR, "software_versions.txt")
    conda:
        "envs/base.yml"
    shell:
        r"""
        set -euo pipefail
        cat > "{output}" << 'EOF'
snakemake: $(snakemake --version 2>/dev/null || echo 'unknown')
tracy: $(tracy --version 2>/dev/null || echo 'unknown')
bcftools: $(bcftools --version 2>/dev/null | head -n1 || echo 'unavailable')
python: $(python --version 2>/dev/null || echo 'unknown')
EOF
        """