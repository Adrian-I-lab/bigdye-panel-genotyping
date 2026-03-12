
# Tracy Sanger Variant Pipeline (HFE H63D / C282Y)

A lightweight pipeline to call **HFE** variants (H63D, C282Y) from **Sanger AB1** traces using:
- `tracy decompose` for per-locus BCFs (and optional VCFs)
- an embedded Python summarizer (cyvcf2 + pandas) that reports genotype/zygosity per locus

The pipeline:
1. **Sanitizes** raw AB1 filenames to avoid spaces/special chars.
2. **Basecalls** to FASTQ (optional, for QC).
3. **Decomposes** each AB1 against single-contig locus FASTA(s).
4. Produces per-sample **BCF** (and optional VCF via bcftools).
5. Writes a summary CSV: `tracy_out/hfe_results.csv`.

---

## Quick start

### Option A: Conda/Micromamba (local)
```bash
micromamba create -f environment.yml -n tracy-pipeline
micromamba activate tracy-pipeline

# Put input AB1 files under ./samples and references under ./refs
bash run_tracy_pipeline.sh   --ab1-dir samples   --outdir tracy_out   --loci C282Y,H63D   --ref-H63D refs/HFE_H63D.clean.fa   --ref-C282Y refs/HFE_C282Y.clean.fa
```

### Option B: Docker (most portable)
```bash
# Build
docker build -t tracy-pipeline:latest .

# Run (mount your local folders into the container)
docker run --rm   -v "$(pwd)/samples:/pipeline/samples"   -v "$(pwd)/refs:/pipeline/refs"   -v "$(pwd)/tracy_out:/pipeline/tracy_out"   tracy-pipeline:latest   --ab1-dir samples   --outdir tracy_out   --loci C282Y,H63D   --ref-H63D refs/HFE_H63D.clean.fa   --ref-C282Y refs/HFE_C282Y.clean.fa
```

> **Outputs:** `tracy_out/fastq`, `tracy_out/vcf`, `tracy_out/logs`, and `tracy_out/hfe_results.csv`.

---

## Directory structure

```
.
├── run_tracy_pipeline.sh
├── environment.yml
├── Dockerfile
├── samples/               # your *.ab1 / *.AB1
├── refs/
│   ├── HFE_H63D.clean.fa  # >HFE_H63D_hg38_chr6_26090750_26091252
│   └── HFE_C282Y.clean.fa # >HFE_C282Y_hg38_chr6_26092607_26093122
└── tracy_out/             # created at runtime
```

> Each reference FASTA must be **single-contig** and its header must match the summarizer metadata:
> - H63D  → `>HFE_H63D_hg38_chr6_26090750_26091252`
> - C282Y → `>HFE_C282Y_hg38_chr6_26092607_26093122`

---

## CLI options

```text
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
  --help, -h             Show help and exit
```

---

## Results

- Per-locus BCF at: `tracy_out/vcf/{sample}_{locus}.bcf`
- Optional VCF if bcftools available (or not disabled)
- Logs: `tracy_out/logs/`
- Summary CSV: `tracy_out/hfe_results.csv` with per-locus genotype/zygosity and control QC flags

---

## Troubleshooting

- **No AB1s processed?** Ensure files exist under `samples/`. Uppercase `.AB1` is supported.
- **Empty calls?** Confirm the FASTA contig headers match the Python metadata (above).
- **Missing cyvcf2/pandas?** Use the provided `environment.yml` or Docker image.
- **HPC without Docker?** Build Singularity from Docker:  
  `singularity build tracy-pipeline.sif docker-daemon://tracy-pipeline:latest`

---

## License / citation
- `tracy` by Rausch et al. (Bioconda distribution)
- `cyvcf2` (Bioconda / Conda-Forge)
- `bcftools` / `htslib` (Bioconda)
