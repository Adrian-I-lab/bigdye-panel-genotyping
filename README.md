
# BigDye Sanger Panel (Tracy + Snakemake)

Reproducible, config-driven pipeline for calling targeted Sanger (BigDye) variants with **Tracy** and summarising per sample/locus. Built for easy use on **Windows/macOS** with a single **Docker** container.

## Quick start (90 seconds)

### Option B — Docker (Windows/macOS with Docker Desktop)
1. Build the container (one-time):
   ```bash
   docker pull adriani25/bigdye-panel:0.1.0
   ```
2. Run via the wrapper (no local installs needed):
   ```bash
   docker run --rm \
      -v "Path to the repository folder"/bigdye-panel-genotyping:/work \
      -v "Path to the repository folder"/bigdye-panel-genotyping/refs:/work/refs \
      -w /work \
      -e AB1_DIR="example_data/ab1" \
      -e OUTDIR="example_data/out" \
      -e PANEL="panel.yaml" \
      bigdye-panel:latest \
      snakemake -j 4
   ```

## Inputs & outputs
- **Inputs**: 
  - `*.ab1` Sanger trace files (.ab1 files) in folders samples/ for your sample files and qc samples in the qc/ folder; 
  - A **panel.yaml** describing each locus/amplicon.
  - reference fasta files in the **refs/** folder for each of the mentioned locus/amplicon in **panel.yaml**
- **Outputs** (in `outdir`):
  - `vcf/*.bcf` (and `.vcf` if bcftools available)
  - `fastq/*.fastq` (basecalls for QC)
  - `logs/*_decompose.log`
  - `hfe_results.csv` (summary across selected loci)
  - `ab1_manifest.csv` (summary of all files and any renamed sample files)

## Adding loci / genes (no code changes)
1. Put a **single-contig FASTA** in `refs/` with a header you control (the `contig`).
2. Add an entry to `panel.yaml` in following format:
   ```yaml
   - locus_id: S65C
     reference_fasta: refs/HFE_S65C.clean.fa
     contig: HFE_S65C_hg38_chr6_...
     rsid: rsXYZ         # or omit and use local_pos
     local_pos: 123      # 1-based within the contig
     expected_genotype: null
   ```
3. Re-run `snakemake`.

## Selecting a subset of loci
Any loci mentioned in the **panel.yaml** file will be searched. 
Its recommemded to only do one at a time if you have multiple samples unless there are multiple mutations covered in the amplified sequence.

## Optional sample sheet
Provide a CSV with columns: `well,sample_id,control_type,expected_genotypes`

## What’s inside
- **Tracy** for basecalling and variant decompose (with an auto-retry `--trim 0`).
- **Snakemake** orchestration (parallel, reproducible, with conda per-rule envs).
- **Config-first** design (`panel.yaml`, `config.yaml`).
- **Summariser** (`scripts/summarize.py`) producing `hfe_results.csv` with per-locus genotype, zygosity, mutation flag, and locus status.

## Versioning & provenance
The pipeline writes `software_versions.txt` capturing Snakemake/Tracy/bcftools/Python versions. The environment is pinned in `envs/base.yml`.

## License
MIT — see `LICENSE`.
