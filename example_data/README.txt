Place a your .ab1 files here.

Set up docker container:

# have docker downloaded and account made first
docker pull adriani25/bigdye-panel:0.1.0

Run this in terminal: 

docker run --rm \
  -v "Path to the repository folder"/bigdye-panel-genotyping:/work \
  -v "Path to the repository folder"/bigdye-panel-genotyping/refs:/work/refs \
  -w /work \
  -e AB1_DIR="example_data/ab1" \
  -e OUTDIR="example_data/out" \
  -e PANEL="panel.yaml" \
  bigdye-panel:latest \
  snakemake -j 4