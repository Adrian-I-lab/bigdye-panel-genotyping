# Fast, reproducible Conda environment via Micromamba
FROM mambaorg/micromamba:1.5.8

# Workdir inside the container
WORKDIR /pipeline

# Copy and create environment
COPY environment.yml /pipeline/environment.yml
RUN micromamba create -y -n tracy-pipeline -f /pipeline/environment.yml && \
    micromamba clean --all --yes

# Activate environment for subsequent RUN/CMD (build time)
ENV MAMBA_DOCKERFILE_ACTIVATE=1
SHELL ["micromamba", "run", "-n", "tracy-pipeline", "/bin/bash", "-c"]

# Copy pipeline script
COPY run_tracy_pipeline.sh /pipeline/run_tracy_pipeline.sh

# Create expected dirs (optional)
RUN mkdir -p /pipeline/samples /pipeline/refs /pipeline/tracy_out

# IMPORTANT: ensure runtime uses the env
ENTRYPOINT ["micromamba","run","-n","tracy-pipeline","/bin/bash","/pipeline/run_tracy_pipeline.sh"]
