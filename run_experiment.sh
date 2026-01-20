#!/bin/bash
#SBATCH --job-name=MPRA
#SBATCH --output=/home/itg/oleg.vlasovets/jobs/output/slurm_mpra_experiment_%j.txt
#SBATCH --error=/home/itg/oleg.vlasovets/jobs/error/slurm_error_experiment_%j.txt
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=cpu_p
#SBATCH --qos=cpu_preemptible
#SBATCH --time=05:00:00

# =======================
# MPRA Snakemake launcher
# =======================

# Create output dirs if missing
mkdir -p "$HOME/jobs/output"
mkdir -p "$HOME/.cache/snakemake/snakemake"

# Set a job-scoped temporary directory (good for Apptainer)
export TMPDIR=/lustre/scratch/users/$USER/$SLURM_JOB_ID
mkdir -p "$TMPDIR"

# 1) Activate your conda 'snakemake' environment
#    (Adjust path if your Miniconda/Mamba is elsewhere)
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate snakemake_env

# 2) Optional: where Snakemake stores per-rule conda envs
export SNAKEMAKE_CONDA_PREFIX="$HOME/.conda/envs/smk_envs"
mkdir -p "$SNAKEMAKE_CONDA_PREFIX"

# 3) Sanity checks (apptainer + snakemake)
command -v apptainer >/dev/null 2>&1 || { echo "ERROR: apptainer not found in PATH"; exit 1; }
command -v snakemake >/dev/null 2>&1 || { echo "ERROR: snakemake not found (conda env?)"; exit 1; }

# 3.5) Change to the working directory where results should be created
cd /lustre/groups/itg/teams/zeggini/projects/GO2/MPRA/mpra_test/

snakemake --unlock --snakefile /home/itg/oleg.vlasovets/projects/MPRAsnakeflow/workflow/Snakefile --configfile /lustre/groups/itg/teams/zeggini/projects/GO2/MPRA/mpra_test/config_experiment.yaml

# 4) Run assignment with your Slurm profile
snakemake \
  --profile /home/itg/oleg.vlasovets/.snakemake_profile \
  --snakefile  /home/itg/oleg.vlasovets/projects/MPRAsnakeflow/workflow/Snakefile \
  --configfile /lustre/groups/itg/teams/zeggini/projects/GO2/MPRA/mpra_test/config_experiment.yaml \
  --software-deployment-method conda apptainer \
  --apptainer-args "-B $TMPDIR -B $HOME/.cache/snakemake/snakemake/ -B /lustre/groups/itg/teams/zeggini/projects/GO2/MPRA/mpra_test/" \
  --apptainer-prefix "$HOME/apptainer_images" \
  --jobs 100 \
  --cores 16