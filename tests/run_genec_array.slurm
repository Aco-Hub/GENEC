#!/bin/bash -l
#SBATCH --job-name=genec_array
#SBATCH --partition=public-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=2-00:00:00
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=quentin.gallien@epfl.ch

module purge
module load GCC/13.2.0
module load AstroPython/3.11.5   # si nécessaire

source /srv/scratch/qgallien/PDM/codes/environments/venv1quentinPy3P11/bin/activate || true

export PATH="/srv/scratch/qgallien/PDM/codes/GENEC/bin:$PATH"
export GENEC_DEFAULT_PROGRAM="/home/aiaiso/GENEC/bin/genec"
export GENEC_EMAIL_ADDRESS="test@example.com"

cd "$SLURM_SUBMIT_DIR"

LIST_FILE="${1:-/srv/scratch/qgallien/PDM/codes/Stars/folders_names.txt}"

if [[ ! -f "$LIST_FILE" ]]; then
  echo "Liste introuvable: $LIST_FILE" >&2
  exit 2
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "SLURM_ARRAY_TASK_ID non défini (lance en job array)" >&2
  exit 3
fi

# Récupère la ligne correspondant à la tâche courante
LINE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST_FILE" | sed 's/[[:space:]]*$//')"
STAR_DIR="$(echo "$LINE" | cut -d ';' -f 1)"
STAR_NAME="$(echo "$LINE" | cut -d ';' -f 2)"

if [[ -z "$STAR_DIR" ]]; then
  echo "Ligne ${SLURM_ARRAY_TASK_ID} vide dans $LIST_FILE" >&2
  exit 4
fi
if [[ ! -d "$STAR_DIR" ]]; then
  echo "Dossier inexistant: $STAR_DIR (depuis $LIST_FILE ligne ${SLURM_ARRAY_TASK_ID})" >&2
  exit 5
fi

echo "[$(date)] Tâche ${SLURM_ARRAY_TASK_ID} → $STAR_NAME ($STAR_DIR) depuis $LIST_FILE"

# (Optionnel) Rediriger les logs *dans* le dossier étoile
# exec > >(tee -a "$STAR_DIR/slurm_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out")
# exec 2> >(tee -a "$STAR_DIR/slurm_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" >&2)

# 2) Aller dans le dossier étoile
cd "$STAR_DIR"

# 3) Choix du fichier .rot/.com à utiliser
INI_FILE=""
if [[ -f "ini_${STAR_NAME}.rot" ]]; then
  INI_FILE="ini_${STAR_NAME}.rot"
elif ls -1 *.rot >/dev/null 2>&1; then
  INI_FILE="$(ls -1 *.rot | head -n1)"
elif ls -1 *.com >/dev/null 2>&1; then
  INI_FILE="$(ls -1 *.com | head -n1)"
fi

if [[ -z "${INI_FILE:-}" || ! -f "$INI_FILE" ]]; then
  echo "Aucun fichier .rot/.com trouvé dans $STAR_DIR"; exit 6
fi

echo "Fichier d'entrée: $INI_FILE"

# 4) Lancer GENEC: Benchmark des différents modes (mesure de vitesse)
BENCHMARK_FILE="$SLURM_SUBMIT_DIR/benchmark_results_${SLURM_ARRAY_TASK_ID}.csv"
echo "Star,Mode,Time_seconds,LogT_core,LogRho_core,Mu_core,HII_core" > "$BENCHMARK_FILE"

# List of modes to benchmark: (name, backend, precision)
MODES=(
  "fortran eager fp64"
  "cpu eager fp64"
  "gpu_eager eager fp64"
  "gpu_compile compile fp64"
  "gpu_triton triton fp64"
  "gpu_fp32 eager fp32"
)

for m in "${MODES[@]}"; do
  read -r name backend precision <<< "$m"

  case "$name" in
    fortran) run_mode="fortran" ;;
    cpu)     run_mode="cpu" ;;
    *)       run_mode="gpu" ;;
  esac

  echo "=== Lancement test de vitesse: $STAR_NAME mode=$name (backend=$backend, prec=$precision) ==="

  CALC_DIR="${STAR_DIR}/calc_${name}"
  mkdir -p "$CALC_DIR"

  start_time=$(date +%s.%N)

  python /home/aiaiso/GENEC/tools/GENEC_launch.py \
    "$STAR_NAME" -i "$INI_FILE" -c "$CALC_DIR" \
    --mode "$run_mode" --backend "$backend" --precision "$precision"

  end_time=$(date +%s.%N)
  elapsed=$(awk "BEGIN {printf \"%.4f\", $end_time - $start_time}")

  # Extract physical values
  PHYS_VALS=$(python /home/aiaiso/GENEC/tools/extract_physics.py "$STAR_NAME" "$CALC_DIR")

  echo "Terminé: $STAR_NAME mode=$name en $elapsed secondes. Physics: $PHYS_VALS"
  echo "$STAR_NAME,$name,$elapsed,$PHYS_VALS" >> "$BENCHMARK_FILE"
done

echo "=== Benchmark complet pour $STAR_NAME ==="
