# run with the command : srun --partition=public-cpu python gen_stars_with_GENEC.py


import argparse
import csv
import shutil
import subprocess
import sys
from decimal import Decimal, InvalidOperation
from pathlib import Path

import matplotlib

# le but est de générer les étoiles que je veux via un fichier .csv


# === adapte ces colonnes à TON CSV ===
COL_NAME = "name_star"
COL_SAVE = "saving_folder_path"
COL_MASS = "Mass_M_"  # ou "mass"
COL_Z = "Metalicity_Z_"  # ou "Metalicity"
COL_ROT = "Rotation_S_"  # ou "Rotation"
COL_ENDPHASE = "End_Phase"  # ou "EndPhase"
COL_DOVPH = "dovph"  # ou "DoVph"


def find_in_parents(
    start_path: Path, relative_target: str, max_levels: int = 10
) -> Path:
    start_path = Path(start_path).resolve()

    current = start_path
    for _ in range(max_levels + 1):
        candidate = current / relative_target
        if candidate.exists():
            return candidate.resolve()

        if current.parent == current:
            break
        current = current.parent

    raise FileNotFoundError(
        f"impossible to find 'GENEC' folder by going up from '{start_path}' within {max_levels} levels."
    )


def run_genec_makeini(
    star_dir: Path,
    name_star: str,
    mass: float,
    z: float,
    rotation: float,
    end_phase=None,
    dovph=None,
    max_parent_search: int = 7,
):
    star_dir = Path(star_dir).resolve()
    star_dir.mkdir(parents=True, exist_ok=True)
    exe = find_in_parents(
        star_dir, "GENEC/bin/genec_makeini", max_levels=max_parent_search
    )
    if not exe.is_file():
        raise FileNotFoundError(f"Executable introuvable: {exe}")
    answers_list = [str(name_star), f"{mass} {z}", str(rotation), "0"]

    if end_phase is not None:
        answers_list.extend(["1", "2", str(end_phase), "0"])

    if dovph is not None:
        answers_list.extend(["7", "3", str(dovph), "0"])

    answers_list.append("0")
    answers = "\n".join(answers_list)
    print("=== INPUT GENEC ===")
    print(answers)
    print("===================")

    res = subprocess.run(
        [str(exe)],
        cwd=star_dir,
        input=answers,
        text=True,
        capture_output=True,
        check=True,
    )

    print(res.stdout)


def clean_cell(value):
    # Convertit une cellule CSV en string propre.
    if value is None:
        return ""
    return str(value).strip()


def parse_decimal(value, column_name):
    # Convertit une valeur CSV en Decimal. Accepte les virgules ou points comme séparateur décimal.
    s = clean_cell(value).replace(",", ".")
    if s == "":
        raise ValueError(f"Valeur vide dans la colonne '{column_name}'.")
    try:
        return Decimal(s)
    except InvalidOperation as e:
        raise ValueError(
            f"Valeur invalide dans la colonne '{column_name}': {value!r}"
        ) from e


def decimal_parts_as_strings(value):
    # Retourne (partie_entiere, partie_decimale) sous forme de strings, sans signe et sans perte des zéros significatifs déjà présents.
    s = format(value, "f")
    if "." in s:
        int_part, frac_part = s.split(".", 1)
    else:
        int_part, frac_part = s, ""
    return int_part.lstrip("-"), frac_part.rstrip("0")


def format_mass(mass):
    if mass < 0:
        raise ValueError(f"La masse doit être positive, reçu: {mass}")
    int_part, frac_part = decimal_parts_as_strings(mass)
    int_part = int_part.zfill(2)
    frac_part = frac_part if frac_part else "00"
    if len(frac_part) < 2:
        frac_part = frac_part.ljust(2, "0")

    return f"M{int_part}P{frac_part}"


def format_metallicity(z):
    if z < 0 or z >= 1:
        raise ValueError(f"La métallicité doit être dans [0, 1), reçu: {z}")
    _, frac_part = decimal_parts_as_strings(z)
    frac_part = frac_part if frac_part else "0"
    if len(frac_part) < 3:
        frac_part = frac_part.ljust(3, "0")

    return f"Z{frac_part}"


def format_rotation(rot):
    if rot <= -1:
        raise ValueError(
            f"La rotation négative doit être strictement supérieure à -1, reçu: {rot}"
        )

    if -1 < rot < 0:
        _, frac_part = decimal_parts_as_strings(abs(rot))
        frac_part = frac_part if frac_part else "00"
        if len(frac_part) < 2:
            frac_part = frac_part.ljust(2, "0")
        return f"S{frac_part}"

    if 0 <= rot < 1:
        _, frac_part = decimal_parts_as_strings(rot)
        frac_part = frac_part if frac_part else "00"
        if len(frac_part) < 2:
            frac_part = frac_part.ljust(2, "0")
        return f"s{frac_part}"

    return f"R{int(rot)}"


def format_dovph(dovph):
    if dovph < 0 or dovph >= 1:
        raise ValueError(f"Dovph doit être dans [0, 1), reçu: {dovph}")

    _, frac_part = decimal_parts_as_strings(dovph)
    frac_part = frac_part if frac_part else "0"
    if len(frac_part) < 3:
        frac_part = frac_part.ljust(3, "0")

    return f"D{frac_part}"


def build_star_name(mass, z, rot, dovph=None):
    name = f"{format_mass(mass)}{format_metallicity(z)}{format_rotation(rot)}"
    if dovph is not None:
        name += format_dovph(dovph)
    return name


def read_and_prepare_stars(csv_file: Path):
    csv_file = Path(csv_file).resolve()
    if not csv_file.is_file():
        raise FileNotFoundError(f"CSV introuvable : {csv_file}")

    stars = []
    updated_rows = []

    with csv_file.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter=";")
        required_cols = [COL_NAME, COL_SAVE, COL_MASS, COL_Z, COL_ROT]

        missing = [col for col in required_cols if col not in reader.fieldnames]
        if missing:
            raise KeyError(f"Colonnes manquantes dans le CSV : {missing}")

        fieldnames = reader.fieldnames

        for line_number, row in enumerate(reader, start=2):
            try:
                raw_name = clean_cell(row.get(COL_NAME))
                raw_save = clean_cell(row.get(COL_SAVE))

                if raw_save == "":
                    raise ValueError(f"Ligne {line_number}: chemin de sauvegarde vide.")

                saving_path = Path(raw_save).expanduser().resolve()
                saving_path.mkdir(parents=True, exist_ok=True)

                mass = parse_decimal(row.get(COL_MASS), COL_MASS)
                z = parse_decimal(row.get(COL_Z), COL_Z)
                rot = parse_decimal(row.get(COL_ROT), COL_ROT)

                raw_dovph = clean_cell(row.get(COL_DOVPH, ""))
                dovph = None if raw_dovph == "" else parse_decimal(raw_dovph, COL_DOVPH)

                end_phase_str = clean_cell(row.get(COL_ENDPHASE, ""))
                end_phase = None if end_phase_str == "" else int(end_phase_str)

                name_star = (
                    raw_name if raw_name else build_star_name(mass, z, rot, dovph)
                )

                # Remplit la cellule vide dans la ligne en mémoire
                if raw_name == "":
                    row[COL_NAME] = name_star

                star_dir = saving_path / name_star
                star_dir.mkdir(parents=True, exist_ok=True)
                stars.append(
                    (
                        star_dir,
                        name_star,
                        float(mass),
                        float(z),
                        float(rot),
                        None if end_phase is None else int(end_phase),
                        None if dovph is None else float(dovph),
                    )
                )

                updated_rows.append(row)
                print(f"folder ready : {star_dir}")

            except Exception as e:
                raise ValueError(f"Erreur à la ligne {line_number} du CSV : {e}") from e

    # Réécriture du CSV avec la colonne name_star complétée
    with csv_file.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter=";")
        writer.writeheader()
        writer.writerows(updated_rows)

    return stars


def write_stars_list(stars, output_file: Path):
    output_file = Path(output_file).resolve()
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with output_file.open("w", encoding="utf-8") as f:
        for star_dir, name_star, mass, z, rot, end_phase, dovph in stars:
            f.write(f"{star_dir};{name_star}\n")

    print(f"Liste des étoiles écrite dans : {output_file}")


def submit_genec_launch_array(list_file: Path, slurm_script: Path, max_parallel: int):
    list_file = Path(list_file).resolve()
    slurm_script = Path(slurm_script).resolve()

    with list_file.open("r", encoding="utf-8") as f:
        nstars = sum(1 for _ in f)

    if nstars == 0:
        raise ValueError("Aucune étoile à lancer.")

    cmd = [
        "sbatch",
        "--wait",
        f"--array=1-{nstars}%{max_parallel}",
        str(slurm_script),
        str(list_file),
    ]

    print("Commande de soumission :", " ".join(cmd))
    try:
        res = subprocess.run(cmd, text=True, capture_output=True, check=True)
        print(res.stdout)
    except FileNotFoundError:
        print("sbatch command not found, executing tasks locally sequentially...")
        import os

        for task_id in range(1, nstars + 1):
            env = os.environ.copy()
            env["SLURM_ARRAY_TASK_ID"] = str(task_id)
            env["SLURM_SUBMIT_DIR"] = str(Path(".").resolve())
            env["GENEC_DEFAULT_PROGRAM"] = "/home/aiaiso/GENEC/bin/genec"
            env["GENEC_EMAIL_ADDRESS"] = "test@example.com"
            print(f"--- Running task {task_id}/{nstars} locally ---")
            subprocess.run(["bash", str(slurm_script), str(list_file)], env=env)


def aggregate_and_plot(list_file: Path):
    import matplotlib.pyplot as plt
    import numpy as np

    stars = []
    with list_file.open("r", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                parts = line.strip().split(";")
                stars.append(parts[1])

    # Modes to track
    modes = ["fortran", "cpu", "gpu_eager", "gpu_compile", "gpu_triton", "gpu_fp32"]
    data = {m: [] for m in modes}

    for i, star_name in enumerate(stars, start=1):
        csv_file = Path(f"benchmark_results_{i}.csv")
        row_data = {m: 0.0 for m in modes}
        if csv_file.exists():
            with csv_file.open("r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    m = row["Mode"]
                    if m == "gpu": m = "gpu_eager"  # Compatibility with old files
                    if m in row_data:
                        row_data[m] = float(row["Time_seconds"])
        for m in modes:
            data[m].append(row_data[m])

    x = np.arange(len(stars))
    n_modes = len(modes)
    total_width = 0.8
    width = total_width / n_modes

    fig, ax = plt.subplots(figsize=(16, 8))

    colors = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e", "#e377c2"]
    for i, m in enumerate(modes):
        offset = (i - (n_modes - 1) / 2) * width
        ax.bar(x + offset, data[m], width, label=m, color=colors[i % len(colors)])

    ax.set_ylabel("Execution Time (seconds)", fontsize=12)
    ax.set_title(
        "Execution Time per Star Model and Backend", fontsize=14, fontweight="bold"
    )
    ax.set_xticks(x)
    ax.set_xticklabels(stars, rotation=45, ha="right", fontsize=10)
    ax.legend(fontsize=11, loc='upper left', bbox_to_anchor=(1, 1))
    ax.grid(True, axis="y", alpha=0.3)

    plt.tight_layout()
    output_png = Path("benchmark_comparison.png").resolve()
    plt.savefig(output_png, dpi=150)
    print(f"\nGraphics generated successfully! Saved to: {output_png}")


def main():
    csv_file = Path("./starParam.csv")
    stars = read_and_prepare_stars(csv_file)

    for star_dir, name_star, mass, z, rot, end_phase, dovph in stars:
        print(f"\n=== MAKEINI : {name_star} ===")
        run_genec_makeini(star_dir, name_star, mass, z, rot, end_phase, dovph)

    list_file = Path("./stars_ready.txt")
    write_stars_list(stars, list_file)

    submit_genec_launch_array(
        list_file=list_file,
        slurm_script=Path("./run_genec_array.slurm"),
        max_parallel=10,
    )

    aggregate_and_plot(list_file)


if __name__ == "__main__":
    main()
