import glob
import gzip
import os
import sys


def extract_physics(star_name, calc_dir):
    # Find latest StrucData file
    files = glob.glob(os.path.join(calc_dir, f"{star_name}_StrucData_*.dat*"))
    if not files:
        return "N/A,N/A,N/A,N/A"

    latest = sorted(files)[-1]

    # Use gzip if needed
    opener = gzip.open if latest.endswith(".gz") else open
    mode = "rt"  # read text

    try:
        with opener(latest, mode) as f:
            lines = f.readlines()
            # Data columns: n, log_r, M_int, log_T, log_rho, log_P, Cv, epsilon, X_H1, X_He4, mu, HII, HeII, HeIII
            for line in lines:
                parts = line.split()
                if len(parts) >= 11 and parts[0].isdigit():
                    # This is a data line. Zone 1 is usually the first one (core).
                    logT = parts[3]
                    logRho = parts[4]
                    mu = parts[10]
                    HII = parts[11]
                    return f"{logT},{logRho},{mu},{HII}"
    except Exception as e:
        return f"Err,Err,Err,Err"

    return "N/A,N/A,N/A,N/A"


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("N/A,N/A,N/A,N/A")
    else:
        print(extract_physics(sys.argv[1], sys.argv[2]))
