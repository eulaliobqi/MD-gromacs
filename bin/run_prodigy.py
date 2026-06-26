#!/usr/bin/env python3
"""
run_prodigy.py — Calcula ΔG de ligação via PRODIGY para todos os sistemas MD.

Estratégia:
  GORE4 + SKTI → PRODIGY (proteína-proteína/peptídeo)
  BEN           → Vina scores (já disponíveis, impressos no resumo)

Para cada sistema, extrai 5 snapshots da trajetória (30, 40, 50, 60, 70 ns),
adiciona chain IDs (receptor=A, ligante=B) e calcula ΔG médio via PRODIGY.

Requisito: pip install prodigy MDAnalysis
Uso:
  mamba run -n md-gromacs python bin/run_prodigy.py
  mamba run -n prodigy-env python bin/run_prodigy.py
"""

import os, csv, tempfile, subprocess, statistics
import MDAnalysis as mda

# ── Sistemas MD disponíveis ────────────────────────────────────────────────────
SISTEMAS = [
    # (label, prod_dir, ndx_file)
    ("ACR157-GORE4",
     "ACR157-GORE4_NEW/MD/acr157-gore4-c1/prod",
     "ACR157-GORE4_NEW/MD/acr157-gore4-c1/analise/lig.ndx"),
    ("QCL936-GORE4",
     "QCL936-GORE4_NEW/MD/qcl936-gore4-c3/prod",
     "QCL936-GORE4_NEW/MD/qcl936-gore4-c3/analise/lig.ndx"),
    ("XP273-GORE4",
     "XP273-GORE4_NEW/MD/xp273-gore4-c1n/prod",
     "XP273-GORE4_NEW/MD/xp273-gore4-c1n/analise/lig.ndx"),
    ("ACR157-SKTI",
     "ACR157-SKTI_NEW/MD/acr157-skti-c2/prod",
     "ACR157-SKTI_NEW/MD/acr157-skti-c2/analise/lig.ndx"),
    ("QCL936-SKTI",
     "QCL936-SKTI-NEW/MD/qcl936-skti-c2/prod",
     "QCL936-SKTI-NEW/MD/qcl936-skti-c2/analise/lig.ndx"),
    ("XP273-SKTI",
     "XP273-SKTI_NEW/MD/xp273-skti-c2/prod",
     "XP273-SKTI_NEW/MD/xp273-skti-c2/analise/lig.ndx"),
]

# Frames a extrair: 30, 40, 50, 60, 70 ns (10 frames/ns → índices 300,400,500,600,700)
FRAMES_NS = [30, 40, 50, 60, 70]
FRAMES_PER_NS = 10


def parse_ndx_groups(ndx_path):
    """Retorna dict {grupo: [atoms indices]} para Receptor e Ligante."""
    groups = {}
    current = None
    with open(ndx_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("["):
                current = line.strip("[] ").strip()
                groups[current] = []
            elif current and line:
                groups[current].extend(int(x) - 1 for x in line.split())
    return groups


def build_complex_pdb(gro, xtc, ndx, frame_idx, outpdb):
    """
    Extrai um frame da trajetória, separa receptor (chain A) e ligante (chain B),
    escreve PDB combinado com chain IDs para PRODIGY.
    """
    import warnings
    warnings.filterwarnings("ignore")

    u = mda.Universe(gro, xtc)
    groups = parse_ndx_groups(ndx)

    rec_ids = groups.get("Receptor", [])
    lig_ids = groups.get("Ligante", [])
    if not rec_ids or not lig_ids:
        raise ValueError(f"Grupos Receptor/Ligante não encontrados em {ndx}")

    u.trajectory[frame_idx]

    rec = u.atoms[rec_ids]
    lig = u.atoms[lig_ids]

    # Assign chain IDs
    rec.chainIDs = "A"
    lig.chainIDs = "B"

    combined = rec + lig
    combined.write(outpdb)
    return True


def run_prodigy(pdb_path):
    """
    Chama PRODIGY e extrai ΔG (kcal/mol) e Kd.
    Retorna (dg_float, kd_str) ou (None, None) em caso de erro.
    """
    try:
        result = subprocess.run(
            ["prodigy", "--quiet", pdb_path],
            capture_output=True, text=True, timeout=30
        )
        output = result.stdout + result.stderr
        dg = None
        kd = "N/A"
        for line in output.splitlines():
            ll = line.lower()
            if "predicted binding affinity" in ll or "binding affinity" in ll:
                try:
                    dg = float(line.split(":")[-1].strip().split()[0])
                except (ValueError, IndexError):
                    pass
            if "dissociation" in ll or "kd" in ll:
                parts = line.split(":")
                if len(parts) > 1:
                    kd = parts[-1].strip().split()[0]
        return dg, kd
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        print(f"  [ERRO PRODIGY] {e}")
        return None, None


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    results = []

    print("=" * 60)
    print("PRODIGY — Predição de ΔG de ligação")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as tmpdir:
        for label, prod_dir, ndx in SISTEMAS:
            gro = os.path.join(prod_dir, "md.gro")
            xtc = os.path.join(prod_dir, "md_fit.xtc")

            missing = [f for f in [gro, xtc, ndx] if not os.path.exists(f)]
            if missing:
                print(f"\n[AUSENTE] {label}: {missing}")
                results.append({"sistema": label, "dG_mean": "N/A", "dG_std": "N/A",
                                 "Kd": "N/A", "n_frames": 0, "nota": "arquivos ausentes"})
                continue

            print(f"\n--- {label} ---")
            dg_values = []
            kd_last = "N/A"

            for ns in FRAMES_NS:
                frame_idx = ns * FRAMES_PER_NS
                pdb_out = os.path.join(tmpdir, f"{label}_{ns}ns.pdb")

                try:
                    build_complex_pdb(gro, xtc, ndx, frame_idx, pdb_out)
                except Exception as e:
                    print(f"  {ns} ns: ERRO extração — {e}")
                    continue

                dg, kd = run_prodigy(pdb_out)
                if dg is not None:
                    dg_values.append(dg)
                    kd_last = kd
                    print(f"  {ns} ns: ΔG = {dg:+.2f} kcal/mol  Kd = {kd}")
                else:
                    print(f"  {ns} ns: PRODIGY sem output")

            if dg_values:
                dg_mean = statistics.mean(dg_values)
                dg_std  = statistics.stdev(dg_values) if len(dg_values) > 1 else 0.0
                print(f"  → ΔG médio = {dg_mean:+.2f} ± {dg_std:.2f} kcal/mol  (n={len(dg_values)})")
                results.append({
                    "sistema": label,
                    "dG_mean": f"{dg_mean:.2f}",
                    "dG_std":  f"{dg_std:.2f}",
                    "Kd": kd_last,
                    "n_frames": len(dg_values),
                    "nota": ""
                })
            else:
                results.append({"sistema": label, "dG_mean": "N/A", "dG_std": "N/A",
                                 "Kd": "N/A", "n_frames": 0, "nota": "PRODIGY falhou"})

    # ── BEN (Vina scores) ──────────────────────────────────────────────────────
    ben_scores = {
        "ACR157-BEN": -4.953,
        "QCL936-BEN": -5.733,
        "XP273-BEN":  -5.484,
        "XP352-BEN":  -4.975,
    }
    print("\n--- Série BEN (AutoDock Vina) ---")
    for label, score in ben_scores.items():
        print(f"  {label}: score = {score:+.3f} kcal/mol")
        results.append({"sistema": label, "dG_mean": f"{score:.3f}", "dG_std": "—",
                         "Kd": "N/A (Vina)", "n_frames": 1, "nota": "Vina score modo 1"})

    # ── Salvar CSV ─────────────────────────────────────────────────────────────
    outcsv = "prodigy_results.csv"
    with open(outcsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["sistema","dG_mean","dG_std","Kd","n_frames","nota"])
        w.writeheader()
        w.writerows(results)

    print(f"\n[OK] Resultados salvos em: {outcsv}")
    print("\n=== Resumo ===")
    print(f"{'Sistema':<20} {'ΔG (kcal/mol)':<18} {'Kd':<15} {'Método'}")
    print("-" * 70)
    for r in results:
        dg_str = f"{r['dG_mean']} ± {r['dG_std']}" if r['dG_std'] not in ("N/A", "—") else r['dG_mean']
        metodo = "Vina" if "BEN" in r["sistema"] else "PRODIGY"
        print(f"{r['sistema']:<20} {dg_str:<18} {r['Kd']:<15} {metodo}")


if __name__ == "__main__":
    main()
