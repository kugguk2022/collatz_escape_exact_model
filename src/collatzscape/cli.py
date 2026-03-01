from __future__ import annotations
import argparse, os
import numpy as np
import yaml
import matplotlib.pyplot as plt

from .plot import basin_plot
from .metrics import escape_rate, finite_time_lyapunov_rough, maximal_lyapunov_stats

def _load_config(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def _complex_from_list(v):
    return complex(float(v[0]), float(v[1]))

def cmd_demo(args):
    cfg = _load_config(args.config)
    g = cfg["grid"]
    m = cfg["model"]
    c = _complex_from_list(m["c"])
    out = args.out or "figures"
    os.makedirs(out, exist_ok=True)
    basin_plot(**g, c=c, alpha=float(m.get("alpha", 0.0)), gamma=float(m["gamma"]), direction=int(m["direction"]), steps=int(m["steps"]),
               kick_mode=str(m.get("kick_mode","tanh")), C_mode=str(m.get("C_mode","cosine")), paired=bool(m.get("paired", True)),
               outpath=os.path.join(out, "demo_escape.png"))

def cmd_sweep(args):
    cfg = _load_config(args.config)
    g = cfg["grid"]
    m = cfg["model"]
    s = cfg["sweep"]
    c = _complex_from_list(m["c"])
    steps = int(m["steps"])
    out = args.out or "figures"
    os.makedirs(out, exist_ok=True)

    # sample points for metrics (cheap)
    rng = np.random.default_rng(int(cfg.get("seed", 0)))
    pts = (rng.uniform(g["re_min"], g["re_max"], size=2000)
           + 1j*rng.uniform(g["im_min"], g["im_max"], size=2000))

    rows = []
    base_seed = int(cfg.get("seed", 0))
    ly_warmup = int(m.get("lyapunov_warmup_steps", 32))
    ly_rough_warmup = int(m.get("lyapunov_rough_warmup_steps", 0))
    ly_trials = int(m.get("lyapunov_trials", 8))
    ly_eps = float(m.get("lyapunov_eps", 1e-8))
    start_z = complex(float(m.get("lyapunov_z0_re", 1.0)), float(m.get("lyapunov_z0_im", 0.0)))

    for dir_idx, direction in enumerate(s["directions"]):
        for gamma_idx, gamma in enumerate(s["gammas"]):
            er = escape_rate(
                pts, c=c, alpha=float(m.get("alpha", 0.0)), gamma=float(gamma), direction=int(direction), steps=steps,
                kick_mode=str(m.get("kick_mode", "tanh")), C_mode=str(m.get("C_mode", "cosine")),
                paired=bool(m.get("paired", True))
            )

            ly_rough = finite_time_lyapunov_rough(
                start_z, c=c, alpha=float(m.get("alpha", 0.0)), gamma=float(gamma), direction=int(direction), steps=steps,
                kick_mode=str(m.get("kick_mode", "tanh")), C_mode=str(m.get("C_mode", "cosine")),
                paired=bool(m.get("paired", True)), eps=ly_eps, warmup_steps=ly_rough_warmup,
                seed=base_seed + 20000 + 1000 * dir_idx + gamma_idx,
            )

            mle_stats = maximal_lyapunov_stats(
                start_z,
                c=c,
                alpha=float(m.get("alpha", 0.0)),
                gamma=float(gamma),
                direction=int(direction),
                steps=steps,
                kick_mode=str(m.get("kick_mode", "tanh")),
                C_mode=str(m.get("C_mode", "cosine")),
                paired=bool(m.get("paired", True)),
                eps=ly_eps,
                warmup_steps=ly_warmup,
                n_trials=ly_trials,
                seed=base_seed + 1000 * dir_idx + gamma_idx,
            )
            rows.append(
                (
                    float(gamma),
                    int(direction),
                    float(er),
                    float(ly_rough),
                    float(mle_stats["mle"]),
                    float(mle_stats["mle_std"]),
                    int(mle_stats["valid_trials"]),
                )
            )

    import csv
    csv_path = os.path.join(out, "sweep_metrics.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "gamma",
                "direction",
                "escape_rate",
                "ft_lyapunov_rough",
                "mle",
                "mle_std",
                "mle_valid_trials",
            ]
        )
        w.writerows(rows)

    # plots
    for direction in sorted(set(r[1] for r in rows)):
        dir_rows = sorted((r for r in rows if r[1] == direction), key=lambda x: x[0])
        xs = [r[0] for r in dir_rows]
        ers = [r[2] for r in dir_rows]
        ly_rough = [r[3] for r in dir_rows]
        ly_mle = [r[4] for r in dir_rows]
        ly_mle_std = [r[5] for r in dir_rows]

        plt.figure()
        plt.plot(xs, ers, marker="o")
        plt.xlabel("gamma")
        plt.ylabel("escape_rate")
        plt.title(f"Escape rate vs gamma (direction={direction})")
        plt.savefig(os.path.join(out, f"sweep_escape_rate_dir{direction}.png"), dpi=200, bbox_inches="tight")
        plt.close()

        plt.figure()
        plt.errorbar(xs, ly_mle, yerr=ly_mle_std, marker="o", capsize=3)
        plt.xlabel("gamma")
        plt.ylabel("MLE (Benettin, mean +/- std)")
        plt.title(f"Maximal Lyapunov vs gamma (direction={direction})")
        plt.savefig(os.path.join(out, f"sweep_lyapunov_dir{direction}.png"), dpi=200, bbox_inches="tight")
        plt.close()

        plt.figure()
        plt.plot(xs, ly_rough, marker="o", label="rough (single direction)")
        plt.errorbar(xs, ly_mle, yerr=ly_mle_std, marker="o", capsize=3, label="accurate MLE")
        plt.xlabel("gamma")
        plt.ylabel("Lyapunov exponent (per step)")
        plt.title(f"Lyapunov estimator comparison (direction={direction})")
        plt.legend()
        plt.savefig(os.path.join(out, f"sweep_lyapunov_compare_dir{direction}.png"), dpi=200, bbox_inches="tight")
        plt.close()



def cmd_fatou(args):
    cfg = _load_config(args.config)
    g = cfg["grid"]
    m = cfg["model"]
    c = _complex_from_list(m["c"])
    out = args.out or "figures"
    os.makedirs(out, exist_ok=True)

    from .fatou import find_attractors, fatou_label_grid, local_stability, basin_measure
    from .plot import fatou_components_plot

    attractors = find_attractors(
        c=c, gamma=float(m["gamma"]), direction=int(m["direction"]), steps=int(m["steps"]),
        re_min=g["re_min"], re_max=g["re_max"], im_min=g["im_min"], im_max=g["im_max"],
        n_samples=int(args.samples), tol=float(args.tol), seed=int(cfg.get("seed", 0)),
    )
    attractors = attractors[: int(args.max_attractors)]

    # stability summary
    rows = []
    for idx, a in enumerate(attractors):
        st = local_stability(a, c=c, gamma=float(m["gamma"]), direction=int(m["direction"]))
        rows.append([idx, a.real, a.imag, st["spectral_radius"], st["eig1"], st["eig2"]])

    import csv
    with open(os.path.join(out, "fatou_attractors.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["id", "re", "im", "spectral_radius", "eig1", "eig2"])
        w.writerows(rows)

    labels = fatou_label_grid(**g, attractors=attractors, c=c, gamma=float(m["gamma"]),
                              direction=int(m["direction"]), steps=int(m["steps"]))
    fatou_components_plot(labels, re_min=g["re_min"], re_max=g["re_max"], im_min=g["im_min"], im_max=g["im_max"],
                          outpath=os.path.join(out, "fatou_components.png"),
                          title=f"Fatou components proxy (gamma={m['gamma']}, direction={m['direction']})")

    meas = basin_measure(labels)
    with open(os.path.join(out, "fatou_basin_measure.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["name", "fraction"])
        for k, v in meas.items():
            w.writerow([k, v])

def main():
    ap = argparse.ArgumentParser(prog="collatzscape")
    sub = ap.add_subparsers(dest="cmd", required=True)

    ap_demo = sub.add_parser("demo", help="Generate a demo escape-set figure")
    ap_demo.add_argument("--config", required=True)
    ap_demo.add_argument("--out", default="figures")
    ap_demo.set_defaults(func=cmd_demo)

    ap_sw = sub.add_parser("sweep", help="Sweep gamma and directions, write metrics+plots")
    ap_sw.add_argument("--config", required=True)
    ap_sw.add_argument("--out", default="figures")
    ap_sw.set_defaults(func=cmd_sweep)

    ap_ft = sub.add_parser("fatou", help="Estimate Fatou components/basins (proxy) and attractor stability")
    ap_ft.add_argument("--config", required=True)
    ap_ft.add_argument("--out", default="figures")
    ap_ft.add_argument("--samples", type=int, default=5000)
    ap_ft.add_argument("--tol", type=float, default=1e-3)
    ap_ft.add_argument("--max-attractors", type=int, default=8)
    ap_ft.set_defaults(func=cmd_fatou)

    args = ap.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
