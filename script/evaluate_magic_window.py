"""
Evaluate funnel vs cylinder simulation results for main.tex.
Tests whether the data support a persistent geometric bias (extracting order from heat noise).
"""
import os
import numpy as np

RESULTS_DIR = "results"
# Steady-state: use last 50_000 steps (paper uses "final 10,000-step window"; we use more for stability)
STEADY_STEPS = 50_000


def load_dat(path):
    """Load ave/time file; return step, n_na_l, n_cl_l, n_na_r, n_cl_r."""
    data = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            data.append([float(x) for x in parts])
    if not data:
        return None
    arr = np.array(data)
    step = arr[:, 0]
    n_na_l, n_cl_l, n_na_r, n_cl_r = arr[:, 1], arr[:, 2], arr[:, 3], arr[:, 4]
    return step, n_na_l, n_cl_l, n_na_r, n_cl_r


def steady_state_means(step, n_na_l, n_cl_l, n_na_r, n_cl_r):
    """Return mean net charge over the last STEADY_STEPS."""
    net_left = n_na_l - n_cl_l
    net_right = n_na_r - n_cl_r
    # polarization: positive = more cations on right
    polarization = net_right - net_left
    last_t = step[-1]
    mask = step >= (last_t - STEADY_STEPS)
    if mask.sum() == 0:
        mask = np.ones(len(step), dtype=bool)
    return {
        "net_left": float(np.mean(net_left[mask])),
        "net_right": float(np.mean(net_right[mask])),
        "polarization": float(np.mean(polarization[mask])),
        "n_pts": int(mask.sum()),
    }


def collect_runs(prefix, max_seeds=20):
    """Load all available seeds for a given prefix (funnel_window or cylinder_control)."""
    out = []
    for i in range(max_seeds):
        path = os.path.join(RESULTS_DIR, f"data.{prefix}_s{i}.dat")
        if not os.path.isfile(path):
            break
        loaded = load_dat(path)
        if loaded is None:
            break
        step, n_na_l, n_cl_l, n_na_r, n_cl_r = loaded
        ss = steady_state_means(step, n_na_l, n_cl_l, n_na_r, n_cl_r)
        ss["seed_idx"] = i
        ss["steps_max"] = float(step[-1])
        out.append(ss)
    return out


def main():
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    funnel = collect_runs("funnel_window")
    cylinder = collect_runs("cylinder_control")
    n_f = len(funnel)
    n_c = len(cylinder)
    print("--- Data summary ---")
    print(f"Funnel runs: {n_f}, Cylinder runs: {n_c}")
    if n_f == 0 or n_c == 0:
        print("Insufficient data. Run run_magic_window.py first.")
        return

    # Use common seed set (min of the two)
    n = min(n_f, n_c)
    F = funnel[:n]
    C = cylinder[:n]

    # Observable: same as paper (net charge separation). Use right-side net charge
    # as "charge at tip" proxy, and polarization = net_right - net_left.
    mu_f_right = np.mean([r["net_right"] for r in F])
    sigma_f_right = np.std([r["net_right"] for r in F], ddof=1) if n > 1 else 0.0
    mu_c_right = np.mean([r["net_right"] for r in C])
    sigma_c_right = np.std([r["net_right"] for r in C], ddof=1) if n > 1 else 0.0

    mu_f_pol = np.mean([r["polarization"] for r in F])
    sigma_f_pol = np.std([r["polarization"] for r in F], ddof=1) if n > 1 else 0.0
    mu_c_pol = np.mean([r["polarization"] for r in C])
    sigma_c_pol = np.std([r["polarization"] for r in C], ddof=1) if n > 1 else 0.0

    delta_right = mu_f_right - mu_c_right
    delta_pol = mu_f_pol - mu_c_pol
    # Combined uncertainty (standard error of the difference)
    se_right = np.sqrt(sigma_f_right**2 / n + sigma_c_right**2 / n) if n > 1 else 0.0
    se_pol = np.sqrt(sigma_f_pol**2 / n + sigma_c_pol**2 / n) if n > 1 else 0.0

    print("\n--- Steady-state net charge (last {} steps) ---".format(STEADY_STEPS))
    print("Right-side net charge (Na+ - Cl- on right reservoir):")
    print("  Funnel  mean = {:.3f}  std = {:.3f}  (n={})".format(mu_f_right, sigma_f_right, n))
    print("  Cylinder mean = {:.3f}  std = {:.3f}  (n={})".format(mu_c_right, sigma_c_right, n))
    print("  Delta (Funnel - Cylinder) = {:.3f}  (SE of diff = {:.3f})".format(delta_right, se_right))
    print("Polarization (net_right - net_left):")
    print("  Funnel  mean = {:.3f}  std = {:.3f}".format(mu_f_pol, sigma_f_pol))
    print("  Cylinder mean = {:.3f}  std = {:.3f}".format(mu_c_pol, sigma_c_pol))
    print("  Delta (Funnel - Cylinder) = {:.3f}  (SE of diff = {:.3f})".format(delta_pol, se_pol))

    # Significance: is |delta| > 2*SE (rough 95% CI)?
    sig_right = abs(delta_right) > 2 * se_right if se_right > 0 else False
    sig_pol = abs(delta_pol) > 2 * se_pol if se_pol > 0 else False

    print("\n--- Evaluation for main.tex (energy from heat noise) ---")
    if n < 4:
        print("Warning: only {} seeds; paper uses N=20 for claims.".format(n))
    if sig_right or sig_pol:
        print("Result: The funnel shows a statistically significant geometric bias vs cylinder")
        print("        (difference larger than ~2*SE). This is consistent with the paper's")
        print("        claim that spatial asymmetry can bias the charge distribution in an")
        print("        isothermal system (rectification of thermal fluctuations).")
    else:
        print("Result: The difference between Funnel and Cylinder is NOT statistically")
        print("        significant (|Delta| < ~2*SE). With the current seed set and run length,")
        print("        the simulation does not clearly support a persistent geometric bias")
        print("        for extracting order from heat noise; the effect may be smaller than")
        print("        the noise or require more/longer runs.")
    print("\n(Conclusion is based on steady-state mean comparison over {} seeds.)".format(n))


if __name__ == "__main__":
    main()
