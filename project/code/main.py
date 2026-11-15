import numpy as np
import matplotlib.pyplot as plt

import graph as g
import simulate as s
import theory as t


ACCENT = "#1F3A5F"      # poster navy
ACCENT_LIGHT = "#5B7FA8"  # lighter navy for the second series
PALE = "#E9EEF4"        # poster block background
GRAY = "#6B7280"

plt.rcParams.update({
    "font.size": 20,                     # keep your current size
    "font.family": "serif",
    "font.serif": ["cmr10"],             # Computer Modern, bundled with matplotlib
    "mathtext.fontset": "cm",            # math in Computer Modern too
    "axes.formatter.use_mathtext": True, # avoids a missing minus-sign glyph with cmr10
    "axes.unicode_minus": False,
    "figure.facecolor": "#E9EEF4",   # pale block colour around the plot
    "savefig.facecolor": "#E9EEF4",  # so the saved PDF keeps it
    "axes.facecolor": "white",       # plotting area stays white
})

SEED = 2026
KS = [3, 4, 7, 8, 15]
N_WALKERS = 10000
V_0 = 0
K_PLOT = 4


#delta settings
EPS = 0.25
R = 30
N_WALKERS_DIST = 10000
N_MAX = 25
N_REPORT = [2, 5, 10, 15, 20] # times reported in the dist table

def t_cross_data(k, N, rng):
    """
    Expected crossing time and N simulated crossing times on the barbell with
    k = k.
    """
    A = g.build_barbell(k)
    P, pi = g.transition_matrix(A)
    g.check_graph(A, P, pi, k)
    nbrs = g.build_neighbor_list(A)
    t_exact = t.crossing_times(P)[V_0]
    T = s.t_cross_sample(nbrs, N, rng, start=V_0)
    return t_exact, T    

def t_cross_analyze(T, t_exact):
    """Mean, SE and difference in standard errors"""
    N = len(T)
    mean = T.mean()
    se = T.std(ddof = 1) / np.sqrt(N)
    n_se = (mean - t_exact) / se  # difference measured in standard errors
    return mean, se, n_se 

def cross_plot(T, t_exact, k):
    """Running mean of T_cross with a +-2 SE region against the theory value"""
    N = len(T)
    m = np.arange(1, N + 1) # number of walkers
    run_mean = np.cumsum(T) / m # computing means for each value of m
    # prevent division by zero for first m = 1 using np.maximum(m - 1, 1)
    run_var = (np.cumsum(T ** 2.0) - m * run_mean ** 2.0) / np.maximum(m - 1, 1)
    run_se = np.sqrt(np.maximum(run_var, 0) / m) # beware of rounding to negatives
    run_se[0] = np.nan # no se for 1 element

    upper = run_mean + 2 * run_se
    lower = run_mean - 2 * run_se

    fig, ax = plt.subplots(figsize=(8.75, 5.625))
    ax.plot(m, run_mean, lw=1.5, label=r"Running mean $\overline{T}_N$")
    ax.fill_between(m, lower, upper, alpha=0.3, label=r"$\pm 2$ SE")
    ax.axhline(t_exact, color="black", ls="--", label=rf"Theory $k^2 = {t_exact:.0f}$")
    ax.set_xscale("log")
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Number of walkers $N$")
    ax.set_ylabel(f"Mean crossing time ($k = {k}$)")
    ax.legend(loc="upper right", fontsize="small")
    fig.tight_layout()
    fig.savefig("../poster/images/cross_plot.pdf")
    plt.close(fig)

def distance_plot(rng):
    """Plot showing D_hat(n) vs D(n) with standard errors"""
    A = g.build_barbell(K_PLOT)
    P, pi = g.transition_matrix(A)
    g.check_graph(A, P, pi, K_PLOT)
    nbrs = g.build_neighbor_list(A)

    Delta = t.exact_distance(P, pi)
    assert len(Delta) > N_MAX
    D_mean, D_se = s.dist_results(pi, nbrs, rng, N_WALKERS_DIST, R, N_MAX, start = V_0)
    lam2 = t.lambda2(P)
    n = np.arange(N_MAX + 1)

    t_mix = t.mixing_time(Delta, eps=EPS)

    fig, ax = plt.subplots(figsize=(8.75, 5.625))
    ax.semilogy(n, Delta[:N_MAX + 1], lw=2, label=r"Exact $\Delta(n)$")
    ax.errorbar(n, D_mean, yerr=2 * D_se, fmt="o", ms=5, capsize=3,
                label=r"Simulated $\hat\Delta(n) \pm 2\,$SE")
    ax.semilogy(n, 0.5 * lam2**n, "--", color="gray",
                label=rf"$\frac{{1}}{{2}}\lambda_2^{{\,n}}$,  $\lambda_2 = {lam2:.3f}$")
    ax.set_xlabel("Step $n$")
    ax.set_ylabel("Distance to equilibrium")
    ax.set_yticks([0.05, 0.1, 0.25, 0.5, 1])
    ax.set_yticklabels(["0.05", "0.1", "0.25", "0.5", "1"])
    ax.annotate(rf"Mixing time $t_{{\mathrm{{mix}}}} = {t_mix}$",
                xy=(t_mix, Delta[t_mix]), xytext=(t_mix + 3, 0.5),
                arrowprops=dict(arrowstyle="->", color="black"), color="black")
    ax.legend(loc="lower left", fontsize="small")
    fig.tight_layout()
    fig.savefig("../poster/images/dist_plot.pdf")
    plt.close(fig)
   
    report = sorted(set(N_REPORT) | {t_mix})
    report = [el for el in report if el <= N_MAX]

    f = open("../poster/tables/table_dist.tex", "w")
    for m in report:
        tag = ""
        n_se = (D_mean[m] - Delta[m]) / D_se[m]
        if m == t_mix:
            tag = r" $= t_{\mathrm{mix}}(\tfrac14)$"
            print(f"n={m:2d} (t_mix({EPS})): exact {Delta[m]:.4f}  simulated {D_mean[m]:.4f} ± {D_se[m]:.4f}  n_se={n_se:+.2f}")
        else:
            print(f"n={m:2d}: exact {Delta[m]:.4f}  simulated {D_mean[m]:.4f} ± {D_se[m]:.4f}  n_se={n_se:+.2f}")
        f.write(f"{m}{tag} & {Delta[m]:.4f} & ${D_mean[m]:.4f} \\pm {D_se[m]:.4f}$ & ${n_se:+.1f}$ \\\\\n")
    f.close()

def main():
    rng = np.random.default_rng(SEED)
    data = {k : t_cross_data(k, N_WALKERS, rng) for k in KS}

    f = open("../poster/tables/table_cross.tex", "w") # prepare a LaTex table
    for k in KS:
        t_exact, T = data[k]
        mean, se, n_se = t_cross_analyze(T, t_exact)
        print(f"k={k}: {mean:.2f} ± {se:.2f}   theory {t_exact:.0f}   n_se={n_se:+.2f}")
        # writing a LaTex row
        f.write(f"{k} & {N_WALKERS} & ${mean:.2f} \\pm {se:.2f}$ & {t_exact:.0f} & ${n_se:+.1f}$ \\\\\n")
    f.close()

    t_exact, T = data[K_PLOT]
    cross_plot(T, t_exact, K_PLOT)

    distance_plot(rng)
    
if __name__ == "__main__":
    main()
