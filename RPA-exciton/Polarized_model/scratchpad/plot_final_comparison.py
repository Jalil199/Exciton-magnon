import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "serif", "font.size": 13,
                      "axes.labelsize": 14, "axes.titlesize": 13, "legend.fontsize": 10})

theta_eq0 = float(np.load("_theta_eq0.npy"))
data = {tx: np.load(f"_dyn_damped_{tx}.npz") for tx in (260, 500, 800)}
colors = {260: "#2166ac", 500: "#762a83", 800: "#b2182b"}

fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))

# left: full 800ps window, envelope visible
for tx in (260, 500, 800):
    t = data[tx]["t"]; th = data[tx]["theta"]
    ax[0].plot(t, th, "-", lw=0.7, color=colors[tx], label=fr"$\tau_X$={tx}ps")
ax[0].axhline(theta_eq0, color="k", ls="--", lw=1, label=r"bare eq. $\theta_0$")
ax[0].set_xlabel(r"$\Delta t$ (ps)")
ax[0].set_ylabel(r"$\theta(t)$ (deg)")
ax[0].set_title("Full run (800 ps)")
ax[0].legend(frameon=False, loc="upper right")
ax[0].grid(alpha=0.25)

# right: zoom on first 300ps so individual oscillation cycles are visible
for tx in (260, 500, 800):
    t = data[tx]["t"]; th = data[tx]["theta"]
    mask = t <= 300
    ax[1].plot(t[mask], th[mask], "-", lw=1.4, color=colors[tx], label=fr"$\tau_X$={tx}ps")
ax[1].axhline(theta_eq0, color="k", ls="--", lw=1)
ax[1].set_xlabel(r"$\Delta t$ (ps)")
ax[1].set_ylabel(r"$\theta(t)$ (deg)")
ax[1].set_title("Zoom: first 300 ps (individual cycles)")
ax[1].legend(frameon=False, loc="upper right")
ax[1].grid(alpha=0.25)

fig.suptitle(r"With added local Gilbert-like damping $\eta$ — clear decay, still no sawtooth", y=1.01)
plt.tight_layout()
outpath = "final_comparison_pretty.png"
plt.savefig(outpath, dpi=150, bbox_inches="tight")
print("saved", outpath)
