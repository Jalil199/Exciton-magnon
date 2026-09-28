import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import csv, os

here = os.path.dirname(os.path.abspath(__file__))
rows = list(csv.DictReader(open(os.path.join(here, "step10_damping_data.csv"))))
theta_pi = np.array([float(r["theta_over_pi"]) for r in rows])
tau1_meV = np.array([1e3*float(r["tau1_eV"]) for r in rows])
eta260 = np.array([float(r["eta_260ps"]) for r in rows])
eta500 = np.array([float(r["eta_500ps"]) for r in rows])
eta800 = np.array([float(r["eta_800ps"]) for r in rows])

fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))

ax[0].plot(theta_pi, tau1_meV, "o-", color="purple", lw=2, ms=4)
ax[0].axhline(0, color="gray", lw=0.6)
ax[0].set_xlabel(r"$\theta/\pi$  (FM$\to$AFM)")
ax[0].set_ylabel(r"$\tau_1 = N_{eq}\cdot\partial_\theta E_X$  (meV)")
ax[0].set_title("Reactive torque (Step 9)")
ax[0].grid(alpha=0.25)

ax[1].plot(theta_pi, eta260, "o-", color="#2166ac", ms=4, label=r"$\tau_X=260$ ps")
ax[1].plot(theta_pi, eta500, "s-", color="#762a83", ms=4, label=r"$\tau_X=500$ ps")
ax[1].plot(theta_pi, eta800, "^-", color="#b2182b", ms=4, label=r"$\tau_X=800$ ps")
ax[1].axhline(0, color="gray", lw=0.6)
ax[1].set_xlabel(r"$\theta/\pi$  (FM$\to$AFM)")
ax[1].set_ylabel(r"$\eta_\theta(\theta)$  (dimensionless)")
ax[1].set_title("Population-lag damping (Step 10, standalone)")
ax[1].legend(frameon=False)
ax[1].grid(alpha=0.25)

plt.tight_layout()
outpath = os.path.join(here, "step10_torque_damping.png")
plt.savefig(outpath, dpi=150)
print("saved", outpath)
