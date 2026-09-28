"""
Scan tau_X and pump-density scale factor (proxy for fluence, since N_bright ~ n_c^2 with an
~universal theta-shape, Step 9b) in the NON-adiabatic 2-macrospin LLG, to find where the
amplitude ratio (300-400ps / 0-100ps) crosses 1 -- i.e. where net damping flips to net growth
(self-sustained oscillation, the thing we're chasing to match Brennan Fig 3a).
"""
import numpy as np
from scipy.integrate import solve_ivp
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
hbar_eVs = 6.582119569e-16
muB_eVperT = 5.7883818060e-5
ex = np.array([1.0,0.0,0.0]); ey = np.array([0.0,1.0,0.0]); ez = np.array([0.0,0.0,1.0])

J_AF = 6e-6; Kx = 14e-6; Kz = 58e-6
gmuB_B = 2.0 * muB_eVperT * 0.2

def crossmat(v): return np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
def Bcl(M, Mbar): return -J_AF*Mbar + gmuB_B*ey - 2*Kz*np.dot(M,ez)*ez + 2*Kx*np.dot(M,ex)*ex

cos_a_eq = gmuB_B/(2*(J_AF+Kx)); a_eq = np.arccos(np.clip(cos_a_eq,-1,1)); theta_eq0 = 2*a_eq

rows = list(csv.DictReader(open(os.path.join(here, "step10_damping_data.csv"))))
theta_arr = np.array([float(r["theta_over_pi"]) for r in rows]) * np.pi
Neq_arr0 = np.array([float(r["Neq"]) for r in rows])
dEXdtheta_arr = np.array([float(r["dEXdtheta_eV"]) for r in rows])

def dEXdtheta(theta): return np.interp(theta, theta_arr, dEXdtheta_arr)

kick_deg = 30.0
theta_kick = theta_eq0 + np.radians(kick_deg); a_kick = theta_kick/2
M1_0 = np.array([np.sin(a_kick), np.cos(a_kick), 0.0])
M2_0 = np.array([-np.sin(a_kick), np.cos(a_kick), 0.0])

def run(tauX_ps, density_scale, t_max_ps=400.0, n_out=6000):
    Neq_arr = Neq_arr0 * density_scale
    def Neq(theta): return np.interp(theta, theta_arr, Neq_arr)
    N0 = Neq(theta_kick)
    y0 = np.concatenate([M1_0, M2_0, [N0]])
    tauX_nat = (tauX_ps*1e-12)/hbar_eVs
    def rhs(t, y):
        M1=y[0:3]; M2=y[3:6]; N=y[6]
        M1=M1/np.linalg.norm(M1); M2=M2/np.linalg.norm(M2)
        costheta=np.clip(np.dot(M1,M2),-1,1); theta=np.arccos(costheta)
        jth = N*dEXdtheta(theta)/np.sin(theta)
        B1 = Bcl(M1,M2) + jth*M2; B2 = Bcl(M2,M1) + jth*M1
        Mdot = np.concatenate([np.cross(M1,B1), np.cross(M2,B2)])
        Ndot = -(N - Neq(theta))/tauX_nat
        return np.concatenate([Mdot, [Ndot]])
    t_max_nat = (t_max_ps*1e-12)/hbar_eVs
    sol = solve_ivp(rhs, [0,t_max_nat], y0, t_eval=np.linspace(0,t_max_nat,n_out),
                     method="RK45", rtol=1e-9, atol=1e-12)
    M1=sol.y[0:3,:]; M2=sol.y[3:6,:]
    costh=np.clip(np.sum(M1*M2,axis=0),-1,1); theta_t=np.degrees(np.arccos(costh))
    t_ps = sol.t*hbar_eVs*1e12
    n=len(theta_t); q=n//4
    amp_first=theta_t[:q].max()-theta_t[:q].min(); amp_last=theta_t[-q:].max()-theta_t[-q:].min()
    return theta_t, t_ps, amp_last/amp_first, sol.success

print("=== scan 1: tau_X (fixed density_scale=1, i.e. s=0.7 as before) ===")
tauX_list = [260, 500, 800, 1500, 3000, 6000, 12000, 25000, 50000]
ratios_tau = []
for tx in tauX_list:
    _,_,ratio,ok = run(tx, 1.0)
    ratios_tau.append(ratio)
    print(f"  tau_X={tx:>7.0f}ps  ratio={ratio:.4f}  {'OK' if ok else 'FAILED'}")

print("\n=== scan 2: density scale factor (fixed tau_X=800ps, Brennan's upper end) ===")
dens_list = [1, 2, 4, 8, 16, 32, 64, 128]
ratios_dens = []
for ds in dens_list:
    _,_,ratio,ok = run(800, ds)
    ratios_dens.append(ratio)
    print(f"  density_scale={ds:>4d}x  ratio={ratio:.4f}  {'OK' if ok else 'FAILED'}")

fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
ax[0].semilogx(tauX_list, ratios_tau, "o-", color="#2166ac")
ax[0].axhline(1.0, color="r", ls="--", lw=1, label="growth threshold")
ax[0].set_xlabel(r"$\tau_X$ (ps)"); ax[0].set_ylabel("amplitude ratio (300-400ps / 0-100ps)")
ax[0].set_title("Threshold scan: exciton lifetime")
ax[0].legend(frameon=False); ax[0].grid(alpha=0.25, which="both")

ax[1].semilogx(dens_list, ratios_dens, "s-", color="#b2182b")
ax[1].axhline(1.0, color="r", ls="--", lw=1, label="growth threshold")
ax[1].set_xlabel(r"pump density scale ($\times$ s=0.7 level)"); ax[1].set_ylabel("amplitude ratio")
ax[1].set_title(r"Threshold scan: pump density ($\tau_X$=800ps)")
ax[1].legend(frameon=False); ax[1].grid(alpha=0.25, which="both")

plt.tight_layout()
outpath = os.path.join(here, "threshold_scan.png")
plt.savefig(outpath, dpi=150)
print("\nsaved", outpath)
