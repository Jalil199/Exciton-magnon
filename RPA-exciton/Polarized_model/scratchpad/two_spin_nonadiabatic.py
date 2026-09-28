"""
Non-adiabatic version: N(t) is now a DYNAMICAL variable (7th component) obeying the same
relaxation equation dN/dt = -(N - N_eq(theta(t)))/tau_X, WITHOUT the first-order-in-thetadot
adiabatic elimination used to build eta_theta(theta) before. Torque = N(t)*dE_X/dtheta(theta(t))
directly, at every instant -- no separate phenomenological eta term needed; any damping/anti-
damping character now emerges purely from the population lag itself.

Same real literature classical-Hamiltonian parameters as two_spin_llg.py (J_int=6ueV,
A_x=14ueV, A_z=58ueV, B=0.2T). Same large kick (+30deg) as two_spin_bigkick.py, for direct
comparison with the adiabatic (eta_theta) result.
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
Neq_arr = np.array([float(r["Neq"]) for r in rows])
dEXdtheta_arr = np.array([float(r["dEXdtheta_eV"]) for r in rows])

def Neq(theta): return np.interp(theta, theta_arr, Neq_arr)
def dEXdtheta(theta): return np.interp(theta, theta_arr, dEXdtheta_arr)

hbar_eVs_ = hbar_eVs
def rhs_nonadiabatic(t, y, tauX_ps):
    M1=y[0:3]; M2=y[3:6]; N=y[6]
    M1=M1/np.linalg.norm(M1); M2=M2/np.linalg.norm(M2)
    costheta=np.clip(np.dot(M1,M2),-1,1); theta=np.arccos(costheta)
    jth = N*dEXdtheta(theta)/np.sin(theta)     # j_theta(t) = N(t)*dE_X/dtheta / sin(theta), NO eta term
    B1 = Bcl(M1,M2) + jth*M2
    B2 = Bcl(M2,M1) + jth*M1
    Mdot = np.concatenate([np.cross(M1,B1), np.cross(M2,B2)])   # no implicit damping now -> explicit
    tauX_nat = (tauX_ps*1e-12)/hbar_eVs_                        # natural units (1/eV)
    Ndot = -(N - Neq(theta))/tauX_nat
    return np.concatenate([Mdot, [Ndot]])

kick_deg = 30.0
theta_kick = theta_eq0 + np.radians(kick_deg); a_kick = theta_kick/2
M1_0 = np.array([np.sin(a_kick), np.cos(a_kick), 0.0])
M2_0 = np.array([-np.sin(a_kick), np.cos(a_kick), 0.0])
N0 = Neq(theta_kick)
y0 = np.concatenate([M1_0, M2_0, [N0]])
print(f"theta_eq0={np.degrees(theta_eq0):.2f}deg, theta_kick={np.degrees(theta_kick):.2f}deg, N0={N0:.4e}")

t_max_ps = 400.0
t_max_nat = (t_max_ps*1e-12)/hbar_eVs

fig, axes = plt.subplots(1, 2, figsize=(13, 5.0))
colors = {260: "#2166ac", 500: "#762a83", 800: "#b2182b"}

for tauX in (260, 500, 800):
    sol = solve_ivp(lambda t,y: rhs_nonadiabatic(t,y,tauX), [0,t_max_nat], y0,
                     t_eval=np.linspace(0,t_max_nat,8000), method="RK45", rtol=1e-9, atol=1e-12)
    print(f"tau_X={tauX}ps: success={sol.success}, msg={sol.message}")
    M1=sol.y[0:3,:]; M2=sol.y[3:6,:]; N_t=sol.y[6,:]
    costh=np.clip(np.sum(M1*M2,axis=0),-1,1); theta_t=np.degrees(np.arccos(costh))
    t_ps = sol.t*hbar_eVs*1e12
    axes[0].plot(t_ps, theta_t, "-", color=colors[tauX], lw=1.0, label=f"$\\tau_X$={tauX}ps")
    axes[1].plot(t_ps, N_t, "-", color=colors[tauX], lw=1.0, label=f"$\\tau_X$={tauX}ps")
    n=len(theta_t); q=n//4
    amp_first=theta_t[:q].max()-theta_t[:q].min(); amp_last=theta_t[-q:].max()-theta_t[-q:].min()
    print(f"  amp(0-100ps)={amp_first:.3f}deg amp(300-400ps)={amp_last:.3f}deg ratio={amp_last/amp_first:.4f}  "
          f"theta range [{theta_t.min():.1f},{theta_t.max():.1f}]deg")

axes[0].axhline(np.degrees(theta_eq0), color="k", ls="--", lw=1, label=r"bare eq. $\theta_0$")
axes[0].set_xlabel("t (ps)"); axes[0].set_ylabel(r"$\theta(t)$ (deg)")
axes[0].set_title(f"NON-adiabatic: N(t) dynamical, +{kick_deg}$^\\circ$ kick")
axes[0].legend(frameon=False, fontsize=8); axes[0].grid(alpha=0.25)

axes[1].axhline(Neq(theta_kick), color="gray", ls=":", lw=1, label=r"$N_{eq}(\theta_{kick})$")
axes[1].set_xlabel("t (ps)"); axes[1].set_ylabel(r"$N(t)$")
axes[1].set_title("Exciton population lag")
axes[1].legend(frameon=False, fontsize=8); axes[1].grid(alpha=0.25)

plt.tight_layout()
outpath = os.path.join(here, "two_spin_nonadiabatic.png")
plt.savefig(outpath, dpi=150)
print("saved", outpath)
