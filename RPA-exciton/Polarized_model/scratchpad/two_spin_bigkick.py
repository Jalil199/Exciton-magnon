"""
Same 2-macrospin LLG (real literature J_int/A_x/A_z, B=0.2T) + our j_theta(theta), eta_theta(theta),
but now with a LARGE initial kick away from the bare equilibrium (instead of starting exactly at
theta_eq0), to test for amplitude-threshold behavior (cf. Brennan Fig 3a: 10uW decays, 130uW
self-sustains -- an amplitude/fluence-dependent effect). Single panel, theta(t) only, styled like
their Fig 3a (canting angle in degrees vs pump-probe-like delay in ps).
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
print(f"theta_eq0 = {np.degrees(theta_eq0):.3f} deg")

rows = list(csv.DictReader(open(os.path.join(here, "step10_damping_data.csv"))))
theta_arr = np.array([float(r["theta_over_pi"]) for r in rows]) * np.pi
tau1_eV = np.array([float(r["tau1_eV"]) for r in rows])
eta_by_tauX = {260: np.array([float(r["eta_260ps"]) for r in rows]),
               500: np.array([float(r["eta_500ps"]) for r in rows]),
               800: np.array([float(r["eta_800ps"]) for r in rows])}
jtheta_arr = tau1_eV/np.sin(theta_arr)
print(f"Step10 theta range covered: [{theta_arr.min()/np.pi:.3f}, {theta_arr.max()/np.pi:.3f}] pi "
      f"= [{np.degrees(theta_arr.min()):.1f}, {np.degrees(theta_arr.max()):.1f}] deg")

def jtheta(t): return np.interp(t, theta_arr, jtheta_arr)
def make_eta(tauX):
    arr = eta_by_tauX[tauX]
    return lambda t: np.interp(t, theta_arr, arr)

def rhs_full(t, y, eta_fn):
    M1=y[0:3]; M2=y[3:6]; M1=M1/np.linalg.norm(M1); M2=M2/np.linalg.norm(M2)
    costheta=np.clip(np.dot(M1,M2),-1,1); theta=np.arccos(costheta)
    jth = jtheta(theta); eta = eta_fn(theta)
    B1 = Bcl(M1,M2) + jth*M2; B2 = Bcl(M2,M1) + jth*M1
    A = np.eye(6); A[0:3,3:6]=-eta*crossmat(M1); A[3:6,0:3]=-eta*crossmat(M2)
    b = np.concatenate([np.cross(M1,B1), np.cross(M2,B2)])
    return np.linalg.solve(A,b)

# ---------------- LARGE initial kick: theta_0 = theta_eq0 + Delta (Delta ~ Brennan's ~34deg swing) ----
kick_deg = 30.0
theta_kick = theta_eq0 + np.radians(kick_deg)
a_kick = theta_kick/2
print(f"large-kick initial condition: theta_0 = theta_eq0 + {kick_deg} deg = {np.degrees(theta_kick):.2f} deg")

M1_0 = np.array([np.sin(a_kick), np.cos(a_kick), 0.0])
M2_0 = np.array([-np.sin(a_kick), np.cos(a_kick), 0.0])
y0 = np.concatenate([M1_0, M2_0])

t_max_ps = 400.0
t_max_nat = (t_max_ps*1e-12)/hbar_eVs

fig, ax = plt.subplots(figsize=(8.5, 5.2))
colors = {260: "#2166ac", 500: "#762a83", 800: "#b2182b"}

for tauX in (260, 500, 800):
    eta_fn = make_eta(tauX)
    sol = solve_ivp(lambda t,y: rhs_full(t,y,eta_fn), [0,t_max_nat], y0,
                     t_eval=np.linspace(0,t_max_nat,8000), method="RK45", rtol=1e-9, atol=1e-11)
    print(f"tau_X={tauX}ps: success={sol.success}, msg={sol.message}")
    M1=sol.y[0:3,:]; M2=sol.y[3:6,:]
    costh=np.clip(np.sum(M1*M2,axis=0),-1,1); theta_t=np.degrees(np.arccos(costh))
    t_ps = sol.t*hbar_eVs*1e12
    ax.plot(t_ps, theta_t, "-", color=colors[tauX], lw=1.1, label=f"$\\tau_X$={tauX}ps")
    n=len(theta_t); q=n//4
    amp_first=theta_t[:q].max()-theta_t[:q].min(); amp_last=theta_t[-q:].max()-theta_t[-q:].min()
    print(f"  amp(0-100ps)={amp_first:.3f}deg, amp(300-400ps)={amp_last:.3f}deg, ratio={amp_last/amp_first:.4f}")

ax.axhline(np.degrees(theta_eq0), color="k", ls="--", lw=1, label=r"bare eq. $\theta_0$")
ax.axhline(np.degrees(theta_kick), color="gray", ls=":", lw=1, label=f"initial kick (+{kick_deg}$^\\circ$)")
ax.set_xlabel("t (ps)"); ax.set_ylabel(r"$\theta(t)$ (deg)")
ax.set_title(f"Large-kick ($+{kick_deg}^\\circ$) canting dynamics -- cf. Brennan Fig 3a")
ax.legend(frameon=False, fontsize=9)
ax.grid(alpha=0.25)
plt.tight_layout()
outpath = os.path.join(here, "two_spin_bigkick.png")
plt.savefig(outpath, dpi=150)
print("saved", outpath)
