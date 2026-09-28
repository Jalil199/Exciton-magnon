"""
Two-macrospin LLG dynamics (Eq. 1 + classical-Hamiltonian structure of Eq. 22 in the draft
paper), using OUR microscopically-computed reactive torque j_theta(theta) (Step 9) and
population-lag damping eta_theta(theta) (Step 10), to see whether the mechanism alone
reproduces Brennan et al. 2026 Fig 3a-like self-sustained canting-angle dynamics.

PARAMETER FIX (2026-07-24): the draft's literal J_AF=19.5meV, Kx=21meV, Kz=57meV give a bare
(no exciton) precession frequency of ~90-96 meV (~22 THz) -- confirmed by a normal-mode
Jacobian calculation (normal_modes.py). That is off by a factor ~1000 from the real CrSBr
magnon frequency (Scheie et al Adv Sci 2022 / Brennan 2026: 0.102-0.141 meV, 24.7-34.1 GHz).

Rather than an arbitrary rescale, we use the ACTUAL literature-fitted macrospin parameters
from "Tunable magnons in a dual-gated 2D antiferromagnet" (LSWT fit):
    A_x = 14 ueV, A_z = 58 ueV, J_int = 6 ueV
(note A_z=58ueV is numerically the same as the draft's Kz=57meV -- almost certainly a
mu-eV/meV unit slip in the draft, not a different physical model). That paper reports, AT
B_0=0.2T: nu_IP=24.4 GHz, nu_OP=33.8 GHz -- so B=0.2T (not Brennan's 0.7T) is the field value
this specific parameter set is validated against; 0.7T is above the ~0.4T interlayer
spin-flop field for these J_int,A_x and pushes the simple smooth-canting ansatz used here
(cos(theta/2)=... single-valued equilibrium) out of its regime. We therefore use B=0.2T,
flagged explicitly -- not silently substituted for Brennan's value.

Units: natural units (hbar=1) for the ODE; all Hamiltonian terms, j_theta, eta_theta are in
eV, used directly as rates (dM/dt = M x [...]). Time -> real ps only for plotting/printing.

Damping: nonlocal only, eta^12=eta^21=eta_theta(theta), eta^11=eta^22=0 (matches the paper's
"nonlocal damping" framing and how eta_theta was derived, from a two-layer quantity).
"""
import numpy as np
from scipy.integrate import solve_ivp
import csv, os, time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
hbar_eVs = 6.582119569e-16
muB_eVperT = 5.7883818060e-5
ex = np.array([1.0,0.0,0.0]); ey = np.array([0.0,1.0,0.0]); ez = np.array([0.0,0.0,1.0])

# ---------------- REAL literature parameters (dual-gated CrSBr magnon paper, LSWT fit) --------
J_AF = 6e-6           # interlayer AFM exchange J_int (eV)
Kx   = 14e-6          # A_x anisotropy (eV)   -- plays the "easy axis" role in Eq.22's structure
Kz   = 58e-6          # A_z anisotropy (eV)   -- "hard axis"
g_factor = 2.0
B_tesla = 0.2         # field value this exact (A_x,A_z,J_int) set is validated against (24.4/33.8 GHz)
gmuB_B = g_factor * muB_eVperT * B_tesla

def crossmat(v):
    return np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])

def Bcl(M, Mbar):
    return -J_AF*Mbar + gmuB_B*ey - 2*Kz*np.dot(M,ez)*ez + 2*Kx*np.dot(M,ex)*ex

cos_a_eq = gmuB_B/(2*(J_AF+Kx))
print(f"J_AF={J_AF*1e6:.2f}ueV, Kx={Kx*1e6:.2f}ueV, Kz={Kz*1e6:.2f}ueV, gmuB*B={gmuB_B*1e6:.2f}ueV (B={B_tesla}T)")
print(f"cos(a_eq) = {cos_a_eq:.4f}  ({'valid canted state' if abs(cos_a_eq)<1 else 'INVALID -- saturated/above spin-flop'})")
a_eq = np.arccos(np.clip(cos_a_eq,-1,1))
theta_eq0 = 2*a_eq
print(f"equilibrium theta_eq0 = {np.degrees(theta_eq0):.3f} deg = {theta_eq0/np.pi:.4f} pi")

# ---------------- validate: normal-mode Jacobian should now land near 24-34 GHz ----------------
def tangent_basis(M0):
    tmp = np.array([1.0,0.0,0.0]) if abs(M0[0])<0.9 else np.array([0.0,1.0,0.0])
    e_a = np.cross(M0,tmp); e_a/=np.linalg.norm(e_a)
    e_b = np.cross(M0,e_a); e_b/=np.linalg.norm(e_b)
    return e_a,e_b

M1_0 = np.array([np.sin(a_eq), np.cos(a_eq), 0.0])
M2_0 = np.array([-np.sin(a_eq), np.cos(a_eq), 0.0])
e1a,e1b = tangent_basis(M1_0); e2a,e2b = tangent_basis(M2_0)

def rhs0(y):
    M1=y[0:3]; M2=y[3:6]; M1n=M1/np.linalg.norm(M1); M2n=M2/np.linalg.norm(M2)
    return np.concatenate([np.cross(M1n,Bcl(M1n,M2n)), np.cross(M2n,Bcl(M2n,M1n))])

def y_of_x(x):
    M1=M1_0+x[0]*e1a+x[1]*e1b; M2=M2_0+x[2]*e2a+x[3]*e2b
    M1/=np.linalg.norm(M1); M2/=np.linalg.norm(M2)
    return np.concatenate([M1,M2])

def tv(x):
    v=rhs0(y_of_x(x)); v1,v2=v[0:3],v[3:6]
    return np.array([np.dot(v1,e1a),np.dot(v1,e1b),np.dot(v2,e2a),np.dot(v2,e2b)])

h=1e-7; Jac=np.zeros((4,4))
for j in range(4):
    xp=np.zeros(4); xp[j]=h; xm=np.zeros(4); xm[j]=-h
    Jac[:,j]=(tv(xp)-tv(xm))/(2*h)
evals = np.linalg.eigvals(Jac)
omegas_eV = sorted(set(np.round(np.abs(evals.imag[np.abs(evals.imag)>1e-14]),14)))
print("bare normal modes (validation against 24.4/33.8 GHz):")
for w in omegas_eV:
    f_GHz = w/(2*np.pi*hbar_eVs)/1e9
    print(f"  omega={w*1e6:.4f} ueV -> f={f_GHz:.3f} GHz  (period={1e12/(f_GHz*1e9):.3f} ps)")

# ---------------- load Step 9/10 data (j_theta, eta_theta vs theta) ----------------
rows = list(csv.DictReader(open(os.path.join(here, "step10_damping_data.csv"))))
theta_arr = np.array([float(r["theta_over_pi"]) for r in rows]) * np.pi
tau1_eV = np.array([float(r["tau1_eV"]) for r in rows])
eta_by_tauX = {260: np.array([float(r["eta_260ps"]) for r in rows]),
               500: np.array([float(r["eta_500ps"]) for r in rows]),
               800: np.array([float(r["eta_800ps"]) for r in rows])}
jtheta_arr = tau1_eV/np.sin(theta_arr)
print(f"\ntheta_eq0/pi={theta_eq0/np.pi:.4f} vs Step10 data range [0.08,0.92]pi: "
      f"{'INSIDE' if 0.08<=theta_eq0/np.pi<=0.92 else 'OUTSIDE -- will clamp/extrapolate flat'}")

def jtheta(theta): return np.interp(theta, theta_arr, jtheta_arr)
def make_eta(tauX):
    arr = eta_by_tauX[tauX]
    return lambda theta: np.interp(theta, theta_arr, arr)

def rhs_full(t, y, eta_fn):
    M1=y[0:3]; M2=y[3:6]
    M1=M1/np.linalg.norm(M1); M2=M2/np.linalg.norm(M2)
    costheta=np.clip(np.dot(M1,M2),-1,1); theta=np.arccos(costheta)
    jth = jtheta(theta); eta = eta_fn(theta)
    B1 = Bcl(M1,M2) + jth*M2
    B2 = Bcl(M2,M1) + jth*M1
    A = np.eye(6)
    A[0:3,3:6] = -eta*crossmat(M1); A[3:6,0:3] = -eta*crossmat(M2)
    b = np.concatenate([np.cross(M1,B1), np.cross(M2,B2)])
    return np.linalg.solve(A,b)

y0 = np.concatenate([M1_0, M2_0])

# period is now ~1/(25-34GHz) ~ 30-40ps -- run a few tens of periods
t_max_ps = 400.0
t_max_nat = (t_max_ps*1e-12)/hbar_eVs
n_out = 6000

fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
colors = {260: "#2166ac", 500: "#762a83", 800: "#b2182b"}

for tauX in (260, 500, 800):
    eta_fn = make_eta(tauX)
    t0 = time.time()
    sol = solve_ivp(lambda t,y: rhs_full(t,y,eta_fn), [0,t_max_nat], y0,
                     t_eval=np.linspace(0,t_max_nat,n_out), method="RK45", rtol=1e-8, atol=1e-10)
    print(f"tau_X={tauX}ps: wall={time.time()-t0:.1f}s, success={sol.success}, msg={sol.message}")
    M1=sol.y[0:3,:]; M2=sol.y[3:6,:]
    costh=np.clip(np.sum(M1*M2,axis=0),-1,1); theta_t=np.degrees(np.arccos(costh))
    t_ps = sol.t*hbar_eVs*1e12
    axes[0].plot(t_ps, theta_t, "-", color=colors[tauX], lw=1.0, label=f"$\\tau_X$={tauX}ps")
    norm1 = np.linalg.norm(M1,axis=0)
    axes[1].plot(t_ps, norm1-1, "-", color=colors[tauX], lw=0.8, label=f"{tauX}ps")

axes[0].axhline(np.degrees(theta_eq0), color="k", ls="--", lw=1, label=r"bare eq. $\theta_0$")
axes[0].set_xlabel("t (ps)"); axes[0].set_ylabel(r"$\theta(t)$ (deg)")
axes[0].set_title("Canting dynamics, real literature J_AF/Kx/Kz/B (0.2T)")
axes[0].legend(frameon=False, fontsize=8); axes[0].grid(alpha=0.25)

axes[1].set_xlabel("t (ps)"); axes[1].set_ylabel(r"$|M_1|-1$")
axes[1].set_title("Numerical sanity: unit-norm drift")
axes[1].legend(frameon=False, fontsize=7); axes[1].grid(alpha=0.25)

plt.tight_layout()
outpath = os.path.join(here, "two_spin_dynamics.png")
plt.savefig(outpath, dpi=150)
print("saved", outpath)
