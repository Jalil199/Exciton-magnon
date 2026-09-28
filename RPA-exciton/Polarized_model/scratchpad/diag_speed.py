import numpy as np
from scipy.integrate import solve_ivp
import csv, os, time

here = os.path.dirname(os.path.abspath(__file__))

J_AF = 0.0195; Kx = 0.021; Kz = 0.057; gmuB_B = 0.042
hbar_eVs = 6.582119569e-16
ex = np.array([1.0,0.0,0.0]); ey = np.array([0.0,1.0,0.0]); ez = np.array([0.0,0.0,1.0])

rows = list(csv.DictReader(open(os.path.join(here, "step10_damping_data.csv"))))
theta_arr = np.array([float(r["theta_over_pi"]) for r in rows]) * np.pi
tau1_eV = np.array([float(r["tau1_eV"]) for r in rows])
eta500_arr = np.array([float(r["eta_500ps"]) for r in rows])
jtheta_arr = tau1_eV/np.sin(theta_arr)

def jtheta(theta):
    return np.interp(theta, theta_arr, jtheta_arr)
def eta500(theta):
    return np.interp(theta, theta_arr, eta500_arr)

cos_a_eq = gmuB_B/(2*(J_AF+Kx)); a_eq = np.arccos(np.clip(cos_a_eq,-1,1)); theta_eq0 = 2*a_eq
print(f"theta_eq0 = {np.degrees(theta_eq0):.3f} deg = {theta_eq0/np.pi:.4f} pi", flush=True)

def crossmat(v):
    return np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
def Bcl(M, Mbar):
    return -J_AF*Mbar + gmuB_B*ey - 2*Kz*np.dot(M,ez)*ez + 2*Kx*np.dot(M,ex)*ex

def rhs_notorque(t, y):
    # bare classical LLG, NO exciton torque, NO damping -- pure precession, to measure the
    # intrinsic oscillation period implied by Eq.22 alone (sanity check vs literature ~26.7GHz)
    M1 = y[0:3]; M2 = y[3:6]
    M1 = M1/np.linalg.norm(M1); M2 = M2/np.linalg.norm(M2)
    B1 = Bcl(M1, M2); B2 = Bcl(M2, M1)
    return np.concatenate([np.cross(M1,B1), np.cross(M2,B2)])

def rhs_full(t, y, eta_fn):
    M1 = y[0:3]; M2 = y[3:6]
    M1 = M1/np.linalg.norm(M1); M2 = M2/np.linalg.norm(M2)
    costheta = np.clip(np.dot(M1,M2),-1,1); theta = np.arccos(costheta)
    jth = jtheta(theta); eta = eta_fn(theta)
    B1 = Bcl(M1,M2) + jth*M2; B2 = Bcl(M2,M1) + jth*M1
    A = np.eye(6)
    A[0:3,3:6] = -eta*crossmat(M1); A[3:6,0:3] = -eta*crossmat(M2)
    b = np.concatenate([np.cross(M1,B1), np.cross(M2,B2)])
    return np.linalg.solve(A,b)

# perturb slightly off equilibrium (M1,M2 tilted a bit out of plane) to seed a visible precession
M1_0 = np.array([np.sin(a_eq), np.cos(a_eq), 0.02])
M2_0 = np.array([-np.sin(a_eq), np.cos(a_eq), -0.02])
M1_0 /= np.linalg.norm(M1_0); M2_0 /= np.linalg.norm(M2_0)
y0 = np.concatenate([M1_0, M2_0])

print("=== bare precession (no torque/damping), find real oscillation period ===", flush=True)
t_test_ps = 50.0
t_test_nat = (t_test_ps*1e-12)/hbar_eVs
t0 = time.time()
sol = solve_ivp(rhs_notorque, [0, t_test_nat], y0, t_eval=np.linspace(0,t_test_nat,4000),
                 method="RK45", rtol=1e-8, atol=1e-10)
print(f"  wall time: {time.time()-t0:.2f}s, npts={sol.t.size}, success={sol.success}", flush=True)
M1 = sol.y[0:3,:]; M2 = sol.y[3:6,:]
costh = np.clip(np.sum(M1*M2,axis=0),-1,1); theta_t = np.degrees(np.arccos(costh))
t_ps = sol.t*hbar_eVs*1e12
print(f"  theta(t) range over {t_test_ps}ps: min={theta_t.min():.4f} max={theta_t.max():.4f} deg", flush=True)
# crude period estimate: count zero-crossings of (theta - mean)
dev = theta_t - theta_t.mean()
crossings = np.where(np.diff(np.sign(dev)))[0]
if len(crossings) >= 2:
    half_periods = np.diff(t_ps[crossings])
    print(f"  ~{len(crossings)} zero-crossings in {t_test_ps}ps -> est. period ~ {2*np.median(half_periods):.4f} ps "
          f"-> ~{1/(2*np.median(half_periods)*1e-12)/1e9:.2f} GHz", flush=True)
else:
    print("  fewer than 2 crossings found in this window -- period may be longer than test window", flush=True)

print("\n=== timing a short full run (with torque+damping, tau_X=500ps) ===", flush=True)
t0 = time.time()
sol2 = solve_ivp(lambda t,y: rhs_full(t,y,eta500), [0, t_test_nat], y0,
                  t_eval=np.linspace(0,t_test_nat,2000), method="RK45", rtol=1e-7, atol=1e-9)
dt = time.time()-t0
print(f"  wall time for {t_test_ps}ps: {dt:.2f}s, npts={sol2.t.size}, success={sol2.success}", flush=True)
print(f"  extrapolated wall time for 1500ps: ~{dt*(1500/t_test_ps):.1f}s", flush=True)
