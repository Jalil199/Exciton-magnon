"""
Linearize the BARE (no exciton torque/damping) 2-macrospin LLG dynamics around the classical
equilibrium to get the actual AFM resonance frequencies -- BEFORE committing to any long time
integration. This replaces guessing at a timescale with a direct, cheap (no ODE solve) answer.
"""
import numpy as np

J_AF = 0.0195; Kx = 0.021; Kz = 0.057; gmuB_B = 0.042
hbar_eVs = 6.582119569e-16
ex = np.array([1.0,0.0,0.0]); ey = np.array([0.0,1.0,0.0]); ez = np.array([0.0,0.0,1.0])

cos_a_eq = gmuB_B/(2*(J_AF+Kx))
a_eq = np.arccos(np.clip(cos_a_eq,-1,1))
theta_eq0 = 2*a_eq
print(f"theta_eq0 = {np.degrees(theta_eq0):.3f} deg = {theta_eq0/np.pi:.4f} pi")

M1_0 = np.array([np.sin(a_eq), np.cos(a_eq), 0.0])
M2_0 = np.array([-np.sin(a_eq), np.cos(a_eq), 0.0])

def Bcl(M, Mbar):
    return -J_AF*Mbar + gmuB_B*ey - 2*Kz*np.dot(M,ez)*ez + 2*Kx*np.dot(M,ex)*ex

def rhs_notorque(y):
    M1 = y[0:3]; M2 = y[3:6]
    M1n = M1/np.linalg.norm(M1); M2n = M2/np.linalg.norm(M2)
    B1 = Bcl(M1n,M2n); B2 = Bcl(M2n,M1n)
    return np.concatenate([np.cross(M1n,B1), np.cross(M2n,B2)])

# sanity: equilibrium should give (near) zero velocity
y0 = np.concatenate([M1_0, M2_0])
v0 = rhs_notorque(y0)
print(f"velocity at equilibrium (should be ~0): |v0| = {np.linalg.norm(v0):.3e}")

# tangent bases at each site (two vectors perpendicular to M_L0)
def tangent_basis(M0):
    # pick any vector not parallel to M0
    tmp = np.array([1.0,0.0,0.0]) if abs(M0[0]) < 0.9 else np.array([0.0,1.0,0.0])
    e_a = np.cross(M0, tmp); e_a /= np.linalg.norm(e_a)
    e_b = np.cross(M0, e_a); e_b /= np.linalg.norm(e_b)
    return e_a, e_b

e1a, e1b = tangent_basis(M1_0)
e2a, e2b = tangent_basis(M2_0)

def y_of_x(x):
    # x = (x1a, x1b, x2a, x2b) small tangent displacements
    M1 = M1_0 + x[0]*e1a + x[1]*e1b
    M2 = M2_0 + x[2]*e2a + x[3]*e2b
    M1 /= np.linalg.norm(M1); M2 /= np.linalg.norm(M2)
    return np.concatenate([M1, M2])

def tangent_velocity(x):
    y = y_of_x(x)
    v = rhs_notorque(y)
    v1, v2 = v[0:3], v[3:6]
    return np.array([np.dot(v1,e1a), np.dot(v1,e1b), np.dot(v2,e2a), np.dot(v2,e2b)])

# numeric 4x4 Jacobian d(tangent_velocity)/d(x) at x=0
h = 1e-6
Jac = np.zeros((4,4))
for j in range(4):
    xp = np.zeros(4); xp[j] = h
    xm = np.zeros(4); xm[j] = -h
    Jac[:, j] = (tangent_velocity(xp) - tangent_velocity(xm)) / (2*h)

evals, evecs = np.linalg.eig(Jac)
print("\nJacobian eigenvalues (natural units, eV, should be purely imaginary +-i*omega for a conservative system):")
for lam in evals:
    print(f"  {lam:.6e}")

omegas_eV = sorted(set(np.round(abs(evals[np.abs(evals.imag)>1e-12].imag), 10)))
print("\ndistinct |omega| (eV):", omegas_eV)
for w in omegas_eV:
    if w > 0:
        f_Hz = w / (2*np.pi*hbar_eVs)
        period_ps = 1/f_Hz * 1e12
        print(f"  omega={w:.6e} eV  ->  f={f_Hz/1e9:.3f} GHz  ->  period={period_ps:.4f} ps")
