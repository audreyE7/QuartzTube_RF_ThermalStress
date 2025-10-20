"""
Analytical & semi-analytical utilities for a quartz tube heated by an inner RF coil.
- Closed-form steady 1D temperature field T(r) for a hollow cylinder with inner heat flux and outer convection
- Axisymmetric thermoelastic hoop/axial stress via 1D finite-difference (plane strain)
Author: Audrey E. (scaffold generated)
"""

import json
import numpy as np

def load_material(json_path):
    with open(json_path, "r") as f:
        m = json.load(f)
    # Convert to SI
    E = m["E_GPa"] * 1e9
    nu = m["nu"]
    alpha = m["alpha_1e6_perK"] * 1e-6
    k = m["k_W_per_mK"]
    sigf = m["sigma_f_MPa"] * 1e6
    rho = m["rho_kg_per_m3"]
    Cp = m["Cp_J_per_kgK"]
    T_inf = m.get("T_inf_C", 25.0)
    return dict(E=E, nu=nu, alpha=alpha, k=k, sigma_f=sigf, rho=rho, Cp=Cp, T_inf=T_inf)

def steady_T_profile(a, b, q_in, h, T_inf, k, npts=400):
    """
    Solve T(r) = C1 ln r + C2 subject to:
      -k dT/dr|_{r=a} = q_in
      -k dT/dr|_{r=b} = h ( T(b) - T_inf )
    Returns r (m) and T(r) in C (relative to T_inf).
    """
    C1 = - q_in * a / k
    C2 = (-k*C1/(b*h)) - C1*np.log(b)
    r = np.linspace(a, b, npts)
    T = C1*np.log(r) + C2 + T_inf
    return r, T

def thermoelastic_stress_plane_strain(r, T, mat, a, b):
    """
    Axisymmetric thermoelastic solution via numerical 1D finite-difference.
    Returns sigma_r, sigma_t, sigma_z arrays.
    """
    E, nu, alpha = mat["E"], mat["nu"], mat["alpha"]
    n = len(r)
    dr = np.gradient(r)
    lam = E*nu/((1+nu)*(1-2*nu))
    mu = E/(2*(1+nu))
    dT = T - T[0]
    A = np.zeros((n, n))
    rhs = np.zeros(n)

    for i in range(1, n-1):
        ri = r[i]; rim = r[i-1]; rip = r[i+1]
        drm = ri - rim; drp = rip - ri
        K = lam + 2*mu
        A[i, i-1] += K / (drm*(0.5*(drm+drp)))
        A[i, i]   += -K * (1/drm + 1/drp) / (0.5*(drm+drp))
        A[i, i+1] += K / (drp*(0.5*(drm+drp)))
        A[i, i]   += -2*mu / (ri**2)
        coeff = E*alpha/(1-2*nu)
        dTdr = (T[i+1]-T[i-1])/(drm+drp)
        rhs[i] = coeff * dTdr

    # Boundary conditions (σr(a)=0, σr(b)=0)
    i=0; K = lam + 2*mu
    A[i,i] = -K/(r[1]-r[0])
    A[i,1] =  K/(r[1]-r[0])
    A[i,i] += lam / r[i]
    rhs[i] =  (E*alpha/(1-2*nu)) * (T[i]-T[0])
    i=n-1
    A[i,i] =  K/(r[-1]-r[-2])
    A[i,n-2] = -K/(r[-1]-r[-2])
    A[i,i] += lam / r[i]
    rhs[i] =  (E*alpha/(1-2*nu)) * (T[i]-T[0])

    u = np.linalg.solve(A, rhs)
    du = np.gradient(u, r)
    er = du; et = u / r
    Kbulk = E/(3*(1-2*nu))
    sigma_r = 2*mu*er + lam*(er+et) - 3*Kbulk*mat["alpha"]*dT
    sigma_t = 2*mu*et + lam*(er+et) - 3*Kbulk*mat["alpha"]*dT
    sigma_z = lam*(er+et) - 3*Kbulk*mat["alpha"]*dT
    return u, sigma_r, sigma_t, sigma_z

def thermal_shock_parameter(mat):
    """ R'' = k * sigma_f / (E * alpha) """
    return mat["k"] * mat["sigma_f"] / (mat["E"] * mat["alpha"])

def safe_window_example(a=0.0125, t=0.0025, q_list=None, h_list=None, material_json="data/material_quartz.json"):
    if q_list is None: q_list = np.linspace(2e4, 2e5, 8)
    if h_list is None: h_list = np.linspace(10, 200, 7)
    b = a + t
    m = load_material(material_json)
    E, nu, alpha, k = m["E"], m["nu"], m["alpha"], m["k"]
    T_inf = m["T_inf"]
    Smax = np.zeros((len(h_list), len(q_list)))
    for i,h in enumerate(h_list):
        for j,q in enumerate(q_list):
            r, T = steady_T_profile(a, b, q_in=q, h=h, T_inf=T_inf, k=k, npts=600)
            u, sr, st, sz = thermoelastic_stress_plane_strain(r, T, dict(E=E,nu=nu,alpha=alpha), a, b)
            Smax[i,j] = np.max(np.abs(st))/1e6
    Rpp = thermal_shock_parameter(m)/1e6
    QkW = np.array(q_list)/1000.0
    return QkW, np.array(h_list), Smax, Rpp
