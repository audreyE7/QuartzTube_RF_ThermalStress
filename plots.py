import json, os
import numpy as np
import matplotlib.pyplot as plt
from analytic_quartz_cylinder import load_material, steady_T_profile, thermoelastic_stress_plane_strain, safe_window_example

def plot_profile(a=0.0125, t=0.0025, q=1.2e5, h=60.0, material_json="../data/material_quartz.json", outdir="../figs"):
    os.makedirs(outdir, exist_ok=True)
    m = load_material(material_json)
    r, T = steady_T_profile(a, a+t, q_in=q, h=h, T_inf=m["T_inf"], k=m["k"], npts=600)
    u, sr, st, sz = thermoelastic_stress_plane_strain(r, T, dict(E=m["E"],nu=m["nu"],alpha=m["alpha"]), a, a+t)

    plt.figure(); plt.plot(r*1e3, T, linewidth=2)
    plt.xlabel("Radius r (mm)"); plt.ylabel("Temperature (°C)")
    plt.title("Steady Temperature Across Tube Wall"); plt.grid(True, alpha=0.3); plt.tight_layout()
    plt.savefig(os.path.join(outdir, "temperature_profile.png"), dpi=180); plt.close()

    plt.figure(); plt.plot(r*1e3, st/1e6, linewidth=2)
    plt.xlabel("Radius r (mm)"); plt.ylabel("Hoop Stress σθ (MPa)")
    plt.title("Hoop Stress Across Tube Wall (Plane Strain)"); plt.grid(True, alpha=0.3); plt.tight_layout()
    plt.savefig(os.path.join(outdir, "hoop_stress_profile.png"), dpi=180); plt.close()

def plot_safe_map(a=0.0125, t=0.0025, material_json="../data/material_quartz.json", outdir="../figs"):
    os.makedirs(outdir, exist_ok=True)
    QkW, H, Smax, Rpp = safe_window_example(a=a, t=t, material_json=material_json)
    QQ, HH = np.meshgrid(QkW, H)
    plt.figure()
    cs = plt.contourf(QQ, HH, Smax, levels=12)
    plt.colorbar(cs, label="Max Hoop Stress (MPa)")
    plt.contour(QQ, HH, Smax, levels=[Rpp], linewidths=2)
    plt.xlabel("Inner heat flux q'' (kW/m²)"); plt.ylabel("External h (W/m²·K)")
    plt.title("Safe Operating Map (σθ vs q'' and h)"); plt.tight_layout()
    plt.savefig(os.path.join(outdir, "safe_operating_map.png"), dpi=180); plt.close()

if __name__ == "__main__":
    plot_profile()
    plot_safe_map()
