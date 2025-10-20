# ANSYS Workbench Notes — Quartz Tube RF Heating (Axisymmetric)

**Geometry:** 2D axisymmetric hollow cylinder, inner radius `a`, thickness `t`.  
**BCs:**
- Inner wall: Heat flux `q''` (from RF coil).
- Outer wall: Convection `h` to `T_inf`.
- Structural: Symmetry on axis; σ_r(a)=σ_r(b)=0. Plane strain in z.

**Mesh:** Bias toward inner wall; 50–100 elements across thickness.

**Outputs:** T(r), σ_r(r), σ_θ(r), σ_z(r). Compare analytic vs FEA.
