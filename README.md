# navier-stokes-lab

Educational Navier–Stokes laboratory in Python: exact solutions, finite-difference and spectral examples, lid-driven cavity flow (ψ–ω and projection solvers), Taylor–Green vortex, boundary layers. 34 diagrams, including 3D surfaces.

**Author:** Harpreet Singh

> **Scope:** this is an educational project, not production CFD. Grids are small, the cavity solvers are teaching versions, and results have not been verified against full benchmark datasets. Do not use them for engineering decisions without grid-convergence and benchmark checks. The `d3_*` 3D views are synthetic teaching illustrations, and `extra_03` shows an analytic field, not a simulation.

## Contents

| File | Purpose |
|---|---|
| `Navier_Stokes_Computational_Lab_v23.ipynb` | Main notebook: theory, equation library, numerical methods, cavity CFD, comparison tables |
| `Navier_Stokes_Diagram_Gallery_v6.ipynb` | Diagram notebook: exact solutions, Burgers, ψ–ω vs MAC cavity solvers, spectral Taylor–Green, Blasius, validation dashboard |
| `navier_stokes.py` | Exact solutions, projection cavity demo, Taylor–Green decay |
| `cavity_cfd.py` | Streamfunction–vorticity cavity solver (`simulate_cavity`) |
| `plot_vorticity.py` | Regenerate 2D/3D vorticity plots from saved fields |
| `viz.py` | Matplotlib helpers |
| `generate_extra_diagrams.py` | Regenerates the `extra_*` diagrams |
| `test_numerics.py` | Pytest tests |
| `sim_cavity_Re100_fields.npz` | Saved Re = 100 cavity fields |
| `cfd_methods_comparison.csv` | CFD methods comparison table |

Image files use prefixes: `sim_` (Re = 100 cavity simulation), `fig` (gallery notebook), `d3_` (3D teaching views), `extra_` (additional diagrams).

## Quick start

```bash
git clone https://github.com/<your-username>/navier-stokes-lab.git
cd navier-stokes-lab
python -m venv .venv
# Windows: .venv\Scripts\activate     Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
jupyter lab
pytest -q
```

Regenerate figures:

```bash
python plot_vorticity.py
python generate_extra_diagrams.py
```

In the gallery notebook, `QUICK_MODE = True` runs in a few seconds on coarse grids. Set it to `False` for finer grids.

## Diagrams

### Cavity flow simulation (Re = 100, 41×41, streamfunction–vorticity)

**Velocity magnitude with streamlines**

![Velocity magnitude with streamlines](sim_01_cavity_velocity_streamlines.png)

**Streamfunction contours**

![Streamfunction contours](sim_02_cavity_streamfunction.png)

**Vorticity field**

![Vorticity field](sim_03_cavity_vorticity.png)

**3D velocity magnitude surface**

![3D velocity magnitude surface](sim_04_cavity_velocity_3d.png)

**3D vorticity surface**

![3D vorticity surface](sim_08_cavity_vorticity_3d.png)

**Vorticity centreline profiles**

![Vorticity centreline profiles](sim_09_cavity_vorticity_centerlines.png)

**Solver update history**

![Solver update history](sim_05_cavity_update_history.png)

**Centreline velocity profiles**

![Centreline velocity profiles](sim_06_cavity_centerline_profiles.png)

**Divergence diagnostic**

![Divergence diagnostic](sim_07_cavity_divergence_diagnostic.png)

**Speed field with velocity vectors**

![Speed field with velocity vectors](extra_01_cavity_speed_vectors.png)

**Streamlines coloured by speed**

![Streamlines coloured by speed](extra_02_cavity_colored_streamlines.png)


### Exact solutions and boundary layers

**Reynolds numbers and flow regimes**

![Reynolds numbers and flow regimes](fig01_reynolds_regimes.png)

**Poiseuille and Couette–Poiseuille: finite difference vs exact**

![Poiseuille and Couette–Poiseuille: finite difference vs exact](fig02_poiseuille_couette.png)

**Hagen–Poiseuille pipe flow (3D)**

![Hagen–Poiseuille pipe flow (3D)](fig03_hagen_poiseuille_3d.png)

**Stokes first problem (3D)**

![Stokes first problem (3D)](fig04_stokes_first_problem_3d.png)

**Lamb–Oseen vortex (3D)**

![Lamb–Oseen vortex (3D)](fig13_lamb_oseen_3d.png)

**Blasius boundary layer (3D)**

![Blasius boundary layer (3D)](fig14_blasius_boundary_layer_3d.png)


### Numerical analysis

**Truncation error and Burgers convergence**

![Truncation error and Burgers convergence](fig05_truncation_burgers_convergence.png)

**Burgers shock, grid comparison**

![Burgers shock, grid comparison](fig06_burgers_shock.png)

**Method error comparison (3D)**

![Method error comparison (3D)](d3_06_method_error_comparison_3d.png)


### Lid-driven cavity: solver comparison (ψ–ω vs MAC projection)

**Centreline validation**

![Centreline validation](fig07_cavity_centerline_validation.png)

**Speed, streamfunction, vorticity surfaces (Re = 400)**

![Speed, streamfunction, vorticity surfaces (Re = 400)](fig08_cavity_3d_surfaces.png)

**Velocity vectors on speed field**

![Velocity vectors on speed field](fig09_cavity_velocity_vectors.png)

**Residual histories**

![Residual histories](fig10_residual_history.png)

**Secondary corner vortices**

![Secondary corner vortices](fig11_secondary_vortices.png)


### Taylor–Green vortex and 3D teaching views

**Energy decay and dissipation rate**

![Energy decay and dissipation rate](fig12_taylor_green_energy.png)

**3D Taylor–Green vorticity slices (analytic field)**

![3D Taylor–Green vorticity slices (analytic field)](extra_03_taylor_green_3d_slices.png)

**Energy decay surface (3D)**

![Energy decay surface (3D)](d3_07_energy_decay_surface.png)

**Channel velocity vectors (3D)**

![Channel velocity vectors (3D)](d3_01_velocity_vector_field.png)

**Pressure surface (3D)**

![Pressure surface (3D)](d3_02_pressure_surface.png)

**Helical vortex tube (3D)**

![Helical vortex tube (3D)](d3_03_vortex_tube.png)

**Structured 3D mesh**

![Structured 3D mesh](d3_04_structured_mesh.png)

**Velocity profile surface (3D)**

![Velocity profile surface (3D)](d3_05_velocity_profile_surface.png)

**Vorticity surface (3D)**

![Vorticity surface (3D)](d3_08_vorticity_surface.png)


## Checks performed

- Poiseuille finite-difference solution matches the exact profile to within 1e-10.
- Blasius f''(0) = 0.33206 matches the known value.
- MAC projection solver: maximum divergence about 1e-14 (quick mode).
- Gallery validation dashboard: 9 of 9 checks pass.
- 5 pytest tests pass.

## Limitations and roadmap

- 2D laminar cavity only; no turbulence modelling.
- Coarse quick-mode grids (cavity ψ_min about -0.100 at Re = 100 vs Ghia et al. -0.1034).
- Roadmap: full Ghia benchmark comparison, grid-convergence study, 3D solver, CI with GitHub Actions.

## License

MIT. See `LICENSE`.
