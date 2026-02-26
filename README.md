# F-16 Pitch Control: Nominal MPC, Tube MPC, and UC-MPC (Uncertainty Compensation)

This repo contains MATLAB code for constrained control of an F-16 pitch attitude model under matched and unmatched uncertainties, including:

- **Nominal MPC** with constraint tightening (offline nominal plan)
- **Vanilla MPC** (baseline MPC without ancillary feedback / tightening)
- **Tube MPC** (nominal MPC on tightened sets + robust ancillary feedback)
- **UC-MPC** （nominal MPC + robust ancillary feedback + adaptive uncertainty compensation)

The model and example are based on the linearized F-16 pitch dynamics 

---

## Repo structure

### Core scripts (you run these)
- `nominal_MPC_unmatched.m`  
  Builds the **nominal MPC** problem with **tightened constraints** and simulates the nominal closed-loop system:
  - output: `umpchist_nominal.mat` (nominal MPC input trajectory)
  - output: `xhist_nominal.mat` (nominal MPC state trajectory)
  - output: `utotalhist_nominal.mat` (total input `u + Kx x` if used)

- `vanilla_MPC_unmatched.m`
  Runs **baseline MPC** (no ancillary feedback in dynamics) and simulates the uncertain system directly.  
  - output: `MPC_umpc.mat`, `MPC_x.mat`, `MPC_u.mat`

- `tube_MPC_unmatched.m` (Tube MPC script)  
  Constructs:
  - disturbance set `W`
  - disturbance invariant tube set `Z` (finite approximation)
  - tightened constraints `X_bar = X ⊖ Z`, `U_bar = U ⊖ KZ`  
  Then runs nominal MPC on tightened sets and simulates the true system using:
  \[
    u = u_{\text{nom}} + K(x_{\text{real}} - x_{\text{nom}})
  \]
  - output: `TMPC_umpc.mat`, `TMPC_x.mat`, `TMPC_u.mat`

- `main_L1MPC_f16.m` (your long main script)  
  Runs **L1 adaptive augmentation** on top of a precomputed nominal MPC trajectory.
  **Important:** this script expects `umpchist_nominal.mat` and `xhist_nominal.mat` to exist.
  - output: `L1MPC_f.mat`, `L1MPC_w.mat`, `L1MPC_sigma1.mat`, `L1MPC_u.mat`, `L1MPC_ul1.mat`, `L1MPC_x.mat`, `L1MPC_xnom.mat`

---

## Dependency files / functions

These functions are called by the main scripts:
- `sim_L1MPC.m`  
  Fast-time simulation that applies baseline + L1 adaptive control and (optionally) constraint tightening logic.

- `compute_L1_bounds.m`  
  Computes/updates the L1 performance bounds used for constraint tightening, including:
  - `rho_til_x`: bound on `||x - x_nom||`
  - `rho_u`: bound on `||u_L1||`
  - derived `rho_til_u = |Kx| rho_til_x + rho_u`

---

## Requirements

- MATLAB (tested with Control System Toolbox)
- YALMIP installed and on MATLAB path
- MOSEK installed (solver used by YALMIP)

In scripts, paths are currently added via:
```matlab
addpath(genpath('.../YALMIP'))
addpath(genpath('.../Mosek/...'))
