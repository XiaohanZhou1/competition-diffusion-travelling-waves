# Competition–Diffusion Travelling Waves

MATLAB code for the numerical and asymptotic study of travelling waves in a
two-species competition–diffusion model.

The full model is

\[
u_t=u_{xx}+u(1-u)(u-a)-\frac{uv}{\delta},
\qquad
v_t=Dv_{xx}+v(1-v)-\gamma\frac{uv}{\delta}.
\]

The repository contains full PDE simulations, travelling-wave visualisations,
Stefan-type and flux-type reduced models, asymptotic calculations, parameter
sweeps, and comparisons between the full and reduced descriptions.

## Repository structure

```text
code/
├── full_model/            Full two-species PDE solvers
├── travelling_waves/     Profiles, snapshots, and composite visualisations
├── reduced_models/
│   ├── stefan/            Small-diffusion Stefan-type approximation
│   └── flux/              Large-diffusion flux-type approximation
├── comparisons/           Full-versus-reduced model comparisons
└── parameter_studies/     Parameter sweeps, contours, and fitted scalings
```

## Requirements

- MATLAB
- `ode15s`
- `pdepe`
- `bvp4c` for scripts containing boundary-value problems

## Quick start

Clone or download the repository, start MATLAB in the repository root, and run:

```matlab
addpath(genpath('code'));
plot_competition_TW_profiles
```

Other useful entry points include:

- `plot_competition_TW_snapshot` — travelling-wave snapshot.
- `plot_contour_surface` — wave-speed parameter sweep.
- `compare_reduced_models_Lt` — interface-location comparison.
- `compare_c_gamma` — wave-speed comparison against reduced models.
- `plot_stefan_smallc_asymptotics` — small-speed Stefan asymptotics.
- `plot_flux_matched_comparison` — matched flux-model asymptotics.

## Notes on reproducibility

Most scripts define their numerical parameters near the beginning of the file.
Parameter sweeps can be computationally expensive. Generated data and figures
should be saved outside `code/` so that source files remain separate from
outputs.

## Citation

If you use this code in academic work, please cite the associated publication
when it becomes available.
