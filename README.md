# syntheticTurbulenceInlet

A velocity inlet boundary condition for **OpenFOAM-13** that generates
synthetic turbulent fluctuations consistent with a target Reynolds-stress
tensor and integral-length-scale field, using a divergence-free synthetic
eddy method (DFSEM).

![Synthetic swirling turbulent inflow at the annulus inlet — velocity magnitude on the full geometry alongside a cross-section coloured by the z-component of velocity](figures/annulus_finer_combined_umag_clipuz_10s.gif)

The implementation derives from ESI's `turbulentDFSEMInlet` with three
types of improvement. Two equation-level corrections bring the boundary
condition into line with the original DFSEM paper (Poletto, Craft &
Revell, 2013), from which the original DFSEM code had deviated. The
global normalisation constant is re-derived from the eddy-shape
integrals, replacing the per-Γ scaling used in both the paper and the
original code. Applicability is broadened through geometry-aware
sampling, automatic profile mapping, runtime-selectable inlet targets,
swirl, periodic and wall-image ghost eddies, and reproducible RNG. The
combined effect is improved fidelity and a substantially reduced
configuration burden.

See `documentation/doc_syntheticTurbulence.pdf` for the full method description.

## Authors

Hursanay Fyhn and Eirik Holm Fyhn — SINTEF Energy Research
(<hursanay.fyhn@sintef.no>)

## Citation

If you are writing academic publications using `syntheticTurbulenceInlet`, please cite one or more of the following articles:

1. H. Fyhn, A. Gruber, H. T. Nygård, J. Dawson, and K.-J. Nogenmyr,
   "Numerical and Experimental Investigation of a Non-Premixed Flame
   of Partially Decomposed Ammonia at Atmospheric Pressure," *Journal
   of Engineering for Gas Turbines and Power* (accepted, 2026).

2. H. Fyhn, O. H. H. Meyer, Q. Wang, N. Worth, J. Koomen, T. Dammers,
   and A. Gruber, "Numerical and Experimental Investigation of a
   Laboratory-Scale Hydrogen-Fired *FlameSheet*™ Burner Operated at
   Atmospheric Pressure With Full Optical Access," *Journal of
   Engineering for Gas Turbines and Power* (accepted, 2026).

## Requirements

- OpenFOAM-13 (foundation version)

## Build

From the repository root, with your OpenFOAM-13 environment sourced:

```bash
wmake
```

This produces `libsyntheticTurbulenceInlet.so` in `$FOAM_USER_LIBBIN`.

## Use

Load the library in `system/controlDict`:

```
libs ("libsyntheticTurbulenceInlet.so");
```

Minimal patch entry in `0/U`:

```
inlet
{
    type            syntheticTurbulenceInlet;
    geometryMode    annulus;
    massFlowMode    true;
    massFlow        2.53e-4;
}
```

`geometryMode` can be `annulus`, `disk`, `rectangle`, `channel`, or
`arbitrary`. The full list of dictionary entries (with defaults) is in
the parameter table of `documentation/doc_syntheticTurbulence.pdf`.

## Tutorials

Ready-to-run cases under `tutorials/`:

| Tutorial | `geometryMode` | Inlet input | Solver |
|---|---|---|---|
| `annulus/` | `annulus` | `massFlow` + non-zero `swirlNr` | compressible (`fluid`) |
| `disk/` | `disk` | `Umean` | incompressible |
| `rectangle_incompress_U/` | `rectangle` | `Umean` | incompressible |
| `rectangle_compress_mdot/` | `rectangle` | `massFlow` | compressible (`fluid`) |
| `channel_poletto/` | `channel` | `Umean` | incompressible (validation case, Re_τ ≈ 395) |
| `channel_fluctuating_U/` | `channel` | `Umean` from a time-varying `Function1` table | incompressible |
| `arbitrary/` | `arbitrary` | `Umean` (half-disk patch) | incompressible |

Each case has `Allrun` and `Allclean` scripts.

## Repository layout

```
syntheticTurbulenceInletFvPatchVectorField.{H,C}   BC implementation
eddy/                                              Eddy class
URL/                                               Bundled DNS profiles (U, R, L)
Make/                                              Build configuration
tutorials/                                         Example cases
documentation/                                     Method documentation (LaTeX + PDF)
```
