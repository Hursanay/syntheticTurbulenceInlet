# syntheticTurbulenceInlet

A velocity inlet boundary condition for **OpenFOAM-12** that generates
synthetic turbulent fluctuations consistent with a target Reynolds-stress
tensor and integral-length-scale field, using a divergence-free synthetic
eddy method (DFSEM).

The implementation derives from ESI's `turbulentDFSEMInlet` with three
classes of improvement: equation-level corrections that bring the BC into
line with the original DFSEM paper (Poletto, Craft & Revell, 2013),
re-derived normalisation constants that improve the match to DNS
reference data (especially near-wall), and a substantially reduced
configuration burden through geometry-aware sampling, automatic profile
mapping, swirl, ghost eddies, and reproducible RNG.

See `documentation/doc_syntheticTurbulence.pdf` for the full method description.

## Authors

Hursanay Fyhn and Eirik Holm Fyhn — SINTEF Energy Research
(<hursanay.fyhn@sintef.no>)

## Citation

If you are writing academic publications using `syntheticTurbulenceInlet`, please cite one or more of the following articles:

1. H. Fyhn, A. Gruber, H. T. Nygård, J. Dawson, and K.-J. Nogenmyr,
   "Numerical and Experimental Investigation of a Non-Premixed Flame
   of Partially Decomposed Ammonia at Atmospheric Pressure," *Journal
   of Engineering for Gas Turbines and Power* (submitted, 2026).

2. H. Fyhn, O. H. H. Meyer, Q. Wang, N. Worth, J. Koomen, T. Dammers,
   and A. Gruber, "Numerical and Experimental Investigation of a
   Laboratory-Scale Hydrogen-Fired *FlameSheet*™ Burner Operated at
   Atmospheric Pressure With Full Optical Access," *Journal of
   Engineering for Gas Turbines and Power* (submitted, 2026).

## Requirements

- OpenFOAM-12 (foundation version)

## Build

From the repository root, with your OpenFOAM-12 environment sourced:

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
| `annulus/` | `annulus` | `Umean` + non-zero `swirlNr` | incompressible |
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
