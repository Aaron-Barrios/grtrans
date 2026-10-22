# grtrans for APERTURE: how one ray is traced

This guide follows one ray of the Cylco-Synchrotron pipeline through grtrans,
from the camera pixel to the number stored in the image. It is written from
the code (2026-10-15, the build of `history.tex` Stage 39, updated for the
grtrans fixes later that day) and cites it as
`file:line`; paths are in `grtrans/` unless they start with `pipeline/` or
`verification_scripts/`. Line numbers drift: search for the routine name when
one does not match.

The physics behind each step is in `docs/RT_background.tex`; the notebook's
use of it, and the measurements, are in `docs/RT_notebook_overview.tex` (P4).
`jnu_calculator/IMPLEMENTATION.md` is the companion guide for the exact
emissivity code that fills the lookup table. Open problems are in
`docs/muchado.md`, cited below by their ids as of the date given (e.g. L8-F7); fixed items are in `docs/history.tex`.

grtrans (Dexter 2016) is the project's git submodule `grtrans/`: branch
`aperture` of the fork `Aaron-Barrios/grtrans`, which starts from commit
`c76cb11` and holds our changes; the project records which commit it uses.
Build: `make grtrans libgrtrans pgrtrans` in `grtrans/`, then commit there
and in the project (`README.md`, "How grtrans is used").

## Contents

1. [The ray at a glance](#1-the-ray-at-a-glance)
2. [Inputs: what the pipeline writes](#2-inputs-what-the-pipeline-writes)
3. [Camera and pixels](#3-camera-and-pixels)
4. [Where the samples go along a ray](#4-where-the-samples-go-along-a-ray)
5. [The fluid binary](#5-the-fluid-binary)
6. [From a sample to the local plasma](#6-from-a-sample-to-the-local-plasma)
7. [The comoving frame](#7-the-comoving-frame)
8. [Emission and absorption](#8-emission-and-absorption)
9. [Radiative transfer](#9-radiative-transfer)
10. [Adaptive refinement](#10-adaptive-refinement)
11. [Outputs and the extra quantities](#11-outputs-and-the-extra-quantities)
12. [Failures](#12-failures)
13. [Threads](#13-threads)
14. [What d and R_high do](#14-what-d-and-r_high-do)
15. [Our changes to stock grtrans](#15-our-changes-to-stock-grtrans)
16. [Known problems](#16-known-problems)

## 1. The ray at a glance

| Step | What happens | Where |
|---|---|---|
| Start | read `files.in`, then the `inputs.in` namelists | `grtrans_program.f90:7-13`, `read_inputs.f90:59-87` |
| Cameras | frequencies, inclination, one camera per (frequency, mass-rate, inclination, time) | `pgrtrans.f90:87-119`, `:144-146` |
| Fluid | read the APERTURE binary once | `pgrtrans.f90:174` → `fluid.f90:203-204` → `fluid_model_aperture.f90:138-232` |
| Pixels | per-pixel constants of motion, all rays at once | `pgrtrans.f90:187-190` → `geodesics.f90:64-111` → `geokerr_wrapper.f:1` |
| Ray loop | OpenMP over rays | `pgrtrans.f90:215-228` → `grtrans_driver` |
| Geodesic | sample points of this ray | `grtrans_driver.f90:106` → `geodesics.f90` → `GEOKERR` (`geokerr_wrapper.f:357-604`) |
| Fluid at samples | interpolate the binary at every sample | `grtrans_driver.f90:149` → `fluid_model_aperture.f90:331-474` |
| Comoving frame | redshift $g$ and field angle $\theta_B$ per sample | `grtrans_driver.f90:151-153` → `kerr.f90:586-730` |
| To cgs | $n_e$, $B$, $\Theta_e$ | `grtrans_driver.f90:188` → `fluid.f90:999-1036` |
| Per frequency | comoving $\nu$, emissivity and absorption | `grtrans_driver.f90:200-201` → `emis.f90:550-551` → `jnu_table.f90:370-427` |
| Invariants | $j\,g^2$, $\alpha/g$, normalization | `grtrans_driver.f90:259-262`, `emis.f90:859-866` |
| Transfer | `explin` along the ray | `grtrans_driver.f90:268` → `radtrans_integrate.f90:890-925` |
| Refinement | optional, off by default | `grtrans_driver.f90:276-278`, `:590-903` |
| Extras | emission-weighted averages | `grtrans_driver.f90:280-405` |
| Store | pixel column of the camera | `grtrans_driver.f90:443-459`; written by `camera.f90:322-341` |

One geodesic and one fluid lookup serve every frequency of a camera
(`grtrans_driver.f90:132-201`), so several frequencies in one run cost little
more than one (muchado §15 item 3). One geodesic also serves several fluid
binaries traced in one run, one camera each (§5, "Several binaries in one
run").

## 2. Inputs: what the pipeline writes

`pipeline/raytrace_q6.py` (`configure_grtrans`, `:318-325`, which refuses a
camera spin other than the run's and writes the text `render_inputs`,
`:328-354`, produces) rewrites `grtrans/inputs.in` before every trace;
`verification_scripts/grtrans_sandbox.py` (`Sandbox`, `:167-262`) writes its
own copy in a scratch directory (with its own substitutions; it could reuse
`render_inputs`). Both set:

| Namelist value | Meaning | Set from |
|---|---|---|
| `fmin`, `fmax` (`&emisdata`), `nfreq=1` | observing frequency [Hz] | `pc.FREQ` |
| `nn=npix,npix,nup` | pixels and samples per ray | `pc.CAMERA["npix"]`, `pc.NUP` |
| `i1=1`, `i2=npix²` | rays traced (both written by both since 2026-10-15, `raytrace_q6._subs_nn`, `:187-191`; before, `configure_grtrans` wrote only `i2`; review item Q4, fixed) | — |
| `mumin=mumax=cos i`, `nmu=1` | inclination as $\cos i$ | `pc.CAMERA["inc_deg"]` |
| `spin` | $a$ | `run_config.spin()` |
| `gridvals=-f,f,-f,f` | field of view, impact parameters [$r_g$] | `pc.CAMERA["fov_rg"]` |
| `uout` | $1/r_{\rm out}$ of the run's grid: rays start at its outer edge | `run_config.grid().r_out` |
| `mbh` (`&emisdata`) | black-hole mass [$M_\odot$], sets the length scale `LBH` | `pc.M_BH_MSUN` |
| `extra=1` | 21 extra slots (§11) | — |
| `iname` | integrator (`explin`) | sandbox only; `inputs.in` keeps `explin` |
| `refine_on` (and `refine_*`) | adaptive refinement (§10) | `raytrace_q6.REFINE` (in `pipeline/raytrace_q6.py`, not `pipeline_config.py`) |
| `fdfile` (`&harm`) | the fluid binary; read by the APERTURE model (`fluid.f90:148`, `fluid_model_aperture.f90:198-213`); at most 100 characters (`fluid.f90:62`), so the sandbox links it under a short name. A name ending in `.lst` is a list of binaries traced in one run (§5, "Several binaries in one run") | — |

Fixed in `inputs.in` and used: `standard=1` (geokerr camera at "infinity"),
`nrotype=2` (rectangular grid), `fname='APERTURE'`, `ename='JNUTABLE'`,
`nvals=1` (Stokes $I$ only), `use_geokerr=T`, `cflag=0` (binary output),
`nt=1`, `nload=1`, `nmdot=1` (the sandbox sets `nmdot` to the number of
binaries when it traces several in one run, §5).

Inert for APERTURE: `phi0` (the model reads only $r$, $\theta$;
`fluid_model_aperture.f90:358-359`), `uin` (never passed to GEOKERR;
`geodesics.f90:226`), `rcut` (circular grids only), `dt`, `mdotmin/max`, and
`sigcut` in `&fluiddata` (`convert_fluidvars_aperture` never reads it, as its
comment and the one in `inputs.in` say; the APERTURE cut comes from the
binary, §5). The other `&harm` and
`&analytic` values are not read by the APERTURE branch.

## 3. Camera and pixels

- **Grid.** Pixel $i$ (1-based) has $I=\lfloor(i-1)/n_\phi\rfloor+1$,
  $J=((i-1)\bmod n_\phi)+1$, and impact parameters at pixel centers,
  $\alpha=a_1+(a_2-a_1)(I-\tfrac12)/n_{\rm ro}$,
  $\beta=b_1+(b_2-b_1)(J-\tfrac12)/n_\phi$ (`geokerr_wrapper.f:258-262`);
  $\beta$ runs fastest. With `gridvals=±f` the pitch is $2f/n_{\rm pix}$
  and the outermost center sits at $f-$ half a pitch.
- **Units.** $\alpha$, $\beta$ are impact parameters at infinity in $r_g$
  (`geokerr_wrapper.f:373-376`).
- **Constants of motion.** $L=-\alpha\sqrt{1-\mu_0^2}$,
  $Q^2=\beta^2-(a^2-\alpha^2)\mu_0^2$, with $\mu_0=\cos i$
  (`geokerr_wrapper.f:275-285`).
- **"Infinity."** The observer sits at $u_0=\min(10^{-4},
  1/(100\,{\rm ABMAX}))$ (`geokerr_wrapper.f:23,155,163`), $\approx6\times10^{-9}$
  for a $\pm30\,r_g$ field ($r\approx1.6\times10^8\,r_g$); time and angle
  integrals start there, but samples are kept only inside `uout` (§4).

## 4. Where the samples go along a ray

GEOKERR places the samples evenly in $u=1/r$, at the midpoints of `nup`
equal steps (`OFFSET = 0.5`, `geokerr_wrapper.f:45,498,519`), so no sample
lies exactly at `uout`, at a turning point or at the horizon.

- **Rays that fall in** ($u$ reaches $u_+$, the horizon): $\Delta u=(u_+-u_{\rm
  out})/n_{\rm up}$, $n_{\rm up}$ samples from `uout` to just outside the
  horizon (`:468-489`). For $a=0.998$ and the run's `uout` $=1.168\times10^{-2}$,
  $\Delta u=0.929/n_{\rm up}$: at $n_{\rm up}=3200$, $\Delta r\approx2.9\times10^{-4}\,r^2$,
  i.e. $0.03\,r_g$ at $r=10$ and $0.7\,r_g$ at $r=50$.
- **Rays with a radial turning point** ($u_b$): $\Delta u=2(u_b-u_{\rm
  out})/n_{\rm up}$, about half the samples on the way in and half on the way
  out (`:481-484`, `:575-576`). Near the turning point equal steps in $u$ are
  large steps in path and in $\mu=\cos\theta$, so `mufill` replaces the last
  and first `kext` $=\lceil3n_{\rm up}/400\rceil$ points around it by `next`
  $=\lceil60n_{\rm up}/400\rceil$ points evenly in $\mu$-path
  (`geodesics.f90:77-92`, `geokerr_wrapper.f:538-571`). Such a ray carries
  $n_{\rm up}+{\tt next}-2\,{\tt kext}$ samples: 226 at 200, 454 at 400,
  3632 at 3200 (`:571,593`).
- **Rays that never get inside `uout`** carry no samples (`:493,597`;
  `geodesics.f90:294`): the pixel is 0.
- **Per sample** grtrans keeps $r$, $\theta$, $\phi$, $t$ (`t0=0` with
  `nload=1`), the affine parameter $\lambda$ (decreasing along the array: index
  1 is at `uout`, the last sample has $\lambda=0$) and the photon momentum
  $k^\mu$ (`geodesics.f90:333-364`).

Because the samples are fixed in $u$, a ray at $d=1$ crosses 1350–4200 grid
cells inside the field of view, and 3200 samples miss 9–36% of them
(`history.tex`, Stage 35, `ray_spacing_cells.py`). That is why the image
converges in samples only where the light is not concentrated in single
cells (§14).

## 5. The fluid binary

Written by `pipeline/aperture_to_binary.py` from a reconstructed fluid file
(`fluid.NNNNN.h5`, or `fluid.NNNNN.c{d}.h5` for the $d\times d$
coarse-grained fluid of `pipeline/coarsen_fluid.py`); read once by
`read_aperture_data` (`fluid_model_aperture.f90:138-232`).

- **Header** (60 bytes, `<iiidddddd`, `aperture_to_binary.py:105-126`):
  magic, $n_r$, $n_\theta$, spin, $\lambda_u$, mass ratio $N$, and the
  emissivity dispatch thresholds `THETA_HYBRID`, `NTHETA_MAX`, `THETA_SYNCH`
  (0 for none). grtrans stops on a magic other than `APL7`/`APL8` and on
  $\lambda_u\le0$ or $N\le0$ (`fluid_model_aperture.f90:149-172`), and on a
  spin that differs from `inputs.in` by more than $10^{-4}$ (`:466-472`; it
  only warned before 2026-10-15; the pipeline also refuses it before tracing,
  `raytrace_q6.check_camera_spin`).
- **Arrays** (float64, $\theta$ slowest, index $k=i_\theta n_r+i_r$): the axes
  $r[n_r]$, $\theta[n_\theta]$, then 14 fields: $n_e$, $\Theta_e$, $|b|$,
  $u^0..u^3$, $b^0..b^3$, valid, $\beta_p$, $\sigma$
  (`aperture_to_binary.py:272-325`).
- **Axes are cell centers** from `config.toml`:
  $r_i=\exp(\ln r_{\rm lo}+(i+\tfrac12)\Delta\ln r)$,
  $\theta_j=\theta_{\rm lo}+(j+\tfrac12)\Delta\theta$ (`:175-191`).
- **Fields:** $n_e=n_e^{\rm rest}$ and $\Theta_e=P_e/n_e^{\rm rest}$ of the
  electrons in their own frame; $u^\mu=u_e^\mu$, $b^\mu=b_e^\mu$ (Kerr–Schild,
  contravariant); $|b|=\sqrt{b_e^2}$ (not used by grtrans, which recomputes
  it, §6, and skipped when grtrans reads the file, as is the validity
  array); $\beta_p=2P/b^2$ and $\sigma=b^2/(2(n_e+n_p))$ of the pair fluid,
  0 where invalid or non-finite (`:158-172, 306-307`). Electron-free and
  non-finite cells get $n_e=\Theta_e=0$; a negative $\Theta_e$ is floored to 0
  (`:284-299`).
- **Relabeled reading** (`pipeline/relabel.py:70-80`): $\Theta_e\to
  N\Theta_e/R(\beta_p)$ with $R=(R_{\rm high}\beta_p^2+R_{\rm
  low})/(1+\beta_p^2)$; this, and $N$ in the header (§6), are the only places
  where $R_{\rm high}$ and the relabeling enter grtrans.
- **The σ cut.** `APL7` with a cut ("cell" mode): cells with
  $\sigma>\sigma_{\rm cut}$ and $n_e>0$ are written with $n_e=\Theta_e=0$,
  keeping $u$, $b$, $\beta_p$, $\sigma$ (`:311-315`). `APL8` ("point" mode):
  the cells are written uncut and the cut value is appended after the arrays
  (`:317-325`); grtrans then tests it at every sample (§6). No cut: `APL7`.
  Why per sample: §6, "The σ cut: why per sample". One grtrans build reads both: the writer's `--sigma-cut-mode` option (or
  `pc.TRACE_SIGMA_CUT_MODE`) decides which of the two identifiers, `APL7` or
  `APL8`, it puts in the binary's first four bytes, and grtrans reads that
  identifier to know whether to test the cut itself.
- **Geometric masks** (`pipeline_config.POLAR_CUT_DEG`, `R_EMIT_MAX`; writer
  options `--polar-cut-deg`, `--r-emit-max`): cells within the polar cone or
  beyond the radius are written with $n_e=\Theta_e=0$, so they neither emit
  nor absorb (cell-level, in either mode). Defaults since 2026-10-16: no cone
  (with the σ cut at 0.5 no cone of 0.25–4° removed any light) and
  $r_{\max}=80\,r_g$, which leaves out the run's damping layer ($r>79.9$,
  numerical gas; up to 2.4% of a $d=8$ image, more at $d=2$). Store tag
  `_pcX`, `_rmX`. Tests: `docs/RT_notebook_overview.tex`, "The emission cuts
  of GRMHD imaging"; `history.tex`, Stage 59.

### Several binaries in one run

Added 2026-10-15 (`planned_run.md`, step 3b). A pixel's geodesic depends only
on the camera, `nn` and the spin, so one run can trace many binaries along the
same paths, through grtrans's existing loop over model parameters (one
parameter per `mdot`, which APERTURE does not use):

- **Input.** `fdfile` ending in `.lst` names a text file with one binary path
  per line (up to 1023 characters; blank lines and `#` lines skipped;
  `read_binary_list`, `fluid_model_aperture.f90:236-271`), and `nmdot` must
  equal their number (`pgrtrans.f90:178-182` stops otherwise). Any other
  `fdfile` is one binary, as before.
- **Storage.** `read_aperture_data` (`:198-232`) reads every binary into its
  own `aperture_bin` (`:65-95`; the array `abins`, `:103`), with its own grid,
  $\lambda_u$, mass ratio and σ cut (cell or per sample). Every binary's spin
  must match the camera's (`:466-472`), and all must carry the same emissivity
  dispatch (`THETA_HYBRID`, `NTHETA_MAX`, `THETA_SYNCH`; `:219-226`), since that
  is one setting of `jnu_table`. 12 of the 16 arrays are kept (the |b| and
  validity arrays are skipped), $96$ bytes per cell: a $d=1$ binary takes
  about 384 MB of memory (measured), 40 of them about 15 GB.
- **Selection.** `aperture_active` (`:108-109`), `threadprivate`, is the binary
  `aperture_vals` (`:349-356`) and `convert_fluidvars_aperture` (through
  `aperture_lambda_u()`, `aperture_mass_ratio()`, `:286-292`;
  `fluid.f90:1032-1034`) use on the calling thread; `select_aperture_binary`
  (`:275-282`) sets it.
- **Driver.** For model parameter $m$, `grtrans_driver` selects binary $m$
  and redoes the fluid at the samples and `comoving_ortho` along the same
  geodesic (`grtrans_driver.f90:146-147, 164-176`); the rest of the loop
  body (conversion, emissivity, transfer, refinement, extras) is unchanged,
  so it is exactly a single-binary trace of binary $m$. With one binary
  nothing changes (`multibin` is false).
- **Output.** One block per binary in `grtrans.out`, in list order (§11).
- **Check** (`verification_scripts/multi_binary_check.py`,
  `output/multi_binary_check_20261015.txt`): single-binary images are bit for
  bit those of the build before; each block of a 9-binary run (dumps 250, 500,
  1000 × $R_{\rm high}$ 1, 100 at $d=4$; $d=2$; a per-cell binary; the $d=1$
  pair-plasma binary) is bit for bit the separate trace, with extras and with
  refinement, at 1 and 24 threads. Time at 3200 samples: 9 binaries at
  $128^2$ 92 s instead of 163 s (1.8×); 6 binaries at $256^2$ 169 s instead
  of 291 s (1.7×): each extra binary costs about half of a separate trace
  (24 s against 48 s at $256^2$), not the 40% the profile suggested. On 28
  threads with the build of 2026-10-17 a further binary costs 0.8–0.9 of a
  trace (0.46 on one thread: the per-binary local frame allocates
  temporaries that contend between threads; `history.tex`, Stage 62).
  Since 2026-10-19 every grtrans run gets glibc malloc tunables
  (`pipeline_config.GRTRANS_MALLOC_TUNABLES`: a 64 MiB trim threshold and
  top pad, through `pc.grtrans_env()`/`raytrace_q6.run_checked`): glibc no
  longer shrinks and regrows the heap around those temporaries, and a
  further binary costs ~0.3 of a trace on 28 threads, single traces ~1.4×
  faster, output byte-identical (`history.tex`, Stage 71). A
  `GLIBC_TUNABLES` already in the environment is kept.
- **Production.** `run_q6` traces a movie's frames in batches this way since
  2026-10-17 (`raytrace_q6.raytrace_batch`; `BATCH` = 4, `--batch`; with the
  tunables about 1.6–2.2× faster per frame than one frame per run, images
  identical).

## 6. From a sample to the local plasma

`aperture_vals` (`fluid_model_aperture.f90:331-474`) runs for all samples of
the ray:

1. **Defaults:** every sample starts empty ($n_e=\Theta_e=|b|=0$, $u=(1,0,0,0)$,
   $b=0$, $\beta_p=\sigma=0$; `:361-368`).
2. **The grid's edge:** samples beyond the grid, half a cell (half a step in
   $\ln r$ and $\theta$) past the outermost cell centers, stay empty; between
   the outermost centers and that edge (at the axis and at the inner and outer
   radii) the position is clamped to the box of centers for the
   interpolation, so the edge cell's values hold, while the metric uses the
   true position (`aperture_vals`; history.tex, Stage 93).
3. **Cell and weights:** bisection on each axis (`:294-309`), weights linear in
   $\ln r$ and $\theta$ (`:378-388`): bilinear interpolation between the four
   surrounding cell centers.
4. **Interpolated:** $n_e$ (`:391`); $\beta_p$ and $\sigma$ over the corners
   that hold fluid ($\sigma>0$), renormalized, as the cut of step 5 tests
   $\sigma$ (`:392-412`; since 2026-10-19, before which empty corners counted
   as 0 and pulled both low at fluid edges, by up to ~20% flux-weighted at
   $d=1$; `history.tex`, Stage 69), kept even where the sample ends up empty,
   for the extras; if $n_e\le0$ the sample stays empty (`:413`).
5. **Per-sample cut** (`APL8` only): $\sigma$ is re-averaged over the corners
   with $\sigma>0$ only (cells with no pair particles, invalid cells and
   field-free cells are left out); if that average exceeds the cut the whole
   sample stays empty, not only its density (`:414-433`).
6. **Temperature:** $\Theta_e=\max(\sum w\,n_e\Theta_e,0)/n_e$, a
   density-weighted interpolation (`:434-437`).
7. **Velocity and field:** $u^\mu$, $b^\mu$ component by component, Kerr–Schild
   (`:439-447`). There is no validity gate: empty and cut cells contribute 0
   to $n_e$ and $n_e\Theta_e$ but keep their $u$ and $b$ (`:323-328`).
8. **To Boyer–Lindquist** (`uks2ubl`, `kerr.f90:131-144`; `:284-287`), then
   normalized: $u\to u/\sqrt{-u\cdot u}$, $b\to b+(b\cdot u)u$ (so $b\cdot u=0$),
   $|b|=\sqrt{b\cdot b}$ (`:461-473`). A sample whose $u$ is not timelike is
   emptied.

Then `convert_fluidvars_aperture` (`fluid.f90:999-1033`) makes the cgs values:
$n_e[{\rm cm^{-3}}]=\lambda_u n_e$, $B[{\rm G}]=|b|\sqrt{4\pi\lambda_u N
m_ec^2}$ (the $\sqrt N$ of the relabeled reading), and `tcgs` $=\Theta_e$
(dimensionless, **not** kelvin). $\beta_p$ and $\sigma$ go on as `f%rho2` and
`f%Be` for the extras (§11).

### The σ cut: why per sample

GRMHD images leave out gas with $\sigma=b^2/\rho>1$ (Heaviside–Lorentz;
our $\sigma=b^2/(2\rho_0)>0.5$), whose thermodynamics is set by numerical
floors. How the main codes apply it (source read 2026-10-15; `history.tex`,
Stage 43):

| Code | Tested | σ at the point | Switched off |
|---|---|---|---|
| ipole (EHT libraries) | each ray sample (default `USE_GEODESIC_SIGMACUT`) | the per-zone σ, interpolated trilinearly | $n_e$ only, so $j=\alpha=0$ |
| grtrans (`andrew_sigcut`) | each sample | $B^2/(4\pi\rho c^2)$ from the point's values (for its iharm/harmpi models those are one zone's, not interpolated) | $\rho$, $n$; sets $T=10^9$ K, keeps $B$ and non-thermal electrons |
| RAPTOR | each sample | $b^2/\rho$ from interpolated primitives | the whole step; also $\Theta_e>100$ and $r>50$ |
| BHOSS | (closed source) "no emission where $b^2/\rho c^2\ge1$" | — | — |
| igrmonty, ebhlight | per zone | the zone's own | $\Theta_e$ set tiny |

The value of the cut is debated (Chael et al. 2019; Fromm et al. 2022;
Cruz-Osorio et al. 2022; Mizuno et al. 2021 find little change between 1 and
5 at 230 GHz), its form hardly: the EHT papers and the code comparison of
Prather et al. (2023) fix $\sigma_{\rm cut}=1$ with emission and absorption
zeroed but leave the interpolation code-specific. Only Wong et al. (2022,
PATOKA) compare forms: ipole's former per-zone cut drew sharp ridges along
the mask edge that change the long-baseline signal, and testing an
interpolated σ at each step removes most of them, which is why ipole cuts
per sample.

The pipeline therefore uses ipole's form (`TRACE_SIGMA_CUT_MODE = "point"`,
the default since 2026-10-15): the per-cell σ is interpolated to each sample
(over the corners that hold particles: a particle-free cell is sparse
sampling, not vacuum) and the sample is left empty above the cut. It is the
convention of most EHT model images; it needs nothing the binary lacks;
forming $b^2/\rho$ from interpolated fields instead (grtrans, RAPTOR) lets
the dense side of an edge dominate and treats our particle-free cells
almost as vacuum (it would mask the samples beside 11–20% of the emitting
cells at $d=1$; `sigma_cut_empty.py`). The per-cell cut (`"cell"`, the
stores traced before 2026-10-15) is ipole's abandoned form. At
$R_{\rm high}\le5$ the choice moves the 230 GHz flux by at most 4% and the
cut's value (0.25–2.5) by about 10%; at $R_{\rm high}=100$, where the light
sits at the mask edge, by factors (§14).

## 7. The comoving frame

`comoving_ortho` (`grtrans_driver.f90:151-153` → `kerr.f90:550-562, 586-730`)
builds the fluid tetrad ($e_{(t)}=u$, after Kulkarni et al. 2011; `:645-658`)
and projects the photon momentum on it:

- **Redshift** $g=\nu_{\rm obs}/\nu_{\rm emit}=1/k_{(t)}$ (`:722`); the
  comoving frequency is $\nu'=\nu/g$ (`grtrans_driver.f90:200`).
- **Field angle:** $\cos\theta_B=\hat b\cdot\hat k$ in the fluid frame (`kerr.f90:694,
  713`), set to 0.5 where $b\cdot b\le0$ (`kerr.f90:714-718`).
- **Clamp:** $\cos\theta_B$ is limited to $[-0.9999,0.9999]$ (`kerr.f90:618,
  719-720`), so every sample within $0.81^\circ$ of the field is evaluated at
  $0.81^\circ$, below the table's first $\theta_B$ column ($1.15^\circ$),
  where the table clamps anyway. Stock grtrans used $\pm0.99$ ($8.1^\circ$)
  here; lowered 2026-10-16 to match the other two angle routines of
  `kerr.f90`. It moved the relabeled 230 GHz flux by at most 0.04% and no
  $d$-ratio (`verification_scripts/clamp_check.py`; `history.tex`,
  Stage 54).

## 8. Emission and absorption

`ename='JNUTABLE'` selects the APERTURE emissivity (`emis.f90:364-370, 452`):
`calc_emissivity` calls `hybrid_emissivity(nu', n, B, theta_B, Theta_e)`
(`emis.f90:550-551` → `jnu_table.f90:370-427`), which splits the samples into
three branches, and leaves the rest empty (`:378-397`):

| Branch | Where | What |
|---|---|---|
| Cyclotron | $\Theta_e<{\tt THETA\_HYBRID}=5\times10^{-3}$ and $(\nu/\nu_B)\Theta_e\le{\tt NTHETA\_MAX}=0.1$, $n_e>0$, $\Theta_e>0$ | `cyclotron_emis.f90`: a sum of Gaussian lines centered at $n\nu_B[1-(n+\tfrac32)\Theta_e]$ (the relativistic mass shift) of width $n\nu_B\sqrt{\Theta_e\max(\lvert\cos\theta_B\rvert,10^{-4})^2+(n+\tfrac32)\Theta_e^2}$ (Doppler and mass shift; `history.tex`, Stage 81), harmonics up to $\min(\nu/\nu_B+20,\,2/(e^2\Theta_e))$ (`:75-110`) |
| Synchrotron fit | above the table's largest $\Theta_e$ ($4.6\times10^4$; `THETA_SYNCH` unset) | Leung, Gammie & Noble (2011), with the exact $K_2(1/\Theta)$ (`besselk2`, a trapezoid sum of its integral, $\sim10^{-14}$ relative; before 2026-10-15 $2\Theta^2-\tfrac12$, 8% off at $\Theta=1$ and 0 at 0.5) and $\max(\lvert\sin\theta_B\rvert,10^{-6})$ (`jnu_table.f90:316-368`); a `THETA_SYNCH` below 1 is warned about, since the fit is ultrarelativistic |
| Table | everything else with $n_e>0$ and $\Theta_e>0$ | interpolation of $\log_{10}F$: tensor-product Catmull-Rom cubic over the $4\times4\times4$ neighboring entries for $\Theta_e\ge$ `th_cubic` $=1$ (linear in an axis's edge interval; the trilinear value where a neighbor is at the $-300$ floor), trilinear below; the linear lookup was low by 2.0–2.4% at 230 GHz (`history.tex`, Stage 79). The $x$ axis is piecewise log-uniform (443 entries, 144 per decade between 0.1 and 30, from $10^{-6}$; Stage 80): its interval is found by bisection, and the cubic in $x$ is used only inside a uniform run |
| None | $n_e\le0$ or $\Theta_e\le0$ | $j=\alpha=0$ (`:383-397`) |

The thresholds come from the binary header (`set_hybrid_dispatch`,
`jnu_table.f90:299-315`, called at `fluid_model_aperture.f90:231`), written
from `pc.THETA_HYBRID`, `pc.NTHETA_MAX`, `pc.THETA_SYNCH`.

**The table** (`jnu_table.bin`, built by `pipeline/make_jnu_table.py` from
`jnu_calculator/jnu`): header `JNUT`, $n_x$, $n_{\theta_B}$, $n_\Theta$, the
three axes, then $\log_{10}F(x,\theta_B,\Theta_e)$ with $x$ slowest, where $F$ is
$j_\nu$ at $B=1$ G, $n_e=1\,{\rm cm^{-3}}$ and $x=\nu/\nu_B$; the sample gets
$j_\nu=n_eB\,10^{\log F}$ (`jnu_table.f90:243`). The production table is
$36\times13\times44$: $x\in[10^{-3},10^5]$, $\theta_B\in[0.02,\pi/2]$ linear,
$\Theta_e\in[7.17\times10^{-5},4.6\times10^4]$ (uniform in $\log$). Interpolation
is linear in $(\ln x,\theta_B,\ln\Theta_e)$ and assumes uniformly spaced
axes (`:207-241`). Clamps: $\theta_B$ folded into $[0,\pi/2]$; $x<10^{-3}$
uses the first column; above $x=10^5$ the edge value is scaled by the ratio
of the synchrotron shape (`:246-250`). 4095 of the 20592 entries are the
$-300$ floor (no emission below the fundamental at low $\Theta_e$); trilinear
interpolation toward them underestimates (muchado L8-F7). `JNU_TABLE_PATH`
overrides the file (`:107-108`); the sandbox clears it.

**Absorption** is Kirchhoff's, $\alpha_\nu=j_\nu/B_\nu(T_e)$ with
$kT_e=\Theta_em_ec^2$: full Planck in the table and fit branches
(`jnu_table.f90:252-259, 274-284, 343-344`), Rayleigh–Jeans in the cyclotron
branch (`cyclotron_emis.f90:113-120`). A sample with $n_e>0$ and
$\Theta_e=0$ (the writer's floor) goes to no branch and gets $j=\alpha=0$
(`jnu_table.f90:383-397`). Before 2026-10-15 it went to the table, which
clamped it to its smallest $\Theta_e$ and gave $j>0$ with $\alpha=0$
($\log F\le-60$: the relabeled 230 GHz and pair-plasma 1.6 MHz images of
dump 500, $d=4$, are bit-identical either way). A NaN $n_e$ or $\Theta_e$
still reaches the table, so the NaN stays visible (§12).

Only Stokes $I$ is carried (`nvals=1`, `e%neq=1`; `emis.f90:366-370`); with
`nvals=4` grtrans stops with a message (`grtrans_driver.f90:155-160`).

## 9. Radiative transfer

- **Invariants.** `invariant_emis` multiplies $j$ by $g^2$ and divides
  $\alpha$ by $g$ (`emis.f90:859-866`; `grtrans_driver.f90:259`), so the
  integrator works with the observed-frame $I_\nu$ directly. The driver then
  divides $j$ by ${\tt fac}=\sum j\,(\lambda_1-\lambda)$ and multiplies
  $\alpha$ by `LBH` (lengths in cm, $GM/c^2$), and multiplies the result back
  (`:227, 260, 271`). A ray with ${\tt fac}\le0$ (no emission) is skipped
  and stored as 0 (`:235-238, 448-450`).
- **Optical depth** runs from the camera end, a trapezoid sum of
  $\alpha\,d\lambda$ (`emis.f90:810-823`; `grtrans_driver.f90:264-265`).
- **`explin`** (`iname='explin'`, our addition; `radtrans_integrate.f90:890-925`):
  $I=\sum_ie^{-\tau_i}[S_i(w_0-w_1)+S_{i+1}w_1]$ with $S=j/|\alpha|$,
  $D=\tau_{i+1}-\tau_i$, $w_0=1-e^{-D}$, $w_1=(1-e^{-D}-De^{-D})/D$: exact
  for a source function linear in $\tau$ across each interval; the trapezoid
  rule on $j\,e^{-\tau}$ where $D<10^{-4}$; it stops where $\tau>700$. It has
  no tolerance, so its error is the sampling error of §4.
- **The others:** `lsoda` (ODE solver with `MAX_TAU = 10`, tolerances 1e-6 /
  1e-30; a failed ray is zeroed; `:24, 212-319`, `odepack.f90:91-95`),
  `delo` and `formal` (with one Stokes component both are the same
  trapezoid quadrature, `:182-189, 927-931`), `lsodasph` (polarized, "under
  development"). `iname` maps to `r%iflag` in `rad_trans.f90:29-41`; with
  `neq` ≠ 1 `explin` falls back to DELO. `lsoda` does not converge on these
  rays (muchado L14, "Integrator checks left over from the sampling
  study"; `history.tex`, "The lsoda integrator fails in optically thick
  gas").

## 10. Adaptive refinement

Off by default (`raytrace_q6.REFINE = False`); on only with `explin` and one Stokes
component (`grtrans_driver.f90:276-278`). `&refine` sets `refine_rtol`
(1e-2), `refine_atol_frac` (1e-3), `refine_tau_cut` (20),
`refine_max_depth` (6), `refine_max_points` (2000) (`read_inputs.f90:33-35,
72-78`).

- **Triggers** on the original intervals (`grtrans_driver.f90:638-642`): $\Delta\tau>0.5$; a
  change of `harmonic_index` where the harmonics are separate
  ($n\le2/(e^2\Theta_e)$); a jump in $j$ or $\alpha$ by more than 10×.
  `harmonic_index` (`:791-816`) reads `tcgs` as $\Theta_e$ for APERTURE and
  as kelvin otherwise. Before 2026-10-15 it always read kelvin, so for
  APERTURE every crossing of $\lfloor\nu/\nu_B\rfloor$ triggered a split; the
  fix removes 32% of the added points on a relabeled 230 GHz frame (dump
  500, $d=4$, $64^2$, 400 samples; flux $-0.05\%$) and hardly any at
  1.6 MHz.
- **Each level** halves the active intervals, places new points on a cubic in
  $\lambda$ through four neighbors (`:842-858`), evaluates the full fluid →
  frame → emissivity chain there, and keeps splitting where the two halves
  change the interval's contribution by more than `rtol` (`:667-719`).
  Intervals behind $\tau\ge$ `tau_cut` or below `atol_frac` of the total are
  dropped. Over budget, the intervals nearest the camera are kept.
- **What it cannot find: a hot cell between two samples.** Every trigger is
  evaluated from the values at the two end points of an existing interval:
  their $\Delta\tau$, their harmonic numbers, the ratio of their $j$ and of
  their $\alpha$. Suppose an isolated emitting cell (one or a few bright
  cells with dark neighbors, as in the relabeled 230 GHz images at high
  $R_{\rm high}$) lies wholly between samples $k$ and $k+1$. Both end points
  sit in dark gas, so $j_k\approx j_{k+1}\approx0$, $\alpha$ too, and
  $\Delta\tau\approx0$: no trigger fires, the interval is never bisected, and
  its `explin` contribution (from the end-point values) is about zero, which
  also passes the convergence test, since there is nothing to compare it
  with. The cell is lost exactly as in a uniform trace with the same
  samples. Refinement resolves structure the samples already touch (a sample
  inside or at the edge of the bright cell gives the factor-10 jump); it does
  not search for structure no sample touched. At $d=1$ a ray crosses
  1350–4200 cells, and 3200 uniform samples miss 9–36% of them (§4), which
  is why refined traces from 200 samples reached only 0.23–0.73 of the
  converged flux there (`history.tex`, Stages 29, 31) and why production uses
  uniform sampling dense enough to touch the cells (`history.tex`, Stage 48).
  Also, the point budget is $\max($`refine_max_points`$,$ `nn(3)`$)$, so at
  `nn(3)` $\ge2000$ (2000 by default) refinement adds no points at all.
- The image uses the refined points; the extras stay on the original ones
  (`:629-630`). Rays with fewer than four points are not refined (the cubic
  needs four; `:644-646`). Statistics (stdout and `refine.out` in the working
  directory) are initialized and written by `grtrans_main`
  (`pgrtrans.f90:171, 242`), so the f2py path writes them too.

## 11. Outputs and the extra quantities

`grtrans.out` (`camera.f90:322-341`; `cflag=0`) is a sequential unformatted
file, one block per camera: `nx, ny, nvals+nextra`; `nkey=3`; three keys
($n_x$, $n_y$, $\nu$); `pixloc(2, nx·ny)` ($\alpha,\beta$ interleaved);
`pixvals(nvals+nextra, nx·ny)` (float32). The blocks follow the camera
index `(l-1)·nfreq·nparams + (m-1)·nfreq + k` (time $l$ slowest, then model
parameter $m$, frequency $k$ fastest; `grtrans_driver.f90:443-459`); a run
over several binaries (`nfreq=1`) writes block $m$ for the $m$-th binary of
the list. The pipeline's one decoder is
`raytrace_q6.read_grtrans_image` (`pipeline/raytrace_q6.py:509-605`), which
takes the trace's `nvals` (or infers it from the row count: 1, 4, 22 or 25
rows) and returns $I$ from row 0 and the extras from the rows after the
`nvals` Stokes rows, so a 4-Stokes trace is decoded correctly; a row count
that fits neither is an error. It returns the first block unless asked for
another (`block=i`, 0-based) or for all (`block="all"`, a list);
`grtrans_sandbox.Sandbox.read_all` reads one per binary.

With `extra=1` each pixel carries 21 more values (`pgrtrans.f90:111-119`;
7 if `nup=1`, which the reader does not expect). Each is a line-of-sight
average weighted by $w=j\,e^{-\tau}$ (the invariant emissivity, attenuated
toward the camera; `grtrans_driver.f90:318`):
$\langle q\rangle=\int wq\,d\lambda\big/\int w\,d\lambda$. Slot 1 is the total
$\tau$; 7–14 are $\langle r\rangle$, $\langle\theta\rangle$,
$\langle\phi\rangle$, $\langle n_e\rangle$, $\langle\Theta_e\rangle$,
$\langle B\rangle$, $\langle2p/b^2\rangle$, $\langle{\rm sign}\cos\theta\rangle$;
20 and 21 are the pipeline's $\langle\beta_p\rangle$ and $\langle\sigma\rangle$
(`:329-390`). The stores keep slots 10, 11, 12, 20, 21. All of 7–14, 20, 21
are accumulated in one pass over the ray (since 2026-10-17: each is
$\sum_i c_iw_iq_i/\sum_i c_iw_i$ with the trapezoid weights $c_i$); a ray where
the weight integrates to 0 gets 0 in all of them. Slot 2 is 0 with one Stokes component, except on
a ray skipped for a NaN emissivity, where it holds the number of NaN
samples (§12). The extras need no new emissivity evaluations and, since the
one-pass rewrite of 2026-10-17, no measurable time (a $128^2$ frame at 3200
samples: 10.5 s with them, against 52.5 s with the earlier dozen `tsum`
passes per ray and 10.5 s without; byte-identical outputs; `history.tex`,
Stage 62).

## 12. Failures

| Case | Result |
|---|---|
| NaN in $j$ | `fac` is NaN, the ray is skipped: pixel and extras 0, but extra slot 2 = the ray's number of NaN samples; a "NaN in emissivity" line per ray and, after the ray loop, "WARNING: N ray(s) had a NaN emissivity" on stdout (`grtrans_driver.f90:202-207, 240, 429-436`; `pgrtrans.f90:233-235`). The pipeline captures stdout and does not read slot 2 (muchado L8-G4, pipeline side) |
| NaN only in $\alpha$ | `explin` returns NaN, "NaN in integrate!" (`radtrans_integrate.f90:204-205`); the pipeline masks non-finite pixels (`raytrace_q6.py:1230-1260`) |
| No samples / geodesic failure | pixel 0 (`grtrans_driver.f90:574-583`) |
| `lsoda` failure | ray zeroed (`odepack.f90:91-95`) |
| Negative or too-hot pixels | not clipped in Fortran; masked by the pipeline (same lines) |

## 13. Threads

The ray loop is `!$omp parallel do schedule(guided)` (`pgrtrans.f90:217-228`;
`static,1` until 2026-10-17, about 5% slower: rays differ in cost).
Each thread owns the module state `f, r, g, e, EXTRA_QUANTS, LBH, MBH`
(`threadprivate`, `grtrans_driver.f90:45`), the index of the fluid binary it
traces (`aperture_active`, `fluid_model_aperture.f90:108-109`; the binaries
themselves are shared and only read) and the integrator's work arrays;
pixels write their own columns. The table is loaded once in a critical
section (`jnu_table.f90:104-140`). The once-per-run flags `warned_hot` and
`warned_cold` are read atomically and set in a critical section
(`:168-205`); the hot-edge warning is reachable only with `THETA_SYNCH` set
above the table's largest $\Theta_e$. The NaN-ray counter `nan_rays` is
updated atomically (`grtrans_driver.f90:205-206`); it counts (ray, frequency, binary) triples, so a list run can count a ray more than once. `calc_emissivity`
makes four work arrays of $352$ bytes per sample in all (1.1 MB at 3200
samples, 9 MB at 25600; `emis.f90:473`). This gfortran build puts them on the
heap, not the thread stack: traces at 25600, 65536 and 131072 samples (up to
46 MB per call, refinement included) ran on 28 threads although the thread
stack is 8 MB and the pipeline sets no `OMP_STACKSIZE` (2026-10-19,
`history.tex`, Stage 70; review item L8-G7, closed). A compiler or flags that
put such arrays on the stack (`-fstack-arrays`, `-Ofast`) would need that
check again, or `OMP_STACKSIZE` raised (e.g. `64M`).

## 14. What d and R_high do

grtrans does not know about $d$ or $R_{\rm high}$. Both act through the
binary of §5: $d$ is the coarse-graining of the fluid before it is written,
$R_{\rm high}$ only the electron temperature $\Theta_e=N\Theta_e^{\rm
sim}/R(\beta_p)$ in it. Everything below is measured on the relabeled
230 GHz images ($128^2$, `explin`, uniform sampling, $\sigma_{\rm cut}=0.5$)
of dumps 250, 500 and 1000 (`history.tex`, Stage 40;
`verification_scripts/rhigh_sweep.py`, `rhigh_trace.py`; outputs
`rhigh_sweep.txt`, `rhigh_trace_{cell,point}.txt`).

### What d changes

- **The grid.** At $d$ the binary has $(2048/d)^2$ cells, each the average of
  $d\times d$ output cells ($2d\times2d$ simulation cells) with
  $\sqrt{-g}$ weights (`coarsen_fluid.py`, rule "datalib"); the cell centers
  move accordingly, and the no-emission box of §6 grows by $d$.
- **The temperature.** Averaging the moments over a tile turns velocity
  differences between its cells into pressure, so a coarser tile is hotter:
  at dump 500 the share of cells with $\Theta_e>0.05$ inside $r=30$ is 15, 36
  and 61% at $d=1$, 2, 4 (`coarse_theta_shares.py`). The heated electrons'
  median $\Theta_e$ moves only a little (table below); the hot tail grows.
- **Sampling.** A ray crosses about $d$ times fewer cells, so the same
  samples resolve them better, and isolated hot cells merge with their
  neighbors.
- **Run time** hardly depends on $d$; it scales with pixels × samples.

### What R_high changes

- **The temperature.** The torus body has $\beta_p\sim10$–80, so there
  $R\approx R_{\rm high}$; only the funnel walls ($\beta_p\lesssim1$) fall to
  $R_{\rm low}=1$ (`output/rhigh_sweep_rtheta.png`). The heated electrons
  ($\Theta_e^{\rm sim}>10^{-4}$, number-weighted median) sit at
  $\Theta_e\approx0.2$ at $R_{\rm high}=100$, 1 at 20, 3–6 at 5 and 16–26 at 1;
  none exceeds the table's top, so the synchrotron-fit branch is not used.
- **Where the light comes from.** At $R_{\rm high}=100$ the torus is too cold
  for 230 GHz, and the light comes from the few hot cells with
  $\beta_p\lesssim0.3$ near the hole and beside the σ mask: often isolated,
  with one to three macro-particles. Lowering $R_{\rm high}$ lets the torus
  and its sheath emit, and the light moves outward (half-light radius
  $r\approx1.1$–1.3 at $R_{\rm high}=100$, 8–21 at 1, $d=1$).
- **The σ cut** matters much less at low $R_{\rm high}$, because less light
  sits beside the mask: changing $\sigma_{\rm cut}$ from 0.25 to 2.5 moves the
  flux by 0.96–1.02 at $R_{\rm high}=1$ and 0.38–1.78 at 100 ($d=4$), and
  testing it per sample instead of per cell by 0.99–1.00 against 0.33–0.95.
- **The run is not a MAD** ($\phi_{\rm BH}\approx5$ Heaviside–Lorentz, 19
  Gaussian; `pipeline/horizon_flux.py`, notebook P3d), so the MAD-based
  preference for $R_{\rm high}=1$–5 (Mizuno et al. 2021) does not carry over
  directly. $R_{\rm high}$ stays a free parameter of each trace (decided
  2026-10-15, `history.tex`, Stage 48; `pc.R_HIGH` = 100 is only the default).

### Both together

Ranges over the three dumps; per-cell σ cut unless noted. *Θ_e electrons:*
number-weighted median of the heated electrons. *Θ_e light:* $j_\nu dV$-weighted
median. *≤3 particles, isolated (fluid):* shares of the local 230 GHz
luminosity from cells with at most three electron macro-particles, and from
emitting cells with no emitting 4-neighbor. *F/F(d=1):* converged flux
against $d=1$. *800/ref, 3200/ref:* flux at 800 and 3200 samples against
the reference (12800 at $d=1$, 2; 3200 at $d=4$); *p90:* 90th-percentile
pixel error at 3200. *Isolated (traced):* image flux lost when isolated cells
are emptied. *n99:* pixels holding 99% of the flux (of 16384).
*Point/cell:* per-sample against per-cell cut.

| d | R_high | Θ_e electrons | Θ_e light | ≤3 particles | isolated (fluid) | F/F(d=1) | 800/ref | 3200/ref | p90 @3200 | isolated (traced) | n99 | point/cell |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1 | 16–22 | 28–41 | 5–44% | 1–9% | 1 | 0.97–1.00 | 0.995–1.000 | 0.4–1.1% | 0.7–7% | 14168–15541 | 0.99–1.00 |
| 1 | 5 | 3.2–4.7 | 14–21 | 50–86% | 19–30% | 1 | 0.86–0.99 | 0.96–1.00 | 2.6–5.1% | 10–22% | 4770–13568 | 0.97–0.99 |
| 1 | 20 | 0.81–1.2 | 13–14 | 92–96% | 55–59% | 1 | 0.79–0.96 | 0.93–1.00 | 4.5–19% | 22–39% | 3014–15499 | 0.87–0.97 |
| 1 | 100 | 0.16–0.24 | 7.9–9.2 | 84–99% | 72–77% | 1 | 0.72–1.19 | 0.81–1.01 | 7.7–50% | 15–87% | 1507–15738 | 0.33–0.95 |
| 2 | 1 | 17–25 | 27–42 | 4–22% | 1–3% | 1.1–2.0 | 1.00–1.01 | 1.000–1.002 | 0.1–0.2% | 0.4–2% | 12857–15356 | 0.99–1.00 |
| 2 | 5 | 3.5–5.4 | 18–26 | 47–54% | 9–11% | 2.4–5.0 | 0.96–1.00 | 0.996–0.999 | 0.4–0.6% | 6–8% | 3454–9441 | 0.96–0.98 |
| 2 | 20 | 0.87–1.4 | 17–19 | 78–88% | 20–31% | 5.9–15 | 0.88–0.98 | 0.97–0.99 | 0.8–1.8% | 13–22% | 1754–12489 | 0.88–0.93 |
| 2 | 100 | 0.17–0.28 | 11–12 | 91–97% | 36–51% | 7.1–15 | 0.80–0.92 | 0.94–0.99 | 1.3–6.8% | 25–40% | 1048–13707 | 0.78–0.81 |
| 4 | 1 | 17–26 | 27–48 | 2–7% | 0–1% | 1.0–2.6 | 1.00 | (ref) | — | 0.1–0.7% | 12402–15333 | 0.99–1.00 |
| 4 | 5 | 3.5–5.6 | 23–33 | 13–22% | 2–3% | 4.1–11 | 1.00–1.01 | (ref) | — | 0.9–2.3% | 3105–8114 | 0.96–0.97 |
| 4 | 20 | 0.87–1.4 | 21–26 | 23–46% | 4–7% | 14–64 | 0.995–1.01 | (ref) | — | 2–6% | 1518–6703 | 0.88–0.94 |
| 4 | 100 | 0.17–0.28 | 13–16 | 36–74% | 9–18% | 32–110 | 0.96–0.98 | (ref) | — | 6–19% | 861–7931 | 0.70–0.85 |

What the table says:

- **Electron temperature** is set by $R_{\rm high}$ (a factor of about 100
  between 1 and 100) and hardly by $d$; the light, though, always comes from
  gas far hotter than the typical electron ($\Theta_e\approx8$–48), because
  only that gas emits at 230 GHz.
- **Convergence in samples** improves with both: lower $R_{\rm high}$ spreads
  the light over many cells, larger $d$ merges the cells. At $R_{\rm high}=1$
  even $d=1$ needs only 3200 samples (800 give 0.97–1.00); at 100 even 12800
  are marginal at $d=1$.
- **The flux** depends on $d$ mostly through $R_{\rm high}$: at 100 the coarse
  grids are 32–110 times brighter than $d=1$, at 1 only 1–2.6 times. At high
  $R_{\rm high}$ the light comes from the few-particle cells that
  coarse-graining heats; at low $R_{\rm high}$ from the bulk gas, which it
  hardly changes.
- **What stays** at low $R_{\rm high}$: at $d=1$, cells with at most three
  macro-particles still give up to 44% of the light (dumps 500 and 1000);
  $d=2$ or 4 removes most of that.

Samples needed (measured 2026-10-16 over $64^2$–$512^2$, three dumps, the
per-sample cut; `verification_scripts/nup_npix_convergence.py`,
`RT_notebook_overview.tex`, "Convergence in `nn(3)` against image size,
$R_{\rm high}$ and $d$"): image size does not matter. The table gives the
fewest samples per ray that keep both the median flux within 1% of the
12800-sample value and the median pixel error below 1%, at that count and at
every larger one (median over dumps 250, 500, 1000; pixel error over the
pixels holding 99% of the reference flux; valid for all four image sizes;
"> 6400" means no count below the 12800-sample reference qualifies):

| | $R_{\rm high}=1$ | 5 | 20 | 100 |
|---|---|---|---|---|
| $d=4$ | 400 | 800 | 1600 | 6400 |
| $d=2$ | 800 | 3200 | 6400 | > 6400 |
| $d=1$ | 1600 | 3200–6400 | > 6400 | > 6400 |

(12800 samples are within 1–2% of 25600 at $d=1$.) These replace the rougher
settings above (judged at $128^2$, per-cell cut). Blurred to the EHT's
20 µas beam, most images look the same much earlier (cross-correlation within
$10^{-3}$ from 200 samples at $d=4$, 200–1600 at $d=2$, 200–6400 at $d=1$;
`verification_scripts/image_similarity.py`; `RT_notebook_overview.tex`, "Is an
undersampled image qualitatively the same?").

## 15. Our changes to stock grtrans

`git -C grtrans diff c76cb11 aperture` shows all of them (the changed
upstream files and the added ones below).

- **Added:** `fluid_model_aperture.f90` (the binary and its interpolation,
  §5–6), `jnu_table.f90` (table, dispatch, synchrotron fit, §8),
  `cyclotron_emis.f90` (§8), `jnu_table.bin`, `test_aperture_model.f90` (at
  cell centers: the stored $n_e$, $\Theta_e$, $u$, $b$ come back, and
  `uks2ubl` preserves $u\cdot u$, $b\cdot b$, $b\cdot u$ between the
  Kerr–Schild and Boyer–Lindquist metrics; build line in its header),
  `inputs.in.default` (the namelists at the pipeline defaults, refinement
  off; `inputs.in` itself is git-ignored, since every trace rewrites it, and
  `raytrace_q6` copies it from the default when it is missing), `files.in`,
  `plot_grtrans_jnu.ipynb`, and this guide.
- **`fluid.f90`:** the `APERTURE` model, its dispatch and
  `convert_fluidvars_aperture` (§6); `rho2`/`Be` carry $\beta_p$ and $\sigma$.
- **Several binaries in one run** (§5): `fluid_model_aperture.f90` (per-binary
  storage, `.lst` lists, `select_aperture_binary`), `grtrans_driver.f90`
  (binary $m$ per model parameter), `pgrtrans.f90` (the `nmdot` check).
- **`emis.f90`:** the `JNUTABLE` emissivity type; `calc_opt_depth` as a
  running trapezoid.
- **`grtrans_driver.f90`:** refinement (`refine_ray`, `refine_eval_points`,
  `harmonic_index`, statistics); the `explin` hook; slot 1 = total $\tau$,
  slots 20/21; guards in the extras; the NaN-ray count and slot-2 flag; the
  `nvals=4` stop.
- **`radtrans_integrate.f90`, `rad_trans.f90`:** the `explin` integrator
  (`iflag=4`); `oatol` 1e-8 → 1e-30; `nextra` 21.
- **`read_inputs.f90`:** the `&refine` namelist. **`pgrtrans.f90`:**
  `nextra` 19 → 21 (7 still if `nup=1`); refinement statistics and the
  NaN-ray summary (moved from `grtrans.f90` on 2026-10-15). **`grtrans_program.f90`:** longer file names.
- **`phys_constants.f90`:** CODATA 2018 constants, $M_\odot=1.98841\times10^{33}$ g.
- **`fits.f90`:** no-op stubs (no cfitsio); **`Makefile`:** the new objects,
  no `-lcfitsio`.
- **The other `fluid_model_*.f90`:** namelists moved after the declarations
  (gfortran).

## 16. Known problems

In `docs/muchado.md` (ids as of 2026-10-15): the cyclotron line width near
$90^\circ$ (L1-R4); the table's floor (L8-F7, no effect: history.tex,
Stage 88), the no-emission box (L8-F9, Stage 93) and the NaN flag of extra
slot 2 (L8-G4, Stage 92) are closed. (L1-R18, the extra slots averaging empty corners as zero, was fixed on 2026-10-19; L8-G7, stack use, was closed the same day: the arrays are on the heap.)
Fixed on 2026-10-15 (`history.tex`): M4-F1, L8-F2, F4 (documented), F5, F10,
F11, F12, G4 (Fortran side), G6, G8, C8, L1-R7, L1-R17.
