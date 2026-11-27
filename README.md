# Standard FEM for study

Fortran implementations of triangular isoparametric FEM for Poisson and Laplace problems, and SUPG for transient transport. The P1, P2 and P3 examples retain their separate entry points and file formats. The `pre-refactor` tag marks the original development snapshot.

## Reproduce an example

Install GNU Fortran, Bash and Python 3.8 or newer. From this repository, run

```sh
./run_case.sh cases/poisson_p2.cfg
./run_case.sh cases/advec_p3.cfg steps=4 dt=0.0025 sample_step=4
python3 tests/check_cases.py
```

The Poisson example uses 25 nodes, 8 quadratic triangles, the source $-\Delta u=-4$ and the boundary values of $u=x^2+y^2$. Its center value is 0.5. The regression checks the maximum nodal error against $10^{-8}$. Inspect `res/out.vtk` in the printed run directory. The advection example writes the initial state and one VTK snapshot per step. To build without running, use `./build.sh all` or `./build.sh poisson_p2`. Executables go in `build/`. The build accepts `FC`, `FFLAGS`, `LDFLAGS` and, on macOS, `SDKROOT`.

| Cases in `cases/` | Purpose |
| --- | --- |
| `poisson_p1`, `poisson_p2`, `poisson_p3` | Scalar Poisson problems using one assembly and boundary treatment with three local bases |
| `advec_p1`, `advec_p2`, `advec_p3` | SUPG advection with a prescribed nodal velocity and initial field |
| `convecdiff_supg` | P1 convection–diffusion with manufactured time-dependent forcing |
| `laplace_p1`, `laplace_p2`, `laplace_test` | Earlier potential solvers with their original edge-based input and post-processing conventions |

Each name above has a corresponding `.cfg` file. The regression command builds and executes all ten cases, tests mixed Dirichlet/Neumann data and polynomial reproduction, checks triangle quadrature and constant-field transport, and verifies shared-edge connectivity in both mesh converters.

## Case settings and results

Case files use `key=value` lines. CLI assignments override file values. Input paths and `output_root` resolve relative to the configuration file, including when supplied on the CLI. The launcher compiles a fresh binary in a unique run directory and records `effective.cfg`, copies of the mesh and all other inputs, compiler, Git state when available, invocation, log and exit code. A run can be replayed with `./run_case.sh /path/to/run/effective.cfg`. Use `output_root=/path/to/results` to select where the runs are stored.

All cases specify `solver`, `mesh` and `boundary`. Poisson accepts `source=zero`, `quadratic` (constant forcing −4) or `sine` (forcing $2\pi^2\sin(\pi x)\sin(\pi y)$). For an order comparison, select the same source and physical boundary data and supply the matching mesh and node IDs for each degree. Standalone Poisson binaries still read the six paths in `file_FEM.dat`, accept an optional `poisson.nml` and allow `--source NAME` to override it.

The Poisson boundary format starts with the number of Dirichlet nodes followed by rows `record_id node_id value`. An optional `N NEUMANN` section has one row per exterior edge: record ID, both endpoint IDs, internal edge IDs ordered from the first endpoint toward the second, and constant outward normal derivative. At least one prescribed potential is required to fix the additive constant. The historical Laplace programs use their own labeled edge sections.

Transport cases also specify `velocity`, `initial`, `dt` and `steps`. Optional `tolerance` and `max_iterations` control the iterative solve. Advection accepts `sample_step` for a cross section (0 disables it); convection–diffusion accepts `diffusivity`. The mass diagnostic records percent of the initial integral, or the change in integral when the initial mass is zero.

## Layout and mesh tools

`src/` contains triangle shapes, quadrature, Poisson assembly, SUPG operators and field I/O. Equation-specific forcing and the public entry points remain under their problem directories. Historical inputs, preprocessing programs and visualization examples are kept with the corresponding studies.

The mesh converters share edge numbering and preserve the direction of internal nodes on a P3 edge. Their historical directory names begin with the full-width character `ｍ`. To convert a mesh using explicit paths, run

```sh
./build.sh mesh_p3
./build/mesh_p3 poisson_isopara/poisson_p1/dat/gene/mesh/cube_test_2.mes /tmp/example-p3.mes
```

The converter also reads the traditional two-line `file.txt` when invoked without arguments. The shared structured-mesh generator is in `tools/gr90.f90`. Existing case-specific `gr90.f90` entry files refer to this implementation.

The SUPG field generators can be built with `./build.sh fields_advection` or `./build.sh fields_diffusion`. Their `file.dat` lists four paths in order: input mesh, output velocity, output initial field, output boundary. Optional lines select `velocity=down|right|rotation`, `initial=box|cosine|zero` and `boundary=all|open_right`. The command-line options `--config FILE`, `--velocity`, `--initial` and `--boundary` override these settings. The legacy advection preset uses downward velocity and a box, and the diffusion preset uses rightward velocity and zero initial field. Build the original structured-grid generator with `./build.sh grid_generator`.
