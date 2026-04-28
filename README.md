# Spacey.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://glwhart.github.io/Spacey.jl/stable)
[![Runtests](https://github.com/glwhart/Spacey.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/glwhart/Spacey.jl/actions/workflows/Runtests.yml)
[![codecov](https://codecov.io/gh/glwhart/Spacey.jl/branch/main/graph/badge.svg?token=OYAPFZI28I)](https://codecov.io/gh/glwhart/Spacey.jl)

Spacey.jl finds the **point group** of a Bravais lattice and the **space group** of a crystal. Other symmetry packages — [FINDSYM](https://iso.byu.edu/) (part of the ISOTROPY suite), [spglib](https://spglib.readthedocs.io/), [aflow-sym](https://aflow.org/), the [Bilbao Crystallographic Server](https://cryst.ehu.eus/), [PLATON](https://platonsoft.nl/) (with ADDSYM), and the symmetry routines in [cctbx](https://github.com/cctbx/cctbx_project) — additionally provide spacegroup names, Niggli cell reduction, and standard settings. Spacey deliberately does not; it provides the *group operations*.

The primary use of this package is to provide spacegroup/pointgroup symmetry operations for downstream computational or analysis tasks. It was created primarily for coupling with [Enumlib.jl](https://github.com/glwhart/Enumlib.jl) cluster-expansion applications and automated k-point generation (such as [autoGR](https://github.com/msg-byu/autoGR)).

One unique feature of the package is its robustness to finite-precision issues — a constant bane of other spacegroup codes. See the [documentation](https://glwhart.github.io/Spacey.jl/stable) for the discussion of tolerance handling, over-promotion, and the `verify_stable` flag.

Another unique feature is the [snap-to-symmetry function](https://glwhart.github.io/Spacey.jl/stable/how-to/snap-to-symmetry/), which takes a noisy lattice and produces the exact-symmetry version closest to the input.

Please see the [documentation landing page](https://glwhart.github.io/Spacey.jl/stable) for examples, how-to guides, and explanations.

Spacey finds symmetries efficiently* by relying on the input being as compact as possible. Internally it Minkowski-reduces the input via [MinkowskiReduction.jl](https://github.com/glwhart/MinkowskiReduction.jl).

(* Normally, efficiency is not important — finding the symmetries of a single lattice is a quick computation. When the symmetries of tens of thousands of cases are needed in a second or two, as in [autoGR](https://github.com/msg-byu/autoGR), efficiency becomes essential.)
