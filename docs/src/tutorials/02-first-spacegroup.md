# Find the space group of NaCl (and three other crystals)

This tutorial walks through finding the **space group** of a crystal — the set of symmetry operations that combine rotations/reflections with translations to map an atomic structure onto itself. We'll work through four common materials, each illustrating a different way that atoms in a cell affect the symmetry count: rocksalt (NaCl), zincblende (GaAs), an ordered alloy (L1₀ CuAu), and a ferroelectric perovskite (BaTiO₃).

By the end you'll know how to:

- Build a `Crystal` from a lattice + atomic positions + per-atom type labels,
- Call `spacegroup`, and
- Read the resulting `SpacegroupOp` objects.

This tutorial picks up where [Find the point group of a cubic lattice](01-first-pointgroup.md) left off. If you haven't read it, do so first — `pointGroup` is the lattice-only step that `spacegroup` builds on. Plan on about fifteen minutes.

## 1. Set up

```jldoctest tut2
julia> using Spacey, LinearAlgebra
```

(See tutorial 01 for installation help if `using Spacey` errors.)

## 2. NaCl — rocksalt

Sodium chloride is the canonical ionic crystal. Geometrically, NaCl is a face-centered cubic lattice with Na on one FCC sublattice and Cl on a second FCC sublattice offset by `(½, ½, ½)`. There are eight atoms in the conventional cubic cell (4 Na + 4 Cl), but the *primitive* cell — the smallest translational unit — has just two:

```jldoctest tut2
julia> A = [0.0  0.5  0.5;        # FCC primitive lattice vectors (each column = one a_i)
            0.5  0.0  0.5;
            0.5  0.5  0.0];

julia> r = [0.0  0.5;              # Na at (0, 0, 0), Cl at (½, ½, ½) — both in primitive coords
            0.0  0.5;
            0.0  0.5];

julia> c = Crystal(A, r, [:Na, :Cl]; coords=:fractional);
```

A few things to notice in that `Crystal` constructor:

- `r` is a 3×N matrix; each *column* is one atom's position. Same column-as-vector convention as `A`.
- `[:Na, :Cl]` is the per-atom type list — Symbols, in the same order as the columns of `r`.
- `coords=:fractional` is **required**. Spacey makes you say it explicitly to prevent the silent footgun of misinterpreting `[0.5, 0.5, 0.5]` as Cartesian when it was meant as fractional or vice-versa. (See [Construct a crystal](../how-to/construct-a-crystal.md) for the full argument.)

Now find the space group:

```jldoctest tut2
julia> ops = spacegroup(c);

julia> length(ops)
48
```

**Forty-eight ops.** That's the full cubic point group (48) times 1 (no centering translations — they're absorbed into the primitive lattice). NaCl's space group is **Fm3̄m**, the highest-symmetry space group compatible with two distinct atom types on an FCC lattice.

The conventional 8-atom FCC cell of the same crystal returns 192 = 48 × 4 ops, because the four FCC centering translations explicitly multiply through. The *physics* is the same; the *count* depends on which cell you wrote down. We'll stick with primitive cells throughout this tutorial — they're smaller, faster, and the 48-or-not story is cleaner.

## 3. GaAs — zincblende

Gallium arsenide has the **same primitive lattice as NaCl** but a different atom arrangement: Ga at the origin and As at `(¼, ¼, ¼)` instead of NaCl's `(½, ½, ½)`. This shift breaks the inversion symmetry — there's no longer a point in the cell about which `r → -r` maps Ga to itself.

```jldoctest tut2
julia> A = [0.0  0.5  0.5;        # same FCC primitive lattice as NaCl
            0.5  0.0  0.5;
            0.5  0.5  0.0];

julia> r = [0.0  0.25;             # Ga at origin, As at (¼, ¼, ¼)
            0.0  0.25;
            0.0  0.25];

julia> c = Crystal(A, r, [:Ga, :As]; coords=:fractional);

julia> length(spacegroup(c))
24
```

**Twenty-four ops** — exactly half of NaCl's 48. The half that's gone is everything with `det = -1` (inversion, mirrors, rotoinversions); only the proper rotations survive. The space group is **F4̄3m** (T_d point group). This same structure (different elements at the two Wyckoff positions) covers an enormous chunk of the optoelectronic-semiconductor world: GaAs, ZnS (the structure's namesake), AlAs, InSb, CdTe, ZnSe, and most III-V and II-VI compounds.

So: **same lattice, different atomic ordering, half the symmetries.** That's the point of the space group — atoms matter, not just the cell.

## 4. L1₀ CuAu — ordered alloy

Now consider an alloy where the lattice itself is geometrically cubic but the *chemistry* breaks the symmetry. CuAu in its ordered low-temperature phase is FCC, but with Cu and Au alternating in layers along z. Picture a stack of (001) planes: Cu, Au, Cu, Au, …

```jldoctest tut2
julia> a1 = [0.5,  0.5, 0.0];        # face-center vector, +y

julia> a2 = [0.5, -0.5, 0.0];        # face-center vector, -y

julia> a3 = [0.0,  0.0, 1.0];        # along z

julia> r_cart = [0.0  0.0;            # Cu at (0, 0, 0), Au at (0, ½, ½) — Cartesian
                 0.0  0.5;
                 0.0  0.5];

julia> c = Crystal(a1, a2, a3, r_cart, [:Cu, :Au]; coords=:cartesian);

julia> length(spacegroup(c))
16
```

**Sixteen ops.** This is **P4/mmm** (D_4h point group). Two new things in this example:

- We used the **three-vector form** of the constructor (`Crystal(a1, a2, a3, r, types; coords)`) instead of the matrix form. Both are equivalent — pick whichever reads better in context.
- We used `coords=:cartesian`. The constructor converts Cartesian → fractional once and stores the fractional form. We could have written `r` in fractional coordinates directly; doing it Cartesian here matches how the L1₀ structure is usually described in textbooks ("Au at the body center of the conventional cell").

The lattice is *geometrically* still cubic (`|a₃| = 1` matches the in-plane vectors). What breaks the symmetry to tetragonal is the chemistry: a 4-fold rotation about x or y would swap Cu and Au layers, which is no longer a symmetry once they're labeled distinct. Only the 4-fold along z survives, plus the mirrors compatible with it.

This L1₀ pattern shows up in many ordered binary alloys: FePt (used in magnetic recording), TiAl, NiPt, CuAu itself.

## 5. BaTiO₃ — ferroelectric perovskite

Barium titanate is the textbook ferroelectric. Above 120 °C it has the cubic Pm3̄m perovskite structure: Ba at the corner, Ti at the body center, O at the three face centers. Cool below the Curie point and the Ti shifts by ~0.05 Å along z, dragging the oxygens with it. The polarization is the displacement; the symmetry drops from 48 (cubic) to 8 (tetragonal P4mm).

We'll build the room-temperature structure, with Ti shifted by `ε = 0.05` Å:

```jldoctest tut2
julia> A = 4.0 * Matrix{Float64}(I, 3, 3);     # cubic lattice, ~4 Å edge

julia> ε = 0.05;                                # Ti displacement in Å

julia> δ = ε / 4.0;                             # in fractional coords

julia> r = [0.0  0.5    0.5  0.5  0.0;          # Ba   Ti      O    O    O (face centers)
            0.0  0.5    0.5  0.0  0.5;
            0.0  0.5-δ  0.0  0.5  0.5];

julia> c = Crystal(A, r, [:Ba, :Ti, :O, :O, :O]; coords=:fractional);

julia> length(spacegroup(c; pos_tol=1e-4))
8
```

**Eight ops** — the ferroelectric P4mm group (C_4v point group). The Ti displacement is small (12 mÅ in fractional units), so we passed an explicit `pos_tol=1e-4` to make sure the algorithm resolves it. With the default `pos_tol`, Spacey would treat the displaced Ti as if it's still at body center (the displacement is below the default tolerance for this cell volume) and report the cubic 48 ops — the *parent* symmetry of the ferroelectric. Both answers are "correct" in the sense that they correspond to a real physical regime; which one you want depends on the question. See [Handle noisy real-world data](../how-to/handle-noisy-data.md) for how to choose `pos_tol`, and [Detect tolerance-dependent answers](../how-to/detect-tolerance-dependence.md) for the `verify_stable` flag that catches near-boundary cases like this automatically.

To see the parent regime, try the same construction with `ε = 0.0`:

```jldoctest
julia> using Spacey, LinearAlgebra

julia> A = 4.0 * Matrix{Float64}(I, 3, 3);

julia> r = [0.0  0.5    0.5  0.5  0.0;
            0.0  0.5    0.5  0.0  0.5;
            0.0  0.5    0.0  0.5  0.5];

julia> c = Crystal(A, r, [:Ba, :Ti, :O, :O, :O]; coords=:fractional);

julia> length(spacegroup(c))
48
```

A one-character change in the input (`0.5-δ` → `0.5`) flips between two physically distinct space groups, separated by less than 0.1 Å of structural distortion.

## What you learned

- `Crystal(A, r, types; coords)` builds a crystal from a lattice + atomic positions + per-atom type labels.
- `coords` is required: `:fractional` or `:cartesian`. The Cartesian form is converted to fractional once at construction.
- `spacegroup(c)` returns a `Vector{SpacegroupOp}` of every symmetry operation; `length(...)` gives the order.
- The number depends on **both** lattice geometry **and** atomic ordering. Same lattice, different atoms (NaCl ↔ GaAs) can halve the count. Same lattice, layered chemistry (L1₀) can drop it further. Atomic *displacements* (BaTiO₃) can drop it again.

## Next steps

- How-to: [Find a space group](../how-to/find-spacegroup.md) — tolerance tuning, sanity checks, primitive vs. conventional cells.
- How-to: [Compose and apply operations](../how-to/compose-and-apply-ops.md) — what to do with the `SpacegroupOp` objects.
- How-to: [Detect tolerance-dependent answers](../how-to/detect-tolerance-dependence.md) — the `verify_stable` flag for ferroelectric-class near-misses.
- Reference: [`Crystal`](../reference/crystals.md), [`spacegroup`](../reference/space-groups.md), [`SpacegroupOp`](../reference/space-groups.md).
