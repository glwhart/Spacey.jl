# AI Code Explanations

Reference notes where Claude has explained unfamiliar Julia syntax or patterns used in this repo. Written up at the point of use so they can be revisited later.

---

## Parametric types, inner/outer constructors, broadcasting

**Context (2026-04-23):** The `Crystal` struct proposed in `spacegroup_plan.md` §4.2 uses several Julia features that aren't elsewhere in the current Spacey.jl code. This section walks through the syntax piece by piece.

The code under discussion:

```julia
struct Crystal{T}
    A::Matrix{Float64}          # 3×3 lattice, columns = a1, a2, a3
    r::Matrix{Float64}          # 3×N fractional positions, columns = atoms
    types::Vector{T}            # length N, any type comparable with ==
    function Crystal{T}(A::AbstractMatrix{<:Real}, r::AbstractMatrix{<:Real},
                        types::AbstractVector{T}; coords::Symbol) where T
        size(A) == (3, 3) || error("A must be 3×3")
        size(r, 1) == 3 || error("r must have 3 rows")
        size(r, 2) == length(types) || error("r columns must match types length")
        coords ∈ (:fractional, :cartesian) || error("coords must be :fractional or :cartesian")
        A64 = Float64.(A)
        r64 = coords === :cartesian ? inv(A64) * Float64.(r) : Float64.(r)
        new{T}(A64, r64, collect(types))
    end
end
Crystal(A, r, types::AbstractVector{T}; coords) where T = Crystal{T}(A, r, types; coords=coords)
Crystal(a1, a2, a3, r, types; coords) = Crystal(hcat(a1, a2, a3), r, types; coords=coords)
```

### The struct header

```julia
struct Crystal{T}
    A::Matrix{Float64}
    r::Matrix{Float64}
    types::Vector{T}
end
```

`{T}` declares `Crystal` as a **parametric type**. `T` is a type parameter — a placeholder for some *other* type that will be chosen per-instance. When you later write `Crystal{Symbol}(...)`, `T` becomes `Symbol` and `types` becomes a `Vector{Symbol}`. When you write `Crystal{Int}(...)`, it becomes `Vector{Int}`. The struct definition is a *template*; each concrete `Crystal{T}` is a distinct type.

This is the Julia generalisation of C++ templates / Java generics. Docs: [Parametric Types](https://docs.julialang.org/en/v1/manual/types/#Parametric-Types).

Why reach for it here: we want the algorithm to work with whatever atom-label type the user chooses, without the library hard-coding the choice. Parametric types are the mechanism.

### Inner constructor signature

```julia
function Crystal{T}(A::AbstractMatrix{<:Real},
                    r::AbstractMatrix{<:Real},
                    types::AbstractVector{T};
                    coords::Symbol) where T
```

Three distinct things happening:

1. **`Crystal{T}(...)` as the function name.** Because this is *inside* the `struct ... end` block, it's an **inner constructor**. Inner constructors are the only place you can call `new{T}(...)` (the low-level "actually build an instance of this struct" primitive). By defining one, you take control away from Julia's default constructor and impose your own validation/conversion. Docs: [Inner Constructor Methods](https://docs.julialang.org/en/v1/manual/constructors/#man-inner-constructor-methods).

2. **Argument type constraints.** `A::AbstractMatrix{<:Real}` means "any matrix whose element type is a subtype of `Real`." The `<:` is Julia's subtype operator. `AbstractMatrix` is the supertype of `Matrix`, `SparseMatrixCSC`, `Adjoint`, etc., so this accepts any matrix-like thing (including integer matrices, rational matrices, views). Docs: [Type Declarations](https://docs.julialang.org/en/v1/manual/types/#Type-Declarations).

3. **`where T` at the end.** This introduces `T` as a method-level type parameter, letting Julia *infer* `T` from whichever `types::AbstractVector{T}` the user passes. If the user calls `Crystal(A, r, [:Fe, :O]; coords=...)`, Julia sees `AbstractVector{Symbol}`, binds `T = Symbol`, and picks the `Crystal{Symbol}` instantiation. Docs: [Parametric Methods](https://docs.julialang.org/en/v1/manual/methods/#Parametric-Methods).

4. **`; coords::Symbol`.** Everything after the semicolon is a keyword argument. `::Symbol` asserts it must be a `Symbol` (so `:fractional` works; `"fractional"` doesn't). No default value means it's **required** — a `MethodError` will be raised if the user omits it.

### Body

```julia
size(A) == (3, 3) || error("A must be 3×3")
```

Short-circuit `||`. Julia evaluates the left operand first; if it's `true`, it returns `true` and skips the right. If it's `false`, it evaluates the right. So this reads as: "assert `size(A) == (3,3)`, otherwise error." A common Julia idiom for guards. Docs: [Short-Circuit Evaluation](https://docs.julialang.org/en/v1/manual/control-flow/#Short-Circuit-Evaluation).

```julia
coords ∈ (:fractional, :cartesian) || error(...)
```

`∈` is Unicode for `in`. `(:fractional, :cartesian)` is a tuple. Same short-circuit pattern.

```julia
A64 = Float64.(A)
```

**Broadcasting.** The dot after `Float64` turns the constructor into an elementwise operation: "apply `Float64` to every element of `A`, return a new array with the converted elements." If `A` is already `Matrix{Float64}`, this is a no-op conversion; if it's `Matrix{Int}` or `Matrix{Rational}`, you get a `Matrix{Float64}`. Docs: [Dot Syntax for Vectorizing Functions](https://docs.julialang.org/en/v1/manual/functions/#man-vectorized).

```julia
r64 = coords === :cartesian ? inv(A64) * Float64.(r) : Float64.(r)
```

Ternary. Reads: "if `coords` is exactly `:cartesian`, compute `inv(A64) * r`; otherwise just convert." `===` is strict equality — for symbols it's the same as `==`, but it avoids any subtle method-overloading surprises. Using it with symbols is idiomatic.

The Cartesian branch converts to fractional via `r_frac = A⁻¹ · r_cart`. (That's the standard change-of-basis: if `r_cart = A · r_frac`, then `r_frac = inv(A) · r_cart`.)

```julia
new{T}(A64, r64, collect(types))
```

`new` is the primitive that actually creates the struct instance; `{T}` says "instantiate `Crystal` with this particular `T`." `collect(types)` materialises `types` into a `Vector` (in case the user handed in e.g. a range or a view).

### Outer constructors

```julia
Crystal(A, r, types::AbstractVector{T}; coords) where T =
    Crystal{T}(A, r, types; coords=coords)
```

This is a one-line method definition using `=` instead of `function`/`end`. It's an **outer constructor**: defined outside the `struct` block, cannot call `new`, but *can* call other constructors (including the inner one). Its job is to let users write `Crystal(...)` without specifying `{T}` — Julia infers `T` from the `types` argument and forwards. Docs: [Outer Constructor Methods](https://docs.julialang.org/en/v1/manual/constructors/#Outer-Constructor-Methods).

```julia
Crystal(a1, a2, a3, r, types; coords) =
    Crystal(hcat(a1, a2, a3), r, types; coords=coords)
```

Same pattern — another outer constructor. This one provides the three-vector convenience form. `hcat(a1, a2, a3)` horizontally concatenates the three column vectors into a 3×3 matrix. Then it delegates to the matrix form, which delegates to the inner constructor — all validation lives in one place.

### Why two "layers" of constructors

The discipline is:

- **Inner constructors** guard invariants (the struct can only be created via `new`, which only inner constructors can call).
- **Outer constructors** provide ergonomic forms (accepting different argument shapes, filling in defaults) by *delegating* to an inner constructor.

So the inner constructor is the one-and-only gatekeeper; outer constructors are convenience facades. Exactly the pattern here: one inner constructor does all validation and conversion; two outer methods provide the matrix-form and three-vector-form entry points.

### Reading order for the Julia manual

To back this up with docs:

1. [Constructors](https://docs.julialang.org/en/v1/manual/constructors/) — inner vs outer, `new`, parametric constructors. The core of what's going on here.
2. [Parametric Types](https://docs.julialang.org/en/v1/manual/types/#Parametric-Types) — the `{T}` syntax.
3. [Parametric Methods](https://docs.julialang.org/en/v1/manual/methods/#Parametric-Methods) — the `where T` clauses.
4. [Dot Syntax for Vectorizing Functions](https://docs.julialang.org/en/v1/manual/functions/#man-vectorized) — `Float64.(A)`.

The first one alone clears up most of the syntax in the inner constructor.
