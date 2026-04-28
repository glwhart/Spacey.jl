# Canonicalizing τ

A `SpacegroupOp` carries a fractional translation `τ` that — by construction in any standard space group — is one of a small number of rationals: `0`, `½`, `⅓`, `¼`, `⅙`, `1/12`, and their multiples. After enough floating-point composition and basis-transform arithmetic, those exact rationals drift to nearby Float64 values that no longer compare equal. Without a canonicalization step, group closure tests start failing on perfectly valid operation sets — and they did, on the first AFLOW corpus run.

This page explains the canonicalization that fixes the drift, why the first attempt (round-to-digits) was wrong, and what the current snap-to-rational approach actually does.

## Why drift breaks group closure

A space-group operation is `r ↦ R·r + τ` in fractional coordinates. Composing two operations `(R₁, τ₁) ∘ (R₂, τ₂)` gives `(R₁·R₂, R₁·τ₂ + τ₁)`. The rotation part `R₁·R₂` is exact (integer matrices), but the translation part involves a Float64 matrix-vector product and accumulates drift on the order of machine epsilon (`~10⁻¹⁶`) per operation.

Spacey verifies group closure during `spacegroup` — it computes `op₁ * op₂` for every pair and checks the result is in the group. With Float64 drift, `(op₁ * op₂).τ` differs from the matching group member's `τ` by `~10⁻¹⁶ × O(10) ≈ 10⁻¹⁵`, which makes Julia's default `==` (exact field equality) fail. Closure tests then fail too. The fix is to canonicalize `τ` at construction so that drifted-but-equivalent values compare equal.

## The first attempt: round to 10 digits

The original `SpacegroupOp` constructor canonicalized with

```julia
SpacegroupOp(R, τ) = new(R, mod.(round.(τ, digits=10), 1.0))
```

The intent: fold `τ` into `[0, 1)` via `mod.(τ, 1.0)`, then `round(…, digits=10)` clamps drift below the `10⁻¹⁰` scale. This worked for `τ` values like `0`, `½`, `¼` — round-trip-stable as 10-digit decimals.

It silently broke for thirds. The Float64 representation of `1/3` rounded to ten digits is `0.3333333333` (ten threes). The Float64 representation of `2/3` rounded to ten digits is `0.6666666667` (rounded *up*: the eleventh digit was a `6` that pushed the tenth up to `7`).

But composing two ops with `τ_a = τ_b = 1/3` gives `τ_result = 0.3333333333 + 0.3333333333 = 0.6666666666` (ten sixes, no carry). After rounding, this stays at `0.6666666666` — which doesn't match the canonical 10-digit form of `2/3` (`0.6666666667`). So `(a * b).τ` wasn't `==` the existing op with `τ = 2/3`. Closure failed — silently, only on the four trigonal-screw structures in the AFLOW corpus (γ-Se, α-Quartz, Cinnabar, CrCl₃, all P3₁21). Cubic, tetragonal, orthorhombic, and hexagonal corpus prototypes have no `⅓` τ-components, so the bug never surfaced there.

## The fix: snap to small rationals

Replace the round-to-digits with a snap-to-rational that prefers exact rational values over decimal approximations:

```julia
function _canonicalize_τ(τ; tol=1e-6)
    out = similar(τ, Float64)
    for i in eachindex(τ)
        x = mod(Float64(τ[i]), 1.0)
        snapped = x
        for q in 1:12
            p = round(Int, q * x)
            cand = (p == q) ? 0.0 : p / q
            if abs(x - cand) < tol
                snapped = cand
                break
            end
        end
        out[i] = snapped
    end
    return out
end
```

For each `τ` component, scan denominators `q = 1, 2, …, 12`, find the closest fraction `p/q`, and if it's within `tol = 1e-6`, snap to that fraction. The rationals `1/2`, `1/3`, `2/3`, `1/4`, `1/6`, `1/12`, `5/12` all live in this set — and crucially, both `1/3` and `2/3` snap to *exactly* their Float64-best representations, so `(1/3) + (1/3)` snapped equals `2/3` snapped.

The bug that started this whole episode is fixed not by rounding finer (which would just push the inconsistency to the next decimal) but by recognizing that `τ` lives in the rationals, not the reals. The Float64 representations are *approximations* of rationals; canonicalizing means snapping each approximation back to the rational it's nearest to.

## Why `q ≤ 12` and `tol = 1e-6`

The International Tables for Crystallography (ITA) `τ`-values for all 230 space groups have denominators in `{1, 2, 3, 4, 6, 12}`. The largest denominator (12) appears in the `c/12` shifts of `P6₅` and similar high-order screw axes. Going to `q ≤ 24` would catch nothing more in the AFLOW corpus and would risk falsely snapping a generic Float64 value to a nearby rational by accident. Stopping at `q ≤ 12` covers all real cases without inviting accidents.

The tolerance `1e-6` is generous enough for any drift accumulated through `O(20)` matrix-vector products at Float64 precision (`machine eps × 10⁵ ≈ 2 × 10⁻¹¹` is the typical accumulated drift), but tight enough to never confuse a genuinely-non-rational `τ` with a rational. The latter case shouldn't happen for legitimate space groups, but it's a useful guarantee against pathological inputs.

`tol` is hardcoded inside `_canonicalize_τ` and not exposed as a kwarg — the right answer is a single value, and exposing it would only invite mis-use.

## What this gives the user

A `SpacegroupOp` constructed with any drifted-rational `τ` (including the result of composing other ops) gets stored with an exact rational. Two ops that differ by a lattice translation or by accumulated Float64 drift now hash and compare equal as expected. Group closure, `Set{SpacegroupOp}`, `Dict{SpacegroupOp,_}` all work without surprises.

The user-facing consequence is invisible — that's the goal. The canonicalization happens in the constructor; the user doesn't pass anything, doesn't see anything, doesn't have to think about it. The previous round-to-digits version *also* hid its existence, which is what made the trigonal-screw bug so subtle. The fix is to keep the silent canonicalization but make it actually correct.

## See also

- Reference: [`SpacegroupOp`](../reference/space-groups.md), [`spacegroup`](../reference/space-groups.md)
- How-to: [Compose and apply operations](../how-to/compose-and-apply-ops.md)
- Explanation: [Validation strategy](validation-strategy.md) — the AFLOW corpus run that surfaced the bug
