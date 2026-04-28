# Tutorials

Step-by-step guided walkthroughs that produce a concrete result. Read these top-to-bottom on first use; they assume nothing about Spacey but assume you can read Julia code at a basic level.

- [Find the point group of a cubic lattice](01-first-pointgroup.md) — ~5 minutes. Build a lattice by hand, call `pointGroup`, see why a simple cubic lattice has 48 symmetries and a tetragonal one has 16.
- [Find the space group of NaCl (and three other crystals)](02-first-spacegroup.md) — ~15 minutes. Build a `Crystal` (lattice + atoms + types), call `spacegroup`, walk through four common materials: rocksalt, zincblende, an ordered alloy, and a ferroelectric perovskite.

If you already know what you want to do, the [How-to guides](../how-to/index.md) are a better fit. If you want background on *why* Spacey decides what it decides, read an [Explanation](../explanation/index.md).
