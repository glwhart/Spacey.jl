# Generates test/aflow_structures.jl from the AFLOW Part 1 paper appendix
# (Mehl et al., Comput. Mater. Sci. 136 (2017) S1–S828).
#
# Assumes the following pdftotext outputs exist:
#   /tmp/aflow_left.txt   (left column of pp. S638–S828)
#   /tmp/aflow_right.txt  (right column of same)
#
# They can be regenerated with:
#   pdftotext -layout -x 0 -y 0 -W 300 -H 800 -f 638 -l 828 \
#             papers_cst.pdf /tmp/aflow_left.txt
#   pdftotext -layout -x 297 -y 0 -W 300 -H 800 -f 638 -l 828 \
#             papers_cst.pdf /tmp/aflow_right.txt
#
# Usage (from the repo root):
#   julia tools/generate_aflow_tests.jl
#
# This writes test/aflow_structures.jl which is checked into the repo and
# consumed by the AFLOW testset in test/runtests.jl.

using LinearAlgebra

# ────────────────────────────────────────────────────────────────────────
# Point-group order of space group n (the number Spacey should return on a
# primitive-cell input — AFLOW POSCARs are primitive, so centering
# translations are absorbed and we see just the point-group operations).
function point_group_order(sg::Int)
    1 ≤ sg ≤ 230 || error("invalid space group: $sg")
    sg == 1           && return 1
    sg == 2           && return 2
    3 ≤ sg ≤ 9        && return 2     # monoclinic 2, m
    10 ≤ sg ≤ 15      && return 4     # monoclinic 2/m
    16 ≤ sg ≤ 24      && return 4     # orthorhombic 222
    25 ≤ sg ≤ 46      && return 4     # orthorhombic mm2
    47 ≤ sg ≤ 74      && return 8     # orthorhombic mmm
    75 ≤ sg ≤ 82      && return 4     # tetragonal 4, -4
    83 ≤ sg ≤ 88      && return 8     # tetragonal 4/m
    89 ≤ sg ≤ 98      && return 8     # tetragonal 422
    99 ≤ sg ≤ 110     && return 8     # tetragonal 4mm
    111 ≤ sg ≤ 122    && return 8     # tetragonal -42m
    123 ≤ sg ≤ 142    && return 16    # tetragonal 4/mmm
    143 ≤ sg ≤ 146    && return 3     # trigonal 3
    147 ≤ sg ≤ 148    && return 6     # trigonal -3
    149 ≤ sg ≤ 155    && return 6     # trigonal 32
    156 ≤ sg ≤ 161    && return 6     # trigonal 3m
    162 ≤ sg ≤ 167    && return 12    # trigonal -3m
    168 ≤ sg ≤ 173    && return 6     # hexagonal 6
    sg == 174         && return 6     # hexagonal -6
    175 ≤ sg ≤ 176    && return 12    # hexagonal 6/m
    177 ≤ sg ≤ 182    && return 12    # hexagonal 622
    183 ≤ sg ≤ 186    && return 12    # hexagonal 6mm
    187 ≤ sg ≤ 190    && return 12    # hexagonal -6m2
    191 ≤ sg ≤ 194    && return 24    # hexagonal 6/mmm
    195 ≤ sg ≤ 199    && return 12    # cubic 23
    200 ≤ sg ≤ 206    && return 24    # cubic m-3
    207 ≤ sg ≤ 214    && return 24    # cubic 432
    215 ≤ sg ≤ 220    && return 24    # cubic -43m
    221 ≤ sg ≤ 230    && return 48    # cubic m-3m
end

# ────────────────────────────────────────────────────────────────────────
# Replace U+2212 MINUS SIGN with ASCII hyphen-minus so Julia's parse can
# read the number. (pdftotext -raw output preserves digits together but
# uses the unicode minus, which Julia's parse(Float64, ...) rejects.)
fix_broken_negatives(s::AbstractString) = replace(s, '−' => '-')

# ────────────────────────────────────────────────────────────────────────
# Prototype label → space group number.  "AB2_cP12_205_a_c" → 205.
function sg_from_prototype(label::AbstractString)
    parts = split(label, '_')
    length(parts) ≥ 3 || return nothing
    try
        return parse(Int, parts[3])
    catch
        return nothing
    end
end

# ────────────────────────────────────────────────────────────────────────
# Parse a single POSCAR block.  `block` is the lines immediately after the
# "X: Y - POSCAR" header line (and before the next header).  Returns
# nothing if the block can't be parsed.
function parse_poscar_block(block::Vector{String})
    # Find the scale-factor line.  It's usually "1.0000000000000000" after
    # some & / description lines from the multi-line header.
    i = 1
    while i ≤ length(block)
        s = strip(block[i])
        if occursin(r"^[0-9]+\.[0-9]+$", s)
            break
        end
        i += 1
    end
    i ≤ length(block) || return nothing
    scale = try
        parse(Float64, strip(block[i]))
    catch
        return nothing
    end
    i += 1

    # Next 3 non-empty lines: lattice vectors (may have split minus signs)
    lattice_rows = Matrix{Float64}(undef, 3, 3)
    for k in 1:3
        while i ≤ length(block) && isempty(strip(block[i]))
            i += 1
        end
        i ≤ length(block) || return nothing
        fixed = fix_broken_negatives(block[i])
        parts = split(strip(fixed))
        length(parts) ≥ 3 || return nothing
        for (j, p) in enumerate(parts[1:3])
            x = tryparse(Float64, p)
            x === nothing && return nothing
            lattice_rows[k, j] = x
        end
        i += 1
    end

    # Skip the element-names / count / Direct header (formats vary — some
    # P1 structures write "( 1a )" / "12" instead of "Fe S" / "4 8").
    # Just scan forward to "Direct" or "Cartesian".
    while i ≤ length(block)
        s = strip(block[i])
        if s == "Direct" || s == "Cartesian"
            break
        end
        i += 1
    end
    i ≤ length(block) || return nothing
    coord_mode = strip(block[i]) == "Direct" ? :fractional : :cartesian
    i += 1

    # Read atom position lines: "<x> <y> <z> <Element> [( <Wyckoff> )]".
    # Stop at the first line that doesn't parse as three floats.
    positions = Vector{Vector{Float64}}()
    types = Symbol[]
    while i ≤ length(block)
        raw = block[i]
        s = strip(fix_broken_negatives(raw))
        if isempty(s)
            i += 1; continue
        end
        parts = split(s)
        if length(parts) < 4
            break
        end
        xs = [tryparse(Float64, p) for p in parts[1:3]]
        if any(x -> x === nothing, xs)
            break
        end
        # Element is the first token from index 4+ that is letters only.
        elem = nothing
        for p in parts[4:end]
            if !isempty(p) && all(isletter, p)
                elem = p; break
            end
        end
        elem === nothing && break
        push!(positions, Float64[xs[1], xs[2], xs[3]])
        push!(types, Symbol(elem))
        i += 1
    end

    isempty(positions) && return nothing

    r = Matrix{Float64}(undef, 3, length(positions))
    for (j, p) in enumerate(positions)
        r[:, j] = p
    end

    A = collect(transpose(lattice_rows)) .* scale

    return (A=A, r=r, types=types, coords=coord_mode)
end

# ────────────────────────────────────────────────────────────────────────
# Walk a stream of lines, find every POSCAR header, and collect the block
# that follows until the next header (POSCAR or CIF).
function find_poscar_blocks(path::String)
    lines = readlines(path)
    blocks = NamedTuple[]
    i = 1
    header_re = r"^(.+?):\s*([A-Za-z0-9_]+)\s*-\s*POSCAR\s*$"
    stop_re   = r"-\s*(POSCAR|CIF)\s*$"
    while i ≤ length(lines)
        m = match(header_re, strip(lines[i]))
        if m === nothing
            i += 1; continue
        end
        name = strip(m.captures[1])
        proto = strip(m.captures[2])
        start = i + 1
        j = start
        while j ≤ length(lines)
            if match(stop_re, strip(lines[j])) !== nothing
                break
            end
            j += 1
        end
        push!(blocks, (name=name, prototype=proto, block=lines[start:j-1]))
        i = j
    end
    return blocks
end

# ────────────────────────────────────────────────────────────────────────
function main()
    blocks = find_poscar_blocks("/tmp/aflow_raw_clean.txt")
    println("Found $(length(blocks)) POSCAR headers")

    parsed = NamedTuple[]
    failed = String[]
    for b in blocks
        sg = sg_from_prototype(b.prototype)
        sg === nothing && (push!(failed, "sg-parse: $(b.name)"); continue)
        result = parse_poscar_block(b.block)
        result === nothing && (push!(failed, "block-parse: $(b.name) [$(b.prototype)]"); continue)
        order = point_group_order(sg)
        push!(parsed, (name=b.name, prototype=b.prototype, sg=sg, order=order,
                       A=result.A, r=result.r, types=result.types,
                       coords=result.coords))
    end

    println("Parsed: $(length(parsed)) / $(length(blocks))")
    println("Failed: $(length(failed))")
    if !isempty(failed)
        println("First failures:")
        for f in failed[1:min(10, length(failed))]
            println("  - $f")
        end
    end

    # Write out the data file
    outpath = joinpath(@__DIR__, "..", "test", "aflow_structures.jl")
    open(outpath, "w") do f
        println(f, "# Auto-generated by tools/generate_aflow_tests.jl from")
        println(f, "# AFLOW Part 1 (Mehl et al., Comput. Mater. Sci. 136 (2017) S1–S828).")
        println(f, "# Do not edit by hand — re-run the generator if the source paper or parser changes.")
        println(f, "#")
        println(f, "# Each entry: (name, prototype, sg, expected_order, A, r, types, coords)")
        println(f, "")
        println(f, "const AFLOW_PART1_STRUCTURES = [")
        for s in parsed
            println(f, "    (")
            println(f, "        name = ", repr(s.name), ",")
            println(f, "        prototype = ", repr(s.prototype), ",")
            println(f, "        sg = ", s.sg, ",")
            println(f, "        expected_order = ", s.order, ",")
            print(f,   "        A = ")
            show_matrix(f, s.A)
            println(f, ",")
            print(f,   "        r = ")
            show_matrix(f, s.r)
            println(f, ",")
            println(f, "        types = ", repr(s.types), ",")
            println(f, "        coords = ", repr(s.coords), ","),
            println(f, "    ),")
        end
        println(f, "]")
    end
    println("Wrote ", outpath)
end

function show_matrix(io::IO, M::Matrix{Float64})
    nrows, ncols = size(M)
    # Use explicit reshape so single-column matrices don't collapse into
    # vectors — [x; y; z] is a Vector in Julia, not a 3×1 Matrix.
    print(io, "reshape([")
    flat = reshape(M, :)
    for (i, v) in enumerate(flat)
        print(io, v)
        i < length(flat) && print(io, ", ")
    end
    print(io, "], ", nrows, ", ", ncols, ")")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
