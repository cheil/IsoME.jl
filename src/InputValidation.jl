#
# Input validation for the `arguments` struct
#
# The fields of `arguments` are concretely typed (type stability), which means a
# wrong input type aborts with Julia's generic `convert` MethodError - printing the
# full ~70-field constructor signature and never naming the offending field.
# `@checked_kwdef` replaces `Base.@kwdef` and installs a keyword constructor that
# reports *all* offending fields by name, with a hint on how to fix them.
#


"""
    _Required

Placeholder for fields declared without a default value, i.e. keyword arguments the
user has to provide.
"""
struct _Required end
const _required = _Required()


"""
    @checked_kwdef mutable struct ... end

Drop-in replacement for `Base.@kwdef` that generates a keyword constructor with
user-friendly error messages (unknown keyword names, missing required keywords and
wrong field types are all reported by field name, see [`_construct`](@ref)).

Defaults are re-evaluated on every call, so mutable defaults (e.g. `[-1.0]`) are not
shared between instances.

Like `Base.@kwdef`, a docstring written above the `@checked_kwdef` line documents the
struct: the generated struct definition is marked with `Base.@__doc__`, which is what tells
the doc system where in the macro output the docstring belongs.
"""
macro checked_kwdef(expr)
    expr isa Expr && expr.head === :struct || error("@checked_kwdef needs a struct definition")
    ismut, sname, body = expr.args
    Tname = sname isa Symbol ? sname : sname.args[1]

    newbody = Expr(:block)
    fnames = Symbol[]
    fdefaults = Any[]
    for line in body.args
        if line isa Expr && line.head === :(=)          # field with default value
            decl, default = line.args
            push!(newbody.args, decl)
            push!(fnames, decl isa Expr ? decl.args[1] : decl)
            push!(fdefaults, default)
        elseif line isa Expr && line.head === :(::)     # field without default value
            push!(newbody.args, line)
            push!(fnames, line.args[1])
            push!(fdefaults, _required)     # marker: keyword argument has to be supplied
        else                                            # LineNumberNode, comments, ...
            push!(newbody.args, line)
        end
    end

    defaults = Expr(:tuple, (Expr(:(=), n, d) for (n, d) in zip(fnames, fdefaults))...)

    structdef = Expr(:struct, ismut, sname, newbody)

    # `Base.@__doc__` marks the struct as the piece of this expansion a docstring above the
    # macro call attaches to; without it `@doc` sees only the generated block and refuses.
    return esc(quote
        Base.@__doc__ $structdef

        # defaults as a NamedTuple, freshly evaluated on every call
        _fielddefaults(::Type{$Tname}) = (; $defaults...)

        # the keyword arguments are funnelled into a type-independent container, so that
        # `_construct` is compiled only once instead of once per set of keyword names
        $Tname(; kwargs...) = _construct($Tname, Pair{Symbol,Any}[k => v for (k, v) in kwargs])
    end)
end


"""
    _removed

Inputs of an earlier IsoME version that no longer exist, each with the line telling the
user what replaces it. Consulted by [`_unknown_msg`](@ref), so that a script written
against an older version fails with the migration step rather than with a bare "not a
field". Add a row here whenever an input is removed or renamed.
"""
const _removed = Dict{Symbol,String}(
    :shiftcut   => "removed in IsoME 2.0: `encut` is now the single ε-cutoff and bounds χ and Nₑ as well. Drop `shiftcut` and set `encut` to the value you used for it (the default moved from 5000 to 2000 meV accordingly)."
)


"""
    _legacySentinel

Fields whose "not set" sentinel changed from `-1` to `NaN` in IsoME 2.0. An explicit `-1`
used to mean "infer this", and would now be taken at face value as a μ*, a mixing factor or
an energy - a change of results with no error anywhere. These are rejected instead, see
[`_check_legacy_sentinels`](@ref).
"""
const _legacySentinel = (:mu, :muc_AD, :muc_ME, :ef, :efW, :typEl, :mixing_beta)


"""
    _nonNegative

Inputs that have no meaning below zero, each with the name used in the error message. A
negative value here is always a typo or a leftover sentinel, and it would otherwise travel
silently into the equations - a negative `muc_ME` flips the sign of the Coulomb term, a
negative `N_it` skips the iteration loop entirely.

`NaN` is not caught by this (every comparison with `NaN` is false), which is what keeps the
"infer this during the run" sentinel of `mu`, `muc_AD`, `muc_ME` and `mixing_beta` working.
"""
const _nonNegative = (
    mu          = "the Coulomb strength μ = N(ε_F)·W(ε_F,ε_F)",
    muc_AD      = "the pseudopotential μ*_AD",
    muc_ME      = "the pseudopotential μ*_ME",
    mixing_beta = "the linear mixing factor",
    conv_thr    = "the convergence threshold",
    minGap      = "the gap threshold",
    N_it        = "the maximum number of iterations",
    min_it      = "the minimum number of iterations",
)


"""
    _check_negative(fields)

Reject a negative value on any input listed in [`_nonNegative`](@ref). `fields` is an
iterable of `name => value` pairs, so that the same check serves the keyword constructor
and [`checkInput!`](@ref), which re-runs it on the struct at the start of every solve to
catch a field assigned after construction.

Deliberately *not* hooked into `setproperty!`: `calcMucME` and `calcMucAD` assign a μ* that
may come out negative and clamp it on the next line, which is a valid recovery, not a bad
input.
"""
function _check_negative(fields)
    bad = Pair{Symbol,Any}[]
    for (k, v) in fields
        haskey(_nonNegative, k) || continue
        v isa Number && !(v isa Bool) && v < 0 && push!(bad, k => v)
    end
    isempty(bad) || _throw_negative(bad)
    return nothing
end


@noinline function _throw_negative(bad::Vector{Pair{Symbol,Any}})
    io = IOBuffer()
    println(io, "invalid input to `arguments`:\n")
    for (k, v) in bad
        println(io, "  * ", k, " = ", v, "  - ", _nonNegative[k], " cannot be negative")
    end
    if any(first(p) in _legacySentinel for p in bad)
        println(io, "\nTo have a value inferred during the run, leave the field out or pass NaN.")
    end
    throw(ArgumentError(String(take!(io))))
end


"""
    _check_legacy_sentinels(kw)

Reject an explicit `-1` on a field that used it as the "not set" sentinel before IsoME 2.0.
"""
function _check_legacy_sentinels(kw::Vector{Pair{Symbol,Any}})
    hit = Symbol[]
    for (k, v) in kw
        k in _legacySentinel || continue
        v isa Number && !(v isa Bool) && v == -1 && push!(hit, k)
    end
    isempty(hit) || _throw_legacy(hit)
    return nothing
end


@noinline function _throw_legacy(hit::Vector{Symbol})
    io = IOBuffer()
    println(io, "outdated input value:\n")
    for k in hit
        println(io, "  * ", k, " = -1")
    end
    println(io, "\nBefore IsoME 2.0, -1 meant \"not set, infer this during the run\". The sentinel")
    println(io, "is now NaN, so a -1 would be used as an actual value here.")
    println(io, "\nTo have the value inferred, leave the field out (or pass NaN). If you really")
    println(io, "mean -1, pass -1.0000001 or set the field after construction.")
    throw(ArgumentError(String(take!(io))))
end


"""
    _construct(T, kw)

Build a `T` from the user supplied keyword arguments `kw`, collecting *all* problems
(unknown keywords, missing required keywords, wrong types) before throwing a single
`ArgumentError` that names every offending field.
"""
function _construct(::Type{T}, kw::Vector{Pair{Symbol,Any}}) where {T}
    names = fieldnames(T)

    # unknown keyword arguments (e.g. fields renamed/removed since an older IsoME version)
    unknown = Symbol[k for (k, _) in kw if !(k in names)]
    isempty(unknown) || _throw_unknown(T, unknown)

    # inputs that are still fields but whose sentinel value changed
    _check_legacy_sentinels(kw)

    # values outside the range an input can have (checked before the struct is built, so
    # that the error points at the `arguments(...)` call rather than at the solve)
    _check_negative(kw)

    defaults = _fielddefaults(T)
    vals = Vector{Any}(undef, length(names))
    missing_fields = Symbol[]
    wrong_types = Pair{Symbol,Any}[]

    for (i, name) in enumerate(names)
        idx = findfirst(p -> first(p) === name, kw)
        val = isnothing(idx) ? defaults[name] : last(kw[idx])
        Ft = fieldtype(T, name)

        if val === _required
            push!(missing_fields, name)
            vals[i] = nothing
        elseif val isa Ft
            vals[i] = val
        else
            try
                vals[i] = convert(Ft, val)      # e.g. Vector{Int} -> Vector{Float64}
            catch ex
                ex isa InterruptException && rethrow(ex)
                push!(wrong_types, name => val)
                vals[i] = nothing
            end
        end
    end

    if !isempty(missing_fields) || !isempty(wrong_types)
        _throw_invalid(T, missing_fields, wrong_types)
    end

    return T(vals...)
end


# Error paths are kept out of line and non-specializing: the message formatting below is
# never inferred as part of a caller, which keeps the compile time of `_construct` and
# `setproperty!` low (the messages are built only when an input is actually wrong).
@noinline function _throw_unknown(@nospecialize(T), unknown::Vector{Symbol})
    throw(ArgumentError(_unknown_msg(T, unknown)))
end

@noinline function _throw_invalid(@nospecialize(T), missing_fields::Vector{Symbol},
                                  wrong_types::Vector{Pair{Symbol,Any}})
    throw(ArgumentError(_input_msg(T, missing_fields, wrong_types)))
end


"""
    _input_msg(T, missing_fields, wrong_types)

Assemble the error message for missing keywords and wrong field types.
"""
function _input_msg(@nospecialize(T), missing_fields::Vector{Symbol},
                    wrong_types::Vector{Pair{Symbol,Any}})
    io = IOBuffer()
    println(io, "invalid input to `", nameof(T), "`:\n")

    if !isempty(missing_fields)
        for name in missing_fields
            println(io, "  * ", name, " is required but was not given",
                    "  (expected ::", fieldtype(T, name), ")")
        end
        isempty(wrong_types) || println(io)
    end

    for (name, val) in wrong_types
        println(io, _type_msg(name, fieldtype(T, name), val))
    end

    println(io, "\nAll other fields are fine. Note that some field types have been")
    println(io, "restricted in recent IsoME versions - see `?", nameof(T),
                "` for the type of every field.")
    return String(take!(io))
end


"""
    _type_msg(name, Ft, val)

Describe a single field whose value could not be converted to the declared type,
including a hint on how to fix the most common cases.
"""
function _type_msg(name::Symbol, @nospecialize(Ft), @nospecialize(val))
    io = IOBuffer()
    println(io, "  * ", name, " :: ", Ft, "  - got ", typeof(val), ": ", repr(val))
    hint = _hint(name, Ft, val)
    isempty(hint) || print(io, "      hint: ", hint)
    return String(take!(io))
end


"""
    _hint(name, Ft, val)

Suggestion on how to fix a wrongly typed input, "" if there is nothing sensible to say.
"""
function _hint(name::Symbol, @nospecialize(Ft), @nospecialize(val))
    if Ft <: AbstractVector && val isa Number
        return "pass a vector, e.g. $name = [$(val)]"
    elseif Ft <: AbstractVector && val isa AbstractVector
        return "element type mismatch: got $(eltype(val)) elements, need $(eltype(Ft))"
    elseif Ft <: Number && val isa AbstractVector
        return "pass a single number instead of a vector"
    elseif Ft === Bool
        return "pass true or false"
    elseif Ft <: Integer && val isa AbstractFloat
        return isfinite(val) ? "pass an integer, e.g. $name = $(round(Int, val))" :
                               "pass a finite integer"
    elseif Ft <: AbstractString
        return "pass a string in double quotes, e.g. $name = \"$(val)\""
    elseif Ft <: Number && val isa AbstractString
        return "pass a number without quotes, e.g. $name = $(val)"
    end
    return ""
end


"""
    _unknown_msg(T, unknown)

Error message for keyword arguments that are not fields of `T`, with a "did you mean"
suggestion for fields that were renamed.
"""
function _unknown_msg(@nospecialize(T), unknown)
    names = collect(fieldnames(T))
    io = IOBuffer()
    println(io, "unknown input to `", nameof(T), "`:\n")
    for k in unknown
        # a field that was removed in a past version gets its migration line instead of a
        # guess: the replacement is rarely the nearest name
        if haskey(_removed, k)
            println(io, "  * ", k, " - ", _removed[k])
            continue
        end
        print(io, "  * ", k, " is not a field of `", nameof(T), "`")
        near = _closest(k, names)
        isempty(near) ? println(io) : println(io, " - did you mean ", join(near, ", "), "?")
    end
    println(io, "\nUse `?", nameof(T), "` for the list of all available fields.")
    return String(take!(io))
end


"""
    _closest(name, names)

Field names similar to `name`, best match first. A field containing `name` (or vice
versa) is always a candidate, otherwise a small edit distance is required - the
tolerance scales with the length of `name` so that short names do not match everything.
"""
function _closest(name::Symbol, names::Vector{Symbol})
    s = lowercase(String(name))
    maxdist = clamp(cld(length(s), 3), 1, 3)

    scored = Tuple{Symbol,Int}[]
    for n in names
        t = lowercase(String(n))
        if length(s) >= 3 && (occursin(s, t) || occursin(t, s))
            push!(scored, (n, 0))                       # substring match: best candidate
        else
            d = _editdist(s, t)
            d <= maxdist && push!(scored, (n, d))
        end
    end

    sort!(scored, by = x -> x[2])
    return first.(scored[1:min(3, end)])
end


"""
    _editdist(a, b)

Levenshtein distance between two strings.
"""
function _editdist(a::AbstractString, b::AbstractString)
    m, n = length(a), length(b)
    prev = collect(0:n)
    curr = similar(prev)
    for (i, ca) in enumerate(a)
        curr[1] = i
        for (j, cb) in enumerate(b)
            curr[j+1] = min(prev[j+1] + 1, curr[j] + 1, prev[j] + (ca == cb ? 0 : 1))
        end
        prev, curr = curr, prev
    end
    return prev[n+1]
end
