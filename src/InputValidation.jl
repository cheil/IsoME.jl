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

    return esc(quote
        $(Expr(:struct, ismut, sname, newbody))

        # defaults as a NamedTuple, freshly evaluated on every call
        _fielddefaults(::Type{$Tname}) = (; $defaults...)

        # the keyword arguments are funnelled into a type-independent container, so that
        # `_construct` is compiled only once instead of once per set of keyword names
        $Tname(; kwargs...) = _construct($Tname, Pair{Symbol,Any}[k => v for (k, v) in kwargs])
    end)
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
            catch
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
