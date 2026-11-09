# The Julia side of PracticalBayesianRegressionR: small helpers the R package
# calls. Models themselves are described to Julia as source text (see
# R/julia.R), so everything here is about moving data and draws across.
module PBRR

import PracticalBayesianRegression as PBR
const FlexiChains = PBR.FlexiChains    # not a direct dependency of the session's project

# A data set R wrote with write_data_csv (exact: numbers as %.17g), with the
# column types R had: "double", "integer", "logical" or "string". Factor
# columns arrive as strings; their levels are then sorted, as R's factor()
# default does.
function read_table(path::AbstractString, names, kinds)
    names = _strings(names); kinds = _strings(kinds)
    t = PBR.read_csv(path)
    cols = Any[]
    for (nm, k) in zip(names, kinds)
        j = PBR._colindex(t, nm)
        miss = t.miss[j]
        v = if k == "double"
            _with_missing(t.num[j], miss)
        elseif k == "integer"
            _with_missing([m ? 0 : round(Int, x) for (x, m) in zip(t.num[j], miss)], miss)
        elseif k == "logical"
            _with_missing([r == "TRUE" for r in t.raw[j]], miss)
        else
            _with_missing(copy(t.raw[j]), miss)
        end
        push!(cols, v)
    end
    return NamedTuple{Tuple(Symbol.(names))}(Tuple(cols))
end
_with_missing(v::Vector{T}, miss) where {T} =
    any(miss) ? Union{Missing,T}[m ? missing : x for (x, m) in zip(v, miss)] : v

# JuliaCall turns a length-one character vector into a String
_strings(x::AbstractString) = [String(x)]
_strings(x) = String[string(s) for s in x]

# ---------------------------------------------------------------- draws out

variable_names(fit) = String[n for n in PBR.variables(fit)]

# iterations x chains x variables, constrained, in `names` order
function draws_array(fit, names)
    names = _strings(names)
    niter, nch = size(fit.draws, 2), size(fit.draws, 3)
    out = Array{Float64,3}(undef, niter, nch, length(names))
    for (k, n) in enumerate(names)
        out[:, :, k] .= Matrix{Float64}(fit.chain[FlexiChains.Parameter(n)])
    end
    return out
end

stat_names(fit) = sort!([string(k.name) for k in keys(fit.chain) if k isa FlexiChains.Extra])

function stats_array(fit)
    names = stat_names(fit)
    niter, nch = size(fit.draws, 2), size(fit.draws, 3)
    out = Array{Float64,3}(undef, niter, nch, length(names))
    for (k, n) in enumerate(names)
        out[:, :, k] .= Matrix{Float64}(fit.chain[FlexiChains.Extra(n)])
    end
    return out
end

# unconstrained draws, iterations x chains x parameters: with them a fit can
# be rebuilt in a later session (from_draws)
theta_array(fit) = permutedims(fit.draws, (2, 3, 1))

chain_times(fit) = Float64.(fit.sampling_time)
response(fit::PBR.BrmsFit) = Float64.(fit.density.y)

nobs(fit::PBR.BrmsFit) = fit.density.N
nobs(fit) = maximum(f -> f.density.N, fit.fits)

# ------------------------------------------------------------ draws back in

"""
    from_draws(formula, data; family, prior, theta, statnames, stats, times, iter, warmup, sampler)

Rebuild a BrmsFit from draws R holds (a compiled model's run, or a fit from
an earlier session): the front end rebuilds the model on the same data, and
the draws are attached unchanged.
"""
function from_draws(formula, data; family, prior, theta, statnames, stats, times, iter, warmup,
                    sampler="nuts")
    ld, info = PBR.build_model(formula, data; family=family, prior=prior)
    th = permutedims(Array{Float64,3}(theta), (3, 1, 2))
    sn = _strings(statnames)
    st = Dict{String,Matrix{Float64}}(sn[k] => Matrix{Float64}(stats[:, :, k]) for k in eachindex(sn))
    tm = times isa Number ? [Float64(times)] : Float64.(collect(times))
    fit = PBR.fit_from_draws(ld, info, th, st, tm; sampler=sampler, iter=Int(iter), warmup=Int(warmup))
    merge!(fit.call, Dict{Symbol,Any}(:formula => formula, :data => data, :family => family, :prior => prior,
                                      :iter => Int(iter), :warmup => Int(warmup), :chains => size(th, 3)))
    return fit
end

# ---------------------------------------------------------------- compiling

function compile(formula, data; family, prior, name, store, force, juliac, key, verbose)
    nm = name == "" ? nothing : String(name)
    jc = juliac == "" ? nothing : String(juliac)
    cm = PBR.compile_model(formula, data; family=family, prior=prior, name=nm, store=String(store),
                           force=force, juliac=jc, keys=[String(key)], verbose=verbose)
    return cm.name
end

# ------------------------------------------------------------- summaries

# get_prior's rows as tab-separated text: prior, class, coef, group, resp,
# dpar, lb, ub, source
function prior_rows(t)
    s(x::Float64) = isfinite(x) ? PBR._g(x) : ""
    return String[join((PBR._prior_string(r.prior), string(r.class), r.coef, r.group, r.resp, r.dpar, s(r.lb),
                        s(r.ub), string(r.source)), "\t") for r in t.rows]
end

# LOOResult as plain numbers: elpd, se, p, p_se, then pointwise elpd, pareto k
loo_numbers(r::PBR.LOOResult) = (Float64[r.elpd, r.se, r.p, r.p_se], r.pointwise, r.pareto_k)

end
