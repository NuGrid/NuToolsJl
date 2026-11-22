# stats.jl: the statistics of a Monte Carlo sensitivity study.

"""
    ranks(v) -> Vector{Float64}

Ranks of `v`, ties sharing their mean rank.
"""
function ranks(v::AbstractVector)
    p = sortperm(v)
    r = Vector{Float64}(undef, length(v))
    i = 1
    n = length(v)
    while i <= n
        j = i
        while j < n && v[p[j+1]] == v[p[i]]
            j += 1
        end
        r[p[i:j]] .= (i + j) / 2
        i = j + 1
    end
    r
end

_zscore(v) = (s = std(v; corrected=false); (v .- mean(v)) ./ (s > 0 ? s : 1.0))

"""
    spearman(x, y) -> Float64

Spearman rank correlation.
"""
spearman(x::AbstractVector, y::AbstractVector) = mean(_zscore(ranks(x)) .* _zscore(ranks(y)))

"""
    spearman_all(xr, y) -> Vector{Float64}

Spearman correlation of `y` with every row of `xr` (nrate × n), threaded over rows
(start Julia with `-t auto`).
"""
function spearman_all(xr::AbstractMatrix, y::AbstractVector)
    ry = _zscore(ranks(y))
    rho = zeros(size(xr, 1))
    Threads.@threads for j in 1:size(xr, 1)
        rho[j] = mean(_zscore(ranks(xr[j, :])) .* ry)
    end
    rho
end

"""
    ksg_mi(x, y; k=3, rng=MersenneTwister(0)) -> Float64

Mutual information I(x; y) in nats of two 1-D samples, by the Kraskov, Stögbauer & Grassberger (2004)
estimator 1 (max-norm, `k` neighbours) — the estimator of Ward et al. (2025). Both variables are
standardized and given a 1e-10 jitter to break ties; negative estimates are clipped to 0.
O(n √n)-ish: about 10 ms for n = 10⁴.
"""
function ksg_mi(x::AbstractVector, y::AbstractVector; k::Integer=3, rng=MersenneTwister(0))
    n = length(x)
    x = _zscore(x) .+ 1e-10 .* randn(rng, n)
    y = _zscore(y) .+ 1e-10 .* randn(rng, n)
    p = sortperm(x)
    xs, ys = x[p], y[p]
    eps = zeros(n)
    best = zeros(k)
    for i in 1:n                                      # k-th neighbour distance, scanning outward in sorted x
        fill!(best, Inf)
        for dir in (-1, 1)
            j = i + dir
            while 1 <= j <= n
                dx = abs(xs[j] - xs[i])
                dx >= best[k] && break
                d = max(dx, abs(ys[j] - ys[i]))
                if d < best[k]
                    best[k] = d
                    sort!(best)
                end
                j += dir
            end
        end
        eps[i] = best[k] * (1 - 1e-12)
    end
    sx, sy = sort(xs), sort(ys)
    H = cumsum([0.0; 1.0 ./ (1:n)])                   # digamma(m) = -γ + H[m] for integer m
    ψ(m) = H[m] - Base.MathConstants.eulergamma
    acc = 0.0
    for i in 1:n
        nx = searchsortedlast(sx, xs[i] + eps[i]) - searchsortedfirst(sx, xs[i] - eps[i])   # excludes the point itself
        ny = searchsortedlast(sy, ys[i] + eps[i]) - searchsortedfirst(sy, ys[i] - eps[i])
        acc += ψ(nx + 1) + ψ(ny + 1)
    end
    max(ψ(k) + ψ(n) - acc / n, 0.0)
end

"""
    spread(v) -> (median, p16, p84, factor)

Median, 16th and 84th percentiles of the finite values of `v`, and the uncertainty factor
√(p84/p16) (Ward et al. 2025 "factor1"/"factor2"; 1σ for a lognormal).
"""
function spread(v::AbstractVector)
    f = filter(isfinite, v)
    p16, p50, p84 = quantile(f, 0.16), quantile(f, 0.50), quantile(f, 0.84)
    (median=p50, p16=p16, p84=p84, factor=sqrt(p84 / p16))
end
