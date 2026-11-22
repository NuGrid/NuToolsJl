# ward2025.jl: the comparison values of Ward et al. (2025, ApJ 986, 109) for CO and ONe novae.

"The isotopic ratios of Ward et al. (2025) Tables 3–5."
const WARD_PAIRS = [("C12", "C13"), ("N14", "N15"), ("O16", "O17"), ("O16", "O18"), ("Mg24", "Mg25"),
                    ("Mg24", "Mg26"), ("Al26", "Al27"), ("Si28", "Si29"), ("Si28", "Si30"), ("S32", "S33"),
                    ("S32", "S34")]

"Dilution ξ of Ward et al. (2025) Eq. 4: 9 (CO) and 17 (ONe)."
const WARD_XI = Dict(:co => 9.0, :one => 17.0)

"""
Ward et al. (2025) Table 3 (`:co`) and Table 4 (`:one`):
`pair => (Xa_burn, Xb_burn, (Xa/Xb)_burn, factor1, (Xa/Xb)_mix, factor2)`.
"""
const WARD2025 = Dict(
    :co => Dict(("C12", "C13") => (9.87e-2, 4.71e-2, 2.10, 1.13, 26.1, 1.11),
                ("N14", "N15") => (8.45e-2, 3.57e-3, 23.6, 1.21, 25.2, 1.21),
                ("O16", "O17") => (3.65e-4, 8.55e-3, 4.30e-2, 1.25, 136.0, 1.08),
                ("O16", "O18") => (3.65e-4, 7.31e-6, 50.0, 1.67, 1.03e4, 1.03),
                ("Mg24", "Mg25") => (1.93e-7, 4.31e-5, 4.47e-3, 1.14, 6.93, 1.02),
                ("Mg24", "Mg26") => (1.93e-7, 1.67e-6, 0.113, 1.94, 6.60, 1.00),
                ("Al26", "Al27") => (1.30e-5, 7.91e-5, 0.164, 1.16, 2.63e-2, 1.23),
                ("Si28", "Si29") => (1.27e-3, 1.49e-5, 85.8, 1.17, 22.7, 1.01),
                ("Si28", "Si30") => (1.27e-3, 2.83e-4, 4.54, 3.46, 13.4, 1.88),
                ("S32", "S33") => (5.44e-4, 1.31e-6, 401.0, 5.77, 139.0, 1.09),
                ("S32", "S34") => (5.44e-4, 3.71e-6, 140.0, 13.0, 24.7, 1.05)),
    :one => Dict(("C12", "C13") => (3.99e-3, 1.33e-3, 2.99, 1.25, 46.0, 1.16),
                 ("N14", "N15") => (1.17e-1, 1.56e-5, 7.47e3, 1.10, 3.04e3, 1.04),
                 ("O16", "O17") => (7.84e-5, 5.38e-7, 145.0, 1.11, 6.39e4, 1.01),
                 ("O16", "O18") => (7.84e-5, 4.69e-10, 1.62e5, 1.43, 1.13e4, 1.00),
                 ("Mg24", "Mg25") => (3.70e-7, 1.36e-4, 2.71e-3, 1.17, 3.51, 1.00),
                 ("Mg24", "Mg26") => (3.70e-7, 5.33e-6, 7.04e-2, 1.45, 5.56, 1.00),
                 ("Al26", "Al27") => (3.68e-5, 3.18e-4, 0.116, 1.20, 7.82e-4, 1.20),
                 ("Si28", "Si29") => (1.51e-2, 1.46e-4, 104.0, 1.18, 38.6, 1.07),
                 ("Si28", "Si30") => (1.51e-2, 5.46e-3, 2.77, 7.72, 4.19, 6.20),
                 ("S32", "S33") => (1.01e-1, 4.14e-4, 211.0, 11.8, 203.0, 6.77),
                 ("S32", "S34") => (1.01e-1, 2.91e-4, 307.0, 9.01, 172.0, 3.05)))

"Ward et al. (2025) Table 5: the reaction(s) with the largest MI for each ratio, `pair => \"reaction MI\"`."
const WARD2025_TOP = Dict(
    :co => Dict(("C12", "C13") => "13N(p,g)14O 0.71", ("N14", "N15") => "15N(p,a)12C 1.4",
                ("O16", "O17") => "16O(p,g)17F 0.37", ("O16", "O18") => "18F(p,a)15O 0.78",
                ("Mg24", "Mg25") => "24Mg(p,g)25Al 0.67", ("Mg24", "Mg26") => "26mAl(p,g)27Si 1.2",
                ("Al26", "Al27") => "26gAl(p,g)27Si 0.85", ("Si28", "Si29") => "29Si(p,g)30P 0.43; 28Si(p,g)29P 0.29",
                ("Si28", "Si30") => "30P(p,g)31S 1.6", ("S32", "S33") => "33S(p,g)34Cl 1.5",
                ("S32", "S34") => "34S(p,g)35Cl 1.4"),
    :one => Dict(("C12", "C13") => "13C(p,g)14N 0.78", ("N14", "N15") => "14N(p,g)15O 1.5",
                 ("O16", "O17") => "17O(p,a)14N 0.95", ("O16", "O18") => "18F(p,a)15O 1.0",
                 ("Mg24", "Mg25") => "24Mg(p,g)25Al 0.56", ("Mg24", "Mg26") => "26mAl(p,g)27Si 0.66",
                 ("Al26", "Al27") => "26gAl(p,g)27Si 0.86", ("Si28", "Si29") => "29Si(p,g)30P 0.48; 28Si(p,g)29P 0.25",
                 ("Si28", "Si30") => "30P(p,g)31S 1.9", ("S32", "S33") => "33S(p,g)34Cl 1.7",
                 ("S32", "S34") => "34S(p,g)35Cl 1.2"))

"""
    ward_report(ens; model=:co, iniab, xi=WARD_XI[model], kwargs...)

[`sensitivity_report`](@ref) on Ward et al. (2025)'s ratios with their Table 3 (`:co`) or Table 4
(`:one`) values and Table 5 top reactions alongside. `iniab` is the initial-abundance file of the runs
(for the mixed ratio). Other keywords (`isotopes`, `top`, `candidates`, `outdir`, `io`) pass through.
"""
function ward_report(ens::MCEnsemble; model::Symbol=:co, iniab::AbstractString, xi::Real=WARD_XI[model], kwargs...)
    haskey(WARD2025, model) || throw(ArgumentError("model must be :co or :one"))
    sensitivity_report(ens; pairs=WARD_PAIRS, X0=read_iniab(iniab), xi, reference=WARD2025[model],
                       reference_top=WARD2025_TOP[model], title="Ward $(uppercase(string(model))) model", kwargs...)
end

_run_dir(d) = isfile(joinpath(d, "x-time.dat")) ? d : joinpath(d, "ppn")

"""
    ward_baseline(runs; model=:co, io=stdout) -> Dict(run => Dict(pair => (burn, mix)))

Ward et al. (2025) Table 3/4 quantities of single finished ppn runs (no Monte Carlo): for each run
directory (one holding x-time.dat, or a case directory with a `ppn/` subdirectory), the burned ratios of
its last x-time.dat row and the mixed ratios using the run's own `ini_filename`, printed next to Ward's,
followed by the isotopes' mass fractions. 26Al includes the isomer.
"""
function ward_baseline(runs::AbstractVector{<:AbstractString}; model::Symbol=:co, io::IO=stdout)
    W, xi = WARD2025[model], WARD_XI[model]
    res = Dict{String,Any}()
    X = Dict{String,Dict{String,Float64}}()
    for c in runs
        d = _run_dir(c)
        x, x0 = final_abundances(d), run_iniab(d)
        x["Al26"] = get(x, "Al26", 0.0) + get(x, "Al26m", 0.0)
        X[c] = x
        res[c] = Dict(p => (get(x, p[1], 0.0) / get(x, p[2], NaN),
                            (get(x, p[1], 0.0) + xi * get(x0, p[1], 0.0)) / (get(x, p[2], 0.0) + xi * get(x0, p[2], 0.0)))
                      for p in WARD_PAIRS)
    end
    names = [basename(normpath(c)) for c in runs]
    @printf(io, "%-10s %10s %s %10s %s\n", "ratio", "Ward burn", join([@sprintf("%13s", n) for n in names], ' '),
            "Ward mix", join([@sprintf("%13s", n) for n in names], ' '))
    for p in WARD_PAIRS
        @printf(io, "%-10s %10.3g %s %10.3g %s\n", "$(p[1])/$(p[2])", W[p][3],
                join([@sprintf("%13.3g", res[c][p][1]) for c in runs], ' '), W[p][5],
                join([@sprintf("%13.3g", res[c][p][2]) for c in runs], ' '))
    end
    wx = Dict{String,Float64}()
    for (p, w) in W
        wx[p[1]], wx[p[2]] = w[1], w[2]
    end
    isos = sort(unique(vcat(first.(WARD_PAIRS), last.(WARD_PAIRS))); by=s -> parse(Int, match(r"\d+", s).match))
    @printf(io, "\n%-8s %10s %s\n", "isotope", "Ward X", join([@sprintf("%13s", n) for n in names], ' '))
    for i in vcat(isos, ["H1", "He4", "Na22", "Ne20", "Na23", "F19"])
        @printf(io, "%-8s %10s %s\n", i, haskey(wx, i) ? @sprintf("%10.3g", wx[i]) : "",
                join([@sprintf("%13.3g", get(X[c], i, 0.0)) for c in runs], ' '))
    end
    res
end
