using Test
using NuGridJl
using DataFrames
using LaTeXStrings
import CairoMakie as CM

const DATA = joinpath(@__DIR__, "data")
const NUPPN = joinpath(DATA, "nuppn_data")

@testset "NuGridJl" begin

@testset "elements" begin
    @test proton_number("He") == 2
    @test proton_number("he") == 2
    @test proton_number("Xx") === nothing
    @test element_symbol(6) == "C"
    @test element_symbol(0) == "n"

    @test Isotope("HE  4") == Isotope(2, 4, 0)
    @test Isotope("he4") == Isotope(2, 4, 0)
    @test Isotope("He-4") == Isotope(2, 4, 0)
    @test Isotope("4He") == Isotope(2, 4, 0)
    @test Isotope("26Alg") == Isotope(13, 26, 0)
    @test Isotope("26Alm") == Isotope(13, 26, 1)
    @test Isotope("PROT") == Isotope(1, 1, 0)
    @test Isotope("p") == Isotope(1, 1, 0)
    @test Isotope("NEUT") == Isotope(0, 1, 0)
    @test tryparse(Isotope, "OOOOO") === nothing
    @test_throws ArgumentError Isotope("not-an-isotope!!")

    @test isotope_name(Isotope(0, 1, 0)) == "n"
    @test isotope_name(Isotope(1, 1, 0)) == "p"
    @test isotope_name(Isotope(13, 26, 1)) == "Al-26m"
    @test mass_number(Isotope(13, 26, 0)) == 26
    @test neutron_number(Isotope(13, 26, 0)) == 13
    @test latex_label(Isotope(13, 26, 0)) isa LaTeXString

    @test is_stable(Isotope(6, 12, 0))
    @test !is_stable(Isotope(13, 26, 0))
    @test !is_stable(Isotope(6, 12, 1))  # isomers are never stable
    @test is_magic(20)
    @test !is_magic(21)
end

@testset "io_ppn: abundances" begin
    ab = read_abundances(joinpath(NUPPN, "iso_massf00020.DAT"))
    @test ab.cycle == 20
    @test ab.t9 > 0
    @test ab.rho > 0
    @test length(ab) == length(ab.iso)
    @test ab[Isotope(2, 4, 0)] > 0          # He-4 present
    @test ab["He-4"] == ab[Isotope(2, 4, 0)]
    @test ab[Isotope(999, 999, 0)] == 0.0   # absent species -> 0.0
    @test mass_fraction(ab, "He-4") == ab["He-4"]

    df = DataFrame(ab)
    @test Set(names(df)) == Set(["Z", "A", "N", "isomer", "isotope", "X"])
    @test nrow(df) == length(ab)

    decay = read_abundances(joinpath(NUPPN, "iso_massfdecay.DAT"))
    @test length(decay) > 0
end

@testset "io_ppn: fluxes" begin
    fx = read_fluxes(joinpath(NUPPN, "flux_00020.DAT"))
    @test nrow(fx) > 0
    @test Set(["index", "z_start", "a_start", "n_start", "z_end", "a_end", "n_end",
               "reactant", "product", "flux", "energy", "timescale"]) ⊆ Set(names(fx))
    @test all(fx.flux .>= 0)
end

@testset "io_ppn: xtime" begin
    xt = read_xtime(joinpath(NUPPN, "x-time.dat"))
    @test length(xt) == length(xt.cycle)
    @test length(xt) > 1
    xv, yv = series(xt, "He-4"; x = :time)
    @test length(xv) == length(yv) == length(xt)
    @test xv == xt.time
    xv9, _ = series(xt, "He-4"; x = :t9)
    @test xv9 == xt.t9
    @test_throws ArgumentError series(xt, "He-4"; x = :bogus)
    @test all(==(0.0), xt[Isotope(999, 999, 0)])  # absent species -> zeros
end

@testset "io_ppn: iniab and trajectory" begin
    ab = read_iniab(joinpath(DATA, "iniab_test.dat"))
    @test ab.cycle == -1
    @test isnan(ab.time)
    @test ab["p"] ≈ 7.06706E-01
    @test ab["n"] ≈ 1.00143E-99
    @test ab["He-4"] ≈ 2.75150E-01
    @test ab["C-12"] ≈ 3.44311E-03

    traj = read_trajectory(joinpath(DATA, "trajectory_test.dat"))
    @test traj.time_s == [0.0, 1.0, 2.0]
    @test traj.temperature_T9 == [0.20, 0.15, 0.10]
    @test traj.density_cgs == [100.0, 90.0, 80.0]

    traj_units = read_trajectory(joinpath(DATA, "trajectory_test_units.dat"))
    @test traj_units.time_s[2] ≈ 365.25 * 24 * 3600           # 1 YRS -> s
    @test traj_units.temperature_T9 == [2.0, 1.0]              # T8K -> T9K
    @test traj_units.density_cgs == [100.0, 1000.0]            # LOG -> CGS
end

@testset "io_network" begin
    net = read_network(joinpath(NUPPN, "networksetup.txt"))
    @test length(net.isotopes) > 0
    @test length(net.reactions) > 0
    @test any(r -> r.active, net.reactions)

    r = first(net.reactions)
    @test r.rtype isa String
    @test r.source isa String
    @test label(r) isa String
    @test latex_label(r) isa LaTeXString

    physics = read_namelist(joinpath(NUPPN, "ppn_physics.input"))
    @test physics[:NVCP] == 55
    @test physics[:NRCP] == 117
    @test physics[:ININET] == 0
end

@testset "io_isotopedatabase" begin
    db = read_isotopedatabase(joinpath(NUPPN, "isotopedatabase.txt"))
    @test length(db) > 0
    neutron_row = first(db)
    @test neutron_row.iso == Isotope(0, 1, 0)
    @test neutron_row.active
    @test all(e -> e.stable_a >= 0, db)
end

@testset "PPNRun" begin
    run = PPNRun(NUPPN)
    @test !isempty(run.cycles)
    @test first(run.cycles) == 0
    @test last(run.cycles) == 804

    ab_final = abundances(run, :final)
    @test ab_final.cycle == last(run.cycles)
    ab_initial = abundances(run, :initial)
    @test ab_initial.cycle == first(run.cycles)
    @test abundances(run, 20).cycle == 20
    @test_throws ArgumentError abundances(run, 999999)

    decay = abundances(run, :decay)
    @test length(decay) > 0

    fx = fluxes(run, :final)
    @test nrow(fx) > 0

    xt = xtime(run)
    @test length(xt) > 0

    net = network(run)
    @test length(net.reactions) > 0

    db = isotopedatabase(run)
    @test length(db) > 0

    in_ = inputs(run)
    @test in_.physics[:NVCP] == 55
end

@testset "PPNSweep" begin
    sweep = PPNSweep(joinpath(DATA, "sweep"))
    @test reactions(sweep) == ["13N_pg_14O"]
    @test factors(sweep, "13N_pg_14O") == [0.5, 2.0]

    lo = sweep_run(sweep, "13N_pg_14O", 0.5)
    hi = sweep_run(sweep, "13N_pg_14O", 2.0)
    @test abundances(lo, :final)["He-4"] == 0.2
    @test abundances(hi, :final)["He-4"] == 0.45
    @test abundances(sweep.baseline, :final)["He-4"] ≈ 0.315401

    @test_throws ArgumentError sweep_run(sweep, "13N_pg_14O", 999.0)
    @test_throws ArgumentError PPNSweep(joinpath(DATA, "no-such-sweep-dir"))
end

@testset "PPNEnsemble" begin
    ens = PPNEnsemble(joinpath(DATA, "ensemble"))
    @test length(ens.samples) == 2
    @test [s.id for s in ens.samples] == [1, 2]

    s1 = sample(ens, 1)
    @test s1.manifest["seed"] == 12345
    @test s1.manifest["rates"][1]["factor"] ≈ 0.82
    @test abundances(s1.run, :final)["He-4"] == 0.28

    s2 = sample(ens, 2)
    @test abundances(s2.run, :final)["He-4"] == 0.36

    @test abundances(ens.baseline, :final)["He-4"] ≈ 0.315401
    @test_throws ArgumentError sample(ens, 999)
    @test_throws ArgumentError PPNEnsemble(joinpath(DATA, "no-such-ensemble-dir"))
end

@testset "charts smoke test" begin
    run = PPNRun(NUPPN)

    @test abundance_chart(run, :final) isa CM.Figure
    @test abundance_chart(abundances(run, :final)) isa CM.Figure

    @test flux_chart(run, :final) isa CM.Figure
    @test flux_chart(fluxes(run, :final), abundances(run, :final); show_abundance = false) isa CM.Figure

    ab1, ab2 = abundances(run, :final), abundances(run, :initial)
    @test ratio_chart(ab1, ab2) isa CM.Figure
    changed = changed_isotopes(ab1, ab2)
    @test nrow(changed) > 0
    @test issorted(abs.(changed.log_ratio); rev = true)

    @test abundance_vs_time(run, ["He-4", "C-12"]) isa CM.Figure
    @test abundance_vs_time(run, Isotope(2, 4, 0)) isa CM.Figure
end

@testset "reaction_lookup" begin
    run = PPNRun(NUPPN)
    net = network(run)
    he4_reactions = reactions_for_isotope(net, Isotope(2, 4, 0))
    @test !isempty(he4_reactions)
    @test all(r -> r.active, he4_reactions)
    @test all(r -> Isotope(2, 4, 0) in r.reactants || Isotope(2, 4, 0) in r.products, he4_reactions)

    all_reactions = reactions_for_isotope(net, Isotope(2, 4, 0); active_only = false)
    @test length(all_reactions) >= length(he4_reactions)

    @test describe_rate(first(he4_reactions)) isa String
end

end # testset "NuGridJl"
