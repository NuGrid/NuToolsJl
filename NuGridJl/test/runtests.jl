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

    traj = trajectory(run)
    @test nrow(traj) > 0
    @test traj.time_s[1] == 0.0
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

@testset "npdata" begin
    test_reaction = Reaction(1, true, [Isotope(1, 1, 0), Isotope(6, 12, 0)], [Isotope(7, 13, 0)],
                              "(p,g)", "NACRR", 1.0, 1.0, 1.943, 4)
    npdata_root = joinpath(DATA, "npdata_test")

    @test reaclib_species(Isotope(1, 1, 0)) == "p"
    @test reaclib_species(Isotope(0, 1, 0)) == "n"
    @test reaclib_species(Isotope(6, 12, 0)) == "c12"

    # the decoy file matches the winvn skip-pattern and carries a rate=1 block
    # for the same reaction; if it weren't skipped, this test would see two
    # (label, file) groups instead of one and the wrong rate.
    curves = rate_curve(test_reaction, npdata_root)
    @test Set(curves.file) == Set(["synthetic_reaclib.txt"])
    @test all(curves.rate .≈ exp(10.0))
    @test all(curves.label .== "nacrr")
    @test all(curves.qvalue .== 1.943)

    # explicit `files` bypasses the skip filter
    curves_all = rate_curve(test_reaction, npdata_root; files = ["synthetic_reaclib.txt", "winvn_decoy.dat"])
    @test Set(curves_all.file) == Set(["synthetic_reaclib.txt", "winvn_decoy.dat"])

    no_match = Reaction(2, true, [Isotope(2, 4, 0)], [Isotope(6, 12, 0)], "(a,g)", "X", 1.0, 1.0, 0.0, 3)
    @test_throws ArgumentError rate_curve(no_match, npdata_root)

    fig = rate_plot(test_reaction, npdata_root)
    @test fig isa CM.Figure

    @test reaclib_rate((10.0, 0, 0, 0, 0, 0, 0), 1.0) ≈ exp(10.0)
    @test_throws DomainError reaclib_rate((0.0, 0, 0, 0, 0, 0, 0), -1.0)
end

@testset "tables" begin
    df = DataFrame(reaction = ["a", "b"], ratio = [1.23456, missing])
    md = dataframe_to_markdown(df)
    @test occursin("| reaction | ratio |", md)
    @test occursin("1.2346", md)  # rounded to 4 digits

    html = dataframe_to_html(df)
    @test html isa RenderedHTML
    @test occursin("<table>", html.html)
    @test occursin("1.2346", html.html)

    out_path = joinpath(mktempdir(), "table.csv")
    result = save_table(df, out_path)
    @test result === df
    @test isfile(out_path)
end

@testset "sensitivity_iliadis" begin
    sweep = PPNSweep(joinpath(DATA, "sweep"))

    sens = sensitivity(sweep, "13N_pg_14O", "He-4")
    @test sens[0.5] ≈ 0.2 / 0.315401
    @test sens[2.0] ≈ 0.45 / 0.315401

    table = sensitivity_table(sweep, "He-4")
    @test nrow(table) == 2
    @test Set(table.factor) == Set([0.5, 2.0])
    @test all(==("13N_pg_14O"), table.reaction)
    @test all(==("He-4"), table.isotope)

    wide = iliadis_table(sweep, "He-4")
    @test nrow(wide) == 1
    @test Set(names(wide)) == Set(["reaction", "isotope", "0.5", "2.0"])

    ranked = rank_reactions(table, "He-4")
    @test ranked.reaction == ["13N_pg_14O"]
    @test ranked.score[1] ≈ maximum(abs.(log10.([0.2 / 0.315401, 0.45 / 0.315401])))
    @test_throws ArgumentError rank_reactions(table, "He-4"; metric = :bogus)
    @test_throws ArgumentError rank_reactions(table, "C-12")  # not in this table

    fig = sensitivity_plot(sweep, "13N_pg_14O", "He-4")
    @test fig isa CM.Figure
end

@testset "reaction_report" begin
    run = PPNRun(NUPPN)

    flux_list = flux_reaction_list(run; cycle = :final)
    @test nrow(flux_list) > 0
    @test issorted(flux_list.flux; rev = true)
    @test all(flux_list.flux .>= 1e-60)
    @test Set(names(flux_list)) == Set(["index", "reaction", "source", "rtype", "active", "flux", "rate"])

    tight_list = flux_reaction_list(run; cycle = :final, threshold = 1.0)  # nothing should pass
    @test nrow(tight_list) == 0

    # whole-trajectory form: peak flux across a set of cycles, not just one.
    # run.cycles[1] == 0 has no flux_00000.DAT (dY/dt isn't meaningful before
    # any integration step) -- that must be skipped, not thrown on, so this
    # deliberately includes it to exercise that.
    @test run.cycles[1] == 0
    @test !isfile(joinpath(NUPPN, "flux_00000.DAT"))
    traj_list = flux_reaction_list(run, run.cycles[1:6])
    @test nrow(traj_list) > 0
    @test issorted(traj_list.flux; rev = true)
    @test Set(names(traj_list)) == Set(["index", "reaction", "source", "rtype", "active", "flux", "peak_cycle", "rate"])
    @test all(c -> c in run.cycles[1:6], traj_list.peak_cycle)
    @test all(c -> c != 0, traj_list.peak_cycle)  # cycle 0 was skipped, so it can never be the peak
    # scanning more cycles can only add reactions or raise peak flux, never remove either
    single = flux_reaction_list(run; cycle = run.cycles[2])
    @test Set(single.index) ⊆ Set(traj_list.index)

    sweep = PPNSweep(joinpath(DATA, "sweep"))
    table = sensitivity_table(sweep, "He-4")
    report = sensitivity_reaction_report(table)
    @test nrow(report) == 1
    @test report.reaction[1] == "13N_pg_14O"
    @test report.n_factors[1] == 2
    @test report.n_isotopes[1] == 1
    @test report.worst_isotope[1] == "He-4"
    @test report.sensitive[1] == true  # ratios of ~0.634/~1.427 give |log10| up to ~0.2, above default threshold=0.1

    strict_report = sensitivity_reaction_report(table; threshold = 10.0)
    @test strict_report.sensitive[1] == false
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

    # a blank tile must only ever mean "not tracked" -- an impossibly strict
    # tolerance/threshold puts every tracked isotope below the color/flux
    # floor, but every one of them must still get a tile+label, so none of
    # these should throw the old "nothing above tolerance" errors anymore.
    @test abundance_chart(run, :final; tolerance = 1.0) isa CM.Figure
    @test flux_chart(run, :final; tolerance = 1.0) isa CM.Figure
    @test ratio_chart(ab1, ab2; tolerance = 1.0) isa CM.Figure

    @test plot_trajectory(run) isa CM.Figure
    @test plot_density_temperature(run) isa CM.Figure
    @test plot_trajectory(trajectory(run)) isa CM.Figure
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

    # flux-aware form: structural candidates narrowed to ones that actually
    # carried flux, single-cycle and whole-trajectory
    he4_flux_final = reactions_for_isotope(run, Isotope(2, 4, 0), :final)
    @test nrow(he4_flux_final) > 0
    @test issorted(he4_flux_final.flux; rev = true)
    structural_indices = Set(r.index for r in he4_reactions)
    @test Set(he4_flux_final.index) ⊆ structural_indices

    he4_flux_traj = reactions_for_isotope(run, Isotope(2, 4, 0), run.cycles[1:6])
    @test "peak_cycle" in names(he4_flux_traj)
    @test Set(he4_flux_traj.index) ⊆ structural_indices
end

end # testset "NuGridJl"
