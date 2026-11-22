using MC, Test, Random, Printf, Statistics

# A stand-in for ppn.exe: reads the seed, draws x for three "rates", and writes a 2-row x-time.dat in
# which C12 depends on the first rate only. Seed 3 fails.
const FAKE_PPN = raw"""#!/bin/bash
seed=$(grep -i 'factor_uncertainty_mc_seed' ppn_physics.input | sed 's/.*= *//')
[ "$seed" = 3 ] && { echo "solver failed"; exit 1; }
awk -v seed="$seed" 'BEGIN {
  srand(seed);
  for (i = 1; i <= 3; i++) { u = rand(); v = rand(); x[i] = sqrt(-2 * log(u + 1e-300)) * cos(6.283185307 * v) }
  f = "rate_uncertainty_x.txt";
  print "# header" > f;
  print "#   row  label  reaction                                                x  how" > f;
  printf "%7d  %5s  %-46s%15.6E  %s\n", 11, "STL02", "1 C  12 + 1 PROT  -> 1 N  13 + 0 OOOOO", x[1], "drawn" > f;
  printf "%7d  %5s  %-46s%15.6E  %s\n", 22, "NAC2R", "1 C  13 + 1 PROT  -> 1 N  14 + 0 OOOOO", x[2], "drawn" > f;
  printf "%7d  %5s  %-46s%15.6E  %s\n", 33, "JINAC", "1 AL 26 + 1 PROT  -> 1 SI 27 + 0 OOOOO", x[3], "drawn default_fu" > f;
  g = "x-time.dat";
  print "#|cycle |time         |  1-C  12  |  2-C  13  |  3-AL 26  |  4-AL26M  " > g;
  printf "%7d %12.5E %12.5E %12.5E %12.5E %12.5E\n", 1, 0.0, 0.1, 0.05, 0.0, 0.0 > g;
  printf "%7d %12.5E %12.5E %12.5E %12.5E %12.5E\n", 2, 1.0, 0.1 * 10^(0.5 * x[1]), 0.05, 1e-5, 2e-5 > g;
}'
echo "simulation complete"
"""


@testset "MC" begin

@testset "parsing" begin
    @test fnum("1.5E-03") == 1.5e-3
    @test fnum("2.24366-100") == 2.24366e-100
    @test isnan(fnum("garbage"))
    @test iso_name("38-AL 26") == "Al26"
    @test iso_name("1103-AL26M") == "Al26m"
    @test iso_name("  2-PROT  ") == "H1"
    @test iso_name("cycle") === nothing
    @test parse_seeds("1-3,7, 9-10") == [1, 2, 3, 7, 9, 10]
    txt = "&ppn_physics\n\tfoo = 1\n\tfactor_uncertainty_mc = .false.\n/\n"
    out = physics_with_seed(txt, 42)
    @test occursin("factor_uncertainty_mc_seed = 42", out)
    @test count("factor_uncertainty_mc =", out) == 1
    @test findfirst("seed = 42", out)[1] < findlast("/", out)[1]
    @test_throws ErrorException physics_with_seed("&ppn_physics\n", 1)
end

@testset "statistics" begin
    rng = MersenneTwister(1)
    @test ranks([3.0, 1.0, 3.0, 2.0]) == [3.5, 1.0, 3.5, 2.0]
    x = randn(rng, 200)
    @test spearman(x, exp.(x)) ≈ 1
    @test spearman(x, -x) ≈ -1
    # Gaussian pair with correlation r: I = -log(1 - r^2)/2
    n, r = 5000, 0.8
    a, b = randn(rng, n), randn(rng, n)
    y = r .* a .+ sqrt(1 - r^2) .* b
    @test isapprox(ksg_mi(a, y), -0.5 * log(1 - r^2); atol=0.03)
    @test ksg_mi(a, b) < 0.02
    @test ksg_mi(a, y) ≈ ksg_mi(a, 3 .* y .+ 2) atol = 1e-6     # variables are standardized first
    s = spread(exp.(log(2.0) .* randn(rng, 20000)))
    @test isapprox(s.factor, 2.0; rtol=0.03)
    @test isapprox(s.median, 1.0; rtol=0.03)
    xr = vcat(a', b')
    @test spearman_all(xr, y) ≈ [spearman(a, y), spearman(b, y)]
end


@testset "run, collect, analyze" begin
    root = mktempdir()
    mkpath(joinpath(root, "NPDATA"))
    tpl = joinpath(root, "ppn")
    mkpath(tpl)
    write(joinpath(tpl, "ppn.exe"), FAKE_PPN); chmod(joinpath(tpl, "ppn.exe"), 0o755)
    write(joinpath(tpl, "ppn_frame.input"), "&ppn_frame\n trajectory_fn = 'traj.input'\n ini_filename = 'ini.dat'\n/\n")
    write(joinpath(tpl, "ppn_physics.input"), "&ppn_physics\n\tstarlib_option = 2\n/\n")
    write(joinpath(tpl, "ppn_solver.input"), "&ppn_solver\n/\n")
    write(joinpath(tpl, "traj.input"), "0 1 1\n")
    write(joinpath(tpl, "ini.dat"), "  6 C  12  1.0E-01\n  6 C  13  5.0E-02\n 13 AL 27  1.0E-03\n")

    out = joinpath(root, "out")
    @test run_mc(tpl, out, 1:30; jobs=4, io=devnull) == (29, 1)
    @test readlines(joinpath(out, "failed.txt")) == ["3\texit 1"]
    @test isfile(joinpath(out, "logs", "s000003.log"))
    @test !isdir(joinpath(out, "work"))
    @test run_mc(tpl, out, 1:30; jobs=2, io=devnull) == (0, 1)     # restart: only the failed seed is retried

    @test collect_mc(out; io=devnull) == 29
    ens = MCEnsemble(out)
    @test nsamples(ens) == 29 && length(ens.rates) == 3
    @test ens.seeds == setdiff(1:30, 3)
    @test ens.rates[3] == (row=33, label="JINAC", reaction="1 AL 26 + 1 PROT -> 1 SI 27 + 0 OOOOO", how="default_fu")
    @test uncertainty_source(ens.rates[3].how) == "default f.u."
    @test all(abundance(ens, "Al26") .≈ 3e-5)                     # ground + isomer
    @test all(abundance(ens, "Al26"; isomers=false) .≈ 1e-5)
    @test all(iszero, abundance(ens, "Fe56"))
    r = ratio(ens, "C12", "C13")
    @test r ≈ 2 .* 10 .^ (0.5 .* ens.xr[1, :]) rtol = 1e-4
    X0 = read_iniab(joinpath(tpl, "ini.dat"))
    @test X0["C13"] == 0.05
    @test ratio(ens, "C12", "C13"; X0, xi=9) ≈ (abundance(ens, "C12") .+ 0.9) ./ (0.05 + 0.45)

    rk = rank_reactions(ens, log10.(r))
    @test rk[1].row == 11 && rk[1].spearman ≈ 1
    @test rk[1].mi > 1
    @test_throws DimensionMismatch rank_reactions(ens, [1.0])

    res = sensitivity_report(ens; pairs=[("C12", "C13")], isotopes=["Al26"], X0, xi=9, io=devnull)
    @test res["C12/C13"][1].row == 11
    for f in ("mc_ratios.txt", "mc_sensitivity.txt", "mc_sensitivity.csv")
        @test isfile(joinpath(out, f))
    end
    @test startswith(readlines(joinpath(out, "mc_sensitivity.csv"))[2], "C12/C13,1,11,STL02,")

    # single-run baseline table on the last x-time.dat row of one worker-style run directory
    run1 = joinpath(root, "single")
    cp(tpl, run1)
    cd(() -> run(pipeline(`./ppn.exe`; stdout=devnull)), run1)
    b = ward_baseline([run1]; model=:co, io=devnull)
    @test b[run1][("Al26", "Al27")][1] |> isnan                       # no Al27 in the fake network
    @test b[run1][("C12", "C13")][2] ≈ (final_abundances(run1)["C12"] + 0.9) / (0.05 + 0.45)
end

end
