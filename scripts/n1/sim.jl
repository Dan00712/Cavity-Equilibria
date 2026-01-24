using DrWatson
HOME = ENV["HOME"]
DrWatson.quickactivate(
	joinpath(HOME, "projects", "bac", "CavityEquilibria"),
	"CavityEquilibria")

using Random

using JLD2
using Plots
using ProgressBars

using CavityEquilibria
using CavityEquilibria.Parameters
using CavityEquilibria.Laser
using CavityEquilibria.RootFinding
using CavityEquilibria.Util

const params = DEFAULT_PARAMS
const Ω = range(1, 400, length = 100) .* 1e3*2π ;

const N = 1
const sz = let
	VALS = 25
	Random.seed!(0)

	#range is [-2, 2]
	logsz = rand(VALS) .* 4 .- 2
	1e-6 .* exp10.(logsz)
end
const vz = Iterators.product([sz for _ in 1:N]...)
const d = 0
const κ = 18e4 * 2π
const vd = [0]
const vϕ = [0]

ω, z = let
	ω = []
	z = []
    for Δ in ProgressBar(Ω)
		for zg in vz
			try
				r = find_roots(collect(zg), vd, vϕ; Δ=Δ, κ=κ, params=params)

				if (length(z) == 0 || !any(v -> isapprox(v, r), z)) && all(abs.(r) .< 1e3)
					push!(z, r)
					push!(ω, Δ)
				end
			catch err
				if !(err isa ConvergenceError)
					rethrow(err)
				end
			end
		end
	end

	z = Iterators.reduce(hcat, z)

	ω, z
end

let
	p = joinpath(
		datadir("P1") |> mkpath,
		"latest.jld2"
	)
	M = (ω=ω, z=z)
	@save p M
end

p = plot(;
	 xlabel="Δ/ 2π kHz",
	 ylabel="z/μm",
	 yscale=:log,

)

scatter!(p,
		 ω,
		 z[1, :]
)
savefig(
	p,
	joinpath(mkpath(plotsdir("n1")), "latest.pdf")
)
savefig(
	p,
	joinpath(mkpath(plotsdir("n1")), "$(now_nodots()).png")
)

if isinteractive()
    gui(p)
end

