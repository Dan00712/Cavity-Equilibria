#!/usr/bin/env julia
using DrWatson
HOME = ENV["HOME"]
DrWatson.quickactivate(
	joinpath(HOME, "projects", "bac", "CavityEquilibria"),
	"CavityEquilibria")

using Random

using JLD2
using Plots
if isinteractive()
    plotlyjs()
end
using ProgressBars

using CavityEquilibria
using CavityEquilibria.Parameters
using CavityEquilibria.Laser
using CavityEquilibria.RootFinding
using CavityEquilibria.Util

const params = DEFAULT_PARAMS
const Ω = range(1, 400, length = 100) .* 1e3*2π ;

const N = 2
const sz = let
	VALS = 20
	Random.seed!(0)

	#range is [-2, 2]
	logsz = rand(VALS) .* 4 .- 2
    1e-6 .* vcat(exp10.(logsz), -1 .* exp10.(logsz))
end
const vz = Iterators.product([sz for _ in 1:N]...)
const d = params.R * 100
const κ = 18e4 * 2π
const vd = [-d/2, d/2]
const d = Dict()
l = ReentrantLock()

function foo(vϕ, p, style)
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

    lock(l) do
        d[vϕ[2]] = (Δ=ω, z=z)
        scatter!(p,
		    z[1, :], # x
            z[2, :], #y
		    ω;
            style...
        )
    end

	ω, z
end

const PREFIX="n2-phi"
let
	p = joinpath(
		datadir(PREFIX) |> mkpath,
		"latest.jld2"
	)
	@save p d
end

p = plot(;
	 zlabel="Δ/ 2π kHz",
	 xlabel="z1/μm",
	 ylabel="z2/μm",
)

for (vϕ, style) in zip(
                       [[0,0],
                        [0, π/2],
                        [0, π],
                        [0, -π/2]],
                       [Dict(:color=>:blue, :label=>"ϕ=0"),
                        Dict(:color=>:red, :label=>"ϕ=π/2"),
                        Dict(:color=>:orange, :label=>"ϕ=π"),
                        Dict(:color=>:green, :label=>"ϕ=-π/2")]
                      )
    foo(vϕ, p, style)
end

savefig(
	p,
	joinpath(mkpath(plotsdir(PREFIX)), "latest.html")
)
savefig(
	p,
	joinpath(mkpath(plotsdir(PREFIX)), "latest.pdf")
)
savefig(
	p,
	joinpath(mkpath(plotsdir(PREFIX)), "$(now_nodots()).png")
)

if isinteractive()
    gui(p)
end

