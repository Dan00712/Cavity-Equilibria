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
	VALS = 15
	Random.seed!(0)

	#range is [-2, 2]
	logsz = rand(VALS) .* 4 .- 2
    1e-6 .* vcat(exp10.(logsz), -1 .* exp10.(logsz))
end
const vz = Iterators.product([sz for _ in 1:N]...)
const d = params.R * 10
const κ = 18e4 * 2π
const vd = [-d/2, d/2]
const vϕ = [0, 0]

ω, z = let
	ω = []
	z = []
    l = ReentrantLock()
    Threads.@threads for Δ in ProgressBar(Ω)
		for zg in vz
			try
				r = find_roots(collect(zg), vd, vϕ; Δ=Δ, κ=κ, params=params)

				if (length(z) == 0 || !any(v -> isapprox(v, r), z)) && all(abs.(r) .< 1)
                    lock(l) do 
					    push!(z, r)
    					push!(ω, Δ)
                    end
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

const PREFIX="n2"
let
	p = joinpath(
		datadir(PREFIX) |> mkpath,
		"latest.jld2"
	)
	M = (ω=ω, z=z)
	@save p M
end

p = plot(;
	 zlabel="Δ/ 2π kHz",
	 xlabel="z1/μm",
	 ylabel="z2/μm",
	 #yscale=:log,
	 #xscale=:log,
)

scatter!(p,
		 z[1, :] .* 1e6, # x
         z[2, :] .* 1e6, # y
		 ω;
         markersize=5
)
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

